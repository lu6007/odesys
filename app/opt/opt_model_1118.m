% A standard model function 
% function model_obj = fyn_gf_model(model_name, varargin)
%     para_name = {'multiple_output', 'best_fit'};
%     default_value = {0, 0};
% model_obj.data = init_data(); 
% model_obj.init_data = @init_data; 
% model_obj.rhs = @rhs; 
% model_obj.output = @output; 

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com); Tongkai Li
% (ltk@pku.edu.cn)
function model_obj = opt_model_1118(model_name, varargin)
%     global optimize_ode_utility_fh;
%     fh = optimize_ode_utility_fh; 

    % model_name = 'opt_model';
    para_name = {'verbose'};
    default_value = {1};
    verbose = parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction opt_model_1118(): model_name = %s\n', model_name);
    end 
    
    ode_model = ode_model_1118('model_1');
    model_obj.ode = ode_model; % data, init_data, rhs, output
    
    % [Ouyang et al 2019 ACS Sens; 350% activation corresponds to 100%]
    model_obj.scale = 148.56;       % nM/au ratio 
    model_obj.index = 3; % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
    
    % Optimization functions and parameters
    model_obj.objective = @objective;
    model_obj.constraint = @constraint;
    model_obj.initial_guess = @initial_guess; 
    model_obj.theta_name = {'kon_4'; 'koff_4'; 'kcaton_7'; 'kdon_7'; ...
        'kcatoff_7'; 'kdoff_7'};
    model_obj.theta_upper_bound = [10; 1; 1; 500; 1; 500]; 

end % function model_obj = fyn_gf_model(model_name)

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
    objective(theta)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 
    
    model_fh = @opt_model_1118;
    % [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh);
    
    % Update the optimized parameters 
    model = model_fh('opt_model_1118', 'verbose', 0);  
    ode_data = model.ode.data;
    theta_name = model.theta_name;
    for i = 1:length(theta_name)
        ode_data.(theta_name{i}) = theta(i);
    end

    if isfield(model, 'index') % fyn_gf_model_nodeg
        index = model.index; 
    else % earlier models. 
        index = 3;
    end
    
    % Load experimental results
    exp_name = 'exp_hela_egf'; 
    scale = model.scale; 
    [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name, 'index', index, 'scale', scale); 
    
    % Calculate simulation results
    [t_ode_temp, y_ode_temp] = batch_fyn_gf(ode_data, 'show_figure', 0, ...
        'test_basal_level', 0, 'rhs_function', model.ode.rhs, 'y0', ode_data.y0, ...
        'output_function', model.ode.output, 'verbose', 0, 'index', index);
    num_exp = size(index, 1);
    t_ode_i = t_ode_temp;
    y_ode_i = y_ode_temp; 
    clear t_ode_temp y_ode_temp; 

    % Reconcile the experiment and simulation
    [t_exp_interp, y_exp_interp, y_ode_interp] = ...
        fh.reconcile_exp_ode(t_exp_i, y_exp_i, t_ode_i, y_ode_i); 

    % Calculate error
    f = 0;
    if num_exp == 1
        temp = sqrt(fh.l2norm_square(t_exp_interp, y_exp_interp)); 
        f = fh.l2norm_square(t_exp_interp, (y_ode_interp - y_exp_interp)/temp); clear temp;
    else % num_exp > 1
        for i = 1:num_exp
            t_i = t_exp_interp{i}; 
            y_i = y_exp_interp{i};
            temp = sqrt(fh.l2norm_square(t_i, y_i)); 
            f = f+ fh.l2norm_square(t_i, (y_ode_interp{i} - y_i)/temp); clear temp;
        end
        f = f/num_exp; 
    end

end

function problem = constraint(problem, theta_var, theta_upper_bound) 
% theta_upper_bound is used in eval()

    % problem.Constraints.thetai_high = (theta_var(i)<= theta_high(i)); 
    % theta(i) <= theta_high(i)
    num_theta = size(theta_var, 1);
    for i = 1:num_theta
        cmd1 = sprintf('problem.Constraints.theta%d_high = ', i);
        cmd2 = sprintf('(theta_var(%d) <= theta_upper_bound(%d));', i, i);
        eval([cmd1 cmd2]); 
    end
end

function theta = initial_guess(x, theta_upper_bound)
theta = x.*theta_upper_bound; % lower_bound = 0; 
end




