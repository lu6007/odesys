% The optmization objective and constraint function for the simple ODE
% 'egfr_huang_v3'

% This is an optimization model for the simple ode defined by 
% the functions fyn_gf_rhs() and fyn_gf_init_data('egfr_huang_v3')
function model_obj = simple_ode_model(model_name, varargin)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 

    para_name = {'verbose'};
    default_value = {1};
    verbose = parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction simple_ode_model(): model_name = %s\n', model_name);
    end 
    
    data = fyn_gf_init_data('egfr_huang_v3'); 
    data.model = model_name;
    ode.data = data; 
    % rhs = @fyn_gf_rhs; % in the function fyn_gd_ode_solve()
    % output_function in the function fyn_gd_ode_solve(); 
    
    % data.scale = 428.7902; 
    model_obj.scale = 428.7902; 
    
    % Optimizatoin functions and parameters
    model_obj.ode = ode; 
    model_obj.objective = @objective;
    model_obj.constraint = fh.constraint;
    model_obj.initial_guess = fh.initial_guess; 
    model_obj.index = 3; % [egf] = 50 ng/ml
    model_obj.theta_name = {'k_caton_4'; 'k_catoff_4'; 'k_on_3'};
    model_obj.theta_upper_bound = [1; 1; 1];
end

% Objective function for the simple ode model 'egfr_huang_v3'
function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective(theta)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 

    % Update the optimized parameters 
    model = simple_ode_model('simple_ode', 'verbose', 0);
    data = model.ode.data;
    theta_name = model.theta_name; 
    % sim_name = 'egfr_huang_v3'; 
    % data = fyn_gf_init_data(sim_name);
    % theta_name = data.opt.theta_name;
    for i = 1:length(theta_name)
        data.(theta_name{i}) = theta(i);
    end

    % Load experimental results
    exp_name = 'exp_hela_egf';
    scale = model.scale; 
    [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name, 'scale', scale); 
    % [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name); 

    % Calculate simulation results
    gf_con = 7.8125; % nM for 50 ng/ml
    [t_ode_i, y_ode_i] = fyn_gf_ode_solve(data, 'gf_1', gf_con, 'show_figure', 0, ...
        'test_basal_level', 0); 

    % Reconcile the experiment and simulation
    [t_exp_interp, y_exp_interp, y_ode_interp] = ...
        fh.reconcile_exp_ode(t_exp_i, y_exp_i, t_ode_i, y_ode_i); 

    % Calculate error
    temp = sqrt(fh.l2norm_square(t_exp_interp, y_exp_interp)); 
    f = fh.l2norm_square(t_exp_interp, (y_ode_interp - y_exp_interp)/temp); clear temp;
end



