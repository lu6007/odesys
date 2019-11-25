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
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 

    % model_name = 'opt_model';
    para_name = {'verbose'};
    default_value = {1};
    verbose = parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction opt_model_1118(): model_name = %s\n', model_name);
    end 
    
    ode_model = ode_model_1118('model_1', 'verbose', 0);
    model_obj.ode = ode_model; % data, init_data, rhs, output
    
    % [Ouyang et al 2019 ACS Sens; 350% activation corresponds to 100%]
    model_obj.scale = 57.1429;       % nM/au ratio = 200 nM /3.5 au; 0.328 * scale = 18.7429
    model_obj.index = 3;             % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
    
    % Optimization functions and parameters
    model_obj.objective = @objective;
    model_obj.constraint = fh.constraint;
    model_obj.initial_guess = fh.initial_guess; 
%     model_obj.theta_name = {'kon_1'; 'koff_1'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; 'kon_4'; 'koff_4'};
%     model_obj.theta_upper_bound = [4.5; 4.5; 38; 10; 1; 2; 100]; 
%     model_obj.theta_name = {'kon_1'; 'koff_1'; 'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
%     model_obj.theta_upper_bound = [4.5; 4.5; 0.1; 12; 1.0; 0.2; 2; 10]; 
    model_obj.theta_name = {'kon_1'; 'koff_1'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
    model_obj.theta_upper_bound = [4.5; 4.5; 38; 10; 1; 1.0; 0.2; 2; 10]; 

end % function model_obj = fyn_gf_model(model_name)

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
    objective(theta)
    
    model_fh = @opt_model_1118;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh);
    
end





