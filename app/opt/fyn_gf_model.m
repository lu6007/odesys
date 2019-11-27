% A standard optimization model function 
% function model_obj = fyn_gf_model(model_name, varargin)
%     para_name = {'multiple_output', 'best_fit', 'verbose', 'fyn_endo'};
%     default_value = {0, 1, 1, 0};

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com); 
function model_obj = fyn_gf_model(model_name, varargin)
    % The global variable optimize_ode_model is necessary since the
    % objective() function cannot take input other than theta. 
    global optimize_ode_utility_fh optimize_ode_model;
    fh = optimize_ode_utility_fh; 

    % model_name = 'complex_ode';
    para_name = {'multiple_output', 'best_fit', 'verbose', 'fyn_endo'};
    default_value = {0, 1, 1, 0};
    [multiple_output, best_fit, verbose, fyn_endo] = parse_parameter(para_name, ...
        default_value, varargin);
    
    if verbose
        fprintf('\nFunction fyn_gf_model(): model_name = %s\n', model_name);
        fprintf('multiple_output = %d, best_fit = %d, fyn_endo = %d\n', ...
            multiple_output, best_fit, fyn_endo); 
    end 
    
    ode = complex_ode(model_name, 'multiple_output', multiple_output, ...
        'best_fit', best_fit, 'fyn_endo', fyn_endo, 'verbose', 0);     
    model_obj.ode = ode; 
    switch fyn_endo
        case 0 
            model_obj.scale = 496.5810; 
        case {1, 2, 3, 4, 5}
            model_obj.scale = 593.2; % 194.2/0.3274
        case {11, 12, 13, 14, 15}
            model_obj.scale = 262.4007; % Initial guess, 85.91/0.3274
    end
    
    % Optimization functions and parameters
    model_obj.objective = @objective;
    model_obj.constraint = fh.constraint;
    model_obj.initial_guess = fh.initial_guess; 
    model_obj.index = 3; % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
    model_obj.theta_name = {'kon_4'; 'koff_4'; 'kcaton_7'; 'kdon_7'; ...
        'kcatoff_7'; 'kdoff_7'};
    % model_obj.theta_upper_bound = [10; 1; 1; 500; 1; 500]; 
    model_obj.theta_bound = fh.get_theta_bound(model_obj, model_obj.theta_name);
    
    optimize_ode_model = model_obj; 

end % function model_obj = fyn_gf_model(model_name)

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
    objective(theta)
    
    model_fh = @fyn_gf_model;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh);
end







