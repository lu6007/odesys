% A standard optimization model function 
% function model_obj = fyn_gf_model(model_name, varargin)
%     para_name = {'multiple_output', 'best_fit', 'verbose', 'fyn_endo'};
%     default_value = {0, 1, 1, 0};

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com); 
function model_obj = fyn_gf_model(model_name, varargin)
    global optimize_ode_utility_fh;
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
    model_obj.scale = 496.5810; 
    
    % Optimization functions and parameters
    model_obj.objective = @objective;
    model_obj.constraint = fh.constraint;
    model_obj.initial_guess = fh.initial_guess; 
    model_obj.index = 3; % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
    model_obj.theta_name = {'kon_4'; 'koff_4'; 'kcaton_7'; 'kdon_7'; ...
        'kcatoff_7'; 'kdoff_7'};
    model_obj.theta_upper_bound = [10; 1; 1; 500; 1; 500]; 

end % function model_obj = fyn_gf_model(model_name)

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
    objective(theta)
    
    model_fh = @fyn_gf_model;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh);
end






