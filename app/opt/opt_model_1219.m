% A standard model function 
% function model = opt_model_1219(model_name, varargin)
%     para_name = {'verbose', 'theta_bound_factor', 'model_id'};
%     default_value = {1, 10, 1};
% model.data = init_data(); 
% model.init_data = @init_data; 
% model.rhs = @rhs; 
% model.output = @output; 

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com)
function model = opt_model_1219(model_name, varargin)
    global optimize_ode_utility_fh; 
    fh = optimize_ode_utility_fh; 

    % model_name = 'opt_model';
    para_name = {'verbose', 'theta_bound_factor', 'model_id'};
    default_value = {1, 10, 1};
    [verbose, theta_bound_factor, model_id] = ...
        parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction opt_model_1118(): model_name = %s\n', model_name);
    end 
    
    model = opt_model_1118(model_name, varargin);
    ode_model = ode_model_1118('model_1', 'verbose', 0, 'ode_id', 2);
    model.ode = ode_model; % data, init_data, rhs, output
    
    % Values inherited from 'model1118'
%     % [Ouyang et al 2019 ACS Sens; 350% activation corresponds to 100%]
%     model.scale = 57.1429;       % nM/au ratio = 200 nM /3.5 au; 0.328 * scale = 18.7429
%     model.index = 3;             % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
%     
%     % Optimization functions and parameters
%     model.objective = fh.objective;
%     % model.constraint = fh.constraint;
%     model.initial_guess = fh.initial_guess; 

    if model_id >=1 
            model.theta_name = {'kon_1'; 'koff_1'; 'kcaton_3'; 'kdon_3'; 'kcatoff_3'; 'kdoff_3';
                'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
    end

    model.theta_bound = fh.get_theta_bound(model, model.theta_name, ...
        'bound_factor', theta_bound_factor);
end % function model = opt_model_1219(model_name, varargin)






