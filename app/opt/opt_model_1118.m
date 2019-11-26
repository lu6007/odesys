% A standard model function 
% function model = opt_model_1118(model_name, varargin)
%     para_name = {'verbose'};
%     default_value = {1};
% model.data = init_data(); 
% model.init_data = @init_data; 
% model.rhs = @rhs; 
% model.output = @output; 

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com)
function model = opt_model_1118(model_name, varargin)
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
    model.ode = ode_model; % data, init_data, rhs, output
    
    % [Ouyang et al 2019 ACS Sens; 350% activation corresponds to 100%]
    model.scale = 57.1429;       % nM/au ratio = 200 nM /3.5 au; 0.328 * scale = 18.7429
    model.index = 3;             % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
    
    % Optimization functions and parameters
    model.objective = @objective;
    model.constraint = fh.constraint;
    model.initial_guess = fh.initial_guess; 
    model.theta_name = {'kon_1'; 'koff_1'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; ...
        'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
    model.theta_upper_bound = [4.5; 45; 1.25; 234.3; 38.02; ...
        0.02; 0.2; 1; 10]; 
    %
    model.theta_fit = [0.26031169 15.35232148 0.131752578 12.56740995 35.81094633 ...
        0.009446086 0.085982789 0.185889161 4.766242707]'; 
    model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
    model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'}; 
    model.theta_upper_bound = fh.get_theta_upper_bound(model, model.theta_name);
    % 
    model.theta_fit = [0.002429011	2.894388505	0.002523534	0.008028778	...
        0.004306559	0.107611631]';
    model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
    model.theta_name = {'kon_4'; 'koff_4'; 'kon_7'; 'koff_7';....
        'kcaton_8';'kdon_8'; 'kcatoff_8'; 'kdoff_8'}; 
    model.theta_upper_bound = fh.get_theta_upper_bound(model, model.theta_name);
end % function model = opt_model_1118(model_name, varargin)

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
    objective(theta)
    
    model_fh = @opt_model_1118;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh);
    
end





