% The optimization model for odes with nodegradation. 
% function model_obj = fyn_gf_model_nodeg()
%     model_obj.ode = ode;  
%     model_obj.objective = @objective;
%     model_obj.constraint = @constraint;
%     model_obj.initial_guess = @initial_guess; 

% Author: 
% Shaoying Kathy Lu (shaoying.lu@gmail.com)
function model_obj = fyn_gf_model_nodeg(model_name, varargin)
    % The global variable optimize_ode_model is necessary since the
    % objective() function cannot take input other than theta. 
    % With this global variable, all the fyn_endo models can be simplified
    global optimize_ode_utility_fh optimize_ode_model;
    fh = optimize_ode_utility_fh; 

    para_name = {'verbose', 'best_fit', 'fyn_endo'};
    default_value = {1, 1, 0};
    [verbose, best_fit, fyn_endo] = parse_parameter(para_name, default_value, varargin);
    if verbose
        fprintf('\nFunction fyn_gf_model_nodeg(): model_name = %s\n', model_name);
        fprintf('multiple_output fixed to %d, best_fit = %d, fyn_endo = %d\n', ...
            0, best_fit, fyn_endo); 
    end 

    % index = 3; for [gf] = 50 ng/ml
    % index = (1:5)'; for all important [gfs];
    model_base = fyn_gf_model('complex_ode', 'multiple_output', 0, 'best_fit', 1, ...
        'verbose', 0, 'fyn_endo', fyn_endo);
    model_obj = model_base; 
    ode = complex_ode_nodeg(model_name, ...
        'best_fit', best_fit, 'fyn_endo', fyn_endo, 'verbose', 0);  
    model_obj.ode = ode; 
    
    % Optimization parameters
    switch best_fit % best_fit = fyn_endo
        case {2, 3}
            model_obj.scale = 597.507106591057; % Smallest error fitting scale
        case 11
            model_obj.scale = 133; % PDGF max 1.5 --> 100% of sensor = 200 nM
    end
            
    
    switch fyn_endo
        case 0 
            model_obj.objective = @objective;
            model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kcaton_7'; 'kdon_7'}; % 
            % model_obj.theta_upper_bound = [10; 1; 0.1; 0.2; 100]; %   
            model_obj.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 1 % fit data.scale
            model_obj.objective = @objective_endo;
%             model_obj.theta_name = {'kon_1'; 'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'; 'kcatoff_7'; 'kdoff_7'; 'scale'}; % 
%             model_obj.theta_upper_bound = [0.18; 5.6; 1; 0.1; 0.49; 10; 0.39; 500; 1000]; % 
            model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'; 'kcatoff_8'; 'kdoff_8'; 'scale'}; % 
            % model_obj.theta_upper_bound = [5.6; 1; 0.1; 0.49; 10; 1; 500; 1000]; % 
            model_obj.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 2 % fit others
            model_obj.objective = @objective_endo2;
%             model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4';'kon_7'; 'koff_7'}; % 
%             model_obj.theta_upper_bound = [5.6; 1; 0.1; 0.49; 10]; % 
%             model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4';'vmaxoff_6'; 'kmoff_6'}; % 
%             model_obj.theta_upper_bound = [5.6; 1; 0.1; 6.62e-1; 500]; % 
            model_obj.theta_name = {'kon_2'; 'koff_2'; 'vmaxoff_6'; 'kmoff_6'}; % 
            % model_obj.theta_upper_bound = [5.6; 1; 6.62e-1; 500]; % 
%             model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4';'vmaxoff_6'; 'kmoff_6'; 'kon_7'; 'koff_7'}; % 
%             model_obj.theta_upper_bound = [5.6; 1; 0.1; 6.62e-1; 500; 0.49; 10]; % 
            %
            model_obj.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 3 % fit others
            model_obj.objective = @objective_endo3;
%             model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'}; % 
%             model_obj.theta_upper_bound = [5.6; 1; 0.1; 0.49; 10]; % 
%             model_obj.theta_name = {'kon_1'; 'koff_1'; 'vmaxoff_6'; 'kmoff_6'; 'kcatoff_8'; 'kdoff_8'}; % 
%             model_obj.theta_upper_bound = [1; 1; 1; 500; 1; 500]; % 
            model_obj.theta_name = {'kon_1'; 'koff_1'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; 'kcatoff_8'; 'kdoff_8'}; % 
            % model_obj.theta_upper_bound = [1; 1; 50; 100; 50; 1; 500]; % 
            model_obj.index = [2;3;5;6]; % [3;4;5]; % 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 11 % fit data.scale
            model_obj.objective = @objective_endo11;
            model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'; 'kcatoff_8'; 'kdoff_8'; 'scale'}; % 
            % model_obj.theta_upper_bound = [0.1; 11.9; 0.1; 0.49; 10; 1; 500; 1000]; % 
            model_obj.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 12 % fit data.scale
            model_obj.objective = @objective_endo12;
            model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'}; % 
            % model_obj.theta_upper_bound = [0.4; 57; 0.011; 1.95; 17.4]; % 
%             model_obj.theta_name = {'kon_2'; 'koff_2'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; 'kon_4'; 'kon_7'; 'koff_7'}; % 
%             model_obj.theta_upper_bound = [0.4; 57; 1.25; 234; 38; 0.011; 1.95; 17.4]; % 
%            model_obj.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
            model_obj.index = 3; % [2;3;5];
    end
    model_obj.theta_bound = fh.get_theta_bound(model_obj, model_obj.theta_name); 
    optimize_ode_model = model_obj; 
end % function model_obj = fyn_gf_model(model_name)

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective(theta)
% Problem: the objective function cannot take any additional parameters

    model_fh = @fyn_gf_model_nodeg;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
        objective_base(theta, model_fh); 
end

function model_obj = fyn_gf_model_nodeg_endo(model_name, varargin)
    para_name = {'verbose', 'best_fit'};
    default_value = {1, 1};
    [verbose, best_fit] = parse_parameter(para_name, default_value, varargin);
    model_obj = fyn_gf_model_nodeg('complex_ode', 'best_fit', best_fit, 'multiple_output', 0, ...
        'verbose', verbose, 'fyn_endo', 1);
    model_obj.data.model = model_name; 
end


function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_endo(theta)

    model_fh = @fyn_gf_model_nodeg_endo;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
        objective_base(theta, model_fh); 
end

function model_obj = fyn_gf_model_nodeg_endo2(model_name, varargin)
    para_name = {'verbose', 'best_fit'};
    default_value = {1, 1};
    [verbose, ~] = parse_parameter(para_name, default_value, varargin);
    model_obj = fyn_gf_model_nodeg('complex_ode', 'best_fit', 2, 'multiple_output', 0, ...
        'verbose', verbose, 'fyn_endo', 2);
    model_obj.data.model = model_name; 
end

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_endo2(theta)

    model_fh = @fyn_gf_model_nodeg_endo2;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
        objective_base(theta, model_fh); 
end

function model_obj = fyn_gf_model_nodeg_endo3(model_name, varargin)
    para_name = {'verbose', 'best_fit'};
    default_value = {1, 1};
    [verbose, ~] = parse_parameter(para_name, default_value, varargin);
    model_obj = fyn_gf_model_nodeg('complex_ode', 'best_fit', 3, 'multiple_output', 0, ...
        'verbose', verbose, 'fyn_endo', 3);
    model_obj.data.model = model_name; 
end

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_endo3(theta)
    model_fh = @fyn_gf_model_nodeg_endo3;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
        objective_base(theta, model_fh); 
end

function model_obj = fyn_gf_model_nodeg_endo11(model_name, varargin)
    para_name = {'verbose', 'best_fit'};
    default_value = {1, 11};
    [verbose, ~] = parse_parameter(para_name, default_value, varargin);
    model_obj = fyn_gf_model_nodeg('complex_ode', 'best_fit', 11, 'multiple_output', 0, ...
        'verbose', verbose, 'fyn_endo', 11);
    model_obj.data.model = model_name; 
end

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_endo11(theta)
    model_fh = @fyn_gf_model_nodeg_endo11;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
        objective_base(theta, model_fh); 
end


function model_obj = fyn_gf_model_nodeg_endo12(model_name, varargin)
    para_name = {'verbose', 'best_fit'};
    default_value = {1, 11};
    [verbose, ~] = parse_parameter(para_name, default_value, varargin);
    model_obj = fyn_gf_model_nodeg('complex_ode', 'best_fit', 11, 'multiple_output', 0, ...
        'verbose', verbose, 'fyn_endo', 12);
    model_obj.data.model = model_name; 
end

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_endo12(theta)
    model_fh = @fyn_gf_model_nodeg_endo12;
    [f, t_exp_interp, y_exp_interp, y_ode_interp] = ...
        objective_base(theta, model_fh); 
end
