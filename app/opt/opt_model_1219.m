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
        model.theta_fit = [0.177946631	3.479475449	0.131169418	125.99938	...
            21.48235746	38.01936756	0.001999939	0.000200841	0.2609356	3.626668782]; 
    end
    
    if model_id >= 2
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_fit = [1.744816973	0.118021428	0.161653954	457.365617	...
            95.20155751	380.0109548	0.01316032 6.21E-06 1.771366636	34.81210749]; 
        % same theta_name as model_id = 1
    end

    if model_id >= 3
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; ...
            'kmoff_6'; 'kon_4'; 'koff_4'; 'kon_5'};
        model.theta_fit = [0.008614418	10.75312475	0.000146671	12.88583716	...
            6.62E-04	0.527207213	0.007604692	6.21E-05	0.000300694];
    end
    
    if model_id >= 4
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        model.theta_fit = [0.077876944	76.8010489	0.010285771	0.000620918	...
            7.780614418	152.2408491];
    end
    
    if model_id >=5 
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_1'; 'koff_1'; 'kcaton_3'; 'kdon_3'; 'kcatoff_3'; 'kdoff_3';
            'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        model.theta_fit = [12.45740524	0.908054327	0.156973969	1098.618859 ...
            951.2502672	1217.450961	0.017120522	0.002715464	77.74608507	492.8449488]; 
    end
    
    if model_id >=6 
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; ...
            'kmoff_6'; 'kon_4'; 'koff_4'; 'kon_5'};
        model.theta_fit = [0.718142337	354.0134501	0.000310683	93.87921434	...
            0.000108461	0.008211286	0.021486414	0.001182615	5.12E-06]; 
    end

    if model_id >=7 
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        % model.theta_fit = []; 
    end

    if model_id >= 11 % concentration dependence
        model.theta_name = {'kon_1'; 'koff_1'; 'kcaton_3'; 'kdon_3'; 'kcatoff_3'; 'kdoff_3';
            'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        % model.theta_fit = []; 
        model.index = [2; 3; 5];
    end
    
    model.theta_bound = fh.get_theta_bound(model, model.theta_name, ...
        'bound_factor', theta_bound_factor);
end % function model = opt_model_1219(model_name, varargin)






