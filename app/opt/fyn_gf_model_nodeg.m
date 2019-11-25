% A standard model function 
% function model_obj = fyn_gf_model_nodeg()
%     model_obj.data = init_data(); 
%     model_obj.init_data = @init_data; 
%     model_obj.rhs = @rhs; 
%     model_obj.output = @output; 
%     model_obj.opt.objective = @objective;
%     model_obj.opt.constraint = @constraint;
%     model_obj.opt.initial_guess = @initial_guess; 


% Authors: 
% Shaoying Kathy Lu (shaoying.lu@gmail.com)
% Tongkai Li (ltk@pku.edu.cn)
function model_obj = fyn_gf_model_nodeg(model_name, varargin)
    para_name = {'verbose', 'best_fit', 'fyn_endo'};
    default_value = {1, 1, 0};
    [verbose, best_fit, fyn_endo] = parse_parameter(para_name, default_value, varargin);
    % index = 3; for [gf] = 50 ng/ml
    % index = (1:5)'; for all important [gfs];

    model_base = fyn_gf_model('complex_ode', 'multiple_output', 0, 'best_fit', 1, ...
        'verbose', 0, 'fyn_endo', fyn_endo);
    data = model_base.data; 
    data.best_fit = best_fit; 
    data.model = model_name; 
    data.fyn_endo = fyn_endo; 
    
    if verbose
        fprintf('\nFunction fyn_gf_model_nodeg(): model_name = %s\n', model_name);
        fprintf('multiple_output fixed to %d, best_fit = %d, fyn_endo = %d\n', ...
            0, best_fit, fyn_endo); 
    end 

    model_obj.data = init_data(data); 
    model_obj.init_data = @init_data; 
    model_obj.rhs = model_base.rhs; % no fyn_endo
    model_obj.output = model_base.output; 
    %
    model_obj.opt = model_base.opt; 
    
    %
    % Optimization parameters
    % model_obj.opt.theta_name = {'kon_4'; 'kcaton_7'; 'kdon_7'; 'kcatoff_7'; 'kdoff_7'};
    % model_obj.opt.theta_upper_bound = [10; 1; 500; 1.0; 500];
%     model_obj.opt.theta_name = {'kon_1'; 'koff_1'}; %{'kon_2'; 'koff_2'; 'kon_4'; 'kcaton_7'; 'kdon_7'}; % 
%     model_obj.opt.theta_upper_bound = [1; 10]; %  [10; 1; 0.1; 0.2; 100]; % 
%     model_obj.opt.index = [3;5;6]; % 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
    switch fyn_endo
        case 0 
            model_obj.opt.objective = @objective;
            model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kcaton_7'; 'kdon_7'}; % 
            model_obj.opt.theta_upper_bound = [10; 1; 0.1; 0.2; 100]; %   
            model_obj.opt.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 1 % fit data.scale
            model_obj.opt.objective = @objective_endo;
%             model_obj.opt.theta_name = {'kon_1'; 'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'; 'kcatoff_7'; 'kdoff_7'; 'scale'}; % 
%             model_obj.opt.theta_upper_bound = [0.18; 5.6; 1; 0.1; 0.49; 10; 0.39; 500; 1000]; % 
            model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'; 'kcatoff_8'; 'kdoff_8'; 'scale'}; % 
            model_obj.opt.theta_upper_bound = [5.6; 1; 0.1; 0.49; 10; 1; 500; 1000]; % 
            model_obj.opt.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 2 % fit others
            model_obj.opt.objective = @objective_endo2;
%             model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4';'kon_7'; 'koff_7'}; % 
%             model_obj.opt.theta_upper_bound = [5.6; 1; 0.1; 0.49; 10]; % 
%             model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4';'vmaxoff_6'; 'kmoff_6'}; % 
%             model_obj.opt.theta_upper_bound = [5.6; 1; 0.1; 6.62e-1; 500]; % 
            model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'vmaxoff_6'; 'kmoff_6'}; % 
            model_obj.opt.theta_upper_bound = [5.6; 1; 6.62e-1; 500]; % 
%             model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4';'vmaxoff_6'; 'kmoff_6'; 'kon_7'; 'koff_7'}; % 
%             model_obj.opt.theta_upper_bound = [5.6; 1; 0.1; 6.62e-1; 500; 0.49; 10]; % 
            %
            model_obj.opt.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 3 % fit others
            model_obj.opt.objective = @objective_endo3;
%             model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'}; % 
%             model_obj.opt.theta_upper_bound = [5.6; 1; 0.1; 0.49; 10]; % 
%             model_obj.opt.theta_name = {'kon_1'; 'koff_1'; 'vmaxoff_6'; 'kmoff_6'; 'kcatoff_8'; 'kdoff_8'}; % 
%             model_obj.opt.theta_upper_bound = [1; 1; 1; 500; 1; 500]; % 
            model_obj.opt.theta_name = {'kon_1'; 'koff_1'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; 'kcatoff_8'; 'kdoff_8'}; % 
            model_obj.opt.theta_upper_bound = [1; 1; 50; 100; 50; 1; 500]; % 
            model_obj.opt.index = [2;3;5;6]; % [3;4;5]; % 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 11 % fit data.scale
            model_obj.opt.objective = @objective_endo11;
            model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'; 'kcatoff_8'; 'kdoff_8'; 'scale'}; % 
            model_obj.opt.theta_upper_bound = [0.1; 11.9; 0.1; 0.49; 10; 1; 500; 1000]; % 
            model_obj.opt.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
        case 12 % fit data.scale
            model_obj.opt.objective = @objective_endo12;
            model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'kon_7'; 'koff_7'}; % 
            model_obj.opt.theta_upper_bound = [0.4; 57; 0.011; 1.95; 17.4]; % 
%             model_obj.opt.theta_name = {'kon_2'; 'koff_2'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; 'kon_4'; 'kon_7'; 'koff_7'}; % 
%             model_obj.opt.theta_upper_bound = [0.4; 57; 1.25; 234; 38; 0.011; 1.95; 17.4]; % 
%            model_obj.opt.index = 3; % (1:6)'; % 6; %(1:6)'; % 3 - 50 ng/ml only
            model_obj.opt.index = 3; % [2;3;5];
    end

    

 
end % function model_obj = fyn_gf_model(model_name)
    
function data = init_data(data) % egfr2
    best_fit = data.best_fit; 
    
    % No degradation model
    data.koff_4 = 0;             % s^(-1) 
    data.kon_5 = 0;              % s^(-1) [Iwamoto et al. 2016 PCB]
    
    switch best_fit
        % default best_fit = 0
        % Use the best model from complex_ode 
        % fitted to {'kon_4'; 'kcaton_7'; 'kdon_7'; 'kcatoff_7'; 'kdoff_7'};
        % data.kon_4 = 0.148; % /10; % for conc. depend.   % s^(-1), to be fitted  
        % data.koff_4 = 0.9205;       % s^(-1) [Iwamoto et al. 2016 PCB]
        % data.kcaton_7 = 0.0143;     % s^(-1) Fyn activation by DGFRp
        % data.kdon_7 = 40.8336;      % nM 
        % data.kcatoff_7 = 0.0015;    % s^(-1) Fyn de-activation by PTP
        % data.kdoff_7 = 170.0807;    % nM   
        case 1 % fit to {'kon_4'; 'kcaton_7'; 'kdon_7'; 'kcatoff_7'; 'kdoff_7'}; 
            % error = 0.002756
            data.kon_4 = 0.004950;          % s^(-1), 
            data.kcaton_7 = 0.01296;        % s^(-1) Fyn activation by DGFRp
            data.kdon_7 = 31.46;            % nM 
            data.kcatoff_7 = 0.003907;      % s^(-1) Fyn de-activation by PTP
            data.kdoff_7 = 	497.8;          %  
            % fit again to {'kon_2'; 'koff_2'; 'kon_4'; 'kcaton_7'; 'kdon_7'};
            % and get the same set of parameters, error = 0.002756
            % data.kon_2 = 0.056;
            % data.koff_2 = 0.001; 
            
        case {2, 3} % after fitting the scale factor and other parameters for case 1
            % fyn_endo = 2 as well
            data.kon_2 = 3.87880195; 										
            data.koff_2 = 0.34335529;
            % data.kon_4 = 4.11E-06;
            data.kon_4 = 0.0057; % best fit
            data.kon_7 = 0.180851762;	
            data.koff_7 = 7.895453269;
            data.kcatoff_8 = 0.07948579;
            data.kdoff_8 = 206.2281065;
            data.scale = 597.507106591057; % Smallest error fitting scale
        
        case 11 % best fit from fyn_endo = 11
            data.kon_2 = 0.040558595;
            data.koff_2 = 5.720696138;
            data.kon_4 = 0.001076275;
            data.kon_7 = 0.195210427;
            data.koff_7 = 1.738250076; 
            data.kcatoff_8 = 0.697117467;
            data.kdoff_8 = 340.8719999;
            % data.scale = 548.4923813; % best fit scale
            data.scale = 133; % PDGF max 1.5 --> 100% of sensor = 200 nM 
    end
    
end

function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective(theta)

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
