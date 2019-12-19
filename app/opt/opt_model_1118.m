% A standard model function 
% function model = opt_model_1118(model_name, varargin)
%     para_name = {'verbose', 'theta_bound_factor', 'model_id'};
%     default_value = {1, 10, 1};
% model.data = init_data(); 
% model.init_data = @init_data; 
% model.rhs = @rhs; 
% model.output = @output; 

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com)
function model = opt_model_1118(model_name, varargin)
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
    
    ode_model = ode_model_1118('model_1', 'verbose', 0);
    model.ode = ode_model; % data, init_data, rhs, output
    
    % [Ouyang et al 2019 ACS Sens; 350% activation corresponds to 100%]
    model.scale = 57.1429;       % nM/au ratio = 200 nM /3.5 au; 0.328 * scale = 18.7429
    model.index = 3;             % [egf] = 50 ng/ml %[3; 5; 6] % 50, 10, 5
    
    % Optimization functions and parameters
    model.objective = fh.objective;
    % model.constraint = fh.constraint;
    model.initial_guess = fh.initial_guess; 
    
    if model_id >= 1
        model.theta_name = {'kon_1'; 'koff_1'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'; ...
            'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        model.theta_fit = [0.26031169 15.35232148 0.131752578 12.56740995 35.81094633 ...
            0.009446086 0.085982789 0.185889161 4.766242707]'; 
    end
    
    if model_id >= 2
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'}; 
        model.theta_fit = [0.002429011	2.894388505	0.002523534	0.008028778	...
            0.004306559	0.107611631]';
    end
    
    if model_id >= 3
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_4'; 'koff_4'; 'kon_7'; 'koff_7';....
            'kcaton_8';'kdon_8'; 'kcatoff_8'; 'kdoff_8'}; 
        model.theta_fit = [0.003861532	0.019056514	0.000675076	0.016057544	...
            0.035225733	48.26366535	0.984973073	208.7484753]';
    end 
    
    if model_id >= 4
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; 'kmoff_6'; ...
            'kon_4'; 'koff_4'; 'kon_5'}; 
        model.theta_fit = [0.004650675	6.293209761	0.008951284	191.4474109	0.000616027	...
            16.57132624	0.002065168	0.102335935	0.09548778]';
    end

    if model_id >= 5
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; 'kmoff_6'; ...
            'kon_4'; 'koff_4'}; 
        model.theta_fit = [0.011413924	9.60380283	0.013044219	211.1324765	0.004019118 ...
            4.195132281	0.000941996	0.290217319]';
    end
    
    if model_id >= 6 % Fit single concentration
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        model.theta_fit = [0.013780574	9.591754187	1.39E-05 0.30794755	0.00053241	0.012744546]';
    end

    if model_id >= 7 && model_id <=9 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_1'; 'koff_1'; 'kon_2'; 'koff_2'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'};
        model.index = [2; 3; 5; 6];
        model.theta_fit = [0.419148398	2.719934631	0.000877458	88.72167814	0.769892656	62.89582402 ...
            139.2952237]';
    end

    if model_id >= 8 && model_id <=9  % Fit concentration dependence, same parameters as model_id = 5
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; 'kmoff_6'; ...
            'kon_4'; 'koff_4'}; 
        % model.index = [2; 3; 5; 6]; 
        model.theta_fit = [0.007515766	529.671426	0.063154203	1799.131895	0.027426714	...
            30.39885476	6.93E-05 1.50340628]';
    end

    if model_id >= 9 && model_id <=9 % Fit concentration dependence, same parameters as model_id = 6
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'}; 
        % model.index = [2; 3; 5; 6]; 
        model.theta_fit = [0.009876822 879.6740399 0.000311732 7.468174309 0.000737206 0.017578657]'; 
    end

    if model_id >= 10 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_1'; 'koff_1'; 'kon_2'; 'koff_2'; 'kon_3'; 'kcatoff_3'; 'kdoff_3'};
        model.index = [2; 3; 5];
        if model_id <= 12
            model.theta_fit = [2.57637294   2.855642128 0.022710822 35.12714557 0.003079185 ...
               46.72058452  279.8101348]'; % error = 0.0464
            
%             % global optimization, error = 0.039
%             model.theta_fit = [0.6602    0.1650    0.0127    7.7818    0.0488   13.2180 ...
%               0.8429]';
        elseif model_id >= 13
            model.theta_fit = [0.535876586	0.489827576	0.003363463	72.03986762	0.039343271	...
                19.66399713	89.77301467]'; 
        end
    end
    
    if model_id >= 11 && model_id <= 12 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; 'kmoff_6'; ...
            'kon_4'; 'koff_4'};
        % model.index = [2; 3; 5];
        model.theta_fit = [0.144434405	121.1713412	0.017469763	225.3910039	0.040067134	...
                36.70762217	4.96E-05 1.614352233]';
    end

    if model_id >= 12 && model_id <=12 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        % model.index = [2; 3; 5];
        model.theta_fit = [0.155118662	121.3401238	0.000238357	8.082622621	...
            0.000511455	0.01223265]';
    end
    
    if model_id >= 13 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; 'kmoff_6'; ...
            'kon_4'; 'koff_4'; 'kon_5'};
        % model.index = [2; 3; 5];
        model.theta_fit = [0.018155344	248.2778007	0.069520447	1800.074546	0.039514857 ...
            3.70E+01 5.13E-05 1.368242319 0.311534566]';
    end
    
    if model_id >= 14 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kcaton_6'; 'kdon_6'; 'vmaxoff_6'; 'kmoff_6'; ...
            'kon_4'; 'koff_4'};
        % model.index = [2; 3; 5];
        model.theta_fit = [0.103543691	1416.217306	0.409069887	12084.979 ...
            0.364703436	360.5672346	0.000249098	8.922180513]';
    end

    if model_id >= 15 % Fit concentration dependence
        model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
        model.theta_name = {'kon_2'; 'koff_2'; 'kon_4'; 'koff_4'; 'kon_7'; 'koff_7'};
        % model.index = [2; 3; 5];
        model.theta_fit = []';
    end

    %
    model.theta_bound = fh.get_theta_bound(model, model.theta_name, ...
        'bound_factor', theta_bound_factor);
end % function model = opt_model_1118(model_name, varargin)






