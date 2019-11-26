% The complex ode model with no degradation
% function model_obj = complex_nodeg()

% Authors: 
% Shaoying Kathy Lu (shaoying.lu@gmail.com)
% Tongkai Li (ltk@pku.edu.cn)
function ode = complex_ode_nodeg(model_name, varargin)

    para_name = {'verbose', 'best_fit', 'fyn_endo'};
    default_value = {1, 1, 0};
    [verbose, best_fit, fyn_endo] = parse_parameter(para_name, default_value, varargin);
    if verbose
        fprintf('\nFunction complex_ode_nodeg(): model_name = %s\n', model_name);
        fprintf('multiple_output fixed to %d, best_fit = %d, fyn_endo = %d\n', ...
            0, best_fit, fyn_endo); 
    end 

    ode_base = complex_ode(model_name, ...
        'best_fit', best_fit, 'fyn_endo', fyn_endo, 'verbose', verbose);
    data = ode_base.data; 
    ode.data = init_data(data); 
    ode.init_data = @init_data; 
    ode.rhs = ode_base.rhs; % no fyn_endo
    ode.output = ode_base.output; 
end % ode = complex_ode_nodeg(model_name, varargin)
    
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
            % data.scale = 597.507106591057; % Smallest error fitting scale
        
        case 11 % best fit from fyn_endo = 11
            data.kon_2 = 0.040558595;
            data.koff_2 = 5.720696138;
            data.kon_4 = 0.001076275;
            data.kon_7 = 0.195210427;
            data.koff_7 = 1.738250076; 
            data.kcatoff_8 = 0.697117467;
            data.kdoff_8 = 340.8719999;
            % data.scale = 548.4923813; % best fit scale
            % data.scale = 133; % PDGF max 1.5 --> 100% of sensor = 200 nM 
    end
    
end

