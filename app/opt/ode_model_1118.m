% A standard model function 
% function model_obj = ode_model_1118(model_name, varargin)
%     para_name = {'multiple_output', 'best_fit'};
%     default_value = {0, 0};
% model_obj.data = init_data(); 
% model_obj.init_data = @init_data; 
% model_obj.rhs = @rhs; 
% model_obj.output = @output; 

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com); Tongkai Li
% (ltk@pku.edu.cn)
function ode = ode_model_1118(model_name, varargin)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 

    para_name = {'multiple_output', 'best_fit', 'verbose', 'fyn_endo'};
    default_value = {0, 1, 1, 0};
    [multiple_output, best_fit, verbose, fyn_endo] = parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction ode_model_1118(): model_name = %s\n', model_name);
        fprintf('multiple_output = %d, best_fit = %d, fyn_endo = %d\n', ...
            multiple_output, best_fit, fyn_endo); 
    end 
    
    data.model = model_name;
    data.multiple_output = multiple_output;
    data.best_fit = best_fit; 
    data.fyn_endo = fyn_endo; 
        
    %
  	ode.data = init_data(data); 
    ode.init_data = @init_data;
    ode.rhs = @rhs; 
    ode.output = fh.output; 

end % function model_obj = fyn_gf_model(model_name)

% Complex ode with multiple output and best_fit results
function data = init_data(data) 
    multiple_output = data.multiple_output;
    model = data.model; 
%     best_fit = data.best_fit; 
%     fyn_endo = data.fyn_endo; 

    data = init_data_original(data);
    data.model = model; % needed for the function batch_fyn_gf() 
    if multiple_output
        data.output_index = [1:4, 6:8]; % 8; % Sensor %[1:4, 6:8];
    end
end
    
function data = init_data_original(data) % complex ode with single output
    % Potentially can separate  this an ode model and an optimization model 
    
    % Based on 'egfr_huang_v3' and 'rhs_2'
    % Concentration Parameters
    data.fyn_total = 174.9847;  % nM, based on cbl kinase [Shin et al. 2018 PCB]
    data.ptp_total = 296.58; % nM, [Shin et al. 2018 PCB, 296.58 nM]
    data.sensor_total = 200; % nM, estimated
    data.base_gfr_total_0 = 398.10; % nM, Shin et al. 2018 PCB
    data.gf_0 = 0;          % [nM] gf at basal level 
    % [nM] gf added at time 0 -> 50 ng/ml
    data.gf_1 = 7.8125; % 50 ng/ml, [78.125, 15.625, 7.8125, 3.90625, 1.5625, 0.78125, 0.15625]; 

    % Reaction Parameters
    data.kon_1 = 0.4509;        % (nM s)^(-1) [Iwamoto et al. 2016 PCB] -> EGF-EGFR
    data.koff_1 = 0.4509*10;       %   % s^(-1) kd = 1 nM [Freed 2017 Cell], kd = 175 nM [Dawson JP et al 2005] 
    data.kon_2 = 9.917e-4;      % (nM s)^-1 [Dawson et al 2005; Coban et al. 2015] 
    data.koff_2 = 1.19;         % s^-1 [Coban et al. 2015 Biophys J]
    data.kon_3 = 0.125;         % s^-1 [Kim Y et al. 2012 Biochmistry]
    data.kcatoff_3 = 23.43;     % (nM s)^-1 [Shin et al. 2018 PCB; Iwamoto et al. 2016 PCB]
    data.kdoff_3 = 3.802;       % nM [Shin et al. 2018 PCB; Iwamoto et al. 2016 PCB]
    %%%
    data.kon_4 = 0.02/100;       % s^(-1), to be fitted  
    data.koff_4 = 0.02;         % s^(-1) [Iwamoto et al. 2016 PCB]
    data.kon_5 = 0.03;       % s^(-1) [Iwamoto et al. 2016 PCB]
    %%%  6: Negative regulation and concentration dependence
    data.kcaton_6 = 1.16e-3*5;     % (nM s)^(-1) [Shin et al. 2018 PCB]
    data.kdon_6 = 46.45;          % nM [Shin et al. 2018 PCB]
%     data.kcaton_6 = 1.16e-3;    % (nM s)^(-1) [Shin et al. 2018 PCB]
%     data.kdon_6 = 46.45;        % nM [Shin et al. 2018 PCB]
    %%% 
    data.vmaxoff_6 = 6.62e-5;    % s^(-1) [Shin et al. 2018 PCB]
    data.kmoff_6 = 52.72;        % nM [Shin et al. 2018 PCB]
    %
    data.kon_7 = 0.04;       % Ref? 
    data.koff_7 = 0.9356; 
    data.kcaton_8 = 0.02; 
    data.kdon_8 = 23.39;
    data.kcatoff_8 = 0.1; 
    data.kdoff_8 = 23.39; 
    
    % Parameters below are needed if run with test_basal_level = 0
    data.fyn_act_0 = 0;
    data.gf_gfr_0 = 0;
    data.dgf_gfr_0 = 0;
    data.gfr_total = data.base_gfr_total_0; 
    
    % Initial values of y at t = 0
    data.y0 = zeros(10, 1); 
    
    % Output species
    data.species_name = {'GF-GFR', 'DGFR-GFR', 'DGFR-GFRp', 'ENDO','DEG',...
        'PTP Act', 'Fyn Act', 'Sensor Act','Fyn Endo'};
    data.output_index = 8; % [1:4, 6:8]; % 8; % Sensor %[1:4, 6:8];
    
end

function ydot = rhs(t, y, data)

    gfgfr = y(1);
    dgfgfr = y(2);
    dgfgfrp = y(3);
    endo = y(4);
    deg = y(5);
    ptp_act = y(6);
    fyn_act = y(7);
    sensor_act = y(8);
    fyn_endo = y(9); 
    fyn_deg = y(10); 

    % Concentration Parameters
    if t <0
        gf = data.gf_0;
    elseif t>=0
        % Add growth factor at a certain concentration
        gf = data.gf_0+data.gf_1;
    end
    gfr_total = data.gfr_total;  
    fyn_total = data.fyn_total;
    ptp_total = data.ptp_total; 
    sensor_total = data.sensor_total;

    gfr=gfr_total-gfgfr-2*(dgfgfr+dgfgfrp+endo+deg);
    fyn=fyn_total-fyn_act-fyn_endo-fyn_deg; 
    ptp=ptp_total-ptp_act;
    sensor=sensor_total-sensor_act;

    kon1 = data.kon_1;
    koff1 = data.koff_1;
    kon2 = data.kon_2;
    koff2 = data.koff_2;
    kon3 = data.kon_3;
    kcatoff3 = data.kcatoff_3;
    kdoff3 = data.kdoff_3;
    kon4 = data.kon_4;
    koff4 = data.koff_4;
    kon5 = data.kon_5;
    kcaton6 = data.kcaton_6;
    kdon6 = data.kdon_6;
    vmaxoff6 = data.vmaxoff_6;
    kmoff6 = data.kmoff_6;
    kon7 = data.kon_7;
    koff7 = data.koff_7;
    kcaton8 = data.kcaton_8;
    kdon8 = data.kdon_8;
    kcatoff8 = data.kcatoff_8;
    kdoff8 = data.kdoff_8; 

    v1=kon1*gf*gfr-koff1*gfgfr;
    v2=kon2*gfgfr^2-2*koff2*dgfgfr;
    v3=kon3*dgfgfr-(kcatoff3*ptp_act)*dgfgfrp/(kdoff3+dgfgfrp);
    v4=kon4*dgfgfrp-koff4*endo;
    v5=kon5*endo;
    v6=kcaton6*dgfgfrp*ptp/(kdon6+ptp)-vmaxoff6*ptp_act/(kmoff6+ptp_act);
    v7=kon7*dgfgfrp*fyn-koff7*fyn_act;
    % v8=kcaton8*(fyn_act)*sensor/(kdon8+sensor)-vmaxoff6*fyn_act/(kmoff6+fyn_act);
    v8=kcaton8*(fyn_act)*sensor/(kdon8+sensor)-kcatoff8*ptp_act*sensor_act/(kdoff8+sensor_act);
    v9 = kon4 * fyn_act - koff4 * fyn_endo; 
    v10 = kon5 * fyn_endo; % fyn_endo deg

    % dy/dt
    ydot = zeros(size(y));
    ydot(1)=v1-v2;      % gfgfr
    ydot(2)=0.5*v2-v3;  % dgfgfr
    ydot(3)=v3-v4;      % dgfgfrp
    ydot(4)=v4-v5;      % endo
    ydot(5)=v5;         % deg
    ydot(6)=v6;         % PTP
    ydot(7)=v7-v9;      % fyn_act
    ydot(8)=v8;         % sensor_act
    ydot(9) = v9-v10;   % fyn_endo
    ydot(10) = v10;     % fyn_deg
    
end


