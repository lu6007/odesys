% The complex ode model 
% function ode_obj = complex_ode(model_name, varargin)

% Authors: Shaoying Kathy Lu (shaoying.lu@gmail.com); Tongkai Li
% (ltk@pku.edu.cn)
function ode_obj = complex_ode(model_name, varargin)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 

    % model_name = 'complex_ode';
    para_name = {'multiple_output', 'best_fit', 'verbose', 'fyn_endo'};
    default_value = {0, 1, 1, 0};
    [multiple_output, best_fit, verbose, fyn_endo] = parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction complex_ode(): model_name = %s\n', model_name);
        fprintf('multiple_output = %d, best_fit = %d, fyn_endo = %d\n', ...
            multiple_output, best_fit, fyn_endo); 
    end 
    
    data.model = model_name;
    data.multiple_output = multiple_output;
    data.best_fit = best_fit; 
    data.fyn_endo = fyn_endo; 
    
    %
    ode_obj.data = init_data(data); 
    ode_obj.init_data = @init_data;
    ode_obj.rhs = @rhs; 
    ode_obj.output = fh.output; 

end % function model_obj = fyn_gf_model(model_name)

% Complex ode with multiple output and best_fit results
function data = init_data(data) 
    multiple_output = data.multiple_output;
    best_fit = data.best_fit; 
    model = data.model; 
    fyn_endo = data.fyn_endo; 

    data = init_data_original(data);
    data.model = model; % needed for the function batch_fyn_gf() 
    if multiple_output
        data.output_index = [1:4, 6:8]; % 8; % Sensor %[1:4, 6:8];
    end
    if best_fit % Best for kon_4, koff_4, kcaton_7, kdon_7, kcatoff_7, kd_off_7
        data.kon_4 = 0.148; % /10;  % for conc. depend.   % s^(-1), to be fitted  
        data.koff_4 = 0.9205;       % s^(-1) [Iwamoto et al. 2016 PCB]
        data.kcaton_7 = 0.0143;     % s^(-1) Fyn activation by DGFRp
        data.kdon_7 = 40.8336;      % nM 
        data.kcatoff_7 = 0.0015;    % s^(-1) Fyn de-activation by PTP
        data.kdoff_7 = 170.0807;    % nM   
    end
    
    switch fyn_endo 
        case {1, 2, 3, 4, 5} % 1 or 2 11/10/2019
            data.kon_1 = 0.4509/10; % 0.4509/250;    % 0.0015/4;     % kon_1/300, kd = 10 nM [Schlesinger J 2002 Cell]
            data.koff_1 = 0.0154*330;     % kd = 100 nM
            data.kon_3 = 5;             % s^(-1) corrected
            data.kon_7 = 0.04;          % Huang et al 2011 Plos one
            data.koff_7 = 0.9356;       % Huang et al 2011 Plos one
            data.y0 = zeros(9,1);
            % data.scale = 593.2; % 194.2/0.3274
        case {11, 12, 13, 14, 15} % 11/17/2019
            data.kon_1 = 0.4509; % (nM s)^-1 ???
            data.koff_1 = 0.4509; % kd = 1 nM Freed 2017 Cell . kd = 175 nM [Dawson JP et al 2005]
            data.kon_2 = 9.917e-4; % (nM s)^-1 [Dawson et al 2005; Coban et al. 2015]
            data.koff_2 = 1.19; % s^-1 [Coban et al. 2015 Biophys J]
            data.kon_3 = 0.125; % s^-1 [Kim Y et al. 2012 Biochmistry]
            % data.kcaton_3 = 0.125; % s^-1 [Kim Y et al. 2012 Biochemistry]
            % data.kdon_3 = 1000/1000; % nM [Kim Y et al. 2012 Biochemistry]
            % from case {1-5}
            data.kon_7 = 0.04;
            data.koff_7 = 0.9356; 
            data.y0 = zeros(9,1); 
            % data.scale = 262.4007; % Inital guess, 85.91/0.3274  
    end
end
    

% In the future, need to separate init_data_original() and rhs() functions
% and merge complex_ode and complex_ode_nodeg
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
    data.kon_1 = 0.4509;         % (nM s)^(-1) [Iwamoto et al. 2016 PCB] -> EGF-EGFR
    data.koff_1 = 0.0154;        %   % s^(-1) [Iwamoto et al. 2016 PCB; Shin et al. 2018 PCB] 
    data.kon_2 = 0.056;          % (nM s)^(-1) [Iwamoto et al. 2016 PCB] -> DEGF-EGFR ? 
    data.koff_2 = 0.001;         % s^(-1) [Iwamoto et al. 2016 PCB]
    data.kon_3 = 7825;           %%% should be 5 (s)^(-1) [Iwamoto et al. 2016 PCB; Shin et al. 2018 PCB] -> DGFGFRp
    data.kcatoff_3 = 23.43;      % (nM s)^-1 [Shin et al. 2018 PCB; Iwamoto et al. 2016 PCB]
    data.kdoff_3 = 3.802;        % nM [Shin et al. 2018 PCB; Iwamoto et al. 2016 PCB]
    data.kon_4 = 0.01;           % s^(-1), to be fitted  
    data.koff_4 = 0.02;          % s^(-1) [Iwamoto et al. 2016 PCB]
    data.kon_5 = 0.03;           % s^(-1) [Iwamoto et al. 2016 PCB]
    data.kcaton_6 = 1.16e-3;     % (nM s)^(-1) [Shin et al. 2018 PCB]
    data.kdon_6 = 46.45;         % nM [Shin et al. 2018 PCB]
    data.vmaxoff_6 = 6.62e-5;    % s^(-1) [Shin et al. 2018 PCB]
    data.kmoff_6 = 52.72;        % nM [Shin et al. 2018 PCB]
    data.kcaton_7 = 0.02;        %%% should be changed to kon_7 and koff_7 % (nM s)^(-1) Fyn activation by DGFRp --> almost like a kon_7
    data.kdon_7 = 23.39;         % nM kd = kr/kf = 0.94/40 * 1000 nM = 23.39 nM by degf-egfr-src;
                                  % [Ref: Huang et al. 2011 PLOS ONE? ]
    data.kcatoff_7 = 0.02;       % s^(-1) Fyn de-activation by PTP
    data.kdoff_7 = 23.39;        
    data.kcaton_8 = 0.02;
    data.kdon_8 = 23.39;
    data.kcatoff_8 = 0.1; % 0.1;
    data.kdoff_8 = 23.39; 
    
    % Parameters below are needed if run with test_basal_level = 0
    data.fyn_act_0 = 0;
    data.gf_gfr_0 = 0;
    data.dgf_gfr_0 = 0;
    data.gfr_total = data.base_gfr_total_0; 
    
    % Initial values of y at t = 0
    data.y0 = zeros(8, 1); 
    
    % Output species
    data.species_name = {'GF-GFR', 'DGFR-GFR', 'DGFR-GFRp', 'ENDO','DEG', 'PTP Act', 'Fyn Act', 'Sensor Act'};
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
    fyn=fyn_total-fyn_act;
    if isfield(data, 'fyn_endo') && data.fyn_endo
        fyn_endo = y(9); 
        fyn = fyn - fyn_endo;
    end
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
    kcaton7 = data.kcaton_7;
    kdon7 = data.kdon_7;
    kcatoff7 = data.kcatoff_7;
    kdoff7 = data.kdoff_7;
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
    % PTP deactivate Fun with kcat_off_8 = 0.1; 
    v7=kcaton7*dgfgfrp*fyn/(kdon7+fyn)-kcatoff7*ptp_act*fyn_act/(kdoff7+fyn_act);
%     % PTP does not deactivate fyn with kcat_off_7 = 0.2; 
%     v7=kcaton7*dgfgfrp*fyn/(kdon7+fyn)-kcatoff7*fyn_act/(kdoff7+fyn_act);
    v8=kcaton8*(fyn_act)*sensor/(kdon8+sensor)-kcatoff8*ptp_act*sensor_act/(kdoff8+sensor_act);

    % dy/dt
    ydot = zeros(size(y));
    ydot(1)=v1-v2;      % gfgfr
    ydot(2)=0.5*v2-v3;  % dgfgfr
    ydot(3)=v3-v4;      % dgfgfrp
    ydot(4)=v4-v5;      % endo
    ydot(5)=v5;         % deg
    ydot(6)=v6;         % PTP
    ydot(7)=v7;         % fyn_act
    ydot(8)=v8;         % sensor_act
    
    if isfield(data, 'fyn_endo') 
        switch data.fyn_endo 
            case {1, 2, 3, 4, 5, 11, 12, 13, 14, 15}
                % make corrections on v7 and ydot(7)
                kon7 = data.kon_7; 
                koff7 = data.koff_7;
                % v7 = kon7*dgfgfrp*fyn-koff7*fyn_act-kcatoff7*ptp_act*fyn_act/(kdoff7+fyn_act);
                v7 = kon7*dgfgfrp*fyn-koff7*fyn_act;
                % Add the fyn_endo reaction
                v9 = kon4*fyn_act;
                ydot(7) = v7 - v9; 
                ydot(9) = v9; 
        end % switch data.fyn_endo                
    end % if isfield(data, 'fyn_endo') 
end



