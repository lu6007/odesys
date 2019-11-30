% fyn signal under growth factor stimulation 
% Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Copyright (2015) The Regents of the University of California
% All Rights Reserved

function data = fyn_gf_init_data( model )
fh = opt_utility();
root = fh.get_root();
%
data.model = model;
data.k_don_4 = 23.39; %nm kd = kr/kf = 0.94/40 *1000 nM = 23.39 nM by degf-egfr-src
data.k_doff_4 = data.k_don_4;
switch(model)

    % Assume that the biosensor readout is Fyn activity
    % Assume that all the growth factors are active in dimers and inactive in
    % monomer; and that only active growth factor dimers can be internalized. 
    % Case number corresponds to the concentration of gfr for each
    % experiment
      
   case 'egfr_huang' % Corrected according to Huang paper.  
        % Huang et al. 2011 PLOS ONE. 
        % Simulating EGFR-ERK paper parameters
        % Concentration Parameters
        data.fyn_total =518;       % [nM]
        data.ptp = 1000;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 0.8065; % 80.6500; % 16.1300; % 1.6130; % 4.0325; % 8.065;            % [nM] gf added at time 0 -> 50 ng/ml
        
        % Reaction Parameters
        data.k_caton_4 = 0.5838;            % [1/(Sec)] activation rate of Fyn by DGFR
        data.k_catoff_4 = 0.2661;            % [1/(Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_on_1 = 0.1;           % [1/(nM Sec)] association rate between gf and gfr
        data.k_off_1 = 0.0038; % 0.0038;        % [1/sec] dis-association rate of gf_gfr
        data.k_on_2 = 0.01;          % [1/(nM sec)] association rate of the growth factor dimers
        data.k_off_2 = 0.02;          % [1/sec] dis-association rate of the growth factor dimers
        data.k_on_3 = 0.000007;         % [1/sec] internalization rate of the growth factor dimers
        
        % Concentration for testing basal level
        data.base_gfr_total_0 = 300; % [nM]

   case 'egfr_huang_kathy' % Corrected according to Huang paper.  
        % Huang et al. 2011 PLOS ONE. 
        % Simulating EGFR-ERK paper parameters
        % Concentration Parameters
        data.fyn_total = 50;        % [nM]
        data.ptp = 100;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 0.8065; % 80.6500; % 16.1300; % 1.6130; % 4.0325; % 8.065;            % [nM] gf added at time 0 -> 5 ng/ml
        
        % Reaction Parameters
        data.k_caton_4 = 0.005;            % [1/(Sec)] activation rate of Fyn by DGFR
        data.k_catoff_4 = 0.01;            % [1/(Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_on_1 = 0.01;           % [1/(nM Sec)] association rate between gf and gfr
        data.k_off_1 = 0.0038; % 0.0038;        % [1/sec] dis-association rate of gf_gfr
        data.k_on_2 = 0.01;          % [1/(nM sec)] association rate of the growth factor dimers
        data.k_off_2 = 0.02;          % [1/sec] dis-association rate of the growth factor dimers
        data.k_on_3 = 0.002;         % [1/sec] internalization rate of the growth factor dimers
        % data.k_on_1 and data.k_on_3 are important
        % pdgf signal can be simulated by making
        % data.k_on_3 100x smaller
        
        % Concentration for testing basal level
        data.base_gfr_total_0 = 300; % [nM]


    case 'egfr_huang_v2' % Chris' version 3/30/2016 
        % Huang et al. 2011 PLOS ONE. 
        % Simulating EGFR-ERK paper parameters
        % Concentration Parameters
        data.fyn_total = 1;       % [nM]
        data.ptp = 0.2;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        % [nM] gf added at time 0 -> 50 ng/ml
        data.gf_1 = 0.78125; % [78.125, 15.625, 7.8125, 3.90625, 1.5625, 0.78125, 0.15625]; 
        
        % Reaction Parameters
        data.k_caton_4 = 1;            % [1/(nM Sec)] activation rate of Fyn by DGFR
        data.k_catoff_4 = 2;            % [1/(nM Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_on_1 = 0.01;           % [1/(nM Sec)] association rate between gf and gfr
        data.k_off_1 = 0.003; % 0.0038;        % [1/sec] dis-association rate of gf_gfr
        data.k_on_2 = 0.02;          % [1/(nM sec)] association rate of the growth factor dimers
        data.k_off_2 = 0.052;          % [1/sec] dis-association rate of the growth factor dimers
        data.k_on_3 = 0.007;         % [1/sec] internalization rate of the growth factor dimers
        
        % Concentration for testing basal level
        data.base_gfr_total_0 = 0.4; % [nM]
        
        % Parameters below are needed if run with test_basal_level = 0
        data.fyn_act_0 = 0;
        data.gf_gfr_0 = 0;
        data.dgf_gfr_0 = 0;
        data.gfr_total_0 = data.base_gfr_total_0; 

  case 'egfr_huang_v3' % Kathy's version 11/2/2019 
        % Based on 'egfr_huang_v2'
        % Concentration Parameters
        data.fyn_total = 174.9847;  % nM, based on cbl kinase [Shin et al. 2018 PCB]
        data.ptp = 296.58/10;          % nM, [Shin et al. 2018 PCB, 296.58 nM]
        data.gf_0 = 0; % 0.0166;         % [nM] gf at basal level 
        % [nM] gf added at time 0 -> 50 ng/ml
        data.gf_1 = 7.8125; % 50 ng/ml, [78.125, 15.625, 7.8125, 3.90625, 1.5625, 0.78125, 0.15625]; 
        
        % Reaction Parameters
        % fit k_caton_4, k_catoff_4, k_on_3, scale
        data.k_caton_4 = 0.01;            % [1/(Sec)] activation rate of Fyn by DGFR
        data.k_catoff_4 = 0.02;            % [1/(Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_on_1 = 0.45;           % (nM s)^(-1) [Iwamoto et al. 2016 PCB]
        data.k_off_1 = 0.015;         % s^(-1) [Shin et al. 2018 PCB; Iwamoto et al. 2016 PCB]
        data.k_on_2 = 3.6e-5;         % [1/(nM sec)] [Iwamoto et al. 2016 PCB]
        data.k_off_2 = 0.001;         % [1/sec] [Iwamoto et al. 2016 PCB]
        data.k_on_3 = 2e-3;           % [1/sec] internalization rate of the growth factor dimers
        
        % Concentration for testing basal level
        data.base_gfr_total_0 = 398.10; % nM, Shin et al. 2018 PCB
        
        % Parameters below are needed if run with test_basal_level = 0
        data.fyn_act_0 = 0;
        data.gf_gfr_0 = 0;
        data.dgf_gfr_0 = 0;
        data.gfr_total_0 = data.base_gfr_total_0; 
        
   case 'egfr_blinov'
        % Blinov et al. 2006 Biosystems. A Network Model paper parameters
        % Concentration Parameters
        data.fyn_total = 0.5;       % [nM]
        data.ptp = 0.2;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 8.065;            % [nM] gf added at time 0 - 50 ng/ml
        
        % Reaction Parameters
        data.k_caton_4 = 1;            % [1/(nM Sec)] activation rate of Fyn by DGFR
        data.k_catoff_4 = 2;            % [1/(nM Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_on_1 = 0.003;         % [1/(nM Sec)] association rate between gf and gfr
        data.k_off_1 = 0.06;          % [1/sec] dis-association rate of gf_gfr
        data.k_on_2 = 0.01;          % [1/(nM Sec)] association rate of the growth factor dimers
        data.k_off_2 = 0.1;           % [1/sec] dis-association rate of the growth factor dimers
        data.k_on_3 = 6e-3;         % [1/sec] internalization rate of the growth factor dimers 
        
        % Concentration for testing basal level
        data.base_gfr_total_0 = 0.4; % [nM]


     case 'pdgfr_huang'
        % PDGF/PDGFR parameters adapted from the model 'egf_huang'
        % Concentration Parameters
        data.fyn_total = 0.5;       % [nM]
        data.ptp = 0.2;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 2.03;             % [nM]  gf added at time 0 - 25 ng/ml
        % Reaction Parameters
        data.k_caton_4 = 1;            % [1/(nM Sec)] activation rate of Fyn by DGFR
        data.k_catoff_4 = 2;            % [1/(nM Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_on_1 = 0.1;           % [1/(nM Sec)] association rate between gf and gfr
        data.k_off_1 = 0.0038;        % [1/sec] dis-association rate of gf_gfr
        data.k_on_2 = 0.02;          % [1/(nM sec)] association rate of the growth factor dimers
        data.k_off_2 = 0.02;          % [1/sec] dis-association rate of the growth factor dimers
        data.k_on_3 = 0.0006;       % [1/sec] internalization rate of the growth factor dimers
%         % Initial Concentration
%         data.fyn_act_0 = 0.00223;  % [nM]
%         data.gf_gfr_0 = 0.04293;     % [nM]
%         data.dgf_gfr_0 = 0.001792;   % [nM]
%         data.gfr_total_0 = 0.1457;    % [nM]

        % Concentration for testing basal level
        data.base_gfr_total_0 = 0.15; % [nM]
        
        % figure; plot(t/60, fyn_act*15.7+0.145); title('[FYN ACT]'); xlabel('Time (min)');
        % hold on; plot(fyn_act_exp(:,1)-5, fyn_act_exp(:,2),'ro');
       
    case 'exp_hela_egf'
        data.file = [root, '0124_2018/b.mat'];     
    case 'exp_mef_pdgf'
        data.file = [root, '1019_2019/pdgf/input/exp.mat'];
end
return;

