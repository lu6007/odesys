% fyn signal under growth factor stimulation 
% Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Copyright (2015) The Regents of the University of California
% All Rights Reserved

function data = init_data( model )
switch(model),

    % Assume that the biosensor readout is Fyn activity
    % Assume that all the growth factors are active in dimers and inactive in
    % monomer; and that only active growth factor dimers can be internalized. 
    case 'egfr_huang',
        % Huang et al. 2011 PLOS ONE. 
        % Simulating EGFR-ERK paper parameters
        % Concentration Parameters
        data.fyn_total = 1;       % [nM]
        data.ptp = 0.2;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 8.065;            % [nM] gf added at time 0 -> 50 ng/ml
        
        % Reaction Parameters
        data.k_catf = 1;            % [1/(nM Sec)] activation rate of Fyn by DGFR
        data.k_catb = 2;            % [1/(nM Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_gff = 0.01;           % [1/(nM Sec)] association rate between gf and gfr
        data.k_gfb = 0.0038;        % [1/sec] dis-association rate of gf_gfr
        data.k_dgf = 0.01;          % [1/(nM sec)] association rate of the growth factor dimers
        data.k_dgb = 0.02;          % [1/sec] dis-association rate of the growth factor dimers
        data.k_enf = 6e-3;         % [1/sec] internalization rate of the growth factor dimers
        
        % Initial Concentration
        data.fyn_act_0 = 0.0030;    % [nM]
        data.gf_gfr_0 = 0.079;     % [nM]
        data.dgf_gfr_0 = 0.0024;    % [nM]
        data.gfr_total_0 = 0.2769;  % [nM]
        
        % Concentration for testing basal level
        data.base_gfr_total_0 = 0.4; % [nM]
        
     case 'pdgfr_huang',
        % PDGF/PDGFR parameters adapted from the model 'egf_huang'
        % Concentration Parameters
        data.fyn_total = 0.5;       % [nM]
        data.ptp = 0.2;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 2.03;             % [nM]  gf added at time 0 - 25 ng/ml
        % Reaction Parameters
        data.k_catf = 1;            % [1/(nM Sec)] activation rate of Fyn by DGFR
        data.k_catb = 2;            % [1/(nM Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_gff = 0.1;           % [1/(nM Sec)] association rate between gf and gfr
        data.k_gfb = 0.0038;        % [1/sec] dis-association rate of gf_gfr
        data.k_dgf = 0.02;          % [1/(nM sec)] association rate of the growth factor dimers
        data.k_dgb = 0.01;          % [1/sec] dis-association rate of the growth factor dimers
        data.k_enf = 6e-4;       % [1/sec] internalization rate of the growth factor dimers
        % Initial Concentration
        data.fyn_act_0 = 0.0037;  % [nM]
        data.gf_gfr_0 = 0.0399;     % [nM]
        data.dgf_gfr_0 = 0.0030;   % [nM]
        data.gfr_total_0 = 0.1386;    % [nM]

        % Concentration for testing basal level
        data.base_gfr_total_0 = 0.15; % [nM]

    case 'egfr_blinov',
        % Blinov et al. 2006 Biosystems. A Network Model paper parameters
        % Concentration Parameters
        data.fyn_total = 0.5;       % [nM]
        data.ptp = 0.2;             % [nM]
        data.gf_0 = 0.0166;         % [nM] gf at basal level 
        data.gf_1 = 8.065;            % [nM] gf added at time 0 - 50 ng/ml
        
        % Reaction Parameters
        data.k_catf = 1;            % [1/(nM Sec)] activation rate of Fyn by DGFR
        data.k_catb = 2;            % [1/(nM Sec)] de-activation rate of active Fyn by Phosphatase 
        data.k_gff = 0.003;         % [1/(nM Sec)] association rate between gf and gfr
        data.k_gfb = 0.06;          % [1/sec] dis-association rate of gf_gfr
        data.k_dgf = 0.01;          % [1/(nM Sec)] association rate of the growth factor dimers
        data.k_dgb = 0.1;           % [1/sec] dis-association rate of the growth factor dimers
        data.k_enf = 6e-3;         % [1/sec] internalization rate of the growth factor dimers 
        
        % Initial Concentration
        data.fyn_act_0 = 0;         % [nM]
        data.gf_gfr_0 = 0.0003;     % [nM]
        data.dgf_gfr_0 = 0;         % [nM]
        data.gfr_total_0 = 0.4;     % [nM]
        % Concentration for testing basal level
        data.base_gfr_total_0 = 0.4; % [nM]

end;

return;

