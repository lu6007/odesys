% fyn signal under growth factor stimulation 
% Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Christopher Blackburn
% Copyright (2015) The Regents of the University of California
% All Rights Reserved

% Right hand side function
function ydot = fyn_gf_rhs(t, y, data)

    ydot = egfr_huang_rhs(t, y, data);

return; 

function ydot = egfr_huang_rhs(t, y, data) 
% Initiliaze variables
fyn_act = y(1); % fyn kinase activity
gf_gfr = y(2);
% dgf_gfr - concentration of dimers of growth-factor-bound growth factor
% receptors. 
dgf_gfr = y(3); 
gfr_total = y(4); % Concentration of all growth factors

% Concentration Parameters
fyn_total = data.fyn_total;
ptp = data.ptp;
if t <0
    gf = data.gf_0;
elseif t>=0
    % Add growth factor at a certain concentration
    gf = data.gf_0+data.gf_1;
end

% Reaction parameters
k_caton_4 = data.k_caton_4;
k_don_4 = data.k_don_4;
k_catoff_4 = data.k_catoff_4;
k_doff_4 =data.k_doff_4;
k_on_1 = data.k_on_1;
k_off_1 = data.k_off_1;
k_on_2 = data.k_on_2;
k_off_2 = data.k_off_2;
k_on_3 = data.k_on_3;

% Reaction flux
gfr = gfr_total-gf_gfr-2*dgf_gfr;
fyn = fyn_total-fyn_act;
v1 = k_on_1*gf*gfr-k_off_1*gf_gfr; 
v2 = k_on_2*gf_gfr*gf_gfr - k_off_2*dgf_gfr;
v3 = k_on_3*dgf_gfr;
v4 = k_caton_4*dgf_gfr*fyn/(k_don_4+fyn) - k_catoff_4*ptp*fyn_act/(k_doff_4+fyn_act);

% Right Hand Side Equations
fyn_act_dot = v4;
gf_gfr_dot = v1 - 2*v2; 
dgf_gfr_dot = v2-v3;
gfr_total_dot = -2*v3;
if gfr_total_dot >0
    gfr_total_dot = 0;
end

% dy/dt
ydot = zeros(size(y));
ydot(1) = fyn_act_dot;
ydot(2) = gf_gfr_dot;
ydot(3) = dgf_gfr_dot;
ydot(4) = gfr_total_dot;
return;

