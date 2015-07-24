% fyn signal under growth factor stimulation 
% Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Chris Blackburn
% Copyright (2015) The Regents of the University of California
% All Rights Reserved

% Right hand side function
function ydot = rhs(t, y, data)
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
if t <0,
    gf = data.gf_0;
elseif t>=0,
    gf = data.gf_0+data.gf_1;
end;

% Reaction parameters
k_catf = data.k_catf;
k_catb = data.k_catb;
k_gff = data.k_gff;
k_gfb = data.k_gfb;
k_dgf = data.k_dgf;
k_dgb = data.k_dgb;
k_enf = data.k_enf;

% Rihgt Hand Side Equations
fyn_act_dot = k_catf*dgf_gfr*(fyn_total-fyn_act)-k_catb*ptp*fyn_act;
gf_gfr_dot = k_gff*gf*(gfr_total-gf_gfr-2*dgf_gfr)-k_gfb*gf_gfr...
           + 2*k_dgb*dgf_gfr-2*k_dgf*(gf_gfr)^2; 
dgf_gfr_dot = k_dgf*(gf_gfr)^2-k_dgb*dgf_gfr-k_enf*dgf_gfr;
gfr_total_dot = -k_enf*2*dgf_gfr;
if gfr_total_dot >0
    gfr_total_dot = 0;
end;

% dy/dt
ydot = zeros(size(y));
ydot(1) = fyn_act_dot;
ydot(2) = gf_gfr_dot;
ydot(3) = dgf_gfr_dot;
ydot(4) = gfr_total_dot;
return;

