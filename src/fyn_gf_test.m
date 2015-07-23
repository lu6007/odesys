function fyn_gf_test
% fyn signal under growth factor stimulation
% Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Copyright (2015) The Regents of the University of California
% All Rights Reserved
%

% Parameters
% ----- Signaling model parameters -------
% a_1d = p(1);
% fyn_total = p(2); % Total Fyn concentration = active + inactive Fyn = constant. 
% b_1p = p(3);
% ptp = p(4); % Phosphatase concentration = constant
% a_2g = p(5);
% % Growth factor concentration is a step function of t
% if t<0,
%     gf = 0;
% elseif t>=0,
%     gf = p(6);
% end;
% b_2d = p(7); % Dis-association rate
% b_3d = p(8);

% % EGF/EGFR parameters 
% % Assume that the biosensor readout is Fyn activity
% % Assume that all the growth factors are active in dimers and inactive in
% % monomer; and that only active growth factor dimers can be internalized. 
% fyn_total = 0.5; % [uM]
% p(1) = 1;       % a_1d        [1/(uM Sec)] activation rate of Fyn by DGFR
% p(2) = fyn_total;     
% p(3) = 2;       % b_1p        [1/(uM Sec)] de-activation rate of active Fyn by Phosphatase 
% p(4) = 0.2;     % PTP         [uM]
% p(5) = 1;       % a_2g        [1/(uM Sec)] activation rate of the growth factor reccpetor by GF
% p(6) = 0.1; % 1.0;     % gf at basal level         [uM]
% p(7) = 0.1;%0.001;       % gf added at time 0         [uM] --- ask Mingxing
% p(8) = 1;    % b_2d        [1/sec] dis-association rate of the growth factor dimers
% p(9) = 0.001; % b_3d        [1/sec] internalization rate of the growth factor dimers
% % Initial conditions 
% fyn_0  = 0.0102;       % [uM]
% dgfr_0 = 0.0083;       % [uM]
% gfr_total_0 = 0.0998;     % [uM]    concentration of growth factor receptors 
% y0     = [fyn_0;dgfr_0;gfr_total_0];    

% PDGF/PDGFR parameters 
% Assume that the biosensor readout is Fyn activity
% Assume that all the growth factors are active in dimers and inactive in
% monomer; and that only active growth factor dimers can be internalized. 
fyn_total = 0.5; % [uM]
p(1) = 1;       % a_1d        [1/(uM Sec)] activation rate of Fyn by DGFR
p(2) = fyn_total;     
p(3) = 2;       % b_1p        [1/(uM Sec)] de-activation rate of active Fyn by Phosphatase 
p(4) = 0.2;     % PTP         [uM]
p(5) = 1;       % a_2g        [1/(uM Sec)] activation rate of the growth factor reccpetor by GF
p(6) = 0.1; % 1.0;     % gf at basal level         [uM]
p(7) = 0.1;       % gf added at time 0         [uM] --- ask Mingxing
p(8) = 0.1;     % b_2d        [1/sec] dis-association rate of the growth factor dimers
p(9) = 0.001; %0.00001; % b_3d        [1/sec] internalization rate of the growth factor dimers
% Initial conditions 
fyn_0  = 0.01195;       % [uM]
dgfr_0 = 0.0098;       % [uM]
gfr_total_0 = 0.0294;     % [uM]    concentration of growth factor receptors
y0     = [fyn_0;dgfr_0;gfr_total_0];    


% Options
tspan = [-300;3000];
options = odeset('RelTol',1e-5,'MaxStep',1.0,'Stats','on'); 

% Run simulation
[t,y] = ode15s(@f,tspan,y0,options,p);

yfinal = y(end,:)

fyn = y(:,1);
dgfr = y(:,2);
gfr_total = y(:,3);

figure; plot(t, fyn); title('[FYN]'); xlabel('Time (sec)');
figure; plot(t, dgfr); title('[DGFR]'); xlabel('Time (sec)');
figure; plot(t, gfr_total); title('[GFR_Total]'); xlabel('Time (sec)');

return;

% Right hand side function
function ydot = f(t, y, p)
% Initiliaze variables
fyn = y(1); % fyn kinase activity
% dgfr - concentration of dimer growth factor receptors, assuming they are active
%        and ligand bound
dgfr = y(2); 
gfr_total = y(3); % Concentration of all growth factors
% Parameters
a_1d = p(1);
fyn_total = p(2); % Total Fyn concentration = active + inactive Fyn = constant. 
b_1p = p(3);
ptp = p(4); % Phosphatase concentration = constant
a_2g = p(5);
% Growth factor concentration is a step function of t
if t<0,
    gf = p(6);
elseif t>=0,
    gf = p(6)+p(7);
end;
b_2d = p(8); % Dissociation rate
b_3d = p(9);

% Equations
fyn_dot = a_1d*dgfr*(fyn_total-fyn)-b_1p*ptp*fyn;
% mgfr - concentration of monomer growth factors 
mgfr = gfr_total-2*dgfr;
dgfr_dot = a_2g*gf* mgfr - b_2d*dgfr;
gfr_total_dot = -b_3d*dgfr;

% dy/dt
ydot = zeros(size(y));
ydot(1) = fyn_dot;
ydot(2) = dgfr_dot;
ydot(3) = gfr_total_dot;

return;

