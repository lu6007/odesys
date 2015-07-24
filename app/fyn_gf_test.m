% fyn signal under growth factor stimulation 
% Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Chris Blackburn
% Copyright (2015) The Regents of the University of California
% All Rights Reserved

function [t, fyn_act] = fyn_gf_test(varargin)
parameter_name = {'model', 'fyn_act_exp', 'test_basal_level','gf_1', 'show_figure'};
default_value = {'egfr_huang', [], 0, [],1};
[model, fyn_act_exp, test_basal_level, gf_1, show_figure] = ...
    parse_parameter(parameter_name, default_value, varargin);

% Initialize parameters
data = init_data(model);
if test_basal_level,
    data.gf_1 = 0;
    y0 = [0; 0; 0; data.base_gfr_total_0];
else
    y0 = [data.fyn_act_0
        data.gf_gfr_0
        data.dgf_gfr_0
        data.gfr_total_0];
end;

if ~isempty(gf_1),
    data.gf_1 = gf_1;
end;

% Run simulation
tspan = [-300;3000];
options = odeset('RelTol',1e-5,'MaxStep',1.0,'Stats','on'); 
[t,y] = ode15s(@rhs,tspan,y0,options,data);
yfinal = y(end,:)

fyn_act = y(:,1);
gf_gfr = y(:,2);
dgf_gfr = y(:,3);
gfr_total = y(:,4);

% Plot results
if show_figure,
    figure; plot(t, fyn_act); title('[FYN ACT]'); xlabel('Time (sec)');
    figure; plot(t/60, fyn_act*10+0.14); title('[FYN ACT]'); xlabel('Time (min)');
    if ~isempty(fyn_act_exp),
        hold on; plot(fyn_act_exp(:,1)-5, fyn_act_exp(:,2),'ro');
    end;
    figure; plot(t, gf_gfr); title('[GF GFR]'); xlabel('Time (sec)');
    figure; plot(t, dgf_gfr); title('[DGF GFR]'); xlabel('Time (sec)');
    figure; plot(t, gfr_total); title('[GFR Total]'); xlabel('Time (sec)');
end
return;


