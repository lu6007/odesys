% fyn signal under growth factor stimulation 
% function [t, fyn_act] = fyn_gf_ode_solve(data, varargin)
% parameter_name = {'fyn_act_exp', 'test_basal_level','gf_1', 'show_figure', ...
%     'rhs_function','y0', 'output_function'};
% default_value = {[], 1, data.gf_1,1, @fyn_gf_rhs, [], []};
%
% Example:
% >> data = init_data('egfr_huang');
% >> [t, fyn_act] = fyn_gf_ode_solve(data, 'model', 'egfr_huang_v2');

% Authors: Shaoying Lu Shirley Wu Mingxing Ouyang Yingxiao Wang
% Chris Blackburn
% Copyright (2015) The Regents of the University of California
% All Rights Reserved

function [t, output] = fyn_gf_ode_solve(data, varargin)
parameter_name = {'fyn_act_exp', 'test_basal_level','gf_1', 'show_figure', ...
    'rhs_function','y0', 'output_function'};
default_value = {[], 1, data.gf_1,1, @fyn_gf_rhs, [], []};
[~, test_basal_level, gf_1, show_figure, rhs_fh, y0, output_fh] = ...
    parse_parameter(parameter_name, default_value, varargin);

if ~isempty(gf_1)
    data.gf_1 = gf_1;
end

% tspan = [-300;3000]; % seconds
% tspan = [-300; 2000];
tspan = [-10; 2000]; % seconds

if isempty(y0)
    % Initialize parameters
    if test_basal_level
        gf_1_save = data.gf_1;
        data.gf_1 = 0;
        % Initial values for fyn_act; gf_gfr; dgf_gfr; gfr_total
        y0 = [0; 0; 0; data.base_gfr_total_0];
        % Run simulation
        options = odeset('RelTol',1e-5,'MaxStep',1.0,'Stats','off'); 
        [~,y] = ode15s(rhs_fh,tspan,y0,options,data);
        y0 = y(end,:);
        clear y options 
        data.gf_1 = gf_1_save;
    else
        y0 = [data.fyn_act_0
            data.gf_gfr_0
            data.dgf_gfr_0
            data.gfr_total_0];
    end
end

% Run simulation
options = odeset('RelTol',1e-5,'MaxStep',1.0,'Stats','off'); 
% ode23s also has lots of warnings sometime. 
[t,y] = ode15s(rhs_fh,tspan,y0,options,data);
% yfinal = y(end,:);

% Plot results
if isempty(output_fh) 
    fyn_act = y(:,1); 
    gf_gfr = y(:,2);
    dgf_gfr = y(:,3);
    gfr_total = y(:,4);
    
    if show_figure
        lw = 1.5; 
        my_figure; plot(t, fyn_act, 'LineWidth', lw); title('[FYN ACT]'); xlabel('Time (sec)');
        my_figure; plot(t, gf_gfr, 'LineWidth', lw); title('[GF GFR]'); xlabel('Time (sec)');
        my_figure; plot(t, dgf_gfr, 'LineWidth', lw); title('[DGF GFR]'); xlabel('Time (sec)');
        my_figure; plot(t, gfr_total, 'LineWidth', lw); title('[GFR Total]'); xlabel('Time (sec)');
    end
    output = fyn_act; 
elseif ~isempty(output_fh)
    output = output_fh(t, y, data, 'show_figure', show_figure); 
end
return;

% [t, output] = fyn_gf_output(t, y, show_figure)
% return
