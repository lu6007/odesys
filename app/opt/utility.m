% function fh = utility()
%     % Useful functions for optimizatoin
%     fh.get_exp_data = @get_exp_data;
%     fh.reconcile_exp_ode = @reconcile_exp_ode;
%     fh.plot_curve = @plot_curve; 
%     % ODE functions
%     fh.output = @output; 
%     % General functions
%     fh.get_subcell = @get_subcell; 
%     fh.l2norm_square = @l2norm_square;
%
% Example: 
%    global optimize_ode_utility_fh;
%    fh = optimize_ode_utility_fh; 
%    % Load experimental results
%     exp_name = 'exp_hela_egf'; 
%     scale = data.scale; 
%     [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name, 'index', index, 'scale', scale); 
% 

function fh = utility()
    % Useful functions for optimizatoin
    fh.get_exp_data = @get_exp_data;
    fh.reconcile_exp_ode = @reconcile_exp_ode;
    fh.plot_curve = @plot_curve; 
    % ODE functions
    fh.output = @output; 
    % General functions
    fh.get_subcell = @get_subcell; 
    fh.l2norm_square = @l2norm_square;
end

% function [t, y] = get_exp_data(exp_name, varargin)
%     data = fyn_gf_init_data(exp_name);
%     para_name = {'scale', 'index'};
%     para_default = {data.scale, 3};
% Read the experimental data and allow it to change by a value scale
% Output: t and y are vectors if index is a number
%         t and y are cells of vectors if index is a vector
function [t, y] = get_exp_data(exp_name, varargin)
    data = fyn_gf_init_data(exp_name);
    para_name = {'scale', 'index'};
    para_default = {1, 3};
    % index = 3 % egf 50 ng/ml, 4-25 ng/ml
    [scale, index] = parse_parameter(para_name, para_default, varargin);
    data.scale = scale; 

    [t_exp, y_exp] = batch_fyn_gf(data, 'show_figure', 0);  
    t = get_subcell(t_exp, index);
    y = get_subcell(y_exp, index); 
end

function [t_exp_interp, y_exp_interp, y_ode_interp] = ...
    reconcile_exp_ode(t_exp, y_exp, t_ode, y_ode)
    if ~iscell(t_exp)
        t_min = max(t_exp(1), t_ode(1)); 
        t_max = min(t_exp(end), t_ode(end));
        ii = (t_exp >= t_min & t_exp <= t_max); 
        t_exp_interp = t_exp(ii);
        y_exp_interp = y_exp(ii);
        y_ode_interp = interp1(t_ode, y_ode, t_exp_interp);
    else % t_exp is a cell of arrays
        num_exp = length(t_exp);
        t_exp_interp = cell(num_exp, 1);
        y_exp_interp = cell(num_exp, 1);
        y_ode_interp = cell(num_exp, 1);
        for i = 1:num_exp 
            t_min = max(t_exp{i}(1), t_ode{i}(1)); 
            t_max = min(t_exp{i}(end), t_ode{i}(end));
            ii = (t_exp{i} >= t_min & t_exp{i} <= t_max); 
            t_exp_interp{i} = t_exp{i}(ii);
            y_exp_interp{i} = y_exp{i}(ii);
            y_ode_interp{i} = interp1(t_ode{i}, y_ode{i}, t_exp_interp{i});
        end
    end
end

% Plot experimental results with the ode solutions from the initial guess
% and the optimized parameters theta. 
function plot_curve(sol0, sol, varargin)
    para_name = {'xy_axis'};
    default_value = {[0 1 0 1]};
    xy_axis = parse_parameter(para_name, default_value, varargin); 
    
    t_old = sol0.t;
    y_ode_old = sol0.y_ode;
    y_exp_old = sol0.y_exp;
    t = sol.t;
    y_ode = sol.y_ode;
    % y_exp = sol.y_exp; 
    
    lw = 1.5; % line_width
    if ~iscell(t_old)
        my_figure; hold on; 
        plot(t_old, y_exp_old, 'r-', 'LineWidth', lw);
        plot(t_old, y_ode_old, 'b--', 'LineWidth', lw); 
        plot(t, y_ode, 'b', 'LineWidth', 1.5);
        axis(xy_axis);
        legend('Experiment', 'Sim-Initial Guess', 'Simlation');
        xlabel('Time (s)'); ylabel('[Fyn\_active] (nM)');
    else % t_old and other variables are cells
        num_cell = length(t_old);
        my_figure; 
        for i = 1:num_cell
            subplot(2, 3, i); hold on; 
            plot(t_old{i}, y_exp_old{i}, 'r-', 'LineWidth', lw);
            plot(t_old{i}, y_ode_old{i}, 'b--', 'LineWidth', lw);
            plot(t{i}, y_ode{i}, 'b-', 'LineWidth', lw); 
            axis(xy_axis);
            if i == 1
                legend('Experiment', 'Sim-Initial Guess','Simulation');
            elseif i ==2
                xlabel('Time (s)'); ylabel('[Fyn\_active] (nM)'); 
            end
        end
    end  
end

% function output = get_subcell(input, index)
% Input : input is a cell, index is a vector
% Output: output is the entry of input{index}, if index is a number. 
%         output is a subcell of input if index is a vector
function output = get_subcell(input, index)
    num_cell = size(index, 1); 
    if num_cell ==1
        output = input{index}; 
    else % n>1
        output = cell(num_cell, 1);
        for i = 1:num_cell
            output{i} = input{index(i)};
        end
    end
end

% L2 norm square function
% y - nM, t - sec, norm_square - nM^2 * sec
function norm_square = l2norm_square(t, y)
    delta_t = diff(t); 
    n = size(t, 1); 
    mid_y = 0.5 *(y(1:n-1) + y(2:n)); 
    norm_square = sum(delta_t.*mid_y.*mid_y); 
end

% plot ode output in a figure
function yy = output(t, y, data, varargin)
    para_name = {'show_figure'};
    default_value = {1};
    show_figure = parse_parameter(para_name, default_value, varargin);
    index = data.output_index;
    yy = y(:, index); 
    if show_figure
        my_figure; plot(t, yy, 'LineWidth', 1.5); 
        legend(data.species_name{index});
        xlabel('Time (sec)'); ylabel('Concentration (nM)'); 
    end
    
end
