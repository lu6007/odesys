% function fh = opt_utility()
%     % Useful functions for optimizatoin
%     fh.objective = @objective; 
%     fh.constraint = @constraint;
%     fh.initial_guess = @initial_guess; 
%     fh.get_exp_data = @get_exp_data;
%     fh.reconcile_exp_ode = @reconcile_exp_ode;
%     fh.plot_curve = @plot_curve; 
%     fh.set_model_theta = @set_model_theta; 
%     fh.get_theta = @get_theta; 
%     fh.get_theta_bound = @get_theta_bound; 
%     fh.get_theta_error = @get_theta_error; 
%     % ODE functions
%     fh.output = @output; 
%     % General functions
%     fh.get_subcell = @get_subcell; 
%     fh.normalize_l2norm_square_difference = @normalize_l2norm_square_difference; 
%     fh.l2norm_square = @l2norm_square;
%     fh.get_latex = @get_latex; 
%
% Example: 
%    global optimize_ode_utility_fh;
%    fh = optimize_ode_utility_fh; 
%    % Load experimental results
%     exp_name = 'exp_hela_egf'; 
%     scale = data.scale; 
%     [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name, 'index', index, 'scale', scale); 
% 

% Author: Shaoying Lu , shaoying.lu@gmail.com
function fh = opt_utility()
    % Useful functions for optimizatoin
    fh.objective = @objective; 
    fh.constraint = @constraint;
    fh.initial_guess = @initial_guess; 
    fh.get_exp_data = @get_exp_data;
    fh.reconcile_exp_ode = @reconcile_exp_ode;
    fh.plot_curve = @plot_curve; 
    fh.set_model_theta = @set_model_theta; 
    fh.get_theta = @get_theta; 
    fh.get_theta_bound = @get_theta_bound; 
    fh.get_theta_error = @get_theta_error; 
    % ODE functions
    fh.output = @output; 
    % General functions
    fh.get_root = @get_root; 
    fh.get_subcell = @get_subcell; 
    fh.normalize_l2norm_square_difference = @normalize_l2norm_square_difference; 
    fh.l2norm_square = @l2norm_square;
    fh.get_latex = @get_latex; 
end

% Calculate the objective function
function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective(theta)
    
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 
    
    global optimize_ode_model; 
    model = optimize_ode_model; 
    
    % Update the optimized parameters 
%     model = model_fh('complex_ode', 'best_fit', 1, 'multiple_output', 0, ...
%         'verbose', 0);  
    model = fh.set_model_theta(model, model.theta_name, theta); 
    data = model.ode.data; 

    if isfield(model, 'index') % fyn_gf_model_nodeg
        index = model.index; 
    else % earlier models. 
        index = 3;
    end
    
    % Load experimental results
    exp_name = 'exp_hela_egf'; 
    scale = model.scale; 
    [t_exp, y_exp] = fh.get_exp_data(exp_name, 'index', index, 'scale', scale); 
    
    % Calculate simulation results
    [t_ode, y_ode] = batch_fyn_gf(data, 'show_figure', 0, ...
        'test_basal_level', 0, 'rhs_function', model.ode.rhs, 'y0', data.y0, ...
        'output_function', model.ode.output, 'verbose', 0, 'index', index);

    % Reconcile the experiment and simulation
    [t_exp_interp, y_exp_interp, y_ode_interp] = ...
        fh.reconcile_exp_ode(t_exp, y_exp, t_ode, y_ode); 

    % Calculate error
    f = fh.normalize_l2norm_square_difference(t_exp_interp, y_exp_interp, y_ode_interp);  
end

% Constraint function for the simple ode model
function problem = constraint(problem, theta_var, theta_upper_bound) 
% theta_upper_bound is used in eval()
% The function is not used in the package yet

    % problem.Constraints.thetai_high = (theta_var(i)<= theta_high(i)); 
    % theta(i) <= theta_high(i)
    num_theta = size(theta_var, 1);
    for i = 1:num_theta
        cmd1 = sprintf('problem.Constraints.theta%d_high = ', i);
        cmd2 = sprintf('(theta_var(%d) <= theta_upper_bound(%d));', i, i);
        eval([cmd1 cmd2]); 
    end
end

% Constraint initial guess function for the simple ode model
function theta = initial_guess(x, theta_bound)
    theta = theta_bound(:,1) + x.* theta_bound(:,2);
    % theta = lower_bound + x.*upper_bound;
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
    
    fs = 18; % font_size
    lw = 2.0; % line_width
    if ~iscell(t_old)
        my_figure('font_size', fs, 'linewidth', lw); hold on; 
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

% function model = set_model_theta(model, theta_name, theta)
% copy parameter values from theta to ode.data.theta_name and
% opt_model.scale
function model = set_model_theta(model, theta_name, theta)
    data = model.ode.data; 
    for i = 1:length(theta_name)
        switch theta_name{i}
            case 'scale'
                model.scale = theta(i);
            otherwise
                data.(theta_name{i}) = theta(i);
        end
    end
    model.ode.data = data; 
end

% function theta_bound = get_theta_bound(model, theta_name)
% returns the lower (column 1) and upper bound (column 2) of theta values based on their current
% values. 
function theta = get_theta(model, theta_name)
    data = model.ode.data; 
    num_theta = length(theta_name);
    theta = zeros(num_theta, 1); 
    for i = 1:length(theta_name)
        switch theta_name{i}
            case 'scale'
                theta(i) = model.scale;
            otherwise
                theta(i) = data.(theta_name{i}); 
        end
    end
end


% function theta_bound = get_theta_bound(model, theta_name)
% returns the lower (column 1) and upper bound (column 2) of theta values based on their current
% values. 
function theta_bound = get_theta_bound(model, theta_name, varargin)
    para_name = {'bound_factor'};
    para_default = {10};
    bound_factor = parse_parameter(para_name, para_default, varargin);
    theta = get_theta(model, theta_name); 
    theta_bound = [theta/bound_factor/bound_factor theta*bound_factor];  

end

function temp = get_theta_error(sol0, sol)
    num_guess = length(sol0);
    num_col = 2* size(sol0{1}.theta, 1)+2; 
    temp = zeros(num_guess, num_col);
    for i = 1: num_guess
        temp(i,:) = [sol0{i}.theta' sol0{i}.error sol{i}.theta' sol{i}.error];
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

% The l2 norm squared of (y2-y1) normalized by the squared l2 norm of y1
% Also normalized by the number of groups if multiple groups if the input
% contains multiple groups of data in cell form. 
function diff = normalize_l2norm_square_difference(t, y1, y2)
    diff = 0; 
    if ~iscell(t)
        temp = sqrt(l2norm_square(t, y1)); 
        diff = l2norm_square(t, (y2 - y1)/temp); clear temp;
    else % iscell(t) 
        num_exp = length(t); 
        for i = 1:num_exp
            t_i = t{i}; 
            y_i = y1{i};
            temp = sqrt(l2norm_square(t_i, y_i)); 
            diff = diff + l2norm_square(t_i, (y2{i} - y_i)/temp); clear temp;
        end
        diff = diff/num_exp; 
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

% get the string for display '_' in a figure
function latex_str = get_latex(str)
    latex_str = regexprep(str, '\_', '\\_'); 
end

% get the root direction for data
function root = get_root()
    system_name = computer;
    switch  system_name
        case {'MACI64', 'MACI32'}
            root = '/Users/kathylu/Documents/sof/odesys/data/';
        case {'PCWIN64', 'PCWIN32'}
            root = 'C:/Users/kalu/Documents/sof/odesys/data/';
        otherwise
            fprintf('\nFunction opt_utililty.get_root(): \n');
            fprintf('Error: system_name = %s, root not defined! \n', ...
                system_name);
    end
    % fprintf('root = %s\n', root); 
end
