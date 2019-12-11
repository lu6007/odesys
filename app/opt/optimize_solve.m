% function [sol0, sol] = optimize_solve(varargin)
% global test_optimize_ode_utility_fh;
% test_optimize_ode_utility_fh = utility();
% para_name = {'num_guess', 'model_name', 'model_id', 'global_optimize'};
% default_value = {1, 'simple_ode', 4, 0}; 
% num_guess = 0 - use initial guess in data. 
% num_guesss = N > 0 - generate N latin hypercube sample

function [sol0, sol, model] = optimize_solve(varargin)
global optimize_ode_utility_fh optimize_ode_model;
optimize_ode_utility_fh = opt_utility();
fh = optimize_ode_utility_fh; 

para_name = {'num_guess', 'model_name', 'model_id', 'global_optimize'};
default_value = {1, 'simple_ode', 4, 0}; 
% num_guess = 0 - use initial guess in data. 
% num_guesss = N > 0 - generate N latin hypercube sample
[num_guess, model_name, model_id, global_optimize] = ...
    parse_parameter(para_name, default_value, varargin);

verbose = 1;
disp('Function optimize_solve():');

switch model_name
    case 'simple_ode'
        model = simple_ode_model(model_name, 'verbose', verbose);
        xy_axis = [-200 2000 0 180];
    case 'complex_ode' 
        % best_fit only makes a difference if num_guess = 0
        model = fyn_gf_model(model_name, 'multiple_output', 0, ...
        'best_fit', 1, 'fyn_endo', 0, 'verbose', 1);     
        xy_axis = [-200 2000 0 200];
        
    case 'complex_ode_nodeg'
        model = fyn_gf_model_nodeg(model_name, 'best_fit', 1, 'verbose', verbose);
        xy_axis = [-200 2000 0 200];
        
    case 'complex_ode_nodeg_endo1'
        set = 1; 
        xy_axis = [-200 2000 0 250];
        model = fyn_gf_model_nodeg(model_name, 'best_fit', set, ...
            'verbose', verbose, 'fyn_endo', set);
    case 'complex_ode_nodeg_endo2'
        set = 2; 
        xy_axis = [-200 2000 0 200];
        model = fyn_gf_model_nodeg(model_name, 'best_fit', set, ...
            'verbose', verbose, 'fyn_endo', set);
    case 'complex_ode_nodeg_endo3'
        set = 3; 
        xy_axis = [-200 2000 0 200];
        model = fyn_gf_model_nodeg(model_name, 'best_fit', set, ...
            'verbose', verbose, 'fyn_endo', set);
    
    case 'complex_ode_nodeg_endo11'
        best_fit = 1;
        fyn_endo = 11;
        xy_axis = [-200 2000 0 200];
        model = fyn_gf_model_nodeg(model_name, 'best_fit', best_fit, ...
            'verbose', verbose, 'fyn_endo', fyn_endo);

    case 'complex_ode_nodeg_endo12'
        best_fit = 11;
        fyn_endo = 12;
        xy_axis = [-200 2000 0 200];
        model = fyn_gf_model_nodeg(model_name, 'best_fit', best_fit, ...
            'verbose', verbose, 'fyn_endo', fyn_endo);
        
    case 'model_1118'
        xy_axis = [-200 2000 0 30];
        model = opt_model_1118(model_name, 'model_id', model_id);
end
disp(model); 
disp(model.ode.data);
optimize_ode_model = model; 

objective_fun = model.objective; 
% get_constraint_fun = model.constraint; 
constraint_initial_guess_fun = model.initial_guess; 

% Initial guess of theta
theta_name = model.theta_name;
num_theta = size(theta_name, 1);
theta = fh.get_theta(model, theta_name); 

% Set up the optimization problem
fprintf('function optimize_solver(): optimize parameters within in the [value*0.01 value*10] interval.\n');
theta_lower_bound = model.theta_bound(:, 1);
theta_upper_bound = model.theta_bound(:, 2); 

% rng(14, 'twister'); 
% uncomment to remove randomness and get the same results. 

% option = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
option = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
% display(option);
time_start = tic; 

if global_optimize
    % global optimization
    problem = createOptimProblem('fmincon', 'objective', objective_fun, ...
        'x0', theta, 'options', option, ...
        'lb', theta_lower_bound, 'ub', theta_upper_bound); 
    go_option = GlobalSearch; % 3 trials
    % go_option = MultiStart; % 50 trials
    [xming, fming, ~, outptm, manyminsm] = run(go_option, problem);
    sol0 = [];
    sol.xming = xming;
    sol.fming = fming;
    sol.outptm = outptm;
    sol.manyminsm = manyminsm; 
else % if global_optimize
    % Local optimization
    theta_var = optimvar('theta', num_theta, 'LowerBound', theta_lower_bound', ...
        'UpperBound', theta_upper_bound'); 
    object_express = fcn2optimexpr(objective_fun, theta_var); % expression
    problem = optimproblem('Objective', object_express); % problem.Objective = object_express; 

    % Latin hypercube sample
    fprintf('\n num_guess = %d\n', num_guess); 
    if num_guess == 0
        use_latin_hypercube_sample = false; 
        % theta = [0.040558595	5.720696138	0.001076275	0.195210427	1.738250076	...
        %           0.697117467	340.8719999	548.4923813]';
        num_guess = 1;
        fprintf('\t set num_guess = 1\n');
    else
        use_latin_hypercube_sample = true;
        X = lhsdesign(num_guess, num_theta);
    end
    warning('off', 'all');
    fprintf('\t use_latin_hypercube_sample = %d\n', use_latin_hypercube_sample);
    fprintf('\t Turned warnings off\n'); 
    sol0 = cell(num_guess, 1);
    sol = cell(num_guess, 1); 
    for i = 1:num_guess
        if use_latin_hypercube_sample
            theta = constraint_initial_guess_fun(X(i, :)', model.theta_bound);
        end

        % Initial guess of the variables
        sol0{i}.theta = theta; 
        [error0, t_old, y_exp_old, y_ode_old] = objective_fun(sol0{i}.theta); 
        sol0{i}.error = error0; 
        sol0{i}.t = t_old; 
        sol0{i}.y_exp = y_exp_old;
        sol0{i}.y_ode = y_ode_old;
        fprintf('\n%d. Initial Guess:\n ', i); 
        disp(theta_name');
        disp((sol0{i}.theta)');
        fprintf('error = %f\n', sol0{i}.error); 

        % Solve the problem with MATLAB optimizer fmincon - local optimization
        [temp, error] = solve(problem, sol0{i}, 'options', option);
        sol{i}.theta = temp.theta; clear temp; 
        sol{i}.error = error; 

        disp(theta_name');
        disp((sol{i}.theta)');
        fprintf('error = %f\n', sol{i}.error); 

        % Visualize the results
        [~, t, y_exp, y_ode] = objective_fun(sol{i}.theta); 
        sol{i}.t = t;
        sol{i}.y_exp = y_exp;
        sol{i}.y_ode = y_ode; 

        % fh = optimize_ode_utility_fh;
        fh.plot_curve(sol0{i}, sol{i}, 'xy_axis', xy_axis); 

    end % for i = 1:num_guess
end % if global_optimize
time_used = toc(time_start); 
fprintf('time_used = %f (sec) \n', time_used);

end % function

