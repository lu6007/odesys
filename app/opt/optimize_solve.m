% function [sol0, sol] = optimize_solve(varargin)
% global test_optimize_ode_utility_fh;
% test_optimize_ode_utility_fh = utility();
% para_name = {'num_guess', 'model_name'};
% default_value = {1, 'simple_ode'}; 
% num_guess = 0 - use initial guess in data. 
% num_guesss = N > 0 - generate N latin hypercube sample

function [sol0, sol] = optimize_solve(varargin)

global optimize_ode_utility_fh;
optimize_ode_utility_fh = opt_utility();

para_name = {'num_guess', 'model_name'};
default_value = {1, 'simple_ode'}; 
% num_guess = 0 - use initial guess in data. 
% num_guesss = N > 0 - generate N latin hypercube sample
[num_guess, model_name] = parse_parameter(para_name, default_value, varargin);

verbose = 1;
switch model_name
    case 'simple_ode'
        model = simple_ode_model(model_name, 'verbose', verbose);
        xy_axis = [-200 2000 0 180];
        
    case 'complex_ode' 
        % best_fit only makes a difference if num_guess = 0
        model = fyn_gf_model(model_name, 'best_fit', 1, 'verbose', verbose);
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
        multiple_output = 0;
        best_fit = 0;
        model = opt_model_1118(model_name, 'multiple_output', multiple_output, ...
            'best_fit', best_fit);
%         model.ode = ode_model_1118(model_name, 'multiple_output', multiple_output, ...
%             'best_fit', best_fit);      
end
disp('Function optimize_solve():');
disp(model); 
disp(model.ode.data);
% disp(model.opt);

data = model.ode.data; 
objective_fun = model.objective; 
get_constraint_fun = model.constraint; 
constraint_initial_guess_fun = model.initial_guess; 

% Initial guess of theta
theta_name = model.theta_name;
num_theta = size(theta_name, 1);
theta = zeros(num_theta, 1);
for i = 1:num_theta
    theta(i) = data.(theta_name{i});
end

% Set up the optimization problem
theta_var = optimvar('theta', num_theta, 'LowerBound', 0.0);
object_express = fcn2optimexpr(objective_fun, theta_var); % expression
% show(theta_var); show(object_express); 
problem = optimproblem('Objective', object_express); % problem.Objective = object_express; 
problem = get_constraint_fun(problem, theta_var, model.theta_upper_bound);

% Latin hypercube sample
time_start = tic; 
fprintf('\nFunction optimize_solve(): num_guess = %d\n', num_guess); 
if num_guess == 0
    use_latin_hypercube_sample = false; 
    % theta = [0.040558595	5.720696138	0.001076275	0.195210427	1.738250076	0.697117467	340.8719999	548.4923813]';
    num_guess = 1;
    fprintf('\t set num_guess = 1\n');
else
    use_latin_hypercube_sample = true;
    X = lhsdesign(num_guess, num_theta);
end
warning('off', 'all');
fprintf('\nFunction optimize_solve(): use_latin_hypercube_sample = %d\n', use_latin_hypercube_sample);
fprintf('\t num_guess = %d\n', num_guess); 
fprintf('\t Turned warnings off\n'); 
sol0 = cell(num_guess, 1);
sol = cell(num_guess, 1); 
for i = 1:num_guess
    if use_latin_hypercube_sample
        theta = constraint_initial_guess_fun(X(i, :)', ...
            model.theta_upper_bound);
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

%     show(problem);

    % Solve the problem with MATLAB optimizer fmincon
    % option = optimoptions('fmincon', 'Display', 'iter');
    option = optimoptions('fmincon', 'Display', 'off');
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
    
    fh = optimize_ode_utility_fh;
    fh.plot_curve(sol0{i}, sol{i}, 'xy_axis', xy_axis); 
    
end % for i = 1:num_guess

time_used = toc(time_start); 
fprintf('Function test_optimize_ode(): time_used = %f (sec) \n', time_used);

end % function
