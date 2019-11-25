% The optmization objective and constraint function for the simple ODE
% 'egfr_huang_v3'
function model_obj = simple_ode_model(model_name, varargin)
    para_name = {'verbose'};
    default_value = {1};
    verbose = parse_parameter(para_name, default_value, varargin);
    
    if verbose
        fprintf('\nFunction simple_ode_model(): model_name = %s\n', model_name);
    end 
    
    data = fyn_gf_init_data('egfr_huang_v3'); 
    data.model = model_name;
    data.scale = 428.7902; 
    
    % Optimizatoin functions and parameters
    model_obj.data = data; 
    model_obj.opt.objective = @objective;
    model_obj.opt.constraint = @constraint;
    model_obj.opt.initial_guess = @initial_guess; 
    model_obj.opt.index = 3; % [egf] = 50 ng/ml
    model_obj.opt.theta_name = {'k_caton_4'; 'k_catoff_4'; 'k_on_3'};
    model_obj.opt.theta_upper_bound = [1; 1; 1];
    
end

% Objective function for the simple ode model 'egfr_huang_v3'
function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective(theta)
global optimize_ode_utility_fh;
fh = optimize_ode_utility_fh; 

% Update the optimized parameters 
model = simple_ode_model('simple_ode', 'verbose', 0);
data = model.data;
theta_name = model.opt.theta_name; 
% sim_name = 'egfr_huang_v3'; 
% data = fyn_gf_init_data(sim_name);
% theta_name = data.opt.theta_name;
for i = 1:length(theta_name)
    data.(theta_name{i}) = theta(i);
end

% Load experimental results
exp_name = 'exp_hela_egf';
scale = data.scale; 
[t_exp_i, y_exp_i] = fh.get_exp_data(exp_name, 'scale', scale); 
% [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name); 

% Calculate simulation results
gf_con = 7.8125; % nM for 50 ng/ml
[t_ode_i, y_ode_i] = fyn_gf_ode_solve(data, 'gf_1', gf_con, 'show_figure', 0, ...
    'test_basal_level', 0); 

% Reconcile the experiment and simulation
[t_exp_interp, y_exp_interp, y_ode_interp] = ...
    fh.reconcile_exp_ode(t_exp_i, y_exp_i, t_ode_i, y_ode_i); 

% Calculate error
temp = sqrt(fh.l2norm_square(t_exp_interp, y_exp_interp)); 
f = fh.l2norm_square(t_exp_interp, (y_ode_interp - y_exp_interp)/temp); clear temp;
end

% Constraint function for the simple ode model
function problem = constraint(problem, theta_var, theta_upper_bound) 
% theta_upper_bound is used in eval()

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
function theta = initial_guess(x, theta_upper_bound)
theta = x.* theta_upper_bound;
% theta = [0; 0; 0; 400] + x.*[0.1; 0.1; 0.1; 100];
end

