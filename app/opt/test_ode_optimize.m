% test_optimize_ode.m
% 
% Example:
% >> test_optimize_ode; 
% 
% To profile:
% >> profile on;
% >> test_ode_optimize;
% >> profile viewer; % the batch_fyn_gf function uses very little time
% >> profile off; 
%% Run simple simulation and load experimental results
% cd /Users/kathylu/Documents/sof/odesys/app/opt

global optimize_ode_utility_fh; % optimize_ode_model;
optimize_ode_utility_fh = opt_utility(); 
fh = optimize_ode_utility_fh; 

%% Optimize for the a simple model  
model_name = 'simple_ode'; % 'complex_ode_nodeg'; % 'simple_ode'
[sol0, sol] = optimize_solve('num_guess', 0, 'model_name', model_name); 

%% Optimize the complex mode with degradation
num_guess = 0; % 0 - draw the best fit solution; N - optimize for N initial guesses
model_name = 'complex_ode'; % 'complex_ode_nodeg'; %'simple_ode'; % 'complex_ode'; 
[sol0, sol] = optimize_solve('num_guess', num_guess, 'model_name', model_name); 
title('Complex ODE with Degradation');

%% Optmize the complex model with no degradation 
num_guess = 0;
model_name = 'complex_ode_nodeg'; % fyn_endo = 0
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name); 
title('Complex ODE with no Degradation'); 
% 
num_col = 2* size(sol0{1}.theta, 1)+2; 
temp = zeros(num_guess, num_col);
for i = 1: num_guess
    temp(i,:) = [sol0{i}.theta' sol0{i}.error sol{i}.theta' sol{i}.error];
end
% Copy temp from MABLAB variables to excel

%% Run the best fit complex-nodeg model with fyn_act ednocytosis
model = complex_ode_nodeg('complex_ode_nodeg', 'best_fit', 1, 'fyn_endo', 1);
batch_fyn_gf(model.data, 'rhs_function', model.rhs, 'y0', model.data.y0, ...
'output_function', model.output);

%% Optmize the complex model with no degradation and endocytosis of Fyn kinase
num_guess = 0;
model_name = 'complex_ode_nodeg_endo1'; % fyn_endo
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name); 
title('Complex ODE with no Degradation and Active Fyn Endocytosis'); 

%% Optmize the complex model with no degradation and endocytosis of Fyn kinase
num_guess = 0;
model_name = 'complex_ode_nodeg_endo2'; % fyn_endo
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name); 
title('Complex ODE with no Degradation and Active Fyn Endocytosis'); 

%% Optmize the complex model with no degradation and endocytosis of Fyn kinase
% fit concentration dependence
num_guess = 0;
model_name = 'complex_ode_nodeg_endo3'; % fyn_endo
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name); 
title('Complex ODE with no Degradation and Active Fyn Endocytosis'); 

%% Optmize the complex model with no degradation and endocytosis of Fyn kinase
num_guess = 0;
endo_str = '11';
model_name = ['complex_ode_nodeg_endo', endo_str]; % fyn_endo
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name); 
title('Complex ODE with no Degradation and Active Fyn Endocytosis'); 

%% Optmize the complex model with no degradation and endocytosis of Fyn kinase
num_guess = 0;
endo_str = '12';
model_name = ['complex_ode_nodeg_endo', endo_str]; % fyn_endo
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name); 
title('Complex ODE with no Degradation and Active Fyn Endocytosis'); 

%% Copy the results to temp
clear temp;
temp = fh.get_theta_error(sol0, sol); 
% Copy temp from MABLAB variables to excel

%% Optmize ode_model_1118.m
num_guess = 0;
model_name = 'model_1118'; 
model_id = 6; % 1; 2; 3; 4; 5; 6; 
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name, ...
    'model_id', model_id); 
title('Model 1118'); 

%% Test the best fitted model for model_1118_6
model_name = 'model_1118'; 
model = opt_model_1118(model_name, 'model_id', 6);
model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
ode = model.ode; 
batch_fyn_gf(ode.data, 'multiple_output', 0, 'best_fit', 0, 'verbose', 1, ...
    'rhs_function', ode.rhs, 'y0', ode.data.y0, ...
'output_function', ode.output);

%% Optmize ode_model_1118 for concentration dependence
num_guess = 0;
model_name = 'model_1118'; 
model_id = 9; % 7; 8; 9;  
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name, ...
    'model_id', model_id); 

%% Optmize ode_model_1118 for concentration dependence
num_guess = 0;
model_name = 'model_1118'; 
model_id = 12; % 10, 13, 14, 15 
% rng(1); 
[sol0, sol] = optimize_solve('num_guess',num_guess, 'model_name', model_name, ...
    'model_id', model_id, 'global_optimization', 0); 


%% Run the best fit model_1118 with different parameter values in batch
model_id = 12; % Initialize to the best fit parameters for model 9
model = opt_model_1118('model_1118', 'model_id', model_id);
% The control group
model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
ode = model.ode; 
batch_fyn_gf(ode.data, 'rhs_function', ode.rhs, 'y0', ode.data.y0, ...
'output_function', ode.output);

% The probed group
% field_name = {'kon_2', 'koff_2', 'kon_4', 'koff_4', 'kon_7', 'koff_7'};
% field_value = {1e+1, 1e-1, 1e-6, 1e6, 1e1, 1e-1}; 
% field_name = {'kon_1', 'koff_1', 'kon_3', 'kcatoff_3', 'kdoff_3'};
% field_value = {1e1, 1e-1, 1e1, 1e-2, 1e-2}; 
% field_name = {'kon_5', 'vmaxoff_6', 'kmoff_6'};
% field_value = {1e1, 1e-2, 0.5}; 
field_name = {'koff_2', 'koff_2', 'gfr_total'};
field_value = {1e+2, 1e+4, 1e3}; 
for i = 1:length(field_name)
    model = fh.set_model_theta(model, model.theta_name, model.theta_fit);
    ode = model.ode; 
    batch_fyn_gf(ode.data, 'rhs_function', ode.rhs, 'y0', ode.data.y0, ...
    'output_function', ode.output, 'probe_field', field_name{i}, ...
    'probe_factor', field_value{i});
end
