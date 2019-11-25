% test_odesys.m
global optimize_ode_utility_fh;
optimize_ode_utility_fh = opt_utility();
fh = optimize_ode_utility_fh; 

%% Run simple simulation and load experimental results
% cd /Users/kathylu/Documents/sof/odesys/app
% Simple ODE models
data = fyn_gf_init_data('egfr_huang_v2');
batch_fyn_gf(data);
%
data = fyn_gf_init_data('exp_hela_egf');
data.scale = 1;
batch_fyn_gf(data);
%
data = fyn_gf_init_data('egfr_huang_v3');
fyn_gf_ode_solve(data); 
[t, output] = batch_fyn_gf(data);
%
data = fyn_gf_init_data('exp_hela_egf');
data.scale = 1;
batch_fyn_gf(data);

%% Run simulation of a complex model
name = 'complex_ode'; 
model= complex_ode(name, 'best_fit', 0); 
fyn_gf_ode_solve(model.data, 'test_basal_level', 0, 'show_figure', 1, ...
'rhs_function', model.rhs, 'y0', model.data.y0, 'output_function', model.output);
title('Complex ODE');
%
model = complex_ode(name, 'best_fit', 1); 
batch_fyn_gf(model.data, 'rhs_function', model.rhs, 'y0', model.data.y0, ...
'output_function', model.output);
% 
model= complex_ode(name, 'multiple_output', 1, 'best_fit', 1, 'verbose', 1);
fyn_gf_ode_solve(model.data, 'test_basal_level', 0, 'show_figure', 1, ...
'rhs_function', model.rhs, 'y0', model.data.y0, 'output_function', model.output);
title('Complex ODE Best Fit Multiple Output'); 

%% Run the complex model with no degradation
name = 'complex_ode_nodeg'; 
model = complex_ode_nodeg(name, 'best_fit', 1);
[t, output]= fyn_gf_ode_solve(model.data, 'test_basal_level', 0, 'show_figure', 1, ...
'rhs_function', model.rhs, 'y0', model.data.y0, 'output_function', model.output);
%
batch_fyn_gf(model.data, 'rhs_function', model.rhs, 'y0', model.data.y0, ...
'output_function', model.output); 

%% Run the best fit complex-nodeg model with different parameter values in batch
model = complex_ode_nodeg('complex_ode_nodeg', 'best_fit', 1);
field_name = {'kon_2', 'koff_2', 'kon_4', 'kcaton_7', 'kdon_7', 'kcatoff_7', 'kdoff_7'};
for i = 1:length(field_name)
    batch_fyn_gf(model.data, 'rhs_function', model.rhs, 'y0', model.data.y0, ...
    'output_function', model.output, 'field_name', field_name{i});
end

%% Run the best fit complex-nodeg model and probe the effect of parameters on 
% concentration dependence. 
model = complex_ode_nodeg('complex_ode_nodeg', 'best_fit', 1, 'fyn_endo', 1);
probe_field = {'kon_1', 'koff_1'};
probe_factor = [1e-3; 1e3];
for i = 1:length(probe_field)
    model.data.output_index = [3;4; 7;8]; % dgfgfr, fyn, sensor
    batch_fyn_gf(model.data, 'rhs_function', model.rhs, 'y0', model.data.y0, ...
    'output_function', model.output, 'probe_field', probe_field{i}, ...
    'probe_factor', probe_factor(i));
end

%% Run the best fit complex-nodeg model with fyn_act ednocytosis
best_fit = 11; fyn_endo = 11; 
model = complex_ode_nodeg('complex_ode_nodeg', 'best_fit', best_fit, 'fyn_endo', fyn_endo);
batch_fyn_gf(model.data, 'rhs_function', model.rhs, 'y0', model.data.y0, ...
'output_function', model.output);

%% Run the best fit complex-nodeg model with fyn_act ednocytosis
best_fit = 11; fyn_endo = 12; 
model = complex_ode_nodeg('complex_ode_nodeg', 'best_fit', best_fit, 'fyn_endo', fyn_endo);
data = model.data; 
batch_fyn_gf(data, 'rhs_function', model.rhs, 'y0', data.y0, ...
'output_function', model.output);

%% Run simulation of ode_model_1118
name = 'ode_model_1118'; 
% Single output
ode_single_output = ode_model_1118(name, 'best_fit', 0); 
ode = ode_single_output; 
fyn_gf_ode_solve(ode.data, 'test_basal_level', 0, 'show_figure', 1, ...
'rhs_function', ode.rhs, 'y0', ode.data.y0, 'output_function', ode.output);
title(fh.get_latex(['Model', name]));

% Multiple output
ode= ode_model_1118(name, 'multiple_output', 1, 'best_fit', 0, 'verbose', 1);
fyn_gf_ode_solve(ode.data, 'test_basal_level', 0, 'show_figure', 1, ...
'rhs_function', ode.rhs, 'y0', ode.data.y0, 'output_function', ode.output);
title(fh.get_latex(['Model ', name, ' Multiple Output'])); 
%
ode = ode_single_output; 
batch_fyn_gf(ode.data, 'multiple_output', 0, 'best_fit', 0, 'verbose', 1, ...
    'rhs_function', ode.rhs, 'y0', ode.data.y0, ...
'output_function', ode.output);