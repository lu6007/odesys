function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 
    
    % Update the optimized parameters 
    model = model_fh('complex_ode', 'best_fit', 1, 'multiple_output', 0, ...
        'verbose', 0);  
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
end % function