function [f, t_exp_interp, y_exp_interp, y_ode_interp] = objective_base(theta, model_fh)
    global optimize_ode_utility_fh;
    fh = optimize_ode_utility_fh; 
    
    % Update the optimized parameters 
    model = model_fh('complex_ode', 'best_fit', 1, 'multiple_output', 0, ...
        'verbose', 0);  
    data = model.data;
    theta_name = model.opt.theta_name;
    for i = 1:length(theta_name)
        data.(theta_name{i}) = theta(i);
    end

    if isfield(model.opt, 'index') % fyn_gf_model_nodeg
        index = model.opt.index; 
    else % earlier models. 
        index = 3;
    end
    
    % Load experimental results
    exp_name = 'exp_hela_egf'; 
    scale = data.scale; 
    [t_exp_i, y_exp_i] = fh.get_exp_data(exp_name, 'index', index, 'scale', scale); 
    
    % Calculate simulation results
    [t_ode_temp, y_ode_temp] = batch_fyn_gf(data, 'show_figure', 0, ...
        'test_basal_level', 0, 'rhs_function', model.rhs, 'y0', data.y0, ...
        'output_function', model.output, 'verbose', 0, 'index', index);
    num_exp = size(index, 1);
    t_ode_i = t_ode_temp;
    y_ode_i = y_ode_temp; 
    clear t_ode_temp y_ode_temp; 

    % Reconcile the experiment and simulation
    [t_exp_interp, y_exp_interp, y_ode_interp] = ...
        fh.reconcile_exp_ode(t_exp_i, y_exp_i, t_ode_i, y_ode_i); 

    % Calculate error
    f = 0;
    if num_exp == 1
        temp = sqrt(fh.l2norm_square(t_exp_interp, y_exp_interp)); 
        f = fh.l2norm_square(t_exp_interp, (y_ode_interp - y_exp_interp)/temp); clear temp;
    else % num_exp > 1
        for i = 1:num_exp
            t_i = t_exp_interp{i}; 
            y_i = y_exp_interp{i};
            temp = sqrt(fh.l2norm_square(t_i, y_i)); 
            f = f+ fh.l2norm_square(t_i, (y_ode_interp{i} - y_i)/temp); clear temp;
        end
        f = f/num_exp; 
    end
end % function