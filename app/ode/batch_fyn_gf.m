% >> model = 'egfr_huang'; %'egfr_huang_v2'; % 'egfr_blinov'; 
% %'egfr_huang_v3'; 
% >> data = fyn_gf_init_data(model);
% >> [t, y] = batch_fyn_gf(data);
function [t, output] = batch_fyn_gf(data, varargin)
para_name = {'show_figure', 'rhs_function', 'y0', 'output_function', 'field_name', 'verbose', ...
    'probe_field', 'probe_factor', 'index'};
default_value = {1, [], [], [], 'gf_1', 1, '', 1, []};
[show_figure, rhs_fh, y0, output_fh, field_name, verbose, probe_field, probe_factor, index] = ...
    parse_parameter(para_name, default_value, varargin);

% concentration string
conc_str = {'500 ng/ml', '100 ng/ml', '50 ng/ml', '25 ng/ml', '10 ng/ml', '5 ng/ml', '1 ng/ml'};

switch data.model
    case {'exp_hela_egf', 'exp_hela_egf_2', 'exp_hela_egf_3'} 
    % Load experimental training data
        res = load(data.file);
        n = length(res.sheets);
        t = cell(n,1);
        y = cell(n,1);
        for i = 1:n
            t{i} = res.times{n+1-i}*60; % convert to seconds
            y{i} = (res.means{n+1-i}-1)*data.scale; % convert to nM
            % Non linear functions could have better fit with smaller error
            % y{i} = sqrt(res.means{n+1-i}-1)*data.scale; % convert to nM nonlinearly
        end
        output = y; 
        clear res;
    otherwise % 'simple_ode', 'complex_ode', 'complex_ode_nodeg'
        % Run simulations with different concentration
        switch field_name
            case 'gf_1' % conc_str already defined above
                % EGF molecular weight = 2600/mol
                % 1 ng/ml EGF = 1/6400 nmol/ ml = 1000/6400 nM = 0.15625 nM
                conc_value = [78.125, 15.625, 7.8125, 3.90625, 1.5625, 0.78125, 0.15625];
                % conc_value = [78.125, 15.625, 7.8125, 3.90625, 1.5625, 0.78125, 0.15625]; % nM
%             case {'kon_2'}
%                 conc_str = {'1', '10', '1e2', '1e3', '1e4', '1e5', '1e6'};
%                 conc_value = data.(field_name) * (10.^(0:6))';
%             case {'koff_2'}
%                 conc_str = {'1e3', '1e4', '1e5', '1e6', '1e7', '1e8','1e9'};
%                 conc_value = data.(field_name) * (10.^(3:9))';
            otherwise %case {'kon_4', 'kcaton_7', 'kdon_7', 'kcatoff_7', 'kdoff_7'}
                conc_str = {'1e3', '1e2', '1e1', '1', '1e-1', '1e-2', '1e-3'};
                conc_value = data.(field_name) * (10.^(3:-1:-3))';
        end
        
        % Only run a portion of experiments chosen by the variable index
        if ~isempty(index)
            num_exp = size(index, 1);
            temp_str = conc_str; temp_value = conc_value; 
            clear conc_str conc_value;
            conc_str = cell(num_exp, 1); conc_value = zeros(num_exp, 1); 
            for i = 1:num_exp
                conc_str{i} = temp_str{index(i)};
                conc_value(i) = temp_value(index(i));
            end
            clear temp_str temp_value;
        end
                 
        % probe parameter field
        if ~isempty(probe_field)
            value = data.(probe_field); 
            if verbose
                fprintf('data.%s = %10f * %6e\n', probe_field, value, probe_factor); 
            end
            data.(probe_field) = value * probe_factor;
        end % if ~isempty(probe_field)
        % data.ptp_total = data.ptp_total*10; 

        
        %  
        num_exp = length(conc_str); 
        t = cell(num_exp,1);
        output = cell(num_exp, 1); 
        % Run different GF concentration
        y = cell(num_exp,1);
        for i = 1:num_exp
            % batch parameter field
            data.(field_name) = conc_value(i); 
            if verbose
                fprintf('data.%s = %f\n', field_name, conc_value(i));
            end
            if isempty(rhs_fh)
                [t{i}, y{i}] = fyn_gf_ode_solve(data, 'show_figure', 0, ...
                    'test_basal_level', 0); 
                output = y; 
            else
                [t{i}, output{i}] = fyn_gf_ode_solve(data, 'test_basal_level', 0, ...
                    'show_figure', 0, 'rhs_function', rhs_fh, 'y0', y0, ...
                    'output_function', output_fh);
            end 
        end % for i = 1:num_exp
end % switch data.model

% Make the plot
num_exp = length(conc_str); 
line_type = {'r-', 'g-', 'b-', 'k-', 'ro', 'go', 'bo'};
fs = 18; lw = 2.0; 
if show_figure 
    if isempty(output_fh)
        y_label = 'Active Fyn Kianse (nM)';
        my_figure('font_size', fs, 'line_width', lw); hold on;
        for i = 1:num_exp
            plot(t{i}/60, output{i}, line_type{i}, 'LineWidth', lw);
            xlabel('Time (min)'); ylabel(y_label);
        end
        legend(conc_str);
        title(regexprep(data.model, '_','\\_'));
    else % has output_fh
        index = data.output_index;
        num_species = length(index); 
        for j = 1:num_species
            my_figure('font_size', fs, 'line_width', lw); hold on;
            for i = 1:num_exp
                plot(t{i}/60, output{i}(:,j), line_type{i}, 'LineWidth', lw);
                xlabel('Time (min)'); 
                % ylabel(data.species_name{index(j)});
                ylabel([data.species_name{index(j)}, ' (nM)']);
            end
            legend(conc_str);
            title(regexprep(field_name, '_','\\_'));
        end
    end
end % if show_figure

% backward compatible
if num_exp == 1
    temp_t = t; temp_output = output; clear t output; 
    t = temp_t{1}; output = temp_output{1};
end

return;