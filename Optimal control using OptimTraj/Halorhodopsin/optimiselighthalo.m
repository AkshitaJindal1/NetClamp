close all;
clear;
clc;

%==========================================================================
% Purpose:
% Determine what types of synaptic currents (parameterised by rise and
% decay) are achievable given the parameters of the opsin.
%==========================================================================

% Define rise and decay values to match the target image
rise_values = [0.01, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, ...
               0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, ...
               0.8, 0.85, 0.9, 1.0];

decay_values = [1.0, 0.95, 0.9, 0.8, ...  
                0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, ...
                0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, ...
                0.075, 0.05, 0.025, 0.01];

%==========================================================================
% Initialize
%==========================================================================

heatmap_data = NaN(length(decay_values), length(rise_values));

addpath(genpath('OptimTraj'));

save_folder = 'figures_new';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Counters
nan_count = 0;
successful_count = 0;
total_count = numel(heatmap_data);

%==========================================================================
% Loop over rise and decay combinations
%==========================================================================

for i = 1:length(decay_values)

    for j = 1:length(rise_values)

        rise  = rise_values(j);
        decay = decay_values(i);

        fprintf('\n=============================================\n');
        fprintf('Rise = %.4g | Decay = %.4g\n', rise, decay);
        fprintf('=============================================\n');

        % Skip invalid values
        if rise == 0 || decay == 0
            heatmap_data(i,j) = NaN;
            nan_count = nan_count + 1;
            continue;
        end

        %==================================================================
        % Target synaptic current
        %==================================================================

        amp = -750;
        duration = 200;

        % Handle rise == decay separately
        if abs(rise - decay) < 1e-12

            syn_signal = @(t) amp .* rise .* t .* exp(-rise .* t);

        else

            syn_signal = @(t) ...
                amp ./ (1/rise - 1/decay) .* ...
                (exp(-rise.*t) - exp(-decay.*t));

        end

        %==================================================================
        % Plot target
        %==================================================================

        dt = 0.1;
        t = (0:dt:duration)';

        fig = figure('Visible','off');
        ax = axes(fig);
        hold(ax,'on');

        p = plot(ax,t,syn_signal(t), ...
            'LineWidth',2.0, ...
            'HandleVisibility','off');

        plot(ax,t(1:10:end),syn_signal(t(1:10:end)), ...
            'LineStyle','none', ...
            'Marker','.', ...
            'MarkerSize',20, ...
            'Color',get(p,'Color'), ...
            'DisplayName','Target');

        xlabel(ax,'Time (ms)');
        ylabel(ax,'Current (pA)');
        set(ax,'FontSize',20);
        legend(ax);

        %==================================================================
        % Set up model
        %==================================================================

        [modelpars,varpars,sysPars,fun] = ...
            setDefaultParameters_opsinmodel_ODE();

        % Hard-coded parameters
        G = 6.4096082676895785;
        gamma = 0.05;
        V = -80;
        E = 0.0;

        I_opsin = @(x) ...
            (G .* (x(1,:) + gamma.*x(3,:)) .* (V-E));

        problem.func.dynamics = ...
            @(t,x,u) opsinmodel_ODE_ctrl(t,x,varpars,u);

        problem.func.pathObj = ...
            @(t,x,u) (I_opsin(x)-syn_signal(t)).^2;

        %==================================================================
        % Initial condition
        %==================================================================

        init_cond = [0.0; 0.8; 0.1];

        try
            x0 = fsolve( ...
                @(x) problem.func.dynamics(0,x,0), ...
                init_cond);

        catch
            warning('fsolve failed for rise = %.4g, decay = %.4g', ...
                rise,decay);

            heatmap_data(i,j) = NaN;
            nan_count = nan_count + 1;

            close(fig);
            continue;
        end

        %==================================================================
        % Bounds
        %==================================================================

        problem.bounds.initialTime.low = 0;
        problem.bounds.initialTime.upp = 0;

        problem.bounds.finalTime.low = duration;
        problem.bounds.finalTime.upp = duration;

        problem.bounds.initialState.low = x0;
        problem.bounds.initialState.upp = x0;

        problem.bounds.finalState.low = zeros(3,1);
        problem.bounds.finalState.upp = ones(3,1);

        problem.bounds.state.low = zeros(3,1);
        problem.bounds.state.upp = ones(3,1);

        problem.bounds.control.low = 0.0;
        problem.bounds.control.upp = 1.0;

        %==================================================================
        % Initial guess
        %==================================================================

        problem.guess.time = [0,duration];
        problem.guess.state = repmat(x0,1,2);
        problem.guess.control = [0.1,0];

        %==================================================================
        % Solver options
        %==================================================================

        problem.options.method = 'hermiteSimpson';

        problem.options.nlpOpt = optimset( ...
            'Display','off', ...
            'MaxFunEvals',1e5);

        %==================================================================
        % Solve
        %==================================================================

        try
            soln = optimTraj(problem);

        catch ME

            warning(['Optimization failed for rise = ', ...
                     num2str(rise), ...
                     ' and decay = ', ...
                     num2str(decay)]);

            heatmap_data(i,j) = NaN;
            nan_count = nan_count + 1;

            close(fig);
            continue;
        end

        %==================================================================
        % Interpolate solution
        %==================================================================

        t_sol = linspace( ...
            soln.grid.time(1), ...
            soln.grid.time(end), ...
            150);

        z = soln.interp.state(t_sol);
        u = soln.interp.control(t_sol);

        plot(ax,t_sol,I_opsin(z), ...
            'LineWidth',2.0, ...
            'DisplayName','Actual');

        yyaxis(ax,'right');

        plot(ax,t_sol,10*u, ...
            'LineWidth',2.0, ...
            'HandleVisibility','off');

        ylabel(ax,'470 nm LED voltage (V)');

        %==================================================================
        % Check solution with ODE
        %==================================================================

        try

            [t_ode,y] = ode45( ...
                @opsinmodel_ODE, ...
                [0,duration], ...
                x0, ...
                [], ...
                varpars, ...
                @(s) interp1(t_sol,u,s));

            I = I_opsin(y');

        catch

            warning('ODE simulation failed.');

            heatmap_data(i,j) = NaN;
            nan_count = nan_count + 1;

            close(fig);
            continue;
        end

        %==================================================================
        % Peak current
        %==================================================================

        if abs(rise-decay) < 1e-12

            % Peak of alpha function
            t_peak_target = 1/rise;

        else

            t_peak_target = ...
                log(decay/rise)/(decay-rise);

        end

        peak_target_I = syn_signal(t_peak_target);

        % Fine time grid around the beginning
        t_init = linspace(0,10,10000);
        z_init = soln.interp.state(t_init);

        I_init = I_opsin(z_init);

        [peak_actual_I,ind] = min(I_init);

        %==================================================================
        % Error calculations
        %==================================================================

        peak_diff = peak_actual_I - peak_target_I;

        peak_percent_error = ...
            -100 * peak_diff / peak_target_I;

        t_peak_actual = t_init(ind);

        latency_diff = ...
            t_peak_actual - t_peak_target;

        latency_percent_error = ...
            100 * latency_diff / t_peak_target;

        %==================================================================
        % Check for NaNs
        %==================================================================

        if any(isnan([ ...
                t_peak_target, ...
                peak_target_I, ...
                peak_actual_I, ...
                peak_diff, ...
                peak_percent_error, ...
                latency_diff, ...
                latency_percent_error]))

            heatmap_data(i,j) = NaN;
            nan_count = nan_count + 1;

            fprintf('RESULT = NaN\n');

        else

            % Store peak percentage error
            heatmap_data(i,j) = peak_percent_error;

            successful_count = successful_count + 1;

            fprintf('Peak difference = %.2f pA\n',peak_diff);
            fprintf('Peak percentage error = %.2f%%\n', ...
                peak_percent_error);
            fprintf('Latency difference = %.2f ms\n',latency_diff);
            fprintf('Latency percentage error = %.2f%%\n', ...
                latency_percent_error);

        end

        %==================================================================
        % Save individual figure
        %==================================================================

        fig_filename = fullfile( ...
            save_folder, ...
            sprintf('figure_%d_%d.png',i,j));

        saveas(fig,fig_filename);

        close(fig);

    end
end

%==========================================================================
% NaN statistics
%==========================================================================

fprintf('\n\n');
fprintf('====================================================\n');
fprintf('                 FINAL RESULTS\n');
fprintf('====================================================\n');

fprintf('Total combinations       : %d\n',total_count);
fprintf('Successful combinations  : %d\n',successful_count);
fprintf('NaN combinations         : %d\n',nan_count);
fprintf('Valid percentage         : %.2f%%\n', ...
    100*successful_count/total_count);
fprintf('NaN percentage           : %.2f%%\n', ...
    100*nan_count/total_count);

%==========================================================================
% Save numerical heatmap data
%==========================================================================

writematrix(heatmap_data, ...
    fullfile(save_folder,'peak_percentage_error_heatmap_data.csv'));

%==========================================================================
% Common heatmap settings
%==========================================================================

white_color = [1 1 1];

%==========================================================================
% HEATMAP
% Styled identically to the requested image
%==========================================================================

fig2 = figure('Color','w', 'Position', [100, 100, 1000, 600]);

h2 = heatmap( ...
    rise_values, ...
    decay_values, ...
    heatmap_data, ...
    'Colormap',parula, ...
    'MissingDataColor',white_color);

h2.ColorLimits = [0 60];

h2.MissingDataLabel = 'NaN';

h2.GridVisible = 'off';

% Hide values inside cells
h2.CellLabelColor = 'none';

h2.XLabel = 'Rise';
h2.YLabel = 'Decay';
h2.Title = 'Peak Percentage Error Heatmap';
h2.FontSize = 14; 

% Save complete heatmap data as CSV
heatmap_table = array2table( ...
    heatmap_data, ...
    'VariableNames', ...
    matlab.lang.makeValidName( ...
        compose('Rise_%g', rise_values)), ...
    'RowNames', ...
    compose('Decay_%g', decay_values));

% Add decay as the first column
heatmap_table = addvars( ...
    heatmap_table, ...
    decay_values(:), ...
    'Before', 1, ...
    'NewVariableNames', 'Decay');

% Save CSV
heatmap_csv = fullfile( ...
    save_folder, ...
    'Peak_Percentage_Error_Heatmap_Data.csv');

writetable( ...
    heatmap_table, ...
    heatmap_csv, ...
    'WriteRowNames', false);

fprintf('\nHeatmap data saved to:\n%s\n', heatmap_csv);
