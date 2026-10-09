
    clc;
    clear;
    close all;
    
    % =========================================================
    % READ DATA
    % =========================================================
    % Reading .dat file
    filename = 'state_history_first_test.dat';
    M = readmatrix(filename);
    t = M(:,1);
    Model_V = M(:,2);
    n = M(:,3);
    f = M(:,4);
    Ca = M(:,5);
    coupling = M(:,6);
    coupling_gradient = M(:,7);
    Data_V = M(:,8);
    
    % Adjust time scaling (5x slower oscillations)
    t = t * 5;
    
    % =========================================================
    % LOAD I_app
    % =========================================================
    filename = 'Iapp_data_new.txt';
    M2 = readmatrix(filename);
    I_app = M2(:,1);
    
    % =========================================================
    % ENSURE SAME LENGTH
    % =========================================================
    min_len = min([length(t), length(I_app)]);
    t_plot       = t(1:min_len);
    Model_V_plot = Model_V(1:min_len);
    Data_V_plot  = Data_V(1:min_len);
    I_app_plot   = I_app(1:min_len);
    
    % Convert time to ms
    t_ms = t_plot / 100;
    
    % =========================================================
    % CUT SIGNAL AT 1736.55 ms (using user's 1712.9)
    % =========================================================
    cut_time = 1712.9;   % ms
    
    % Keep data only before cutoff
    cut_idx = find(t_ms <= cut_time);
    t_cut       = t_ms(cut_idx);
    Model_V_cut = Model_V_plot(cut_idx);
    Data_V_cut  = Data_V_plot(cut_idx);
    I_app_cut   = I_app_plot(cut_idx);
    
    % =========================================================
    % COMPUTE RMSE
    % =========================================================
    % Calculate RMSE between the Cut Model and Data Voltages
    squared_errors = (Model_V_cut - Data_V_cut).^2;
    rmse_val = sqrt(mean(squared_errors));
    
    % Print the result to the MATLAB Command Window
    fprintf('\n==================================\n');
    fprintf('Calculated RMSE: %.4f mV\n', rmse_val);
    fprintf('==================================\n\n');
    
    % =========================================================
    % STRETCH TIME AXIS TO 2000 ms
    % =========================================================
    final_time = 2000;  % desired end time
    t_stretch = linspace(0, final_time, length(t_cut));
    
    % =========================================================
    % PLOTTING
    % =========================================================
    fig = figure('Color', 'w');
    
    % ---------------------------------------------------------
    % LEFT AXIS : VOLTAGE
    % ---------------------------------------------------------
    yyaxis left
    p1 = plot(t_stretch, Model_V_cut, ...
        'r', ...
        'LineWidth', 1.5);
    hold on;
    p2 = plot(t_stretch, Data_V_cut, ...
        'k--', ...
        'LineWidth', 1.7);
    ylim([-100 20]);
    ax = gca;
    ax.XColor = 'k';
    ax.YColor = 'k';
    ylabel('Voltage (mV)', ...
        'FontSize', 35, ...
        'Color', 'k');
    % xlabel('Time (ms)', ...
    %     'FontSize', 40, ...
    %     'Color', 'k');
    
    % ---------------------------------------------------------
    % LEGEND
    % ---------------------------------------------------------
    lgd = legend([p2 p1], ...
        {'Data Voltage', 'Model Voltage'}, ...
        'Location', 'northeast', ...
        'Orientation', 'horizontal', ...
        'Box', 'off');
    lgd.FontSize = 35;
    
    % ---------------------------------------------------------
    % RIGHT AXIS : I_app
    % ---------------------------------------------------------
    yyaxis right
    plot(t_stretch, 8 * I_app_cut, ...
        'b', ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    ylabel('I_{app} (pA)', ...
        'FontSize', 35, ...
        'Color', 'b');
    ylim([-1 10]);
    ax.YColor = 'b';
    
    % ---------------------------------------------------------
    % TITLE + FORMATTING
    % ---------------------------------------------------------
    % Added the dynamic RMSE value into the title
    
    
    xlim([0 2000]);
    set(gca, ...
        'FontSize', 35, ...
        'LineWidth', 1.2);
    % ax.XTick = []; % REMOVES X-TICKS
    grid off;
    box off;
    drawnow;
