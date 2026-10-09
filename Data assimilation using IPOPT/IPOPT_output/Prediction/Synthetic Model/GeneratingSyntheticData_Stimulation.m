% % 11.10.2022
% % Akshita Jindal
% %
% % COMPLETE UPDATED VERSION
% %
% % MODIFIED:
% % - Voltage data now loaded from: voltage_data_new.txt
% % - voltage_data_new.txt contains ONLY ONE COLUMN (voltage values)
% % - Protocol1.dat is now used for stimulation
% % - Experimental time vector generated automatically from protocol timing
% %
% % Features:
% % - Uses one common time vector
% % - Experimental data + model + current aligned
% % - Automatic amplitude matching
% % - Adjustable stimulation scaling
% % - Prevents plotting mismatch issues
% 
% clear;
% close all;
% clc;
% 
% %% ================= LOAD MODEL PARAMETERS =================
% 
% [modelPars,varPars,sysPars,fun] = setDefaultParameters_LM();
% 
% %% ================= READ PROTOCOL DATA =================
% 
% filename = 'Protocol1.dat';
% 
% M = readmatrix(filename);
% 
% % Time from protocol
% % Alain's data frequency scaling
% t_protocol = 5 * M(:,1);
% 
% % Current protocol
% I_protocol = M(:,2);
% 
% %% ================= LOAD VOLTAGE DATA =================
% % voltage_data_new.txt contains ONLY voltage values
% 
% V_data = readmatrix('V.txt');
% 
% % Ensure column vector
% V_data = V_data(:);
% 
% %% ================= CREATE EXPERIMENTAL TIME VECTOR =================
% % Since voltage file has no time column,
% % generate time vector using protocol timing
% 
% N_data = length(V_data);
% 
% t_data = linspace(t_protocol(1), ...
%                   t_protocol(end), ...
%                   N_data)';
% 
% %% ================= INITIAL CONDITIONS =================
% 
% % [Voltage, n-gate, f-gate, Ca concentration]
% y0 = [-4.190985e+01;
%       1.347435e-01;
%       3.752549e-01;
%       4.310361e-01];
% 
% %% ================= ODE OPTIONS =================
% 
% odeopts = odeset('AbsTol', 1e-8, ...
%                  'RelTol', 1e-6);
% 
% %% ================= PRE-EQUILIBRATION =================
% % Run with zero current first
% 
% I_app = @(s) 0.0;
% 
% [~, y_pre] = ode45(fun, ...
%                    [0 10000], ...
%                    y0, ...
%                    odeopts, ...
%                    varPars, ...
%                    I_app);
% 
% % Updated initial condition
% y0 = y_pre(end,:);
% 
% %% ================= CREATE COMMON TIME VECTOR =================
% 
% % Use overlapping time range
% t_end = min([t_protocol(end), t_data(end)]);
% 
% % Number of interpolation points
% Npts = 10000;
% 
% % Common time vector
% t_common = linspace(0, t_end, Npts);
% 
% %% ================= INTERPOLATE DATA =================
% 
% % Current interpolated onto common time
% I_common = interp1(t_protocol, ...
%                    I_protocol, ...
%                    t_common, ...
%                    'pchip', ...
%                    'extrap');
% 
% % Experimental voltage interpolated onto common time
% V_common = interp1(t_data, ...
%                    V_data, ...
%                    t_common, ...
%                    'linear', ...
%                    'extrap');
% 
% %% ================= DEFINE APPLIED CURRENT =================
% % IMPORTANT:
% % Increase stim_scale to increase model oscillation amplitude
% 
% stim_scale = 18;
% 
% I_app = @(s) stim_scale * interp1(t_common, ...
%                                   I_common, ...
%                                   s, ...
%                                   'nearest', ...
%                                   'extrap');
% 
% %% ================= SOLVE MODEL =================
% 
% [t_model, y] = ode45(fun, ...
%                      t_common, ...
%                      y0, ...
%                      odeopts, ...
%                      varPars, ...
%                      I_app);
% 
% %% ================= EXTRACT MODEL VOLTAGE =================
% 
% V_model_raw = y(:,1);
% 
% %% ================= AUTOMATIC AMPLITUDE MATCHING =================
% % This rescales model amplitude to match experimental data
% 
% % Remove means
% Vm0 = V_model_raw - mean(V_model_raw);
% Vd0 = V_common - mean(V_common);
% 
% % Compute amplitude scaling factor
% amp_scale = std(Vd0) / std(Vm0);
% 
% % Scale model voltage
% V_model = mean(V_common) + amp_scale * Vm0;
% 
% %% ================= OPTIONAL OFFSET CORRECTION =================
% % Uncomment if baseline mismatch exists
% 
% % offset = mean(V_common) - mean(V_model);
% % V_model = V_model + offset;
% 
% %% ================= PLOT RESULTS =================
% 
% fig = figure('Position', [100 100 1000 450]);
% 
% %% ----- LEFT AXIS: VOLTAGE -----
% 
% yyaxis left
% 
% % Synthetic voltage: black dashed
% p1 = plot(t_common/1000, ...
%           V_common, ...
%           'Color', 'k', ...
%           'LineStyle', '--', ...
%           'LineWidth', 1.7);
% 
% hold on;
% 
% % Model voltage: red solid
% p2 = plot(t_model/1000, ...
%           V_model, ...
%           'Color', 'r', ...
%           'LineStyle', '-', ...
%           'LineWidth', 1.2);
% 
% % xlabel('Time (s)', ...
% %        'FontSize', 35);
% 
% ylabel('Voltage (mV)', ...
%        'FontSize', 35);
% 
% ylim([-120 20]);
% 
% ax = gca;
% 
% ax.XColor = 'k';
% ax.YColor = 'k';
% 
% %% ----- RIGHT AXIS: APPLIED CURRENT -----
% 
% yyaxis right
% 
% p3 = plot(t_common/1000, ...
%           8 * I_common, ...
%           'b', ...
%           'LineWidth', 1.2);
% 
% ylabel('I_{app} (pA)', ...
%        'FontSize', 35, ...
%        'Color', 'b');
% 
% ylim([-1 10]);
% 
% ax.YColor = 'b';
% 
% %% ================= LEGEND =================
% 
% lgd = legend([p1 p2], ...
%              {'Data Voltage', ...
%               'Predicted Voltage'}, ...
%              'Location', 'northeast', ...
%              'NumColumns', 2, ...
%                 'FontSize', 35);
% 
% lgd.Box = 'off';
% 
% %% ================= FIGURE SETTINGS =================
% 
% set(gca, 'FontSize', 35);
% % ax.XTick = []; % REMOVES X-TICKS
% 
% box off;
% grid off;
% 
% drawnow;
% 
% %% ================= DISPLAY FIT INFORMATION =================
% 
% fprintf('\n====================================\n');
% fprintf('Stimulus scaling factor = %.2f\n', stim_scale);
% fprintf('Amplitude scaling factor = %.2f\n', amp_scale);
% fprintf('====================================\n\n');

% 11.10.2022
% Akshita Jindal
%
% COMPLETE UPDATED VERSION
%
% MODIFIED:
% - Voltage data now loaded from: voltage_data_new.txt
% - voltage_data_new.txt contains ONLY ONE COLUMN (voltage values)
% - Protocol1.dat is now used for stimulation
% - Experimental time vector generated automatically from protocol timing
%
% Features:
% - Uses one common time vector
% - Experimental data + model + current aligned
% - Automatic amplitude matching
% - Adjustable stimulation scaling
% - Prevents plotting mismatch issues

clear;
close all;
clc;

%% ================= LOAD MODEL PARAMETERS =================
[modelPars,varPars,sysPars,fun] = setDefaultParameters_LM();

%% ================= READ PROTOCOL DATA =================
filename = 'Protocol1.dat';
M = readmatrix(filename);

% Time from protocol
% Alain's data frequency scaling
t_protocol = 5 * M(:,1);

% Current protocol
I_protocol = M(:,2);

%% ================= LOAD VOLTAGE DATA =================
% voltage_data_new.txt contains ONLY voltage values
V_data = readmatrix('V.txt');

% Ensure column vector
V_data = V_data(:);

%% ================= CREATE EXPERIMENTAL TIME VECTOR =================
% Since voltage file has no time column,
% generate time vector using protocol timing
N_data = length(V_data);
t_data = linspace(t_protocol(1), ...
                  t_protocol(end), ...
                  N_data)';

%% ================= INITIAL CONDITIONS =================
% [Voltage, n-gate, f-gate, Ca concentration]
y0 = [-4.190985e+01;
      1.347435e-01;
      3.752549e-01;
      4.310361e-01];

%% ================= ODE OPTIONS =================
odeopts = odeset('AbsTol', 1e-8, ...
                 'RelTol', 1e-6);

%% ================= PRE-EQUILIBRATION =================
% Run with zero current first
I_app = @(s) 0.0;
[~, y_pre] = ode45(fun, ...
                   [0 10000], ...
                   y0, ...
                   odeopts, ...
                   varPars, ...
                   I_app);

% Updated initial condition
y0 = y_pre(end,:);

%% ================= CREATE COMMON TIME VECTOR =================
% Use overlapping time range
t_end = min([t_protocol(end), t_data(end)]);

% Number of interpolation points
Npts = 10000;

% Common time vector
t_common = linspace(0, t_end, Npts);

%% ================= INTERPOLATE DATA =================
% Current interpolated onto common time
I_common = interp1(t_protocol, ...
                   I_protocol, ...
                   t_common, ...
                   'pchip', ...
                   'extrap');

% Experimental voltage interpolated onto common time
V_common = interp1(t_data, ...
                   V_data, ...
                   t_common, ...
                   'linear', ...
                   'extrap');

%% ================= DEFINE APPLIED CURRENT =================
% IMPORTANT:
% Increase stim_scale to increase model oscillation amplitude
stim_scale = 18;
I_app = @(s) stim_scale * interp1(t_common, ...
                                  I_common, ...
                                  s, ...
                                  'nearest', ...
                                  'extrap');

%% ================= SOLVE MODEL =================
[t_model, y] = ode45(fun, ...
                     t_common, ...
                     y0, ...
                     odeopts, ...
                     varPars, ...
                     I_app);

%% ================= EXTRACT MODEL VOLTAGE =================
V_model_raw = y(:,1);

%% ================= AUTOMATIC AMPLITUDE MATCHING =================
% This rescales model amplitude to match experimental data
% Remove means
Vm0 = V_model_raw - mean(V_model_raw);
Vd0 = V_common - mean(V_common);

% Compute amplitude scaling factor
amp_scale = std(Vd0) / std(Vm0);

% Scale model voltage
V_model = mean(V_common) + amp_scale * Vm0;

%% ================= OPTIONAL OFFSET CORRECTION =================
% Uncomment if baseline mismatch exists
% offset = mean(V_common) - mean(V_model);
% V_model = V_model + offset;

%% ================= PLOT RESULTS =================
fig = figure('Position', [100 100 1000 450]);

%% ----- LEFT AXIS: VOLTAGE -----
yyaxis left

% Synthetic voltage: black dashed
p1 = plot(t_common/1000, ...
          V_common, ...
          'Color', 'k', ...
          'LineStyle', '--', ...
          'LineWidth', 1.7);
hold on;

% Model voltage: red solid
p2 = plot(t_model/1000, ...
          V_model, ...
          'Color', 'r', ...
          'LineStyle', '-', ...
          'LineWidth', 1.2);

% xlabel('Time (s)', ...
%        'FontSize', 35);
ylabel('Voltage (mV)', ...
       'FontSize', 35);
ylim([-120 20]);
ax = gca;
ax.XColor = 'k';
ax.YColor = 'k';

%% ----- RIGHT AXIS: APPLIED CURRENT -----
yyaxis right
p3 = plot(t_common/1000, ...
          8 * I_common, ...
          'b', ...
          'LineWidth', 1.2);

ylabel('I_{app} (pA)', ...
       'FontSize', 35, ...
       'Color', 'b');
ylim([-1 10]);
ax.YColor = 'b';

%% ================= LEGEND =================
lgd = legend([p1 p2], ...
             {'Data Voltage', ...
              'Predicted Voltage'}, ...
             'Location', 'northeast', ...
             'NumColumns', 2, ...
                'FontSize', 35);
lgd.Box = 'off';

%% ================= FIGURE SETTINGS =================
set(gca, 'FontSize', 35);
% ax.XTick = []; % REMOVES X-TICKS
box off;
grid off;
drawnow;

%% ================= COMPUTE RMSE =================
% Calculate Root Mean Square Error between Data and Predicted Voltage
% Using (:) to ensure column vectors for element-wise subtraction
rmse_val = sqrt(mean((V_common(:) - V_model(:)).^2));

%% ================= DISPLAY FIT INFORMATION =================
fprintf('\n====================================\n');
fprintf('Stimulus scaling factor = %.2f\n', stim_scale);
fprintf('Amplitude scaling factor = %.2f\n', amp_scale);
fprintf('RMSE                    = %.4f mV\n', rmse_val);
fprintf('====================================\n\n');
