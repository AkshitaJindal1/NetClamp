% 11.10.2022
% Akshita Jindal
clear;
close all;
clc;
[modelPars,varPars,sysPars,fun] = setDefaultParameters_LM();

%Reading .dat file
filename = 'I_app_protocol.dat';
M = readmatrix ("I_app_protocol.dat");
t = 5*M(:,1);  %In Alain's data, the frequency is 50KHz (for neurons) which is too high for our cells (10KHz). So, to make the oscillations 5times slower, multiply by 5 instead of dividing
I = M(:,2);

%y0 are the initial conditions.
y0 = [-60;0.01;0.01;0.5]; %Voltage, n-activation of potassium current, f-activation of Bk Channels, concentartion of Ca ions

% Anonymous functions - no step current i.e. zero amplitude
% stepCurrent(t, pretime, duration, amplitude)
% Define anonymous function using @
% @t because I_app is time dependent
% Syntax - [t,y] = ode45(odefun, tspan, y0, options)
odeopts = odeset('AbsTol', 1e-8, 'RelTol', 1e-6); %makes the plot run faster/slower
I_app = @(s) 0.0; %Zero amplitude
[~,y] = ode45(fun, [0,10000], y0, odeopts, varPars, I_app);
y0 = y(end,:);
I_app = @(s) 10*interp1(t,I,s,'nearest');
s = linspace(t(1), t(end),10000);
[t,y] = ode45(fun, [0,10000], y0, odeopts, varPars, I_app); %time span is now 20000 because Alain's data has 10,000 rows for t and for our cells, it is 5 times slower. So, 10000*5 = 50000. 20000 is within the range

% Anonymous functions - 10 pA (picoamperes) amplitude*interp1
% I_app = @(s) 10*interp1(t,I,s);
% s = linspace(t(1), t(end),50000);
%[t,y] = ode45(fun, [0,10000], y(end,:)', odeopts, varPars, I_app);

fig = figure('Position', [100 100 800 400]);
% Left axis (Voltage)
yyaxis left
plot(t/1000, y(:,1), 'k', 'LineWidth', 1);
xlabel('Time (s)', 'FontSize', 14, 'Color', 'k');
ylabel('Voltage (mV)', 'FontSize', 14, 'Color', 'k');
xlim([0 10]);
ylim([-160 20]);
ax = gca;
ax.XColor = 'k';   % X-axis black
ax.YColor = 'k';   % Left Y-axis black

% Right axis (Applied current)
yyaxis right
plot(s/1000, I_app(s), 'b', 'LineWidth', 1);
ylabel('I_{app} (pA)', 'FontSize', 14, 'Color', 'b');
ylim([-2 10]);
ax.YColor = 'b';   % Right Y-axis blue
drawnow;

% ========================================================================
% NEW CODE TO SAVE VOLTAGE DATA AS A SINGLE COLUMN TXT
% ========================================================================
% Extract just the voltage data (first column of y)
Voltage = y(:,1);

% Save the voltage array to a plain text file
writematrix(Voltage, 'Generated_Voltage.txt');
disp('Voltage data successfully saved to Generated_Voltage.txt');
