clc;
clear;
close all;

%% ============================================================
% LOAD DATA
%% ============================================================

folder = 'C:\Users\aj609\OneDrive - University of Exeter\Year 4_Thesis_Material\Voltage Imaging\Voltron\Spike_detection';
file = 'Kinetix_data.csv';

data = readtable(fullfile(folder,file), ...
    'VariableNamingRule','preserve');

Frame = data{:,1};
Time  = data{:,2};

names = data.Properties.VariableNames;


%% ============================================================
% AUTO-DETECT CELLS + BACKGROUND
%% ============================================================

cellIdx = contains(names,'Cell','IgnoreCase',true);
bgIdx   = contains(names,'Background','IgnoreCase',true);

signals_raw = data{:,cellIdx};
signalNames = names(cellIdx);

nCells = sum(cellIdx);
nSamples = height(data);

if any(bgIdx)
    Background = data{:,bgIdx};
else
    Background = zeros(nSamples,1);
end

fprintf('Detected %d cells\n',nCells);

%% ============================================================
% BACKGROUND SUBTRACTION
%% ============================================================

signals = signals_raw - Background;
signals(signals <= 0) = eps;

%% ============================================================
% ΔF/F COMPUTATION (YOUR PIPELINE)
%% ============================================================

dFF = zeros(nSamples,nCells);

for c = 1:nCells
    F0 = median(signals(:,c));
    dFF(:,c) = (signals(:,c) - F0) / F0;
end

dFF = dFF * 100;

%% ============================================================
% PARAMETERS (PYTHON-LIKE EMA MODEL)
%% ============================================================

alpha_mean = 0.08;
alpha_var  = 0.02;

thresh = 2;
min_spike_gap = 5;    % 62.5 ms (spike interval) for the given camera speed. 1000/78 (fps) = 12.5 ms

meanF = zeros(nSamples,nCells);
varF  = zeros(nSamples,nCells);

spike_idx = cell(nCells,1);
spike_amp = cell(nCells,1);

overallSTD = zeros(nCells,1);
firstQuarterSTD = zeros(nCells,1);

%% ============================================================
% MAIN LOOP (PER CELL)
%% ============================================================

for c = 1:nCells

    sig = dFF(:,c);

    firstQuarter = floor(nSamples/4);

    firstQuarterSTD(c) = std(sig(1:firstQuarter));
    overallSTD(c) = std(sig);

    meanF(1,c) = mean(sig(1:min(200,nSamples)));
    varF(1,c)  = var(sig(1:min(200,nSamples)));

    last_spike = -Inf;

    for i = 1:nSamples-1

        % EMA mean
        meanF(i+1,c) = alpha_mean*sig(i+1) + ...
            (1-alpha_mean)*meanF(i,c);

        % EMA variance
        varF(i+1,c) = alpha_var*(sig(i+1)-meanF(i+1,c))^2 + ...
            (1-alpha_var)*varF(i,c);

        % refractory
        if (i-last_spike) < min_spike_gap
            continue
        end

        sigma = sqrt(varF(i+1,c));

        if (sig(i+1) - meanF(i+1,c)) < -thresh*sigma

            spike_idx{c}(end+1) = i+1;
            spike_amp{c}(end+1) = sig(i+1) - meanF(i+1,c);

            last_spike = i+1;

        end

    end

end

%% ============================================================
% PRINT STATISTICS
%% ============================================================

fprintf('\n================ SPIKE SUMMARY ================\n');

for c = 1:nCells
    fprintf('\n%s\n',signalNames{c});
    fprintf('First-quarter STD : %.4f\n',firstQuarterSTD(c));
    fprintf('Overall STD       : %.4f\n',overallSTD(c));
    fprintf('Spike count       : %d\n',numel(spike_idx{c}));
end

%% ============================================================
% COLOURS
%% ============================================================

blue = [0.53 0.81 0.98];
pink = [1.00 0.41 0.71];

%% ============================================================
% PLOT 1 — CELL 1
%% ============================================================

figure('Color','w');

subplot(2,1,1); hold on;

plot(Time, dFF(:,1), 'Color', blue, 'LineWidth', 1.5);

if ~isempty(spike_idx{1})
    scatter(Time(spike_idx{1}), dFF(spike_idx{1},1), 40, 'k', 'filled');
end


plot(Time, meanF(:,1)-thresh*sqrt(varF(:,1)),'--k','LineWidth',1.2);

ylabel('\DeltaF/F (%)','FontSize',20 , 'FontWeight','bold');

title('Cell 1','FontSize',20,'FontWeight','bold');

xlim([0 30]);
ylim([-8 8]);

ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
yticks([-8 -4 0  4  8]);

legend({'Signal','Spikes','Threshold'}, ...
       'Location', 'northeast',... 
       'NumColumns', 3,...
       'Box', 'off');

box off;

%% ============================================================
% PLOT 2 — CELL 2
%% ============================================================

subplot(2,1,2); 

hold on;

plot(Time, dFF(:,2), 'Color', pink, 'LineWidth', 1.5);

if ~isempty(spike_idx{2})
    scatter(Time(spike_idx{2}), dFF(spike_idx{2},2), 40, 'r', 'filled');
end

plot(Time, meanF(:,2)-thresh*sqrt(varF(:,2)),'--r','LineWidth',1.2);

title('Cell 2','FontSize',20,'FontWeight','bold');

xlabel('Time (s)','FontSize',20,'FontWeight','bold');
ylabel('\DeltaF/F (%)','FontSize',20 , 'FontWeight','bold');

xlim([0 30]);
ylim([-8 8]);

ax = gca;
ax.FontSize = 20;
ax.LineWidth = 1.2;
ax.TickDir = 'out';
yticks([-8 -4 0  4  8]);

legend({'Signal','Spikes','Threshold'}, ...
       'Location', 'northeast',... 
       'NumColumns', 3,...
       'Box', 'off');

box off;
