
% ==============================
% CELL COMPARISON: Ionic Currents
% ==============================
% Load the synthetic data (File 1) and experimental data (File 2)

file1 = 'state_history_first_test.dat';
file2 = 'state_history_exp_cell.dat';

data1 = readtable(file1);
data2 = readtable(file2);

%% ==============================
% Preprocessing
% ==============================

% Extracting columns:
% Step/time, Membrane voltage (V), potassium gate (n),
% BK gate (f), Intracellular calcium (Ca)

preprocess = @(data) deal( ...
    data{:,1}, data{:,2}, data{:,3}, data{:,4}, data{:,5});

[steps1,V1,n1,f1,Ca1] = preprocess(data1);
[steps2,V2,n2,f2,Ca2] = preprocess(data2);

% Removing any NaNs
clean = @(s,V,n,f,Ca) deal( ...
    s(~isnan(s)&~isnan(V)&~isnan(n)&~isnan(f)&~isnan(Ca)), ...
    V(~isnan(s)&~isnan(V)&~isnan(n)&~isnan(f)&~isnan(Ca)), ...
    n(~isnan(s)&~isnan(V)&~isnan(n)&~isnan(f)&~isnan(Ca)), ...
    f(~isnan(s)&~isnan(V)&~isnan(n)&~isnan(f)&~isnan(Ca)), ...
    Ca(~isnan(s)&~isnan(V)&~isnan(n)&~isnan(f)&~isnan(Ca)) ...
);

[steps1,V1,n1,f1,Ca1] = clean(steps1,V1,n1,f1,Ca1);
[steps2,V2,n2,f2,Ca2] = clean(steps2,V2,n2,f2,Ca2);

%% ==============================
% Model Parameters
% ==============================

% Computing currents using default parameters

gCa = 5;
gk  = 20;
gSK = 10;
gBK = 0.998;

VCa = 120;
Vk  = -80;

Vm = -20;
Sm = 6;

Ks = 1;

compute_currents = @(V,n,f,Ca) deal( ...
    gCa .* (1 ./ (1 + exp((Vm - V) ./ Sm))) .* (V - VCa), ...
    gk  .* n .* (V - Vk), ...
    gSK .* ((Ca.^2) ./ (Ca.^2 + Ks^2)) .* (V - Vk), ...
    gBK .* f .* (V - Vk) ...
);

[ICa1,IK1,ISK1,IBK1] = compute_currents(V1,n1,f1,Ca1);
[ICa2,IK2,ISK2,IBK2] = compute_currents(V2,n2,f2,Ca2);

%% ==============================
% Truncate to match the timescales
% ==============================

L1 = min(32000,length(V1));
L2 = min(32000,length(V2));

ICa1 = ICa1(1:L1);
IK1  = IK1(1:L1);
ISK1 = ISK1(1:L1);
IBK1 = IBK1(1:L1);
Ca1  = Ca1(1:L1);

ICa2 = ICa2(1:L2);
IK2  = IK2(1:L2);
ISK2 = ISK2(1:L2);
IBK2 = IBK2(1:L2);
Ca2  = Ca2(1:L2);

%% ==============================
% WINDOW ANALYSIS
% ==============================

% Window length and shift
W = 500;
step = 250;

% Windowed metrics
compute_window_metric = @(x) ...
    arrayfun(@(k) std(x(k:k+W-1)), ...
    1:step:length(x)-W);

rms_metric = @(x) ...
    arrayfun(@(k) rms(x(k:k+W-1)), ...
    1:step:length(x)-W);

ptp_metric = @(x) ...
    arrayfun(@(k) peak2peak(x(k:k+W-1)), ...
    1:step:length(x)-W);

rough_metric = @(x) ...
    arrayfun(@(k) std(diff(x(k:k+W-1))), ...
    1:step:length(x)-W);

%% ==============================
% CURRENT VARIABILITY
% ==============================

ICav1 = compute_window_metric(ICa1);
ICav2 = compute_window_metric(ICa2);

IKv1 = compute_window_metric(IK1);
IKv2 = compute_window_metric(IK2);

ISKv1 = compute_window_metric(ISK1);
ISKv2 = compute_window_metric(ISK2);

IBKv1 = compute_window_metric(IBK1);
IBKv2 = compute_window_metric(IBK2);

%% ==============================
% TOTAL K+ VARIABILITY
% ==============================

Ktot_v1 = compute_window_metric(IK1 + ISK1 + IBK1);
Ktot_v2 = compute_window_metric(IK2 + ISK2 + IBK2);

%% ==============================
% K+ RMS ACTIVITY
% ==============================

IKrms1 = rms_metric(IK1 + ISK1 + IBK1);
IKrms2 = rms_metric(IK2 + ISK2 + IBK2);

%% ==============================
% PEAK-TO-PEAK
% ==============================

Kptp1 = ptp_metric(IK1 + ISK1 + IBK1);
Kptp2 = ptp_metric(IK2 + ISK2 + IBK2);

%% ==============================
% TEMPORAL ROUGHNESS
% ==============================

rough1 = rough_metric(IK1 + ISK1 + IBK1);
rough2 = rough_metric(IK2 + ISK2 + IBK2);

Ca_rough1 = rough_metric(ICa1);
Ca_rough2 = rough_metric(ICa2);

%% ==============================
% TABLE OUTPUT
% ==============================

names = {
    'Ca current variability'
    'Total K variability'
    'IK variability'
    'ISK variability'
    'IBK variability'
    'K RMS activity'
    'K peak-to-peak'
    'K temporal roughness'
    'Ca temporal roughness'
};

Cell1 = [
    mean(ICav1)
    mean(Ktot_v1)
    mean(IKv1)
    mean(ISKv1)
    mean(IBKv1)
    mean(IKrms1)
    mean(Kptp1)
    mean(rough1)
    mean(Ca_rough1)
];

Cell2 = [
    mean(ICav2)
    mean(Ktot_v2)
    mean(IKv2)
    mean(ISKv2)
    mean(IBKv2)
    mean(IKrms2)
    mean(Kptp2)
    mean(rough2)
    mean(Ca_rough2)
];

Delta = Cell2 - Cell1;

T = table(names,Cell1,Cell2,Delta);

disp(T);

%% ==============================
% PLOT SETTINGS
% ==============================

% Same figure size for every figure
figPos = [100 100 1200 750];

% Same font size everywhere
fontSize = 45;

% Model colours
c1 = [0.85 0.3 0.5];    % Synthetic
c2 = [0.2 0.5 0.9];     % Experimental


%% ============================================================
% FIGURE 1: WINDOWED CURRENT VARIABILITY
% ============================================================

figure('Color','w','Position',figPos);

hold on;

% --------------------------------
% ICa
% --------------------------------

scatter(ones(size(ICav1))*1,ICav1,100,c1,...
    'filled','MarkerFaceAlpha',0.45);

scatter(ones(size(ICav2))*2,ICav2,100,c2,...
    'filled','MarkerFaceAlpha',0.45);

% --------------------------------
% IK
% --------------------------------

scatter(ones(size(IKv1))*4,IKv1,100,c1,...
    'filled','MarkerFaceAlpha',0.45);

scatter(ones(size(IKv2))*5,IKv2,100,c2,...
    'filled','MarkerFaceAlpha',0.45);

% --------------------------------
% ISK
% --------------------------------

scatter(ones(size(ISKv1))*7,ISKv1,100,c1,...
    'filled','MarkerFaceAlpha',0.45);

scatter(ones(size(ISKv2))*8,ISKv2,100,c2,...
    'filled','MarkerFaceAlpha',0.45);

% --------------------------------
% IBK
% --------------------------------

scatter(ones(size(IBKv1))*10,IBKv1,100,c1,...
    'filled','MarkerFaceAlpha',0.45);

scatter(ones(size(IBKv2))*11,IBKv2,100,c2,...
    'filled','MarkerFaceAlpha',0.45);

% --------------------------------
% Mean markers
% --------------------------------

plot([1 2],[mean(ICav1) mean(ICav2)],...
    'k_','LineWidth',7,'MarkerSize',35);

plot([4 5],[mean(IKv1) mean(IKv2)],...
    'k_','LineWidth',7,'MarkerSize',35);

plot([7 8],[mean(ISKv1) mean(ISKv2)],...
    'k_','LineWidth',7,'MarkerSize',35);

plot([10 11],[mean(IBKv1) mean(IBKv2)],...
    'k_','LineWidth',7,'MarkerSize',35);

% --------------------------------
% Axes
% --------------------------------

xlim([0 12]);

xticks([1.5 4.5 7.5 10.5]);

xticklabels({...
    '$I_{\mathrm{Ca}}$',...
    '$I_{\mathrm{K}}$',...
    '$I_{\mathrm{SK}}$',...
    '$I_{\mathrm{BK}}$'});

set(gca,'TickLabelInterpreter','latex');

xlabel('Reconstructed Ionic Current',...
    'FontSize',fontSize);

ylabel('Current variability, $\sigma_{\alpha,j}$ (pA)',...
    'Interpreter','latex',...
    'FontSize',fontSize);

% --------------------------------
% Legend
% --------------------------------

% Create dummy scatter objects for legend
hLeg1 = scatter(nan,nan,100,c1,...
    'filled','MarkerFaceAlpha',0.45);

hLeg2 = scatter(nan,nan,100,c2,...
    'filled','MarkerFaceAlpha',0.45);

% Create legend using the dummy markers
lgd = legend([hLeg1 hLeg2],...
    {'Synthetic model','Experimental model'},...
    'Orientation','horizontal',...
    'Location','northwest',...
    'Box','off');

lgd.FontSize = 38;

% --------------------------------
% Formatting
% --------------------------------

ax = gca;

ax.FontSize = fontSize;
ax.LineWidth = 1.5;

grid off;
box off;

hold off;


%% ============================================================
% FIGURE 2: WINDOWED POTASSIUM RMS ACTIVITY
% ============================================================

figure('Color','w','Position',figPos);

hold on;

% --------------------------------
% Individual analysis windows
% --------------------------------

s1 = scatter(ones(size(IKrms1))*1,...
    IKrms1,...
    100,...
    c1,...
    'filled',...
    'MarkerFaceAlpha',0.45);

s2 = scatter(ones(size(IKrms2))*2,...
    IKrms2,...
    100,...
    c2,...
    'filled',...
    'MarkerFaceAlpha',0.45);

% --------------------------------
% Mean markers
% --------------------------------

plot(1,mean(IKrms1),...
    'k_','MarkerSize',35,'LineWidth',7);

plot(2,mean(IKrms2),...
    'k_','MarkerSize',35,'LineWidth',7);

% --------------------------------
% Axes
% --------------------------------

xlim([0.5 2.5]);

xticks([1 2]);

xticklabels({...
    'Synthetic model',...
    'Experimental model'});

ylabel('Potassium activity, $\mathrm{RMS}_{K,j}$ (pA)',...
    'Interpreter','latex',...
    'FontSize',fontSize);

% --------------------------------
% Formatting
% --------------------------------

ax = gca;

ax.FontSize = fontSize;
ax.LineWidth = 1.5;

grid off;
box off;

hold off;


%% ============================================================
% FIGURE 3: MEAN CHANNEL CONTRIBUTION
% ============================================================

epsilon = 1e-12;

% --------------------------------
% Total instantaneous current
% --------------------------------

Itotal1 = ICa1 + IK1 + ISK1 + IBK1 + epsilon;
Itotal2 = ICa2 + IK2 + ISK2 + IBK2 + epsilon;

% --------------------------------
% Instantaneous normalised contributions
% --------------------------------

pCa1 = ICa1 ./ Itotal1;
pK1  = IK1  ./ Itotal1;
pSK1 = ISK1 ./ Itotal1;
pBK1 = IBK1 ./ Itotal1;

pCa2 = ICa2 ./ Itotal2;
pK2  = IK2  ./ Itotal2;
pSK2 = ISK2 ./ Itotal2;
pBK2 = IBK2 ./ Itotal2;

% --------------------------------
% Time-averaged contributions
% --------------------------------

pbar1 = [
    mean(pCa1)
    mean(pK1)
    mean(pSK1)
    mean(pBK1)
];

pbar2 = [
    mean(pCa2)
    mean(pK2)
    mean(pSK2)
    mean(pBK2)
];

% --------------------------------
% Shannon entropy
% --------------------------------

H1 = -sum(pbar1 .* log(pbar1 + epsilon));
H2 = -sum(pbar2 .* log(pbar2 + epsilon));


%% ============================================================
% PERMUTATION TEST FOR ENTROPY
% ============================================================

Nperm = 5000;

obs_diff = abs(H2-H1);

% Pool instantaneous contributions

poolCa = [pCa1; pCa2];
poolK  = [pK1;  pK2];
poolSK = [pSK1; pSK2];
poolBK = [pBK1; pBK2];

n1 = length(pCa1);
n2 = length(pCa2);

perm_stat = zeros(Nperm,1);

for i = 1:Nperm

    idx = randperm(n1+n2);

    % Random split

    CaA = poolCa(idx(1:n1));
    CaB = poolCa(idx(n1+1:end));

    KA = poolK(idx(1:n1));
    KB = poolK(idx(n1+1:end));

    SKA = poolSK(idx(1:n1));
    SKB = poolSK(idx(n1+1:end));

    BKA = poolBK(idx(1:n1));
    BKB = poolBK(idx(n1+1:end));

    % Mean contribution vectors

    pA = [
        mean(CaA)
        mean(KA)
        mean(SKA)
        mean(BKA)
    ];

    pB = [
        mean(CaB)
        mean(KB)
        mean(SKB)
        mean(BKB)
    ];

    % Renormalise

    pA = pA ./ sum(pA + epsilon);
    pB = pB ./ sum(pB + epsilon);

    % Entropy

    H_A = -sum(pA .* log(pA + epsilon));
    H_B = -sum(pB .* log(pB + epsilon));

    perm_stat(i) = abs(H_B-H_A);

end

% --------------------------------
% Permutation p-value
% --------------------------------

p_value = ...
    (sum(perm_stat >= obs_diff)+1)/(Nperm+1);

fprintf('\n');
fprintf('Observed Delta H = %.4f\n',obs_diff);
fprintf('Permutation p-value = %.5f\n',p_value);


%% ============================================================
% PLOT MEAN CHANNEL CONTRIBUTION
% ============================================================

figure('Color','w','Position',figPos);

hold on;

% --------------------------------
% Bar plot
% --------------------------------

b = bar([pbar1 pbar2]);

b(1).FaceColor = c1;
b(2).FaceColor = c2;

b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';

% --------------------------------
% X-axis
% --------------------------------

set(gca,'XTick',1:4);

ax = gca;

ax.TickLabelInterpreter = 'latex';

xticklabels({...
    '$I_{\mathrm{Ca}}$',...
    '$I_{\mathrm{K}}$',...
    '$I_{\mathrm{SK}}$',...
    '$I_{\mathrm{BK}}$'});

xlim([0.5 4.5]);

% --------------------------------
% Axis labels
% --------------------------------

xlabel('Reconstructed Ionic Currents',...
    'FontSize',fontSize);

ylabel('Mean channel contribution, $\bar{p}_{\alpha}$',...
    'Interpreter','latex',...
    'FontSize',fontSize);

% --------------------------------
% Legend
% --------------------------------

lgd = legend({...
    sprintf('H_{Synthetic} = %.3f',H1),...
    sprintf('H_{Experimental} = %.3f',H2)},...
    'Orientation','vertical',...
    'Location','northwest',...
    'Box','off');

lgd.FontSize = fontSize;


% --------------------------------
% Formatting
% --------------------------------

ax.FontSize = fontSize;
ax.LineWidth = 1.5;

grid off;
box off;

hold off;
