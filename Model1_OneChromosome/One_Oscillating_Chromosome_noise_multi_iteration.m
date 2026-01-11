%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% One_Oscillating_noise_multi_iteration.m
%
% PURPOSE
%   Run many independent “cells” (iterations) of the noisy single
%   sister-chromosome model to estimate population statistics for:
%   - CM (center-of-mass) oscillation amplitude & period
%   - KK (inter-kinetochore) oscillation amplitude & period
%    Chromosome activity metrics (Avg_KE, mean |v_L|, mean |v_R|)
%
% MODEL LINEAGE
%   This script matches the dynamics of the single-run noise model
%   (one_oscillating_chromosome_noise.m), but repeats the simulation for
%   N_iter independent runs and aggregates summary statistics.
%
% OUTPUTS 
%   - Printed summary statistics across iterations for:
%       CM amplitude1/2, CM period, KK amplitude1/2, KK period, Average KE
%   - Optional boxplots for each metrics 
%   - CSV files:
%       CM_oscillation_amplitude.csv   (Amplitude1, Amplitude2)
%       CM_oscillation_period.csv      (Period)
%       KK_oscillation_amplitude.csv   (Amplitude1, Amplitude2)
%       KK_oscillation_period.csv      (Period)
%       Average_KE.csv                 (Avg_KE per iteration)
% 
% FIGURE REFERENCE
%   This script was used to generate simulation data shown in Figure 3D.
%
%   Parameter sets corresponding to experimental conditions:
%     - Control (Ctrl):
%         Alpha scale = 3.20
%         Kct         = 14.8  pN/µm
%
%     - Nocodazole (Noc):
%         Alpha scale = 3.30
%         Kct         = 14.8  pN/µm
%
%     - Trichostatin A (TSA):
%         Alpha scale = 3.30
%         Kct         = 11.0  pN/µm
%   Note:
%     Alpha is defined as Alpha = n_dot * (Alpha scale) / (1 - Beta).
%
% AUTHOR
%   Jiali Zhu, 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;close all; 

% ---------------- Simulation settings ----------------
N_iter = 150;            % number of simulated cells (independent iterations)
dt = 2e-3;               % time step (min)
Nsteps = 5000;           % number of simulation steps
duration = 10;           % Duration (min) used in KE analysis
% Results containers:
% results_osc: one row per iteration with CM/KK oscillation measurements 
% results_ke : one row per iteration with Avg_KE, vL, vR
results_osc = table();  
results_ke  = table();  

% ---------------- Noise parameters ----------------
noise = 0.05;   % magnitude of uniform parameter perturbations for n_dot and Nbar
noise_b = 0.05; % magnitude of common-mode positional jitter added each time step

for iter = 1:N_iter
    % ---------------- Model parameters (per iteration) ----------------
    Kct = 14.8;   % pN/µm  centromere spring constant
    I0 = 2;       % µm     centromere rest length
    Kkt = 1;      % pN/µm  KT–MT spring constant
    Gamma = 0.1;  %        effective drag coefficient
    Beta = 0.7;   %        scaling for MT tip velocity
    Nmax = 25;    % count  max number of KT–MT attachments
    %n_dot and Nbar are perturbed ONCE per iteration to mimic cell-to-cell variability.
    n_dot = 1*(1 + noise*(2*rand(1) - 1));    % baseline attachment gain term
    Nbar = 20*(1 + noise*(2*rand(1) - 1));    % attachment number at steady state (used for initialization)
    Lambda = n_dot / (Nbar);                  % detachment rate constant
    Alpha = n_dot * 3.40/ (1 - Beta);         % velocity-attachment coupling coefficient (scale factor = 3.4)  
    epsilon = 0.1;                            % small symmetry-breaking perturbation (dimensionless)

    % ---------------- Pre-allocate state variables ----------------
    xL = zeros(Nsteps+1,1); xR = zeros(Nsteps+1,1);         % left/right MT tip x-position over time
    NL = zeros(Nsteps+1,1); NR = zeros(Nsteps+1,1);         % left/right KT–MT attachment number over time 
    cL = zeros(Nsteps+1, 2, 1); cR = zeros(Nsteps+1, 2, 1); % left/right chromosome position: (time, [x y], chromosome index)
    vL = zeros(Nsteps, 1); vR = zeros(Nsteps, 1);           % left/right chromosome velocity over time

    % ---------------- Initial conditions ----------------
    % Initial attachment numbers (slightly perturbed to break symmetry)
    NR(1) = Nbar*(1 - epsilon);
    NL(1) = Nbar*(1 + epsilon);
    % Initial chromosome positions (x in µm; y fixed at 0) 
    cL(1,1,1) =-I0/2-epsilon;
    cR(1,1,1) =I0/2+epsilon;
    cL(:,2,1) = 0;
    cR(:,2,1) = 0;
    % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
    xbar = (Kct * I0) / (Nbar * Kkt);
    xL(1) = cL(1,1,1) - 0.3;
    xR(1) = cR(1,1,1) + 0.3;
    % Initial chromosome velocities
    vL(1) = 0; vR(1) = 0;

    % ---------------- Time integration loop (explicit Euler) ----------------
    for t = 1:Nsteps
       % ----- Force calculations -----
        F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));  % force on left chromosome from KT–MT bundle
        F_CT_L = Kct * (cR(t,1,1) - cL(t,1,1) - I0);  % Centromere forces on Left Chromosome
        F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));  % force on right chromosome from KT–MT bundle
        F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0); % Centromere forces on Right Chromosome
        % ----- Velocities (overdamped dynamics) ----- 
        vL(t) = (F_KT_L + F_CT_L)/Gamma;
        vR(t) = (F_KT_R + F_CT_R)/Gamma;
        % ----- Update chromosome positions ----- 
        dx_diff = noise_b * randn(1,1) * sqrt(dt);   % Common-mode positional jitter
        cL(t+1,1,1) = cL(t,1,1) + dt*vL(t) + dx_diff;
        cR(t+1,1,1) = cR(t,1,1) + dt*vR(t) + dx_diff;
        % ----- Update MT tip positions -----
        xL(t+1) = xL(t) + dt * Beta * vL(t);
        xR(t+1) = xR(t) + dt * Beta * vR(t);
       % ----- Update attachment numbers -----
        NL(t+1) = NL(t) + dt * (...
            n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
            - Lambda * NL(t));
        NR(t+1) = NR(t) + dt * (...
            n_dot + (Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
            - Lambda * NR(t));
    end

    % ---------------- Oscillation analysis (per iteration) ----------------
    % oscillation_measurement returns summary rows for CM and KK (if detected).
    % Append summary rows only if values are finite (skip NaN/Inf cases).
    [table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, false);
    if ~isempty(table_summary)
    % Check all numeric entries are finite across ALL rows (CM & KK)
    nums = table2array(table_summary(:, vartype('numeric')));
        if ~isempty(nums) && all(isfinite(nums(:)))
            results_osc = [results_osc; table_summary];   % append both rows
        else
            fprintf('Iter %d: NaN/Inf found in oscillation summary — skipping append.\n', iter);
        end
    else
        fprintf('Iter %d: no oscillation rows — skipping append.\n', iter);
    end
     % ---------------- Activity analysis (per iteration) ----------------
     % Store one row per iteration with Avg_KE and mean |v| metrics.
    activity_table = chromosome_activity_measurement(cL, cR, dt,duration);
    ke_row = table(iter, activity_table.Avg_KE, activity_table.vL, activity_table.vR, ...
    'VariableNames', {'iteration','Avg_KE','vL','vR'});
    results_ke = [results_ke; ke_row];

    % Append to results % Only Oscillating ones 
    fprintf('[%s] Iteration %2d/%2d done\n', datestr(now,'HH:MM:SS'), iter, N_iter);
end

% ---------------- Summary statistics across simulated cells ----------------
% Metrics reported below: CM amp1/amp2/period, KK amp1/amp2/period, Average KE
fprintf('\nSummary Statistics Across Simulated Cells:\n');
fprintf('%-20s %6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
    'Metric','count','Mean','Std','Min','25%','Median','75%','Max');

% Make sure 'metric' is easy to filter
metric_str = string(results_osc.metric);

% Pull vectors for each metric from the correct table
cm_amp1 = results_osc.mean_amplitude1(metric_str=="CM");
cm_amp2 = results_osc.mean_amplitude2(metric_str=="CM");
cm_per = results_osc.mean_period(  metric_str=="CM");
kk_amp1 = results_osc.mean_amplitude1(metric_str=="KK");
kk_amp2 = results_osc.mean_amplitude2(metric_str=="KK");
kk_per = results_osc.mean_period(  metric_str=="KK");
avg_ke = results_ke.Avg_KE;   % one row per iteration

% Define metrics and values
metric_defs = {
    'CM Amplitude1', cm_amp1;
    'CM Amplitude2', cm_amp2;
    'CM Period',    cm_per;
    'KK Amplitude1', kk_amp1;
    'KK Amplitude2', kk_amp2;
    'KK Period',    kk_per;
    'Average KE',   avg_ke;
};

% Loop through metrics and display stats (+ optional boxplot)
for i = 1:size(metric_defs,1)
    label = metric_defs{i,1};
    vals  = metric_defs{i,2};

    n  = numel(vals);
    m  = mean(vals, 'omitnan');
    sd = std(vals,  'omitnan');
    mn = min(vals,  [], 'omitnan');
    q1 = prctile(vals,25);
    md = median(vals, 'omitnan');
    q3 = prctile(vals,75);
    mx = max(vals,  [], 'omitnan');

    fprintf('%-20s %6d %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f\n', ...
        label, n, m, sd, mn, q1, md, q3, mx);

    % Boxplot (optional)
    % figure; boxplot(vals);
    % title(label, 'FontSize', 16);
    % ylabel(label, 'FontSize', 14);
    % set(gca, 'FontSize', 14); set(gcf, 'Color', 'w'); grid on;
end


%% Save results as 4 csv files of the following format: 
% Format:
%   - CM/KK files: iteration, chromosome, Amplitude1, Amplitude2 (or Period)
%   - KE file:     iteration, Measurement  (Avg_KE)

% Masks using the updated results_osc table
maskCM = strcmp(results_osc.metric, 'CM');
maskKK = strcmp(results_osc.metric, 'KK');

% 1) CM Amplitude (both components)
T = results_osc(maskCM, {'iteration','chromosome','mean_amplitude1','mean_amplitude2'});
T = renamevars(T, {'mean_amplitude1','mean_amplitude2'}, {'Amplitude1','Amplitude2'});
writetable(T, 'CM_oscillation_amplitude.csv');

% 2) CM Period
T = results_osc(maskCM, {'iteration','chromosome','mean_period'});
T = renamevars(T, {'mean_period'}, {'Period'});
writetable(T, 'CM_oscillation_period.csv');

% 3) KK Amplitude (both components)
T = results_osc(maskKK, {'iteration','chromosome','mean_amplitude1','mean_amplitude2'});
T = renamevars(T, {'mean_amplitude1','mean_amplitude2'}, {'Amplitude1','Amplitude2'});
writetable(T, 'KK_oscillation_amplitude.csv');

% 4) KK Period
T = results_osc(maskKK, {'iteration','chromosome','mean_period'});
T = renamevars(T, {'mean_period'}, {'Period'});
writetable(T, 'KK_oscillation_period.csv');

% 5) Average KE (from results_ke; per-iteration only)
T = results_ke(:, {'iteration','Avg_KE'});
T.Properties.VariableNames{'Avg_KE'} = 'Measurement';
writetable(T, 'Average_KE.csv');

