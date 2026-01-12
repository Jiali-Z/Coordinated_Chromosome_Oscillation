%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi_chromosomes_noise_iteration.m
%
% PURPOSE
%   Run many independent “cells” (iterations) of a multi-chromosome spindle,
%   where each “cell” contains N independent sister-chromosome pairs moving in
%   1D (x). This is a negative control: because chromosome pairs are NOT
%   mechanically coupled, cross-chromosome coordination metrics (Pearson r,
%   displacement cross-correlation) should be ~0
%   on average (up to finite-sample noise). 
%
%   For each valid iteration (no NaNs), run and pool:
%     1) Oscillation metrics (CM & KK): amplitude1, amplitude2, period
%     2) Activity metrics per chromosome: Avg_KE, |vL|, |vR|
%     3) Pearson’s r across all chromosome-pair CM traces (all i<j)
%     4) Displacement cross-correlation as a function of:
%           (a) neighbor order
%           (b) lag time tau
%        (reported as normalized C_value = C(neighbor)/N(neighbor))
%
% WHAT'S DIFFERENT FROM single-iteration multi_chromosomes_noise scripts
%   - Repeats simulations until `target_iterations` valid runs are collected
%   - Pools per-iteration outputs into master tables/records for statistics
%   - Optional sanity-check CM plot every ~20 valid iterations
%
% OUTPUTS
%   In-memory pooled results:
%     - results_osc       : oscillation summaries (rows labeled metric = CM / KK)
%     - results_MSV       : activity metrics per chromosome per iteration
%     - results_pearsons  : Pearson r for all CM pair combinations per iteration
%     - records_table     : cross-correlation records across tau & neighbor
%
%   Saved CSV files (for downstream plotting / Python):
%     - CM_oscillation_amplitude.csv
%     - CM_oscillation_period.csv
%     - KK_oscillation_amplitude.csv
%     - KK_oscillation_period.csv
%     - Average_MSV.csv
%     - Pearsons_coefficients.csv
%     - Crosscorr_records.csv   (includes tau_seconds)
%
%   Figures (summary across pooled iterations):
%     - Histogram of Pearson r
%     - Avg cross-correlation vs neighbor level (grouped by tau bands)
%     - Avg cross-correlation vs tau for each neighbor level
%     - (Optional) CM traces every 20th valid iteration (sanity check)
%
% DEPENDENCIES (must be on MATLAB path)
%   oscillation_measurement.m
%   chromosome_activity_measurement.m
%   pearsons_correlation_analysis.m
%   crosscorrelation_measurement.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC–Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

% ---------------- Time information ---------------- 
dt = 2e-3;         % time step (min)
Nsteps = 5000;     % number of simulation steps 

% ---------------- Number of chromosome ---------------- 
Nchromosomes=5; %Number of chromosome pairs 

% ---------------- Noise parameters ----------------
noise=0.05;        % magnitude of uniform parameter perturbations for n_dot and Nbar
noise_b=0.05;      % magnitude of common-mode positional jitter added each time step

% ---------------- Model parameters ----------------
Kct = 12.3;       % pN/µm  centromere spring constant  
I0 = 2;            % µm     centromere rest length
Kkt = 1;           % pN/µm  KT–MT spring constant 
Gamma = 0.1;       %        effective drag coefficient 
Beta = 0.7;        %        scaling for MT tip velocity
Nmax = 25;         % count  max number of KT–MT attachments
epsilon =0.1;      % small symmetry-breaking perturbation (dimensionless) 

% ---------------- Number of iterations ---------------- 
target_iterations = 50; % Number of valid rounds of data needed
valid_iterations = 0; % Counter for valid rounds

% ---------------- Pooled outputs across iterations ----------------
results_osc= table();          % iteration, chromosome, amplitude1, amplitude2, period
results_MSV = table();         % iteration, chromosome, Avg_KE, vL, vR
results_pearsons= table();     % iteration, pair, coefficient 
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {}); % Cross-correlation records stored as a struct array, converted to table at end

% ---------------- Analysis parameters ----------------
duration = 10;                    % minutes used by pearsons_correlation_analysis (sampling window)
tau_step = round(5 / (dt * 60));  % steps per 5 seconds (dt is in minutes)
tau_values = tau_step * (1:30);  % 5s..150s in steps    
tau_groups = reshape(tau_values, 10, [])'; % 3 tau-bands, each containing 10 tau values

% ---------------- Parameter noise per chromosome ----------------
% Each chromosome i gets its own kinetic parameters (chromosome-to-chromosome variability).
% These are drawn ONCE (at initialization) and held fixed during the simulation.
while valid_iterations < target_iterations
    Nbar   = zeros(Nchromosomes, 1);  % steady-state attachment number (initialization reference)
    n_dot  = zeros(Nchromosomes, 1);  % baseline gain term for attachments
    Lambda = zeros(Nchromosomes, 1);  % detachment rate constant
    Alpha  = zeros(Nchromosomes, 1);  % velocity-attachment coupling coefficient
    for i = 1:Nchromosomes
        Nbar(i)  = 20 * (1 + noise * (2*rand(1) - 1));
        n_dot(i) =  1 * (1 + noise * (2*rand(1) - 1));
        Lambda(i) = n_dot(i) / Nbar(i);
        Alpha(i)  = n_dot(i) * 6.2 / (1 - Beta);
    end
 
    % ---------------- Pre-allocate state variables ----------------
    cL = zeros(Nsteps+1, 2,Nchromosomes);       % left chromosome position: (time, [x y], chromosome index)
    cR = zeros(Nsteps+1, 2,Nchromosomes);       % right chromosome position: (time, [x y], chromosome index)
    xL = zeros(Nsteps+1, Nchromosomes);         % left MT tip x-position over time
    xR = zeros(Nsteps+1, Nchromosomes);         % right MT tip x-position over time
    NL = zeros(Nsteps+1, Nchromosomes);         % left KT–MT attachment number over time 
    NR = zeros(Nsteps+1, Nchromosomes);         % right KT–MT attachment number over time 
    vL = zeros(Nsteps, Nchromosomes);           % left chromosome velocity over time 
    vR = zeros(Nsteps, Nchromosomes);           % right chromosome velocity over time
    
    % ---------------- Initial conditions ----------------
    % Fixed y positions are assigned
    y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
    % Each Pairs have their own 
    for i = 1: Nchromosomes 
        % Initial attachment numbers (slightly perturbed to break symmetry)
        NR(1,i)=Nbar(i)*(1-epsilon);
        NL(1,i)=Nbar(i)*(1+epsilon); 
         % Initial chromosome positions (x in µm; y fixed) 
        cL(1,:,i) =-I0/2-epsilon*(2*rand(1)-1); % left sister initial x-position
        cR(1,:,i) =I0/2+epsilon*(2*rand(1)-1);  % right sister initial x-position
        cL(:,2,i) =  y_positions(i);            % y fixed
        cR(:,2,i) =  y_positions(i);            % y fixed
        % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
        xL(1,i) = cL(1,1,i)-0.3;
        xR(1,i) = cR(1,1,i)+0.3;
        % Initial chromosome velocities
        vL(1,i)=0;
        vR(1,i)=0;
    end 
    
    % ---------------- Time integration loop (explicit Euler) ----------------
    for t = 1:Nsteps
        for i=1:Nchromosomes 
        % ----- Force calculations ----- 
        % Define centromere and MT forces on the given chromosome
        F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));  % force on left chromosome from KT–MT bundle
        F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0);      % Centromere forces on Left Chromosome
        F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));  % force on right chromosome from KT–MT bundle
        F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0);     % Centromere forces on Right Chromosome
        % ----- Velocities (overdamped dynamics) -----
        vL(t,i) = (F_KT_L + F_CT_L)/Gamma;
        vR(t,i) = (F_KT_R + F_CT_R)/Gamma;
        
        % ----- Update chromosome positions -----
        % Common-mode positional jitter: the same dx_diff is added to both sisters,
        dx_diff=noise_b *(randn(1))*sqrt(dt);
        cL(t+1,1,i) = cL(t,1,i) + dt * vL(t,i)+dx_diff ;
        cR(t+1,1,i) = cR(t,1,i) + dt * vR(t,i)+dx_diff ;
         % Maintain fixed Y positions
        cL(t+1,2,i) = cL(1,2,i);
        cR(t+1,2,i) = cR(1,2,i);
        % ----- Update MT tip positions -----
        xL(t+1,i) = xL(t,i) + dt * Beta * (vL(t,i));
        xR(t+1,i) = xR(t,i) + dt * Beta * (vR(t,i));
        % ----- Update attachment numbers -----
        NL(t+1,i) = NL(t,i) + dt * (...
            n_dot(i) +(-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
            - Lambda(i)* NL(t,i));
        NR(t+1,i) = NR(t,i) + dt * (...
            n_dot(i) + (Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
            - Lambda(i)* NR(t,i));  
        end
    end
    %%%%%%%%%%%%%%%%%% Analysis at the end of each iteration %%%%%%%%%%%%%%
    % Skip failed runs (NaNs), otherwise append results to pooled tables.
    % Reject runs with NaNs in trajectories
    if any(isnan(cL),'all') || any(isnan(cR),'all')
        fprintf('NaN detected, skipping this round and starting a new one\n');
        continue; % Skip this round and start a new iteration
    end
    
    % Optional sanity-check plot every ~20 valid iterations
    if (mod(valid_iterations,20)==0)
        % side-by-side CM positions overtime 
        figure;
        hold on;
        time = (0:Nsteps) * dt;
        for i = 1:Nchromosomes
            CMx = ((cL(:,1,i) + cR(:,1,i)) / 2)+i;
            plot(time, CMx, 'LineWidth', 2);
        end
        xlabel('Time (min)');
        ylabel('Center of Mass X Position (µm)');
        title('Chromosome CM Oscillations');
        legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
        set(gca, 'FontSize', 14);
        grid on
    end
    % ---- (1) Oscillation analysis: (CM + KK) amplitude+period ---- 
     [~, ~, table_summary] = oscillation_measurement(cL(:,:,:), cR(:,:,:), dt, valid_iterations+1, false);
    if ~isempty(table_summary)
        % Check all numeric entries are finite across ALL rows (CM & KK)
        nums = table2array(table_summary(:, vartype('numeric')));
        if ~isempty(nums) && all(isfinite(nums(:)))
            results_osc= [results_osc; table_summary];   % append both rows
        else
            fprintf('Iter %d: NaN/Inf found in oscillation summary — skipping append.\n', valid_iterations+1);
        end
    else
        fprintf('Iter %d: no oscillation rows — skipping append.\n', valid_iterations+1);
    end
    % ---- (2) Activity metrics ----
    activity_all= chromosome_activity_measurement(cL(:,:,:), cR(:,:,:), dt);
    activity_row = table( ...
        repmat(valid_iterations+1, height(activity_all), 1), ...
        activity_all.PairID, ...
        activity_all.Avg_KE, ...
        activity_all.vL, ...
        activity_all.vR, ...
        'VariableNames', {'iteration','chromosome','Avg_MSV','vL','vR'} );
    results_MSV = [results_MSV; activity_row];
    % ---- (3) Pearson correlations between CM traces of all chromosome pairs ---- 
    pearson_table = pearsons_correlation_analysis(cL, cR, dt, duration);
    pearson_row = table( ...
        repmat(valid_iterations+1, height(pearson_table), 1), ...
        pearson_table.pair, ...
        pearson_table.pearson_r, ...
        'VariableNames', {'iteration','pair','pearson_r'});
    results_pearsons = [results_pearsons; pearson_row]; 
    % ---- (4) Displacement cross-correlation vs neighbor & tau ---- 
    CMx = (cL(:,1,:) + cR(:,1,:)) / 2;
    iteration_id=valid_iterations + 1;
    % Compute correlations at selected tau values
    for tau = tau_values
        [C, N] = crosscorrelation_measurement(CMx, tau, Nsteps, Nchromosomes);
        for neighbor_idx = 1:(Nchromosomes - 1)
            if N(neighbor_idx) > 0
                C_norm = C(neighbor_idx) / N(neighbor_idx);
            else
                C_norm = NaN;
            end
            records(end+1) = struct( ...
                'iteration', iteration_id, ...
                'tau', tau, ...
                'neighbor', neighbor_idx, ...
                'C_value', C_norm, ...
                'N_value', N(neighbor_idx) ...
            );
        end
    end
    % If no NaNs, proceed with correlation calculation
    fprintf('Iteration %d complete\n', valid_iterations + 1);
    valid_iterations = valid_iterations + 1; % Increment valid rounds counte
end
records_table = struct2table(records);

% =======================================================================
% VISUALIZATION of output 
% This section summarizes the distributions of:
%   - CM/KK oscillation metrics (amplitudes + periods)
%   - Activity metric (Avg_MSV column; currently populated with Avg_KE in code)
%   - Pearson r across all chromosome-pair CM traces
%   - Displacement cross-correlation vs neighbor and vs lag time tau
% =======================================================================

% ---------------- Summary statistics across simulated “cells” ----------------
% Print a compact stats table for each pooled metric:
%   count, mean, std, min, quartiles, max
fprintf('\nSummary Statistics Across Simulated Cells:\n');
fprintf('%-20s %6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
    'Metric','count','Mean','Std','Min','25%','Median','75%','Max');

% Convert metric labels to string for easy filtering (expects 'CM' and 'KK')
metric_str = string(results_osc.metric);

% Pull pooled values from results_osc (one value per row in results_osc)
cm_amp1 = results_osc.mean_amplitude1(metric_str=="CM");
cm_amp2 = results_osc.mean_amplitude2(metric_str=="CM");
cm_per = results_osc.mean_period(  metric_str=="CM");
kk_amp1 = results_osc.mean_amplitude1(metric_str=="KK");
kk_amp2 = results_osc.mean_amplitude2(metric_str=="KK");
kk_per = results_osc.mean_period(  metric_str=="KK");
MSV = results_MSV.Avg_MSV;   % one row per iteration

% Define metrics and values
metric_defs = {
    'CM Amplitude1', cm_amp1;
    'CM Amplitude2', cm_amp2;
    'CM Period',    cm_per;
    'KK Amplitude1', kk_amp1;
    'KK Amplitude2', kk_amp2;
    'KK Period',    kk_per;
    'Average MSV',  MSV;
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

    % % Boxplot (optional)
    % figure; boxplot(vals);
    % title(label, 'FontSize', 16);
    % ylabel(label, 'FontSize', 14);
    % set(gca, 'FontSize', 14); set(gcf, 'Color', 'w'); grid on;
end

% ---------------- Pearson correlation distribution (pooled) ----------------
figure; 
histogram(results_pearsons.pearson_r, 20); % 20 bins; adjust as desired
xlabel('Pearson''s r (CM vs CM)');
ylabel('Count');
title('Distribution of Pearson''s Correlations ');
grid on; box on; hold on;
xline(0, '--k');                                % reference at r=0
xline(mean(pearson_table.pearson_r,'omitnan'), 'r-', 'LineWidth', 1.5);  % mean r
hold off;

% ---------------- Cross-correlation summary plots (pooled) ----------------
% records_table contains per-iteration cross-correlation values:
%   neighbor = 1..(Nchromosomes-1)
%   tau      = lag in simulation steps (converted to seconds in plotting)
%   C_value  = normalized correlation value for that (iteration, tau, neighbor)
%
% We plot pooled mean ± std across ALL iterations.

% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
ylim([-1,1]);
xticks(1:1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:Nchromosomes-1
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2);
end
ylim([-1,1]);
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for 1st, 2nd, 3rd Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);

%%
% =====================================================================
% EXPORT (CSV)
% Save pooled results in “long format” CSVs for plotting in Python/R.
% =====================================================================
% 1) CM oscillation amplitude (Amplitude1/Amplitude2)
% 2) CM oscillation period
% 3) KK oscillation amplitude (Amplitude1/Amplitude2)
% 4) KK oscillation period
% 5) Average activity metric (currently Avg_KE, exported as avg_MSV)
% 6) Pearson coefficients (pair1, pair2, pearson_r)
% 7) Cross-correlation records (tau, tau_seconds, neighbor, C_value, N_value)
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

% 3) KK Amplitude 
T = results_osc(maskKK, {'iteration','chromosome','mean_amplitude1','mean_amplitude2'});
T = renamevars(T, {'mean_amplitude1','mean_amplitude2'}, {'Amplitude1','Amplitude2'});
writetable(T, 'KK_oscillation_amplitude.csv');

% 4) KK Period
T = results_osc(maskKK, {'iteration','chromosome','mean_period'});
T = renamevars(T, {'mean_period'}, {'Period'});
writetable(T, 'KK_oscillation_period.csv');

% 5) Average KE 
T = results_MSV(:, {'iteration','chromosome','Avg_MSV'});
T.Properties.VariableNames{'Avg_MSV'} = 'avg_MSV';
writetable(T, 'Average_MSV.csv');
% 7) Pearson's coefficients 
P = results_pearsons;
Pair1 = arrayfun(@(r) P.pair(r,1), (1:height(P))');
Pair2 = arrayfun(@(r) P.pair(r,2), (1:height(P)))';
pearsons_export = table(P.iteration, Pair1, Pair2, P.pearson_r, ...
    'VariableNames', {'iteration','pair1','pair2','pearson_r'});
writetable(pearsons_export, 'Pearsons_coefficients.csv');
% 8) Cross-correlation records (pooled across iterations)
%    Add tau in seconds using your tau_step definition (steps per 5 s)
tau_step = round(5 / (dt * 60));               % steps per 5 s on native grid
records_export = records_table;
records_export.tau_seconds = (records_export.tau / tau_step) * 5;  
writetable(records_export, 'Crosscorr_records.csv')

