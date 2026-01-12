%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi_chromosomes_noise_iteration.m
%
% PURPOSE
%   Run a multi-iteration simulation of N independent sister-chromosome pairs
%   with small stochasticity added to:
%     (i) kinetic parameters per chromosome (n_dot, Nbar), and
%     (ii) chromosome x–positions via positional jitter (driving both COM & KK).
%   After each valid iteration, perform four analyses and pool results:
%     1) Oscillation metrics for CM & KK: amplitude1, amplitude2, period
%     2) MSV / activity metrics per chromosome: Avg_MSV, vL, vR
%     3) Pearson’s r between all CM pairs 
%     4) Cross-correlation vs neighbor order and vs lag τ 
%
% WHAT’S NEW VS. SINGLE-RUN SCRIPTS
%   • Iterative loop until `target_iterations` valid runs are collected
%   • Per-iteration parameter noise drawn independently for each chromosome
%   • Per-iteration pooling into master tables:
%       - results_osc      (oscillation rows, includes `iteration`)
%       - results_MSV      (Avg_MSV, vL, vR per chromosome & iteration)
%       - results_pearsons (pairwise Pearson r with `pair` = [i j], per iteration)
%       - records_table    (cross-correlation records across τ & neighbor)
%
% OUTPUTS 
% FIGURES 
%   – Histogram of pooled Pearson’s r (across all iterations)
%   – Cross-correlation:
%       (a) Avg C vs neighbor level, grouped by τ bands
%       (b) Avg C vs τ for each neighbor level
%   – (Occasional) CM traces for sanity check every 20th valid iteration
%
% DEPENDENCIES (on MATLAB path)
%   oscillation_measurement.m
%   chromosome_activity_measurement.m
%   pearsons_correlation_analysis.m   
%   crosscorrelation_measurement.m    
%
%
% AUTHOR
%   Jiali Zhu, 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;
close all;
% Time Information 
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
Nchromosomes=5; %Number of chromosome pairs 

%Parameter Values for chromosome movement
noise=0.05;
noise_b=0.05;
Kct = 14.20; % pN/µm centromee spring constant, Harasymiw et al, 2019
I0 = 2; % µm Rest length of centromere spring
Kkt = 1;% pN/m kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.1; % kg/s Drag coefficient
Beta = 0.7; % Scaling factor
Nmax = 25; % Maximum number of attachments
% n_dot = 1*(1+noise*(2*rand(1)-1)); %Add noise 
% Nbar=20*(1+noise*(2*rand(1)-1)); % Steady state number of MTs when Ch is centered
% Lambda=n_dot/(Nbar);% s^-2 KMT detach rate Akioshi et al. 2010
%Alpha=n_dot*6.3/(1-Beta); 
epsilon =0.1; % Small perturbation factor
%Parameter Values for coupling sping dynamics
l0=5; % natural length of the spring 
koff_0=40; % base disconnection rate
kon_0=20; % base connection rate 
Nspring = 1; % Number of springs per chromosome not used 
k =2; % interchromosome coupling value, make weak?

target_iterations = 50; % Number of valid rounds of data needed
valid_iterations = 0; % Counter for valid rounds

% Create a structure array to store results
results_osc= table();     % iteration, chromosome, amplitude1, amplitude2, period
results_MSV = table();      % iteration, chromosome, Avg_KE, vL, vR
results_pearsons= table(); % iteration. pair, coefficient 
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
duration = 10; 
% Tau setup with 5 s resampling of chromosome movement 
tau_steps   = 1:30;                 % 5 s, 10 s, ..., 150 s (in resampled steps)
tau_seconds = 5 * tau_steps;        % for labels only
% Three bands of 10 lags each: [5–50 s], [55–100 s], [105–150 s]
tau_groups_steps = reshape(tau_steps, 10, []).';   % 3x10 in step units
tau_groups_sec   = reshape(tau_seconds, 10, []).'; % 3x10 in seconds (optional)


while valid_iterations < target_iterations
    %Varying parameters for each chromoosme
    Nbar = zeros(Nchromosomes,1); %The steady-state number of MTs when the Ch is centered
    for i = 1:Nchromosomes
        Nbar(i) = 20*(1+noise*(2*rand(1)-1));
    end
    n_dot = zeros(Nchromosomes,1); %The steady-state number of MTs when the Ch is centered
    for i = 1:Nchromosomes
        n_dot(i) = 1*(1+noise*(2*rand(1)-1));
    end
    Lambda = zeros(Nchromosomes, 1);% s^-2 KMT detach rate Akioshi et al. 2010
    for i = 1:Nchromosomes
        Lambda(i)=n_dot(i)/Nbar(i);
    end
    Alpha = zeros(Nchromosomes, 1);
    for i = 1:Nchromosomes
        Alpha(i)=n_dot(i)*5.20/(1-Beta);
    end
 
    % Set up the vectors for results
    cL = zeros(Nsteps+1, 2,Nchromosomes); % Left chromosome position (time,axis,Nchromosome) 
    cR = zeros(Nsteps+1, 2,Nchromosomes); % Right chromosome position (time,axis,Nchromosome) 
    xL = zeros(Nsteps+1, Nchromosomes); % Left MT tip position 
    xR = zeros(Nsteps+1, Nchromosomes); % Right MT tip position 
    NL = zeros(Nsteps+1, Nchromosomes); % Left MT-KT attachment amount 
    NR = zeros(Nsteps+1, Nchromosomes); % Right MT-KT attachment amount 
    vL = zeros(Nsteps, Nchromosomes); % Left chromosome position 
    vR = zeros(Nsteps, Nchromosomes); % Right chromosome position
    
    % Initial Condition 
    % Fixed y positions
    y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
    % Coupling state
    flag_c=zeros(Nchromosomes,Nchromosomes);  % current state, update at each timestep 
    flag_cp=flag_c; % Previous state 
    flag_c_sum_accumulated = zeros(Nsteps, 1); % To accumulate sum of flag_c over time
    for i = 1: Nchromosomes 
       % Number of attachments
        NR(1,i)=Nbar(i)*(1-epsilon);
        NL(1,i)=Nbar(i)*(1+epsilon); 
        % Ch position 
        cL(1,:,i) =-I0/2-epsilon*(2*rand(1)-1);%[x,y]
        cR(1,:,i) =I0/2+epsilon*(2*rand(1)-1);%[x,y]
        cL(:,2,i) =  y_positions(i);
        cR(:,2,i) =  y_positions(i);
        % MT tip position 
        xL(1,i) = cL(1,1,i)-0.3;
        xR(1,i) = cR(1,1,i)+0.3;
        vL(1,i)=0;
        vR(1,i)=0;
    end 
    
    % ODE solver loop
    for t = 1:Nsteps
    %Update connection state 
    CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);  % N x 1 vector
    YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);  % N x 1 vector
    [flag_c,dist]=spring_connect(CM,YM,flag_cp,Nchromosomes,dt,l0,koff_0,kon_0);
    flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
    flag_cp = flag_c;  % Update prior state
        for i=1:Nchromosomes
            % Force Calculations 
            %Coupling force 
            F_coupling_x=0;
            for j = 1:Nchromosomes
                if j ~= i  % Skip self
                    delta_x = CM(i) - CM(j);
                    F_coupling_x = F_coupling_x - k * delta_x * flag_c(i, j);
                end
            end
            % Random Force 
            F_noise = 0*noise_b *(randn(1,1)) / sqrt(dt); % Brownian Force
            % Define centromere and MT forces on the given chromosome
            F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i)); % MT forces on Left Chromosome
            F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Left Chromosome
            F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i)); % MT forces on Right Chromosome
            F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Right Chromosome
        
            % Calculate velocities 
            vL(t,i) = (F_KT_L + F_CT_L+F_noise+F_coupling_x)/Gamma;
            vR(t,i) = (F_KT_R + F_CT_R+F_noise+F_coupling_x)/Gamma;
            
            % Update Positions 
            % Calculate the current chromosome position 
            dx_diff=noise_b *(randn(1))*sqrt(dt);
            cL(t+1,1,i) = cL(t,1,i) + dt * vL(t,i)+dx_diff ;
            cR(t+1,1,i) = cR(t,1,i) + dt * vR(t,i)+dx_diff ;
            % Maintain fixed Y positions
            cL(t+1,2,i) = cL(1,2,i);
            cR(t+1,2,i) = cR(1,2,i);
            % Calculate the current microtubule tip position in both directions
            xL(t+1,i) = xL(t,i) + dt * Beta * (vL(t,i));
            xR(t+1,i) = xR(t,i) + dt * Beta * (vR(t,i));
            % Update number of attachments 
            NL(t+1,i) = NL(t,i) + dt * (...
                n_dot(i) +(-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
                - Lambda(i)* NL(t,i));
            NR(t+1,i) = NR(t,i) + dt * (...
                n_dot(i) + (Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
                - Lambda(i)* NR(t,i));  
        end
    end
    %%%%%%%%%%%%%%%%%% Analysis at the end of each iteration %%%%%%%%%%%%%%
    % Check for NaNs in cL and cR
    if any(isnan(cL),'all') || any(isnan(cR),'all')
        fprintf('NaN detected, skipping this round and starting a new one\n');
        continue; % Skip this round and start a new iteration
    end
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
    % Ocsillation Analysis 
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
    % Activity analysis 
    activity_all= chromosome_activity_measurement(cL(:,:,:), cR(:,:,:), dt);
    activity_row = table( ...
        repmat(valid_iterations+1, height(activity_all), 1), ...
        activity_all.PairID, ...
        activity_all.Avg_KE, ...
        activity_all.vL, ...
        activity_all.vR, ...
        'VariableNames', {'iteration','chromosome','Avg_MSV','vL','vR'} );
    results_MSV = [results_MSV; activity_row];
    % Pearson's analysis 
    pearson_table = pearsons_correlation_analysis(cL, cR, dt, duration);
    pearson_row = table( ...
        repmat(valid_iterations+1, height(pearson_table), 1), ...
        pearson_table.pair, ...
        pearson_table.pearson_r, ...
        'VariableNames', {'iteration','pair','pearson_r'});
    results_pearsons = [results_pearsons; pearson_row]; 
    % Crosscorrelation Analysis 
    sampling_interval = round((5/60) / dt);           % native steps per 5 s
    total_steps   = min(round(duration / dt), Nsteps);  % cap by duration & Nsteps    
    % Sample L/R x positions and compute CMx (x only)
    cL_sampled   = cL(1:sampling_interval:total_steps, 1, :);                 % [Nsamp x 1 x Nchromosomes]
    cR_sampled   = cR(1:sampling_interval:total_steps, 1, :);
    CMx_sampled  = ((cL_sampled + cR_sampled) / 2)/sqrt(10);     % [Nsamp x 1 x Nchromosomes]
    CMx          = reshape(CMx_sampled, [], Nchromosomes);  % [Nsamp x Nchromosomes]
    Nsteps_res   = size(CMx, 1);                      % resampled length
    
    % Safety: must have enough points for max τ
    if Nsteps_res <= max(tau_steps)
        warning('Not enough resampled time points (%d) for max tau (%d). Skipping iteration.', ...
                Nsteps_res, max(tau_steps));
        continue;
    end
    iteration_id=valid_iterations + 1;
    % Compute correlations at selected tau values
    for tau = tau_steps
        [C, N] = crosscorrelation_measurement(CMx, tau, Nsteps_res, Nchromosomes);
        for neighbor_idx = 1:(Nchromosomes - 1)
            if N(neighbor_idx) > 0
               C_norm = C(neighbor_idx) / N(neighbor_idx);
            else
                C_norm = NaN;
            end
            records(end+1) = struct( ...
                'iteration', iteration_id, ...
                'tau',       tau, ...            % <-- store τ in resampled steps (1..30)
                'neighbor',  neighbor_idx, ...
                'C_value',   C_norm, ...
                'N_value',   N(neighbor_idx) );
        end
    end
    % If no NaNs, proceed with correlation calculation
    fprintf('Iteration %d complete\n', valid_iterations + 1);
    valid_iterations = valid_iterations + 1; % Increment valid rounds counte
end
records_table = struct2table(records);
%Pearson's correlation coefficient 
figure; 
histogram(results_pearsons.pearson_r, 20);                        % 20 bins; adjust if you like
xlabel('Pearson''s r (CM vs CM)');
ylabel('Count');
title('Distribution of Pearson''s Correlations ');
grid on; box on; hold on;
xline(0, '--k');                                % reference at r=0
xline(mean(results_pearsons.pearson_r,'omitnan'), 'r-', 'LineWidth', 1.5);  % mean r
hold off;
%Visualizing the averaged results from multiple iterations 
% Print the average oscillatioin and MSV
% Summary Statistics and Boxplots for All 5 Measurements (

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
%% Save results as 4 csv files of the following format: 
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
records_export = records_table;
records_export.tau_seconds = records_export.tau * 5;   % 1 resampled step = 5 s
writetable(records_export, 'Crosscorr_records.csv');


%% Visualizing the averaged results from multiple iterations 
% Print the average oscillatioin and MSV
% Summary Statistics and Boxplots for All 5 Measurements (

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

%% Pearson's correlation coefficient 
figure; 
histogram(results_pearsons.pearson_r, 20);                        % 20 bins; adjust if you like
xlabel('Pearson''s r (CM vs CM)');
ylabel('Count');
title('Distribution of Pearson''s Correlations ');
grid on; box on; hold on;
xline(0, '--k');                                % reference at r=0
xline(mean(results_pearsons.pearson_r,'omitnan'), 'r-', 'LineWidth', 1.5);  % mean r
hold off;

%%
% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups_steps, 1));
for g = 1:size(tau_groups_steps, 1)
    current_group = tau_groups_steps(g, :);   % in steps
    subset  = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    t_start = current_group(1)  * 5;   % seconds
    t_end   = current_group(end) * 5;  % seconds
    errorbar(grouped.neighbor, grouped.mean_C_value, grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d s', t_start, t_end), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
xticks(1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;


% Plot 2: C vs Time Lag for Neighbor Level                                                                                                                   
figure; hold on;
for n = 1:min(Nchromosomes-1, 3)   % show 1..3 if available
    subset  = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");
    tau_sec = grouped.tau * 5;  % steps -> seconds

    errorbar(tau_sec, grouped.mean_C_value, grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('Neighbor %d', n), 'LineWidth', 2);
end
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for Neighbor Levels', 'FontSize', 16);
legend('Location', 'best'); grid on; box on; set(gca, 'FontSize', 12);

set(gca, 'FontSize', 12);
