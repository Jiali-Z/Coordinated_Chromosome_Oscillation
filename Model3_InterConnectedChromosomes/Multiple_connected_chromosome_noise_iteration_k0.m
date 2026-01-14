%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiple_connected_chromosome_noise_iteration_k.m 
%
%
% PURPOSE
%   Sweep the inter-chromosomal spring stiffness (k; paper notation Kic)
%   in a mechanically-coupled N-chromosome model, and quantify how
%   distance-dependent chromosome coordination changes with k.
%
%
%   For each k value:
%     - Run N_iter independent “cells” (valid iterations; NaN-free)
%     - Simulate N sister-chromosome pairs oscillating in 1D (x)
%     - Mechanically couple chromosome centers-of-mass (CMs) using a stochastic
%       spring network that forms/breaks over time (spring_connect)
%     - Compute displacement cross-correlation as a function of:
%           (a) neighbor order (1..Nchromosomes-1)
%           (b) lag time tau (5–150 s in 5 s steps)
%       and store the normalized correlation:
%           C_norm(neighbor, tau) = C(neighbor, tau) / N(neighbor, tau)
%
%
%HOW THIS SCRIPT IS USED
%   This script is a parameter sweep to test how the *shape* and *magnitude*
%   of the correlation-vs-neighbor curve depend on k (Kic), while keeping the
%   coupling spring effective length l0 (Iic) and kinetics (kon_0, koff_0) fixed.
%
%   Output is cross-correlation table suitable for:
%     - plotting mean±std C_norm vs neighbor for each k (in selected tau
%     bands) (Figure6C)
%     - plotting C_norm vs tau for selected neighbors, stratified by k
%
% OUTPUTS
%   In-memory:
%     - all_records.(k_*) : table of cross-correlation records for each k
%     - combined           : concatenated table with an added l0 column
%
%   Saved CSV files:
%     - crosscorr_records_all_l0.csv
%     - crosscorr_records_l0_<value>.csv   (optional per-k export)
%
%   Figures:
%     - C_norm vs neighbor level for different k values, averaged over τ = 
%       5–50 s, 55–100 s, and 105–150 s
%     - C_norm vs lag time τ (5–105 s)for different k values,for 
%       neighbor orders 1–3 
%
% DEPENDENCIES 
%   spring_connect.m
%   crosscorrelation_measurement.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC–Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% ---------------- Time information ----------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps

% ---------------- Number of chromosome ---------------- 
Nchromosomes = 5;           % Number of chromosome pairs

% ---------------- Inter-chromosome coupling parameters ----------------
k_values = [0 0.2 0.5 1 ];     % coupling strengths to scan
l0        = 3;                 % fixed natural length
koff_0    = 40;
kon_0     = 20;

% ---------------- Number of iterations ---------------- 
N_iter  = 50;               % iterations per k
% ---------------- Noise parameters ----------------
noise   = 0.05;             % magnitude of uniform parameter perturbations for n_dot and Nbar
noise_b = 0.05;             % magnitude of common-mode positional jitter added each time step
epsilon = 0.1;              % small symmetry-breaking perturbation (dimensionless) 

% ---------------- Model parameters ----------------
Kct  = 12.3;                 % pN/µm  centromere spring constant  
I0   = 3;                    % µm     centromere rest length
Kkt  = 1;                    % pN/µm  KT–MT spring constant 
Gamma= 0.1;                  %        effective drag coefficient 
Beta = 0.7;                  %        scaling for MT tip velocity
Nmax = 25;                   % count  max number of KT–MT attachments

% ---------------- Analysis parameters ----------------
duration = 10;              % minutes to analyze (for resampling window)
% Tau setup with 5 s resampling,1 resampled step = 5 seconds
tau_steps   = 1:30;                         % 5 s .. 150 s in 5 s steps
tau_seconds = 5 * tau_steps;                % labels if needed
tau_groups_steps = reshape(tau_steps, 10, []).';   % 3x10 bands: [1..10], [11..20], [21..30]
% tau_groups_sec   = reshape(tau_seconds, 10, []).'; % optional

% ---------------- Pooled outputs across iterations ----------------
% Store results for all k
all_records = struct();

% ----------------------- MAIN SWEEP -----------------------------------
for k = k_values
    fprintf('Running k = %.1f...\n', k);
    valid_iterations = 0;
    records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
    % ---------------- Parameter noise per chromosome ----------------
    % Each chromosome i gets its own kinetic parameters (chromosome-to-chromosome variability).
    while valid_iterations < N_iter
        %Varying parameters for each chromoosme
        Nbar   = zeros(Nchromosomes,1); % steady-state attachment number (initialization reference)
        n_dot  = zeros(Nchromosomes,1); % baseline gain term for attachments
        Lambda = zeros(Nchromosomes,1); % detachment rate constant
        Alpha  = zeros(Nchromosomes,1); % velocity-attachment coupling coefficient
        for i = 1:Nchromosomes
            Nbar(i)   = 20*(1+noise*(2*rand(1)-1));
            n_dot(i)  =  1*(1+noise*(2*rand(1)-1));
            Lambda(i) = n_dot(i)/Nbar(i);
            Alpha(i)  = n_dot(i)*6.2/(1-Beta);
        end

        % ---------------- Pre-allocate state variables ----------------
        cL = zeros(Nsteps+1, 2, Nchromosomes);    % left chromosome position: (time, [x y], chromosome index)
        cR = zeros(Nsteps+1, 2, Nchromosomes);    % right chromosome position: (time, [x y], chromosome index)
        xL = zeros(Nsteps+1, Nchromosomes);       % left MT tip x-position over time
        xR = zeros(Nsteps+1, Nchromosomes);       % right MT tip x-position over time
        NL = zeros(Nsteps+1, Nchromosomes);       % left KT–MT attachment number over time 
        NR = zeros(Nsteps+1, Nchromosomes);       % right KT–MT attachment number over time 
        vL = zeros(Nsteps,   Nchromosomes);       % left chromosome velocity over time 
        vR = zeros(Nsteps,   Nchromosomes);       % right chromosome velocity over time

        % ---------------- Initial conditions ----------------
        % Spring connectivity state (updated each time step)
        flag_c=zeros(Nchromosomes,Nchromosomes);   % current adjacency (1 = connected)
        flag_cp=flag_c;                            % Previous state 
        flag_c_sum_accumulated = zeros(Nsteps, 1); % To accumulate sum of flag_c over time
        % Fixed y positions are assigned
        y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
        for i = 1: Nchromosomes 
            % Initial attachment numbers (slightly perturbed to break symmetry)
            NR(1,i)=Nbar(i)*(1-epsilon);
            NL(1,i)=Nbar(i)*(1+epsilon); 
            % Initial chromosome positions (x in µm; y fixed) 
            cL(1,:,i) =-I0/2-epsilon*(2*rand(1)-1);%[x,y]
            cR(1,:,i) =I0/2+epsilon*(2*rand(1)-1);%[x,y]
            cL(:,2,i) =  y_positions(i);
            cR(:,2,i) =  y_positions(i);
            % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
            xL(1,i) = cL(1,1,i)-0.3;
            xR(1,i) = cR(1,1,i)+0.3;
            % Initial chromosome velocities
            vL(1,i)=0;
            vR(1,i)=0;
        end 


        % ---------------- Time integration loop (explicit Euler) ----------------
        for t = 1:Nsteps
            % ----- Update stochastic spring network  -----
            CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);   % CM x-position (N x 1)
            YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);   % CM fixed y-position (N x 1)
            % spring_connect returns:
            %   flag_c : NxN adjacency matrix (1 = connected spring, 0 = no spring)
            %   dist   : pairwise distance matrix (optional diagnostic)
            [flag_c, ~] = spring_connect(CM, YM, flag_cp, Nchromosomes, dt, l0, koff_0, kon_0);
            % Track total number of unique connections (upper triangle excludes double counting)
            flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
            % Update prior state
            flag_cp = flag_c;

            for i = 1:Nchromosomes
                % ----- Force calculations -----
                % Inter-chromosomal Coupling force 
                F_coupling_x = 0;
                for j = 1:Nchromosomes
                    if j ~= i
                        delta_x = CM(i) - CM(j);
                        F_coupling_x = F_coupling_x - k * delta_x * flag_c(i,j);
                    end
                end 
                 % Define centromere and MT forces on the given chromosome
                F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));    % force on left chromosome from KT–MT bundle
                F_CT_L =  Kct    * (cR(t,1,i) - cL(t,1,i) - I0);    % Centromere forces on Left Chromosome
                F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));    % force on right chromosome from KT–MT bundle
                F_CT_R = -Kct    * (cR(t,1,i) - cL(t,1,i) - I0);    % Centromere forces on Right Chromosome

                % ----- Velocities (overdamped dynamics) -----
                vL(t,i) = (F_KT_L + F_CT_L + F_coupling_x) / Gamma;
                vR(t,i) = (F_KT_R + F_CT_R + F_coupling_x) / Gamma;

                % ----- Update chromosome positions -----
                % Common-mode positional jitter: the same dx_diff is added to both sisters,
                dx_diff=noise_b *(randn(1))*sqrt(dt);
                cL(t+1,1,i) = cL(t,1,i) + dt*vL(t,i) + dx_diff;
                cR(t+1,1,i) = cR(t,1,i) + dt*vR(t,i) + dx_diff;
                % Maintain fixed Y positions
                cL(t+1,2,i) = cL(1,2,i);
                cR(t+1,2,i) = cR(1,2,i);
                % ----- Update MT tip positions -----
                xL(t+1,i) = xL(t,i) + dt*Beta*vL(t,i);
                xR(t+1,i) = xR(t,i) + dt*Beta*vR(t,i);
                % ----- Update attachment numbers -----
                NL(t+1,i) = NL(t,i) + dt * (n_dot(i) ...
                    + (-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
                    - Lambda(i)* NL(t,i));
                NR(t+1,i) = NR(t,i) + dt * (n_dot(i) ...
                    + ( Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
                    - Lambda(i)* NR(t,i));
            end
        end

        %%%%%%%%%%%%%%%%%% Analysis at the end of each iteration %%%%%%%%%%%%%%
        % Skip failed runs (NaNs), otherwise append results to pooled tables.
        % Reject runs with NaNs in trajectories
        if any(isnan(cL), 'all') || any(isnan(cR), 'all')
            fprintf('NaN detected for k = %.1f, skipping iteration\n', k);
            continue;
        end

       % ---- Displacement cross-correlation vs neighbor & tau ---- 
        sampling_interval = round((5/60) / dt);           % native steps per 5 s
        total_steps   = min(round(duration / dt), Nsteps);  % cap by duration & Nsteps    
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
        % --- Cross-correlation across tau = 5..150 s (in 5 s steps) ---
        iteration_id = valid_iterations + 1;
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
                    'tau',       tau, ...
                    'neighbor',  neighbor_idx, ...
                    'C_value',   C_norm, ...
                    'N_value',   N(neighbor_idx));
            end
        end

        valid_iterations = valid_iterations + 1;
        fprintf('Iteration %d complete for k = %.1f\n', iteration_id, k);
    end

    % Store with a clean per-k field name
    k_key = sprintf('k_%.1f', k);           % e.g., 'k_0.5'
    k_key = strrep(k_key, '.', '_');        % 'k_0_5'
    all_records.(k_key) = struct2table(records);
end

% -------- Plot: Avg Cross-Correlation vs Neighbor for each τ Group (overlay k values) --------
group_labels = ["\tau = 5–50 s", "\tau = 55–100 s", "\tau = 105–150 s"];
colors = lines(numel(k_values));

for g = 1:size(tau_groups_steps, 1)
    current_group = tau_groups_steps(g, :);   % resampled steps for this band
    figure('Color','w'); hold on; grid on; box on;

    for ki = 1:numel(k_values)
        k_val = k_values(ki);
        k_key = sprintf('k_%.1f', k_val);
        k_key = strrep(k_key, '.', '_');

        if ~isfield(all_records, k_key) || isempty(all_records.(k_key))
            continue;
        end
        records_table = all_records.(k_key);

        % subset to the current τ band
        subset = records_table(ismember(records_table.tau, current_group), :);
        if isempty(subset), continue; end

        % Aggregate: mean ± std across τ in the band and across iterations
        neighbors = 1:(Nchromosomes-1);
        meanC = nan(size(neighbors));
        stdC  = nan(size(neighbors));
        for n = neighbors
            vals = subset.C_value(subset.neighbor == n);
            meanC(n) = mean(vals, 'omitnan');
            stdC(n)  = std(vals,  'omitnan');
        end

        errorbar(neighbors, meanC, stdC, '-o', ...
            'Color', colors(ki,:), 'LineWidth', 2, ...
            'MarkerSize', 8, 'DisplayName', sprintf('k = %.1f', k_val));
    end

    xlim([0.5, Nchromosomes-0.5]);
    xticks(1:Nchromosomes-1);
    xlabel('Neighbor Level');
    ylabel('Avg Norm Cross-Corr');
    title(sprintf('C vs Neighbor Level (%s)', group_labels(g)));
    legend('Location', 'eastoutside');
    ylim([-0.1, 0.2]);
    set(gca, 'FontSize', 13);
end

% -------- Plot: C vs tau (5–105 s) for neighbors 1,2,3, for each k --------
% ONE figure per neighbor (1,2,3). Within each figure, each curve is a
% different k value. Tau is limited to 5..105 s (tau_steps 1..21).

tau_max_sec   = 105;
tau_max_steps = tau_max_sec / 5;     % 105 s -> 21 steps
tau_keep      = 1:tau_max_steps;     % tau_steps to include

neighbors_to_plot = 1:min(Nchromosomes-1, 3);

for n = neighbors_to_plot
    figure('Color','w'); hold on; grid on; box on;
    colors = lines(numel(k_values));

    for ki = 1:numel(k_values)
        k_val = k_values(ki);

        % Field name used in your k-sweep script:
        %   k_key = sprintf('k_%.1f', k_val); k_key = strrep(k_key,'.','_');
        k_key = sprintf('k_%.1f', k_val);
        k_key = strrep(k_key, '.', '_');  % e.g. k_0_5

        if ~isfield(all_records, k_key) || isempty(all_records.(k_key))
            continue
        end

        T = all_records.(k_key);

        % Subset: this neighbor + tau range
        subset = T(T.neighbor == n & ismember(T.tau, tau_keep), :);
        if isempty(subset), continue; end

        % Mean ± std across iterations at each tau
        grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");
        tau_sec = grouped.tau * 5;

        errorbar(tau_sec, grouped.mean_C_value, grouped.std_C_value, '-o', ...
            'Color', colors(ki,:), 'LineWidth', 1.8, ...
            'DisplayName', sprintf('k = %.1f', k_val));
    end

    xlabel('\tau (seconds)', 'FontSize', 14);
    ylabel('Avg normalized cross-correlation', 'FontSize', 14);
    title(sprintf('C vs \\tau (5–%d s), Neighbor %d', tau_max_sec, n), 'FontSize', 16);
    legend('Location','best');
    set(gca, 'FontSize', 12);
end

%% -------- Output combined CSV (long form) --------
% all_records fields: 'k_0_0', 'k_0_5', ..., each a table with:
% iteration, tau, neighbor, C_value, N_value
raw = table();
fn = fieldnames(all_records);

for f = 1:numel(fn)
    key = fn{f};                 % e.g., 'k_1_0' or 'k_0_5'
    T = all_records.(key);

    % Parse numeric k from field name
    k_str = erase(key, 'k_');    % '1_0' -> '1.0'
    k_str = strrep(k_str, '_', '.');
    k_val = str2double(k_str);

    T.k = repmat(k_val, height(T), 1);

    % Optional per-k file
    writetable(T, sprintf('raw_records_%s.csv', key));

    raw = [raw; T]; %#ok<AGROW>
end

% One combined long CSV for Python
writetable(raw, 'crosscorr_records_all_k.csv');
disp('Wrote effect_of_k_records_all.csv (and raw_records_k_*.csv)');
