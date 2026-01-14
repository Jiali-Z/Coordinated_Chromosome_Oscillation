%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kon_Koff_Sweep_CrossCorrelation.m
%
% PURPOSE
%   Run a stochastic (multi-iteration; noisy) parameter sweep over the
%   inter-chromosomal spring network *kinetics*:
%       (1) kon_0  -> base connection rate (formation) for inter-chromosomal links
%       (2) koff_0 -> base disconnection rate (breakage) for inter-chromosomal links
%   while simulating N sister-chromosome pairs that oscillate in 1D (x) using
%   the same force-balance oscillation dynamics as the single-chromosome model,
%   with stochastic inter-chromosomal coupling updated by spring_connect(...).
%
%   For each (kon_0, koff_0) pair, the script simulates N_iter independent
%   “cells” (noise realizations). For each valid cell (no NaNs / non-degenerate
%   values), chromosome center-of-mass motion is resampled to 5-second
%   resolution and displacement cross-correlation is quantified using
%   crosscorrelation_measurement(...). The primary metric reported is the
%   neighbor-order-1 displacement cross-correlation averaged over short time
%   lags (τ = 5–50 s; i.e., the first τ-group).
%
% NOTE / PAPER NOTATION
%   The paper uses the same kon_0 and koff_0 notation 
% 
%
% HOW THIS SCRIPT IS USED IN THE PAPER
%   - Figure S4A: plot the C_{ij} contour overlays on the
%     connectivity heatmap, enabling identification of (kon_0, koff_0)
%     parameter regimes that match the experimentally measured
%     neighbor-order-1 displacement cross-correlation (τ = 5–50 s)
%
% NOTES / IMPLEMENTATION DETAILS
%   To facilitate reproducibility without rerunning the full sweep, the
%   processed CSV output has been COMMITTED to the repository:Users can 
%   directly load this CSV file to reproduce Figure 4B without re-running
%   this script.
%
%
% OUTPUTS
%   1) results_all table (one row per valid iteration):
%        kon_0, koff_0, iteration, mean_Cij, std_Cij
%   2) CSV export for downstream plotting / panel assembly:
%        kon_koff_meanCij_results.csv
%   3) Phase map grid (Z) constructed from results_all for visualization
%      (rows = koff_0, cols = kon_0) with contours (e.g., Cij=0.01,0.05,0.1).
%
% DEPENDENCIES
%   spring_connect.m
%   crosscorrelation_measurement.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC–Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% ------------------ Sweep definitions ------------------
koff0_list  = 0:5:100;      % base disconnection rate(s)
kon0_list   = 0:5:100;      % base connection rate(s)
nTotal  = numel(kon0_list) * numel(koff0_list) * N_iter;
counter = 0;

% ---------------- Simulation settings ----------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps

% ---------------- Number of chromosome ---------------- 
Nchromosomes = 5;           % Number of chromosome pairs


% ------------------ Stochasticity parameters ------------------
N_iter  = 50;               % Number of cells per (kon_0, koff_0)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (position + force)

% ------------------ Fixed model constants ------------------
Kct = 12.3;                 % pN/µm: centromere spring constant
I0   = 2;                   % µm, rest length centromere spring
Kkt  = 1;                   % pN/µm, kinetochore spring constant
Gamma= 0.1;                 % effective drag coefficient
Beta = 0.7;                 % MT tip motion coupling factor (x tip follows chromosome)
Nmax = 25;                  % max attachments per kinetochore
epsilon = 0.1;              % small perturbation for initial conditions


% ---------------- Inter-chromosome coupling parameters ----------------
l0 = 3;                     % µm     inter-chromosomal spring effective length 
k  = 0.5;                   % pN/µm: inter-chromosomal spring constant 
Nspring = 1;                 


% ------------------ analysis parameter ------------------
duration = 10;                              % (min) simulation duration to analyze
% Tau setup with 5 s resampling of chromosome movement 
tau_steps   = 1:30;                         % 5 s, 10 s, ..., 150 s (in resampled steps)
tau_seconds = 5 * tau_steps;                % for labels only
tau_groups_steps = reshape(tau_steps, 10, []).';   % 3x10 in step units
tau_groups_sec   = reshape(tau_seconds, 10, []).'; % 3x10 in seconds (optional)
% Cij sanity thresholds
Cij_min_abs = 1e-4;   % "very very small" -> effectively zero / degenerate
Cij_max_abs = 10;     % "very very large" -> outside correlation bounds (allow tiny overshoot)


% ------------------ Storage for results ------------------
results_all = table( ...
    'Size',[0 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'kon_0','koff_0','iteration','mean_Cij','std_Cij'});

% ========================================================================
%                         PARAMETER SWEEP LOOP
% =========================================================================
for iKon = 1:numel(kon0_list)
    for iKoff = 1:numel(koff0_list)
        kon_0  = kon0_list(iKon);
        koff_0 = koff0_list(iKoff);
        for iter = 1:N_iter

            % ---- cell-specific parameters with noise ----
            Nbar   = zeros(Nchromosomes, 1);  % steady-state attachment number (initialization reference)
            n_dot  = zeros(Nchromosomes, 1);  % baseline gain term for attachments
            Lambda = zeros(Nchromosomes, 1);  % detachment rate constant
            Alpha  = zeros(Nchromosomes, 1);  % velocity-attachment coupling coefficient
            for i = 1:Nchromosomes
                Nbar(i)  = 20 * (1 + noise * (2*rand(1) - 1));
                n_dot(i) =  1 * (1 + noise * (2*rand(1) - 1));
                Lambda(i) = n_dot(i) / Nbar(i);
                Alpha(i)  = n_dot(i) * 6.2 / (1 - Beta); alpha_effective = mean(Alpha); 
            end 

            % ------------------ Allocate state variables ------------------
            cL = zeros(Nsteps+1, 2, Nchromosomes);       % left chromosome position: (time, [x y], chromosome index)
            cR = zeros(Nsteps+1, 2, Nchromosomes);       % right chromosome position: (time, [x y], chromosome index)
            xL = zeros(Nsteps+1, Nchromosomes);          % left MT tip x-position over time
            xR = zeros(Nsteps+1, Nchromosomes);          % right MT tip x-position over time
            NL = zeros(Nsteps+1, Nchromosomes);          % left KT–MT attachment number over time 
            NR = zeros(Nsteps+1, Nchromosomes);          % right KT–MT attachment number over time 
            vL = zeros(Nsteps, Nchromosomes);            % left chromosome velocity over time 
            vR = zeros(Nsteps, Nchromosomes);            % right chromosome velocity over time

            % ---------------- Initial conditions ----------------
            % Spring connectivity state (updated each time step)
            flag_c = zeros(Nchromosomes);                   % current adjacency (1 = connected)
            flag_cp = flag_c;                               % Previous state 
            flag_c_sum_accumulated = zeros(Nsteps, 1);      % To accumulate sum of flag_c over time
            % Fixed y positions are assigned 
            y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
            % Each Pairs have their own 
            for i = 1:Nchromosomes
                % Initial attachment numbers (slightly perturbed to break symmetry)
                NR(1,i) = Nbar(i)*(1 - epsilon);
                NL(1,i) = Nbar(i)*(1 + epsilon);
                % Initial chromosome positions (x in µm; y fixed) 
                cL(1,:,i) = -I0/2 - epsilon*(2*rand(1)-1);   % left sister initial x-position
                cR(1,:,i) =  I0/2 + epsilon*(2*rand(1)-1);   % right sister initial x-position
                cL(:,2,i) =  y_positions(i);                 % y fixed
                cR(:,2,i) =  y_positions(i);                 % y fixed
                % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
                xL(1,i) = cL(1,1,i) - 0.3;
                xR(1,i) = cR(1,1,i) + 0.3;
                % Initial chromosome velocities
                vL(1,i)=0; vR(1,i)=0;
            end

            % ---------------- Time integration loop (explicit Euler) ----------------
            for t = 1:Nsteps
                % ----- Update stochastic spring network  -----
                CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);
                YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);
                [flag_c, ~] = spring_connect(CM, YM, flag_cp, Nchromosomes, dt, l0, koff_0, kon_0);
                flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
                flag_cp = flag_c;
               % Per-chromosome forces & updates
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
                    F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));     % force on left chromosome from KT–MT bundle
                    F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0);         % Centromere forces on Left Chromosome
                    F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));     % force on right chromosome from KT–MT bundle
                    F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0);        % Centromere forces on Right Chromosome
                     % ----- Velocities (overdamped dynamics) -----
                    vL(t,i) = (F_KT_L + F_CT_L + F_coupling_x) / Gamma;
                    vR(t,i) = (F_KT_R + F_CT_R + F_coupling_x) / Gamma;
                    % ----- Update chromosome positions -----
                    % Common-mode positional jitter: the same dx_diff is added to both sisters,
                    dx_diff = noise_b * randn(1) * sqrt(dt);
                    cL(t+1,1,i) = cL(t,1,i) + dt*vL(t,i) + dx_diff;
                    cR(t+1,1,i) = cR(t,1,i) + dt*vR(t,i) + dx_diff;
                    cL(t+1,2,i) = cL(1,2,i);  % fixed y
                    cR(t+1,2,i) = cR(1,2,i);  % fixed y
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
            % ---------------Displacement cross-correlation vs neighbor & tau-------- 
            skipped = false;
            mean_C = NaN; 
            std_C  = NaN; 

            if any(isnan(cL), 'all') || any(isnan(cR), 'all')
                skipped = true;
            else
                % Resample to every 5 seconds on the native arrays 
                sampling_interval = round((5/60) / dt);          % native steps per 5 s
                total_steps_nat   = min(round(duration / dt), Nsteps); % cap by duration & Nsteps    
                idx_vec = 1:sampling_interval:total_steps_nat;
                % Sample L/R x positions and compute CMx (x only)
                cL_sampled  = cL(idx_vec, 1, :);                        % [Nsamp x 1 x Nchromosomes]
                cR_sampled  = cR(idx_vec, 1, :);
                % tidy scaling equivalent to ((cL+cR)/2)/sqrt(10):
                S = 1/sqrt(40);
                CMx_sampled = (cL_sampled + cR_sampled) * S;            % [Nsamp x 1 x Nchromosomes]
                CMx         = reshape(CMx_sampled, [], Nchromosomes);   % [Nsamp x Nchromosomes]
                Nsteps_res  = size(CMx, 1);                              % resampled length
                % Safety: enough points for max τ
                if Nsteps_res <= max(tau_steps)
                    warning('Not enough resampled time points (%d) for max tau (%d). Skipping iteration.', ...
                            Nsteps_res, max(tau_steps));
                    skipped = true;
                else
                    % Cross-correlation over the FIRST τ group: 5–50 s ----
                    taus   = tau_groups_steps(1, :);     % [1..10] resampled steps (5..50 s)
                    cvals  = nan(1, numel(taus));
                    for ti = 1:numel(taus)
                        tau = taus(ti);                  % resampled steps
                        [C, N] = crosscorrelation_measurement(CMx, tau, Nsteps_res, Nchromosomes);
                        if N(1) > 0
                            cvals(ti) = C(1) / N(1);     % neighbor level 1
                        end
                    end

                    mean_C = mean(cvals, 'omitnan');
                    std_C  = std(cvals,  'omitnan');
                    too_small = ~isnan(mean_C) && (abs(mean_C) < Cij_min_abs);
                    too_large = (abs(mean_C) > Cij_max_abs);
                    if  too_small || too_large
                        skipped = true;
                    else
                       % Append the results 
                        results_all = [results_all; table( ...
                            kon_0, koff_0, double(iter), mean_C, std_C, ...
                            'VariableNames', {'kon_0','koff_0','iteration','mean_Cij','std_Cij'})];
                    end
                end
            end

           % --- --------------------Progress note ---------------
            counter = counter + 1;
            if skipped
                fprintf('Progress: %d/%d (%.1f%%) | kon_0=%.1f, koff_0=%.1f, iter=%d | SKIPPED\n', ...
                    counter, nTotal, 100*counter/nTotal, kon_0, koff_0, iter);
            else
                fprintf('Progress: %d/%d (%.1f%%) | kon_0=%.1f, koff_0=%.1f, iter=%d | mean C_{ij}=%.4f\n', ...
                    counter, nTotal, 100*counter/nTotal, kon_0, koff_0, iter, mean_C);
            end
        end % iter
    end
end

%% ========================================================================
%   FIGURE S4A: HEATMAP  
% =========================================================================

% --- Load sweep results from CSV (exported earlier) ---
csvFile = "Output_Sweep_Bifurcation/kon_koff_meanCij_results.csv";
results_all = readtable(csvFile);


kon_key  = round(results_all.kon_0,  6);
koff_key = round(results_all.koff_0, 6);

unique_kon  = unique(kon_key);     % x-axis
unique_koff = unique(koff_key);    % y-axis

% Preallocate (rows = koff, cols = kon)
Z = nan(numel(unique_koff), numel(unique_kon));

% Fill grid with mean over iterations (omit NaNs)
for ia = 1:numel(unique_kon)
    for ik = 1:numel(unique_koff)
        mask = (kon_key == unique_kon(ia)) & (koff_key == unique_koff(ik));
        if any(mask)
            Z(ik, ia) = mean(results_all.mean_Cij(mask), 'omitnan');
        end
    end
end

% Plot
figure('Color','w');
imagesc(unique_kon, unique_koff, Z);
set(gca, 'YDir','normal');              % low koff at bottom
axis tight;
colormap(parula);
cb = colorbar;
ylabel(cb, 'Mean C_{ij} (neighbor=1, \tau=5–50 s)', 'Interpreter','tex');

xlabel('$k_{\mathrm{on},0}$','Interpreter','latex','FontSize',16);
ylabel('$k_{\mathrm{off},0}$','Interpreter','latex','FontSize',16);
title('Mean C_{ij}','Interpreter','tex','FontSize',16);
hold on;

% Overlay contour lines at selected Cij levels
[KKON, KKOFF] = meshgrid(unique_kon, unique_koff);

hold on;
[C01, h01] = contour(KKON, KKOFF, Z, [0.01 0.01], '--', ...
    'LineWidth', 2, 'Color',[0.6 0.6 0.6]);

[C05, h05] = contour(KKON, KKOFF, Z, [0.05 0.05], '--k', ...
    'LineWidth', 2);

[C10, h10] = contour(KKON, KKOFF, Z, [0.1 0.1], '--', ...
    'LineWidth', 2, 'Color',[0.85 0.2 0.2]);

% ---- Label each contour line ----
clabel(C01, h01, 'FontSize',12, 'Color',[0.6 0.6 0.6], ...
    'Interpreter','tex', 'LabelSpacing',400);

clabel(C05, h05, 'FontSize',12, 'Color','k', ...
    'Interpreter','tex', 'LabelSpacing',400);

clabel(C10, h10, 'FontSize',12, 'Color',[0.85 0.2 0.2], ...
    'Interpreter','tex', 'LabelSpacing',400);

hold off;


%% Export the data for python plots 

writetable(results_all, 'kon_koff_meanCij_results.csv');   % CSV for Python