clc; clear; close all;

% ------------------ Global Settings ------------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
Nchromosomes = 5;           % Number of chromosome pairs
koff0_list  = 0:5:100;      % base disconnection rate(s)
kon0_list   = 0:5:100;      % base connection rate(s)

% Noise / cell-to-cell variability
N_iter  = 50;               % Number of cells per (kon_0, koff_0)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (position + force)
epsilon = 0.1;              % small perturbation for initial conditions

% Fixed (per-run) constants
Kct = 12.3; 
Beta = 0.7;                 % Scaling factor
I0   = 2;                   % µm, rest length centromere spring
Kkt  = 1;                   % pN/µm, kinetochore spring constant
Gamma= 0.1;                 % kg/s, drag
Nmax = 25;                  % Max attachments

% Inter-chromosomal spring parameters
l0 = 3; 
k  = 0.5; 
Nspring = 1; 

duration = 10;              % (min) simulation duration to analyze

% Tau setup with 5 s resampling of chromosome movement 
tau_steps   = 1:30;                         % 5 s, 10 s, ..., 150 s (in resampled steps)
tau_seconds = 5 * tau_steps;                % for labels only
tau_groups_steps = reshape(tau_steps, 10, []).';   % 3x10 in step units
tau_groups_sec   = reshape(tau_seconds, 10, []).'; % 3x10 in seconds (optional)

% Progress counter
nTotal  = numel(kon0_list) * numel(koff0_list) * N_iter;
counter = 0;

% Output table: per-iteration summary for neighbor 1, τ=5–50 s
results_all = table( ...
    'Size',[0 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'kon_0','koff_0','iteration','mean_Cij','std_Cij'});

% Cij sanity thresholds
Cij_min_abs = 1e-4;   % "very very small" -> effectively zero / degenerate
Cij_max_abs = 10;     % "very very large" -> outside correlation bounds (allow tiny overshoot)

% ------------------ Parameter Sweep with Noise ------------------
for iKon = 1:numel(kon0_list)
    for iKoff = 1:numel(koff0_list)
        kon_0  = kon0_list(iKon);
        koff_0 = koff0_list(iKoff);

        for iter = 1:N_iter
            % Initialize chromosome parameters
            Nbar = zeros(Nchromosomes,1);
            for i = 1:Nchromosomes
                Nbar(i) = 20*(1+noise*(2*rand(1)-1));
            end
            n_dot = zeros(Nchromosomes,1);
            for i = 1:Nchromosomes
                n_dot(i) = 1*(1+noise*(2*rand(1)-1));
            end
            Lambda = zeros(Nchromosomes, 1);
            for i = 1:Nchromosomes
                Lambda(i)=n_dot(i)/Nbar(i);
            end
            Alpha = zeros(Nchromosomes, 1);
            for i = 1:Nchromosomes
                Alpha(i)=n_dot(i)*6.2/(1-Beta);    
            end

            % Set up the vectors for results
            cL = zeros(Nsteps+1, 2, Nchromosomes);
            cR = zeros(Nsteps+1, 2, Nchromosomes);
            xL = zeros(Nsteps+1, Nchromosomes);
            xR = zeros(Nsteps+1, Nchromosomes);
            NL = zeros(Nsteps+1, Nchromosomes);
            NR = zeros(Nsteps+1, Nchromosomes);
            vL = zeros(Nsteps, Nchromosomes);
            vR = zeros(Nsteps, Nchromosomes);

            % Coupling state
            flag_c_sum_accumulated = zeros(Nsteps, 1); 
            flag_c = zeros(Nchromosomes);
            flag_cp = flag_c;

            % Initial Condition 
            y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
            for i = 1:Nchromosomes
                NR(1,i) = Nbar(i)*(1 - epsilon);
                NL(1,i) = Nbar(i)*(1 + epsilon);
                cL(1,:,i) = -I0/2 - epsilon*(2*rand(1)-1); % [x,y]
                cR(1,:,i) =  I0/2 + epsilon*(2*rand(1)-1); % [x,y]
                cL(:,2,i) =  y_positions(i);
                cR(:,2,i) =  y_positions(i);
                xL(1,i) = cL(1,1,i) - 0.3;
                xR(1,i) = cR(1,1,i) + 0.3;
                vL(1,i)=0; vR(1,i)=0;
            end

            % ODE solver loop
            for t = 1:Nsteps
                CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);
                YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);
                [flag_c, ~] = spring_connect(CM, YM, flag_cp, Nchromosomes, dt, l0, koff_0, kon_0);
                flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
                flag_cp = flag_c;

                for i = 1:Nchromosomes
                    % Coupling force 
                    F_coupling_x = 0;
                    for j = 1:Nchromosomes
                        if j ~= i
                            delta_x = CM(i) - CM(j);
                            F_coupling_x = F_coupling_x - k * delta_x * flag_c(i,j);
                        end
                    end
                    % Forces
                    F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));
                    F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0);
                    F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));
                    F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0);
                    % Velocities 
                    vL(t,i) = (F_KT_L + F_CT_L + F_coupling_x) / Gamma;
                    vR(t,i) = (F_KT_R + F_CT_R + F_coupling_x) / Gamma;
                    % Positions 
                    dx_diff = noise_b * randn(1) * sqrt(dt);
                    cL(t+1,1,i) = cL(t,1,i) + dt*vL(t,i) + dx_diff;
                    cR(t+1,1,i) = cR(t,1,i) + dt*vR(t,i) + dx_diff;
                    cL(t+1,2,i) = cL(1,2,i);  % fixed y
                    cR(t+1,2,i) = cR(1,2,i);  % fixed y
                    % MT tips
                    xL(t+1,i) = xL(t,i) + dt*Beta*vL(t,i);
                    xR(t+1,i) = xR(t,i) + dt*Beta*vR(t,i);
                    % Attachments 
                    NL(t+1,i) = NL(t,i) + dt * (n_dot(i) ...
                        + (-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
                        - Lambda(i)* NL(t,i));
                    NR(t+1,i) = NR(t,i) + dt * (n_dot(i) ...
                        + ( Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
                        - Lambda(i)* NR(t,i));
                end
            end

            % ---- NaN guard & per-iteration progress reporting ----
            skipped = false;
            mean_C = NaN; 
            std_C  = NaN; 

            if any(isnan(cL), 'all') || any(isnan(cR), 'all')
                skipped = true;
            else
                % ===== Resample to every 5 seconds on the native arrays =====
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
                    % ---- Cross-correlation over the FIRST τ group: 5–50 s ----
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
                        % Keep this iteration
                        results_all = [results_all; table( ...
                            kon_0, koff_0, double(iter), mean_C, std_C, ...
                            'VariableNames', {'kon_0','koff_0','iteration','mean_Cij','std_Cij'})];
                    end
                end
            end

            % ---- progress print AFTER each iteration ----
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

%% Export the data for python plots 

writetable(results_all, 'kon_koff_meanCij_results.csv');   % CSV for Python

%% === Mean Cij phase map (heatmap + Cij=0.05 contour) ===

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

% Overlay contour at Cij = 0.05
[KKON, KKOFF] = meshgrid(unique_kon, unique_koff);
C_level  = 0.05;
contour(KKON, KKOFF, Z, [C_level C_level], 'k', 'LineWidth', 2);
hold off;
