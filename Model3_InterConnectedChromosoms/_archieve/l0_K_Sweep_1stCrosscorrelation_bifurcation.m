clc; clear; close all;

% ------------------ Global Settings ------------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
Nchromosomes = 5;           % Number of chromosome pairs
l0_values = 1:0.5:5;      % Alpha scale factors
k_values = 0:0.1:2;   % Kct values to sweep

% Noise / cell-to-cell variability
N_iter  = 50;              % Number of cells per (alpha, Kct)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;              % Brownian noise strength (position + force)
epsilon = 0.1;              % small perturbation for initial conditions

% Fixed (per-run) constants
Kct =12.3; 
Beta = 0.7;     % Scaling factor
I0   = 2;       % µm, rest length centromere spring
Kkt  = 1;       % pN/µm, kinetochore spring constant
Gamma= 0.1;     % kg/s, drag
Nmax = 25;      % Max attachments

% Inter-chromosomal spring parameters
koff_0=20; 
kon_0=40; 
%l0 = 2; 
%k = 0.5; 
Nspring = 1; 

duration = 10; % (min) simulation duration to analyze
% Cross-correlation parameters
% Tau setup with 5 s resampling of chromosome movement 
tau_steps   = 1:30;                 % 5 s, 10 s, ..., 150 s (in resampled steps)
tau_seconds = 5 * tau_steps;        % for labels only
% Three bands of 10 lags each: [5–50 s], [55–100 s], [105–150 s]
tau_groups_steps = reshape(tau_steps, 10, []).';   % 3x10 in step units
tau_groups_sec   = reshape(tau_seconds, 10, []).'; % 3x10 in seconds (optional)

% Progress counter
nTotal  = numel(l0_values) * numel(k_values) * N_iter;
counter = 0;

% Output table: per-iteration summary for neighbor 1, τ=5–50 s
results_all = table( ...
    'Size',[0 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'l0','k','iteration','mean_Cij','std_Cij'});
% Cij sanity thresholds
Cij_min_abs = 1e-4;   % "very very small" -> effectively zero / degenerate
Cij_max_abs = 10;    % "very very large" -> outside correlation bounds (allow tiny overshoot)

% ------------------ Parameter Sweep with Noise ------------------
for l0 = l0_values
    for k = k_values
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
            alpha_effective = mean(Alpha); 

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
            
            % default so fprintf never sees an undefined variable
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
                CMx_sampled  = ((cL_sampled + cR_sampled) / 2)/sqrt(10);        % [Nsamp x 1 x Nchromosomes]
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
                    std_C  = std( cvals, 'omitnan');

                    too_small = ~isnan(mean_C) && (abs(mean_C) < Cij_min_abs);
                    too_large =  (abs(mean_C) > Cij_max_abs);

                    if  too_small || too_large
                        skipped = true;
                    else
                        % Keep this iteration
                        results_all = [results_all; table( ...
                            l0, k, double(iter), mean_C, std_C, ...
                            'VariableNames', {'l0','k','iteration','mean_Cij','std_Cij'})];
                    end
                end
            end
            
            % ---- progress print AFTER each iteration ----
            counter = counter + 1;
            if skipped
                fprintf('Progress: %d/%d (%.1f%%) | l0=%.1f, k=%.2f, iter=%d | SKIPPED\n', ...
                    counter, nTotal, 100*counter/nTotal, l0, k, iter);
            else
                fprintf('Progress: %d/%d (%.1f%%) | l0=%.1f, k=%.2f, iter=%d | mean C_{ij}=%.4f\n', ...
                    counter, nTotal, 100*counter/nTotal, l0, k, iter, mean_C);
            end

        end % iter
    end
end

%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_all, 'l0_k_pearsonR_results.csv');   % CSV for Python



%% === Mean C_ij phase map (heatmap + contour) ===


% 1) Pull columns and round keys to avoid float-equality issues
l0_key = round(results_all.l0, 6);
k_key  = round(results_all.k,  6);
mcij   = results_all.mean_Cij;

% 2) Unique axis values (sorted)
unique_l0 = unique(l0_key(:));   % x-axis
unique_k  = unique(k_key(:));    % y-axis

% 3) Pivot into Z(k, l0)
Z = nan(numel(unique_k), numel(unique_l0));
for ia = 1:numel(unique_l0)
    for ik = 1:numel(unique_k)
        mask = (l0_key == unique_l0(ia)) & (k_key == unique_k(ik));
        if any(mask)
            Z(ik, ia) = mean(mcij(mask), 'omitnan');
        end
    end
end

% 4) Plot
figure('Color','w');
imagesc(unique_l0, unique_k, Z);      % imagesc(x, y, C)
set(gca,'YDir','normal');             % low k at bottom
axis tight
colormap(parula);
cb = colorbar;
ylabel(cb, 'Mean C_{ij}', 'Interpreter','tex');

xlabel('$l_0$','Interpreter','latex','FontSize',16);
ylabel('$k$ (pN/$\mu$m)','Interpreter','latex','FontSize',16);
title('Mean C_{ij} heatmap','Interpreter','tex','FontSize',16);

set(gca,'XTick', unique_l0, 'YTick', unique_k);

% 5) Optional: draw a contour at a target C_ij level
hold on
[L0, KK] = meshgrid(unique_l0, unique_k);
C_level = 0.05;                       % threshold example
contour(L0, KK, Z, [C_level C_level], 'k', 'LineWidth', 2);
hold off




