clc; clear; close all;

% ------------------ Global Settings ------------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
Nchromosomes = 5;           % Number of chromosome pairs
Alpha_scales = 3:1:15;      % Alpha scale factors
Kct_values   = 10:0.5:20;   % Kct values to sweep

% Noise / cell-to-cell variability
N_iter  = 50;              % Number of cells per (alpha, Kct)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (position + force)
epsilon = 0.1;              % small perturbation for initial conditions

% Fixed (per-run) constants
Beta = 0.7;     % Scaling factor
I0   = 2;       % µm, rest length centromere spring
Kkt  = 1;       % pN/µm, kinetochore spring constant
Gamma= 0.1;     % kg/s, drag
Nmax = 25;      % Max attachments

% Inter-chromosomal spring parameters
koff_0=20; 
kon_0=60; 
l0 = 2; 
k = 0.5; 
Nspring = 0.5; 

% Progress counter
nTotal  = numel(Alpha_scales) * numel(Kct_values) * N_iter;
counter = 0;

% Output table: per-iteration summary for neighbor 1, τ=5–50 s
results_all = table( ...
    'Size',[0 7], ...
    'VariableTypes', {'double','double','double','double','double','double','double'}, ...
    'VariableNames', {'alpha_scale','alpha_effective','Kct','iteration','mean_r','std_r','Npairs'});

% ------------------ Parameter Sweep with Noise ------------------
for a = Alpha_scales
    for Kct = Kct_values
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
                Alpha(i)=n_dot(i)*a/(1-Beta);    %% FIX: scale by 'a'
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
            if any(isnan(cL), 'all') || any(isnan(cR), 'all')
                skipped = true;
            else
                % ---- Cross-correlation: neighbor=1, τ in first group (5–50 s) ----
                CMx  = squeeze((cL(1:Nsteps,1,:) + cR(1:Nsteps,1,:)) / 2);   % [Nsteps x Nchromosomes]
                taus = tau_groups(1,:);                                      % first τ-group
                cvals = nan(1, numel(taus));
                for k = 1:numel(taus)
                    tau = taus(k);
                    [C, N] = crosscorrelation_measurement(CMx, tau, Nsteps, Nchromosomes);
                    if N(1) > 0
                        cvals(k) = C(1) / N(1);
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
                        a, ...
                        alpha_effective, ...
                        Kct, ...
                        double(iter), ...
                        mean_C, ...
                        std_C, ...
                        'VariableNames', {'alpha_scale','alpha_effective','Kct','iteration','mean_Cij','std_Cij'})];
                end
            end
            % ---- progress print AFTER each iteration ----
            counter = counter + 1;
            if skipped
                tag = ' [skipped: NaN]';
            else
                tag = '';
            end
            fprintf('Progress: %d/%d (%.1f%%) | a=%.1f, Kct=%.2f, iter=%d\n', ...
                counter, nTotal, 100*counter/nTotal, a, Kct, iter);
        end % iter
    end
end

%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_all, 'alpha3-10_noise_results_all.csv');   % CSV for Python



%% === Mean Cij phase map (heatmap + Cij=0.05 contour) ===
% Round keys a bit to avoid floating-point mismatches when selecting subsets
alpha_key = round(results_all.alpha_scale, 6);
Kct_key   = round(results_all.Kct,         6);

unique_alpha = unique(alpha_key);     % x-axis
unique_Kct   = unique(Kct_key);       % y-axis

% Preallocate (rows = Kct, cols = alpha)
Z = nan(numel(unique_Kct), numel(unique_alpha));

% Fill grid with mean over iterations (omit NaNs)
for ia = 1:numel(unique_alpha)
    for ik = 1:numel(unique_Kct)
        mask = (alpha_key == unique_alpha(ia)) & (Kct_key == unique_Kct(ik));
        if any(mask)
            Z(ik, ia) = mean(results_all.mean_Cij(mask), 'omitnan');
        end
    end
end

% Plot
figure('Color','w');
imagesc(unique_alpha, unique_Kct, Z);
set(gca, 'YDir','normal');              % so low Kct at bottom, high at top
axis tight;
colormap(parula);
cb = colorbar;
ylabel(cb, 'Mean C_{ij} (neighbor=1, \tau=5–50 s)', 'Interpreter','tex');

xlabel('$\alpha$ scale','Interpreter','latex','FontSize',16);
ylabel('$K_{ct}$ (pN/$\mu$m)','Interpreter','latex','FontSize',16);
title('Mean C_{ij}','Interpreter','tex','FontSize',16);
hold on;

% Overlay contour at Cij = 0.05
% meshgrid returns size [numel(Kct) x numel(alpha)], matching Z
[AA, KK] = meshgrid(unique_alpha, unique_Kct);
C_level  = 0.05;
contour(AA, KK, Z, [C_level C_level], 'k', 'LineWidth', 2);  % black contour at 0.05

% Optional cosmetics (uncomment if you like)
% caxis([-0.2 0.6]);                         % fix color range if helpful
% grid on; box on;                            % add axes box/grid
% set(gca,'TickDir','out','LineWidth',1.2);   % styling




