clc; clear; close all;

% ------------------ Global Settings ------------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
Beta      = 0.7;
duration  = 10;             % Duration (min) used in KE analysis
Alpha_scales = 3:0.2:15;    % Alpha scale factors
Kct_values   = 10:0.1:20;   % Kct values to sweep

% Noise / cell-to-cell variability
N_iter  = 50;              % Number of cells per (alpha, Kct)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (position + force)

% Fixed (per-run) constants
I0   = 2;       % µm, rest length centromere spring
Kkt  = 1;       % pN/µm, kinetochore spring constant
Gamma= 0.1;     % kg/s, drag
Nmax = 25;      % Max attachments

results_all = table( ...
    'Size',[0 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'alpha_scale','alpha','Kct','iteration','mean_KE'});
nTotal  = numel(Alpha_scales) * numel(Kct_values) * N_iter;
counter = 0;

% ------------------ Parameter Sweep with Noise ------------------
for a = Alpha_scales
    for Kct = Kct_values
        for iter = 1:N_iter
            % ---- cell-specific parameters with noise ----
            n_dot  = 1 * (1 + noise*(2*rand - 1));  % production rate
            Nbar   = 20   * (1 + noise*(2*rand - 1));  % steady-state MTs
            Lambda = n_dot / Nbar;                             % detach rate
            Alpha  = n_dot * a / (1 - Beta);                   % scaled Alpha
            epsilon = 0.1;                                     % small perturb.

            % ---- allocate & initialize ----
            % Note: Nsteps+1 because we use t and t+1
            xL = zeros(Nsteps+1,1); xR = zeros(Nsteps+1,1);
            NL = zeros(Nsteps+1,1); NR = zeros(Nsteps+1,1);
            cL = zeros(Nsteps+1,2,1); cR = zeros(Nsteps+1,2,1);
            vL = zeros(Nsteps+1,1);  vR = zeros(Nsteps+1,1);

            NR(1)     = Nbar*(1 - epsilon);
            NL(1)     = Nbar*(1 + epsilon);
            cL(1,1,1) =-I0/2-epsilon;
            cR(1,1,1) =I0/2+epsilon;
            cL(:,2,1) = 0;
            cR(:,2,1) = 0;

            % MT tips start slightly offset from centromeres
            xL(1) = cL(1,1,1) - 0.3;
            xR(1) = cR(1,1,1) + 0.3;
            vL(1) = 0; vR(1) = 0;

            % ---- time stepping ----
            for t = 1:Nsteps
                % Forces on L and R centromeres
                F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));
                F_CT_L =  Kct   * (cR(t,1,1) - cL(t,1,1) - I0);
                F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));
                F_CT_R = -Kct   * (cR(t,1,1) - cL(t,1,1) - I0);

                % Velocities
                vL(t) = (F_KT_L + F_CT_L ) / Gamma;
                vR(t) = (F_KT_R + F_CT_R ) / Gamma;

                % Brownian step (shared to keep the pair coherent; change if desired)
                dx_diff = noise_b * randn(1,1) * sqrt(dt);
                % Update centromeres
                cL(t+1,1,1) = cL(t,1,1) + dt*vL(t) + dx_diff;
                cR(t+1,1,1) = cR(t,1,1) + dt*vR(t) + dx_diff;
                % Update MT tips
                xL(t+1) = xL(t) + dt * Beta * vL(t);
                xR(t+1) = xR(t) + dt * Beta * vR(t);
                % Attachments
                NL(t+1) = NL(t) + dt * ( ...
                    n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
                    - Lambda * NL(t));
                NR(t+1) = NR(t) + dt * ( ...
                    n_dot + ( Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
                    - Lambda * NR(t));
            end

            % ---- Activity measurement (3-pt slope, sampled every 5 s) ----
            activity_table = chromosome_activity_measurement(cL, cR, dt);
            % Aggregate over noisy cells for this (alpha, Kct)
            results_all = [results_all; ...
                table(a, Alpha, Kct, double(iter), activity_table.Avg_KE, ...
                'VariableNames', {'alpha_scale','alpha','Kct','iteration','mean_KE'})];
            % --- Progress note ---
            counter = counter + 1;
            fprintf('Progress: %d/%d (%.1f%%) | a=%.1f, Kct=%.2f, iter=%d\n', ...
                counter, nTotal, 100*counter/nTotal, a, Kct, iter);
        end
    end
end

%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_all, 'MSV_noise_alpha_3_0.2_15.csv');   % CSV for Python



%% Mean kinetic energy phase map
unique_alpha = unique(results_all.alpha_scale);
unique_Kct   = unique(results_all.Kct);
mean_KE_map = zeros(length(unique_alpha), length(unique_Kct));

for iA = 1:length(unique_alpha)
    for iK = 1:length(unique_Kct)
        % Extract subset for this alpha_scale & Kct
        subset = results_all(results_all.alpha_scale == unique_alpha(iA) & ...
                             results_all.Kct == unique_Kct(iK), :);

        if ~isempty(subset)
            % Compute mean kinetic energy
            mean_KE_map(iA, iK) = mean(subset.mean_KE);
        else
            mean_KE_map(iA, iK) = NaN;
        end
    end
end
%
figure('Color','w');

x = unique_alpha(:);    % α scale
y = unique_Kct(:);      % Kct
Z = mean_KE_map;

% --- Make sure Z matches [numel(y) x numel(x)] ---
ny = numel(y); nx = numel(x);
[rz, cz] = size(Z);
if ~(rz==ny && cz==nx)
    if (rz==nx && cz==ny)
        Z = Z.';  % transpose if swapped
    else
        error('mean_KE_map size is %dx%d but expected %dx%d.', rz, cz, ny, nx);
    end
end

imagesc(x, y, Z);
set(gca,'YDir','normal');
colormap(parula); colorbar; 

xlabel('$\alpha$ Scale','Interpreter','latex','FontSize',16);
ylabel('$K_{ct}$ (pN/$\mu$m)','Interpreter','latex','FontSize',16);
title('Mean MSV','Interpreter','latex','FontSize',16);
% --- Overlay contoursVer1 ---
hold on;
[AA, KK] = meshgrid(x, y);

% Levels
level_main   = 5.6101;
level_second = 4.59;
level_third  = 3.49;

% Custom colors
custom_colors = struct( ...
    'main',  '#5E6C82', ...  % dark blue-gray
    'second','#299d8f', ...  % teal
    'third', '#f3a361');     % orange

% Contour at KE = 5.6101
[C1, h1] = contour(AA, KK, Z, [level_main level_main], ...
    'LineColor', custom_colors.main, 'LineWidth', 2);

% Contour at KE = 4.59
[C2, h2] = contour(AA, KK, Z, [level_second level_second], ...
    'LineColor', custom_colors.second, 'LineWidth', 2);

% Contour at KE = 3.49
[C3, h3] = contour(AA, KK, Z, [level_third level_third], ...
    'LineColor', custom_colors.third, 'LineWidth', 2);

% Legend
legend([h1 h2 h3], { ...
    '$\langle \mathrm{KE} \rangle = 5.6101$', ...
    '$\langle \mathrm{KE} \rangle = 4.59$', ...
    '$\langle \mathrm{KE} \rangle = 3.49$' ...
}, 'Interpreter','latex', 'Location','southoutside');

hold off;


% --- Overlay contoursVer2 ---
hold on;
[AA, KK] = meshgrid(x, y);

% Levels
level_main = 5.6101;
levels_dotted = [1.57 17.24];

% Solid contour at KE = 5.6101
[C1, h1] = contour(AA, KK, Z, [level_main level_main], 'k', 'LineWidth', 2);

% Dotted contours at KE = 1.57 and 17.24
[C2, h2] = contour(AA, KK, Z, levels_dotted, ...
                   'LineStyle', '--', 'LineColor', 'k', 'LineWidth', 1.2);

legend([h1 h2], { ...
    '$\langle \mathrm{KE} \rangle = 5.61$', ...
    '$\langle \mathrm{KE} \rangle = 1.57 \ \mathrm{and}\ 17.24$' ...
}, 'Interpreter','latex', 'Location','southoutside');


hold off;



%% Probabilistic phase boundary 
threshold_KE = 5.61;  % Example threshold
unique_alpha = unique(results_all.alpha_scale);
unique_Kct   = unique(results_all.Kct);
P_osc = zeros(length(unique_alpha), length(unique_Kct));
for iA = 1:length(unique_alpha)
    for iK = 1:length(unique_Kct)
        % Extract subset for this alpha_scale & Kct
        subset = results_all(results_all.alpha_scale == unique_alpha(iA) & ...
                             results_all.Kct == unique_Kct(iK), :);

        if ~isempty(subset)
            % Count fraction of oscillating runs
            count_osc = sum(subset.mean_KE > threshold_KE);
            P_osc(iA, iK) = count_osc / height(subset);
        else
            P_osc(iA, iK) = NaN;
        end
    end
end
%
figure('Color','w');

x = unique_alpha(:);         % x = α scale
y = unique_Kct(:);                % y = Kct
Z = P_osc;

% --- Make sure Z matches [numel(y) x numel(x)] ---
ny = numel(y); nx = numel(x);
[rz, cz] = size(Z);
if ~(rz==ny && cz==nx)
    if (rz==nx && cz==ny)
        Z = Z.';  % transpose if swapped
    else
        error('P_osc size is %dx%d but expected %dx%d. Fix the data layout.', rz, cz, ny, nx);
    end
end

imagesc(x, y, Z);
set(gca,'YDir','normal');
colormap(parula); colorbar; 

xlabel('$\alpha$ Scale','Interpreter','latex','FontSize',16);
ylabel('$K_{ct}$ (pN/$\mu$m)','Interpreter','latex','FontSize',16);
title('Probability of Oscillation');

hold on;
[AA, KK] = meshgrid(x, y);

p_star = 0.5;
p_band = [0.2 0.8];

% Solid boundary at p=0.5
[~, h1] = contour(AA, KK, Z, [p_star p_star], 'k', 'LineWidth', 2);

% Dotted bands at p=0.2 and p=0.8
[~, h2] = contour(AA, KK, Z, p_band, 'LineStyle','--', 'LineColor','k', 'LineWidth', 1.2);

legend([h1 h2], 'P_{osc}=0.5', 'P_{osc}=0.2 & 0.8', 'Location','southoutside');
hold off;



%%
% Build vectors of outcomes per (alpha_scale, Kct)
% y = [0/1 outcomes], X = [alpha_index, Kct], or fit per alpha row.

K50 = nan(size(unique_alpha));
for iA = 1:numel(unique_alpha)
    % For this alpha row, expand counts into 0/1 trials
    k  = unique_Kct(:);
    p  = P_osc(iA,:).';                      % probability estimates
    n  = n_trials(iA,:).';                   % # iterations per cell
    y1 = round(p.*n); y0 = n - y1;
    % Use binomial GLM with frequency weights
    tbl = table(k, y1, y0, 'VariableNames', {'Kct','Y1','Y0'});
    mdl = fitglm(tbl, [y1 y0] ~ Kct, 'Distribution','binomial', 'Link','logit');
    % K50 solves logit(0.5) = 0 = b0 + b1*K50 => K50 = -b0/b1
    b = mdl.Coefficients.Estimate;
    K50(iA) = -b(1)/b(2);
end

plot(unique_alpha, K50, 'k-', 'LineWidth', 2)


%% Same plot as the ones for determinate 

summary_stats = groupsummary(results_all, {'alpha_scale','Kct'}, ...
    {'mean','std'}, {'alpha','mean_KE'});
summary_stats=removevars(summary_stats, 'GroupCount'); 
% --- Aggregate KE per (alpha_scale, Kct) ---
% --- Plot (errorbar only) ---
figure('Color','w'); hold on;
colors = turbo(length(Alpha_scales));

alpha_unique = unique(summary_stats.alpha_scale);

for ia = 1:numel(alpha_unique)
    a = alpha_unique(ia);

    % Select rows for this alpha
    rows = summary_stats.alpha_scale == a;
    
    kct_vals = summary_stats.Kct(rows);
    m = summary_stats.mean_mean_KE(rows);
    s = summary_stats.std_mean_KE(rows);

    % Errorbar plot (mean ± std)
    errorbar(kct_vals, m, s, ...
        'LineStyle','-', 'Marker', 'o', ...
        'Color', colors(ia,:), 'LineWidth', 1.25, ...
        'MarkerSize', 4, 'DisplayName', ['$\alpha=$ ' num2str(a)]);
end

xlabel('\bf $K_{ct}$ (pN/$\mu$m)','Interpreter','latex', 'FontName','Arial','FontSize',18);
ylabel('\bf Avg. MSV ($\mu$m$^2$/min$^2$)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
legend('Location','eastoutside','Interpreter','latex');
grid on; set(gca,'FontSize',14);



%% --------------------------------------------------------------
% Plot 2: Avg KE vs Alpha (separate curves for each Kct)
% --------------------------------------------------------------
figure('Color','w'); hold on;

% Unique Kct values and a colormap for them
kct_unique = unique(summary_stats.Kct);
colors = turbo(length(Alpha_scales));

for ik = 1:numel(kct_unique)
    k = kct_unique(ik);

    % Select rows for this Kct and sort by alpha for a clean line
    rows = summary_stats.Kct == k;
    a_vals = summary_stats.alpha_scale(rows);
    [a_vals, ord] = sort(a_vals);

    m = summary_stats.mean_mean_KE(rows);
    s = summary_stats.std_mean_KE(rows);

    m = m(ord);
    s = s(ord);

    % Errorbar plot (mean ± std) across alpha for this fixed Kct
    errorbar(a_vals, m, s, ...
        'LineStyle','-', 'Marker', 'o', ...
        'Color', colors(ik,:), 'LineWidth', 1.25, ...
        'MarkerSize', 4, 'DisplayName', sprintf('$K_{ct} = %.0f$', k));
end

xlabel('$\alpha$ (scale factor)', 'Interpreter','latex', 'FontSize', 18);
ylabel('Avg. Kinetic Energy ($\mu$m$^2$/min$^2$)', 'Interpreter','latex', 'FontSize', 18);
title('Bifurcation: Avg KE vs $\alpha$ (by $K_{ct}$)', 'Interpreter','latex', 'FontSize', 20);
legend('Location','eastoutside','Interpreter','latex');
grid on; set(gca,'FontSize',14);
