%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determinate_Sweep_Activity_bifurcation.m
%
% PURPOSE
%   Run a deterministic (noise-free) parameter sweep of:
%     (1) AlphaScale  -> controls Alpha (velocity-dependent attachment term)
%     (2) Kct         -> centromere stiffness (pN/µm)
%   For each (AlphaScale, Kct) pair, simulate sister-chromosome dynamics and
%   quantify chromosome "activity" using chromosome_activity_measurement
%   function. 
%
% IMPORTANT:
%   Avg_KE (code) ≡ MSV (paper)
%   These terms are used interchangeably to describe the same chromosome
%   activity measure derived from time-resolved chromosome velocities.
%
% HOW THIS WAS USED IN THE PAPER
%   - Figure S2B: Phase diagram of Avg MSV/Avg KE 
%   - Figure S2D: Heatmap visualization of Avg MSV/Avg KE on the same grid
%
% KEY OUTPUTS
%   1) CSV summary (used for Python replotting / panel assembly):
%        MSV_determinate_Alpha_3_0.2_15.csv
%      Columns: AlphaScale, Alpha, Kct, Avg_KE
%   2) MATLAB heatmap figure (imagesc) showing Avg_KE over the sweep grid
%
% NOTES
%   - This script stores ONLY per-parameter summary metrics (Avg_KE).
%     It does NOT store full time series to keep output lightweight.
%   - Alpha is derived from AlphaScale by:
%        Alpha = n_dot * AlphaScale / (1 - Beta)
%   - chromosome_activity_measurement(cL, cR, dt) return a table with
%     variable "Avg_KE" (single-row for a single simulation run).
%
% AUTHOR: 
%   Jiali Zhu, August, 2025, UNC-Chapel Hill 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% ---------------- Simulation settings ----------------
dt = 2e-3;              % Timestep (min)
Nsteps = 5000;          % Number of integration steps
duration = 10;          % Duration (min) used in KE analysis
% Kinetochore attachment baseline / coupling parameters
n_dot = 1;              % Basal attachment gain rate (used in NL/NR dynamics)
Beta = 0.7;             % MT tip motion coupling factor (x tip follows chromosome)

% ------------------ Sweep definitions ------------------
Alpha_scales = 3:0.2:15;       % AlphaScale: dimensionless scale factor;
Kct_values = 10:0.1:20;        % Kct: centromere stiffness parameter to sweep (pN/µm)
nTotal=numel(Alpha_scales)*numel(Kct_values); % Total number of parameter combinations (for progress reporting)

% ------------------ Storage for results ------------------
results = [];
counter = 0;

% ========================================================================
%                         PARAMETER SWEEP LOOP
% =========================================================================
for a = Alpha_scales
    Alpha = n_dot * a / (1 - Beta);  % Recalculate Alpha
    for Kct = Kct_values
        counter = counter + 1;
        % ------------------ Fixed model parameters (not swept) ------------------
        I0 = 2;                  % µm     centromere rest length
        Kkt = 1;                 % pN/µm  KT–MT spring constant
        Gamma = 0.1;             %        effective drag coefficient
        Nmax = 25;               % count  max number of KT–MT attachments
        Nbar = 20;               % attachment number at steady state (used for initialization)
        Lambda = n_dot / Nbar;   % detachment rate constant
        epsilon = 0.1;           % small symmetry-breaking perturbation (dimensionless)

        % ------------------ Allocate state variables ------------------
        xL = zeros(Nsteps+1, 1); xR = zeros(Nsteps+1, 1);       % left/right MT tip x-position over time
        NL = zeros(Nsteps+1, 1); NR = zeros(Nsteps+1, 1);       % left/right KT–MT attachment number over time 
        cL = zeros(Nsteps+1, 2, 1); cR = zeros(Nsteps+1, 2, 1); % left/right chromosome position: (time, [x y], chromosome index)
        vL = zeros(Nsteps, 1); vR = zeros(Nsteps, 1);           % left/right chromosome velocity over time
        % ------------------ Initial conditions ------------------
        % Initial attachment numbers (slightly perturbed to break symmetry)
        NR(1) = Nbar * (1 - epsilon);
        NL(1) = Nbar * (1 + epsilon);
        % Initial chromosome positions (x in µm; y fixed at 0) 
        cL(1,:,1) = [-I0/2 - epsilon, 0];
        cR(1,:,1) = [ I0/2 + epsilon, 0];
        % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
        xL(1) = cL(1,1,1) - 0.3;
        xR(1) = cR(1,1,1) + 0.3;

        % ------------------ Run Simulation ------------------
        for t = 1:Nsteps
            % ----- Force calculations -----
            F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));  % force on left chromosome from KT–MT bundle
            F_CT_L = Kct * (cR(t,1,1) - cL(t,1,1) - I0);  % Centromere forces on Left Chromosome
            F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));  % force on right chromosome from KT–MT bundle
            F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0); % Centromere forces on Right Chromosome
            % ----- Velocities (overdamped dynamics) ----- 
            vL(t) = (F_KT_L + F_CT_L) / Gamma;
            vR(t) = (F_KT_R + F_CT_R) / Gamma;
            % ----- Update chromosome positions ----- 
            cL(t+1,1,1) = cL(t,1,1) + dt * vL(t);
            cR(t+1,1,1) = cR(t,1,1) + dt * vR(t);
            % ----- Update MT tip positions -----
            xL(t+1) = xL(t) + dt * Beta * vL(t);
            xR(t+1) = xR(t) + dt * Beta * vR(t);
            % ----- Update attachment numbers -----
            NL(t+1) = NL(t) + dt * ( ...
                n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
                - Lambda * NL(t) );

            NR(t+1) = NR(t) + dt * ( ...
                n_dot + ( Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
                - Lambda * NR(t) );
        end

        % ------------------ Activity Measurement ------------------
        activity_tbl = chromosome_activity_measurement(cL, cR, dt);
        avg_KE = activity_tbl.Avg_KE;   % assumes single-row output
        % ------------------ Store result ------------------
        results = [results; struct( ...
            'AlphaScale', a, ...
            'Alpha',      Alpha, ...
            'Kct',        Kct, ...
            'Avg_KE',     avg_KE )];
        % ------------------ Progress Print ------------------
        fprintf('Progress: %d / %d (%.1f%%) | Alpha = %.1f, Kct = %.2f\n', ...
            counter, nTotal, 100*counter/nTotal, a, Kct);
    end
end
results_table = struct2table(results); %Convert struct array to table 
%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_table, 'MSV_determinate_Alpha_3_0.2_15.csv');      % per-point CSV

%% ========================================================================
%   FIGURE S2B / S2D: PHASE MAP HEATMAP 
% =========================================================================

% --- Load sweep results from CSV (exported by MATLAB) ---
csvFile = "Output_Sweep_Bifurcation/MSV_determinate_Alpha_3_0.2_15.csv";
results_table = readtable(csvFile);

% --- Build grid: rows = Kct (y), cols = AlphaScale (x) ---
aVals = sort(unique(results_table.AlphaScale));   % AlphaScale (x)
kVals = sort(unique(results_table.Kct));          % Kct (y)

[~, ia] = ismember(results_table.AlphaScale, aVals);
[~, ik] = ismember(results_table.Kct,        kVals);

% Z(y,x) = mean Avg_KE for each (Kct, AlphaScale)
Z = accumarray([ik ia], results_table.Avg_KE, ...
               [numel(kVals) numel(aVals)], @mean, NaN);

% --- Custom sci31 colormap (teal = low, red = high) ---
sci31_hex = ["#C65B3F"; "#DEA091"; "#EDCCC5"; "#D7E6EA"; "#77A5B3"; "#408094"];
sci31_rgb = zeros(numel(sci31_hex), 3);
for i = 1:numel(sci31_hex)
    h = char(sci31_hex(i));
    sci31_rgb(i,:) = [hex2dec(h(2:3)) hex2dec(h(4:5)) hex2dec(h(6:7))] / 255;
end
cmap = flipud(interp1(linspace(0,1,size(sci31_rgb,1)), sci31_rgb, linspace(0,1,256)));

% --- Plot heatmap ---
fig = figure('Color','w');
ax = axes(fig);

imagesc(ax, aVals, kVals, Z);
set(ax, 'YDir','normal');  
colormap(ax, cmap);
axis(ax, 'tight');

% Limits (match the Python view window)
xlim(ax, [3.0 15.0]);
ylim(ax, [10 20]);          

% Ticks similar to Python
xticks(ax, 3:1:15);
yticks(ax, 10:1:20);
ax.XTickLabel = arrayfun(@(v) sprintf('%.1f', v), ax.XTick, 'UniformOutput', false);

set(ax, 'TickDir','out', 'LineWidth',1, 'FontSize',8);
xlabel(ax, '$\alpha$ Scale', 'Interpreter','latex');
ylabel(ax, '$K_{ct}$ (pN/$\mu$m)', 'Interpreter','latex');

% --- Horizontal colorbar at bottom ---
cb = colorbar(ax, 'southoutside');
cb.Label.String = 'Mean MSV';
cb.Ticks = 0:20:80;
cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), cb.Ticks, 'UniformOutput', false);

%% ========================================================================
% OPTIONAL DIAGNOSTIC PLOTS (commented out)
%   - Avg KE vs Kct grouped by AlphaScale
%   - Avg KE vs AlphaScale grouped by Kct
%   These are useful for sanity checks but not required for S2B/S2D panels.
% =========================================================================

% % ------------------ Plot: Avg KE vs Kct (grouped by AlphaScale) ------------------
% figure; hold on;
% alpha_values = unique(results_table.AlphaScale);
% colors = turbo(length(alpha_values));
% for i = 1:length(alpha_values)
%     a_val = alpha_values(i);
%     subset = results_table(results_table.AlphaScale == a_val, :);
%     plot(subset.Kct, subset.Avg_KE, ...
%         'LineStyle','-', 'Marker','o', 'Color', colors(i,:), ...
%         'LineWidth', 2, 'MarkerSize', 5, ...
%         'DisplayName', sprintf('$\\alpha = %.1f$', a_val));
% end
% xlabel('\bf $K_{ct}$ (pN/$\mu$m)','Interpreter','latex', 'FontName','Arial','FontSize',18);
% ylabel('\bf Avg. MSV ($\mu$m$^2$/min$^2$)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
% legend('Location','best', 'Interpreter','latex');
% grid on; set(gca,'FontSize',14);
%
% % ------------------ Plot: Avg KE vs alpha scale (grouped by Kct) ------------------
% figure; hold on;
% kct_values = unique(results_table.Kct);
% colors = turbo(length(kct_values));
% for i = 1:length(kct_values)
%     kct = kct_values(i);
%     subset = results_table(results_table.Kct == kct, :);
%     plot(subset.AlphaScale, subset.Avg_KE, ...
%         'LineStyle','-', 'Marker','o', 'Color', colors(i,:), ...
%         'LineWidth', 2, 'MarkerSize', 5, ...
%         'DisplayName', sprintf('$K_{ct} = %.2f$', kct));
% end
% xlabel('\bf$\alpha$ (scale factor)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
% ylabel('\bf Avg. MSV ($\mu$m$^2$/min$^2$)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
% legend('Location','best', 'Interpreter','latex');
% grid on; set(gca,'FontSize',14);
