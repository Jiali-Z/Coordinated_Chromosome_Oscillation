%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise_Sweep_Activity_bifurcation.m
%
% PURPOSE
%   Run a stochastic (noise-enabled) parameter sweep over:
%     (1) AlphaScale  -> controls Alpha (velocity-dependent attachment term)
%     (2) Kct         -> centromere stiffness (pN/µm)
%   For each (AlphaScale, Kct) pair, simulate N_iter independent “cells”
%   (replicates) with cell-to-cell parameter variability and Brownian motion,
%   then quantify chromosome "activity" using chromosome_activity_measurement().
%
% STOCHASTICITY INCLUDED
%   1) Cell-to-cell variability:
%        n_dot  and Nbar are randomly perturbed each replicate (uniform noise).
%   2) Brownian motion:
%        A diffusive displacement term is added each timestep.
%
% IMPORTANT TERMINOLOGY
%   Avg_KE (code) ≡ MSV (paper)
%   These terms are used interchangeably to describe the same chromosome
%   activity measure derived from time-resolved chromosome velocities.
%
% HOW THIS WAS USED IN THE PAPER
%   - Figure S2D: Heatmap visualizationof Avg KE/ Avg MSV across (AlphaScale, Kct)
%
% KEY OUTPUTS
%   1) CSV summary (used for Python replotting / panel assembly):
%        MSV_noise_alpha_3_0.2_15.csv
%      Columns: alpha_scale, alpha, Kct, iteration, mean_KE
%      NOTE: mean_KE column stores Avg_KE ≡ MSV.
%
% NOTES
%   - This script stores ONLY per-replicate summary metrics (mean_KE) and
%     does NOT store full trajectories to keep outputs lightweight.
%   - Because the full Alpha–Kct sweep with N_iter replicates is
%     computationally expensive,the processed CSV summary output 
%     (MSV_noise_Alpha_3_0.2_15.csv) is committed to the repository 
%     and can be used directly to reproduce Figure S2D and Figure S2D.
%
% AUTHOR:
%   Jiali Zhu, August 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% ---------------- Simulation settings ----------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
Beta      = 0.7;            % MT tip motion coupling factor (x tip follows chromosome)
duration  = 10;             % Duration (min) used in KE analysis
% ------------------ Sweep definitions ------------------
Alpha_scales = 3:0.2:15;    % Alpha scale factors
Kct_values   = 10:0.1:20;   % Kct values to sweep
nTotal  = numel(Alpha_scales) * numel(Kct_values) * N_iter;
counter = 0;

% ------------------ Stochasticity parameters ------------------
N_iter  = 50;               % number of noisy replicates per (AlphaScale, Kct)
noise   = 0.05;             % variability for n_dot and Nbar (cell-to-cell)
noise_b = 0.05;             % Brownian noise strength (diffusive displacement)

% ------------------ Fixed model constants ------------------
I0   = 2;        % µm: centromere rest length
Kkt  = 1;        % pN/µm: KT–MT spring constant
Gamma= 0.1;      % effective drag coefficient
Nmax = 25;       % max attachments per kinetochore
epsilon = 0.1;    % small symmetry-breaking perturbation (dimensionless)                             

% ------------------ Storage for results ------------------
results_all = table( ...
    'Size',[0 5], ...
    'VariableTypes', {'double','double','double','double','double'}, ...
    'VariableNames', {'alpha_scale','alpha','Kct','iteration','mean_KE'});  %mean_KE == Avg_KE == MSV


% ========================================================================
%                         PARAMETER SWEEP LOOP
% =========================================================================
for a = Alpha_scales
    for Kct = Kct_values
        for iter = 1:N_iter
            % ---- cell-specific parameters with noise ----
            n_dot  = 1 * (1 + noise*(2*rand - 1));             % Basal attachment gain rate
            Nbar   = 20   * (1 + noise*(2*rand - 1));          % attachment number at steady state
            Lambda = n_dot / Nbar;                             % detachment rate constant
            Alpha  = n_dot * a / (1 - Beta);                   % velocity-attachment coupling coefficient 

            % ------------------ Allocate state variables ------------------
            xL = zeros(Nsteps+1,1); xR = zeros(Nsteps+1,1);    % left/right MT tip x-position over time
            NL = zeros(Nsteps+1,1); NR = zeros(Nsteps+1,1);    % left/right KT–MT attachment number over time 
            cL = zeros(Nsteps+1,2,1); cR = zeros(Nsteps+1,2,1);% left/right chromosome position: (time, [x y], chromosome index)
            vL = zeros(Nsteps+1,1);  vR = zeros(Nsteps+1,1);   % left/right chromosome velocity over time
            % ------------------ Initial conditions ------------------
            % Initial attachment numbers (slightly perturbed to break symmetry)
            NR(1)     = Nbar*(1 - epsilon);
            NL(1)     = Nbar*(1 + epsilon);
            % Initial chromosome positions (x in µm; y fixed at 0) 
            cL(1,1,1) =-I0/2-epsilon;
            cR(1,1,1) =I0/2+epsilon;
            cL(:,2,1) = 0;
            cR(:,2,1) = 0;
            % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
            xL(1) = cL(1,1,1) - 0.3;
            xR(1) = cR(1,1,1) + 0.3;

            % ------------------ Run Simulation ------------------
            for t = 1:Nsteps
                % ----- Force calculations ----- 
                F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));     % force on left chromosome from KT–MT bundle
                F_CT_L =  Kct   * (cR(t,1,1) - cL(t,1,1) - I0);  % Centromere forces on Left Chromosome
                F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));     % force on right chromosome from KT–MT bundle
                F_CT_R = -Kct   * (cR(t,1,1) - cL(t,1,1) - I0);  % Centromere forces on Right Chromosome
                % ----- Velocities (overdamped dynamics) ----- 
                vL(t) = (F_KT_L + F_CT_L ) / Gamma;
                vR(t) = (F_KT_R + F_CT_R ) / Gamma;
                % ----- Update chromosome positions ----- 
                % Brownian step 
                dx_diff = noise_b * randn(1,1) * sqrt(dt);
                cL(t+1,1,1) = cL(t,1,1) + dt*vL(t) + dx_diff;
                cR(t+1,1,1) = cR(t,1,1) + dt*vR(t) + dx_diff;
                % ----- Update MT tip positions -----
                xL(t+1) = xL(t) + dt * Beta * vL(t);
                xR(t+1) = xR(t) + dt * Beta * vR(t);
                % ----- Update attachment numbers -----
                NL(t+1) = NL(t) + dt * ( ...
                    n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
                    - Lambda * NL(t));
                NR(t+1) = NR(t) + dt * ( ...
                    n_dot + ( Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
                    - Lambda * NR(t));
            end

            % ------------------ Activity Measurement ------------------
            activity_table = chromosome_activity_measurement(cL, cR, dt);
           % ------------------ Store result ------------------
            results_all = [results_all; ...
                table(a, Alpha, Kct, double(iter), activity_table.Avg_KE, ...
                'VariableNames', {'alpha_scale','alpha','Kct','iteration','mean_KE'})];
             % ------------------ Progress Print ------------------
            counter = counter + 1;
            fprintf('Progress: %d/%d (%.1f%%) | a=%.1f, Kct=%.2f, iter=%d\n', ...
                counter, nTotal, 100*counter/nTotal, a, Kct, iter);
        end
    end
end

%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_all, 'MSV_noise_alpha_3_0.2_15.csv');   % CSV for Python

%% ========================================================================
%   FIGURE S2D /3C: HEATMAP  
% =========================================================================

% --- Load sweep results from CSV (exported by MATLAB) ---
csvFile = "Output_Sweep_Bifurcation/MSV_noise_alpha_3_0.2_15.csv";
results_all = readtable(csvFile);

% --- Build grid: rows = Kct (y), cols = alpha_scale (x) ---
aVals = sort(unique(results_all.alpha_scale));   % AlphaScale (x)
kVals = sort(unique(results_all.Kct));           % Kct (y)

[~, ia] = ismember(results_all.alpha_scale, aVals);
[~, ik] = ismember(results_all.Kct,        kVals);

% Z(y,x) = mean Avg_KE for each (Kct, AlphaScale)
Z = accumarray([ik ia], results_all.mean_KE, [numel(kVals) numel(aVals)], @mean, NaN);

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

% Limits
xlim(ax, [3.0 15.0]);
ylim(ax, [10 20]);
xticks(ax, 3:1:15);
yticks(ax, 10:1:20);
ax.XTickLabel = arrayfun(@(v) sprintf('%.1f', v), ax.XTick, 'UniformOutput', false);

set(ax, 'TickDir','out', 'LineWidth',1, 'FontSize',6);
xlabel(ax, '$\alpha$ Scale', 'Interpreter','latex');
ylabel(ax, '$K_{ct}$ (pN/$\mu$m)', 'Interpreter','latex');

% --- Horizontal colorbar (bottom) ---
cb = colorbar(ax, 'southoutside');
cb.Label.String = 'Mean MSV';
cb.Ticks = 0:20:80;
cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), cb.Ticks, 'UniformOutput', false);

% --- Overlay treatment contours (Ctl/Noc/TSA) ---
hold(ax, 'on');
[AA, KK] = meshgrid(aVals, kVals);
levels = struct('Ctl', 5.6101, 'Noc', 4.59, 'TSA', 3.49);
colors = struct('Ctl', '#5E6C82', 'Noc', '#D6CDBE', 'TSA', '#81B3A9');

% Draw contours and keep handles for legend
hCtl = contour(ax, AA, KK, Z, [levels.Ctl levels.Ctl], 'LineWidth',2, 'LineColor', colors.Ctl);
hNoc = contour(ax, AA, KK, Z, [levels.Noc levels.Noc], 'LineWidth',2, 'LineColor', colors.Noc);
hTSA = contour(ax, AA, KK, Z, [levels.TSA levels.TSA], 'LineWidth',2, 'LineColor', colors.TSA);



