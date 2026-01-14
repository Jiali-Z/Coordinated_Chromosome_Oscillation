%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha_Kct_Sweep_Pearson.m
%
% PURPOSE
%   Run a stochastic (multi-iteration; noisy) parameter sweep of:
%     (1) AlphaScale  -> scales Alpha (velocity-dependent attachment term)
%     (2) Kct         -> centromere stiffness (pN/µm)
%   while INCLUDING inter-chromosomal mechanical coupling via a stochastic
%   spring network between chromosome centers-of-mass (CMs).
%
%   For each (AlphaScale, Kct) pair, simulate Nchromosomes sister-chromosome
%   dynamics for N_iter independent “cells” (noise realizations) and quantify
%   chromosome coordination using Pearson’s correlation across CM traces
%   (pearsons_correlation_analysis).
%
% HOW THIS SCRIPT IS USED IN THE PAPER
%   - Figure S3A: Parameter-phase heatmap of mean Pearson’s r on the
%     (AlphaScale, Kct) grid in the presence of inter-chromosomal coupling.
%
% NOTES / IMPLEMENTATION DETAILS
%   To facilitate reproducibility without rerunning the full sweep, the
%   processed CSV output has been COMMITTED to the repository:Users can 
%   directly load this CSV file to reproduce Figure S3A without re-running
%   this script.
%
% KEY OUTPUTS
%   1) CSV summary (for Python replotting / panel assembly):
%      Columns:alpha_scale, alpha_effective, Kct, iteration, mean_pearsonR, std_pearsonR
%
% DEPENDENCIES 
%   spring_connect.m
%   pearsons_correlation_analysis.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC–Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% ---------------- Simulation settings ----------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps

% ---------------- Number of chromosome ---------------- 
Nchromosomes = 5;           % Number of chromosome pairs

% ------------------ Sweep definitions ------------------
Alpha_scales = 3:0.2:15;      % Alpha scale factors
Kct_values   = 10:0.1:20;   % Kct values to sweep
nTotal  = numel(Alpha_scales) * numel(Kct_values) * N_iter;
counter = 0;

% ------------------ Stochasticity parameters ------------------
N_iter  = 50;               % Number of cells per (alpha, Kct)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (diffusive displacement)


% ------------------ Fixed model constants ------------------
I0   = 2;                   % µm: centromere rest length
Kkt  = 1;                   % pN/µm: KT–MT spring constant
Gamma= 0.1;                 % effective drag coefficient
Nmax = 25;                  % max attachments per kinetochore
Beta = 0.7;                 % MT tip motion coupling factor (x tip follows chromosome)
epsilon = 0.1;              % small symmetry-breaking perturbation (dimensionless)   

% ---------------- Inter-chromosome coupling parameters ----------------
koff_0=40;                  % 1/s    base disconnection rate
kon_0=20;                   % 1/s    base disconnection rate 
l0 = 3;                     % µm     centromere rest length 
k = 0.5;                    % pN/µm  coupling strength 
Nspring = 1;                % Number of springs per chromosome not used 


% ------------------ Storage for results ------------------
results_all = table( ...
    'Size',[0 6], ...
    'VariableTypes', {'double','double','double','double','double','double'}, ...
    'VariableNames', {'alpha_scale','alpha_effective','Kct','iteration','mean_pearsonR','std_pearsonR'});

% ========================================================================
%                         PARAMETER SWEEP LOOP
% =========================================================================
for a = Alpha_scales
    for Kct = Kct_values
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
                    F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0);        % Centromere forces on Left Chromosome
                    F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));    % force on right chromosome from KT–MT bundle
                    F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0);       % Centromere forces on Right Chromosome
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
            % Skip failed runs (NaNs), otherwise append results to pooled tables.
            % Reject runs with NaNs in trajectories
            skipped = false;
            if any(isnan(cL), 'all') || any(isnan(cR), 'all')
                skipped = true;
            else
                % ---- Cross-correlation: neighbor=1, τ in first group (5–50 s) ----
                CMx  = squeeze((cL(1:Nsteps,1,:) + cR(1:Nsteps,1,:)) / 2);   % [Nsteps x Nchromosomes]\
                pearson_table = pearsons_correlation_analysis(cL, cR, dt, 10);
                mean_r=mean(pearson_table.pearson_r,'omitnan');
                std_r = std(pearson_table.pearson_r,'omitnan');
                % Append the results 
                results_all = [results_all; table( ...
                    a, alpha_effective, Kct, double(iter), mean_r, std_r, ...
                    'VariableNames', {'alpha_scale','alpha_effective','Kct','iteration','mean_pearsonR','std_pearsonR'})];
            end
          % --- --------------------Progress note ---------------
          counter = counter + 1;
          fprintf('Progress: %d/%d (%.1f%%) | a=%.1f, Kct=%.2f, iter=%d | mean r=%.3f\n', ...
            counter, nTotal, 100*counter/nTotal, a, Kct, iter, mean_r);
        end % iter
    end
end



%% ========================================================================
%   FIGURE S3B: HEATMAP  
% =========================================================================


% --- Load sweep results from CSV (exported by MATLAB) ---
csvFile = "Output_Sweep_Bifurcation/kct_alpha_pearson_results.csv";
results_all = readtable(csvFile);

% --- Aggregate across all iterations and chromosomes per (alpha_scale, Kct) ---
summary_stats = groupsummary(results_all, {'alpha_scale','Kct'}, 'mean', 'mean_pearsonR');
% The table has columns: alpha_scale, Kct, GroupCount, mean_MSV

% Extract unique sorted parameter values
unique_alpha = unique(summary_stats.alpha_scale);
unique_Kct   = unique(summary_stats.Kct);

% Preallocate grid: rows = Kct, columns = alpha_scale
Z = nan(numel(unique_Kct), numel(unique_alpha));

% Fill the grid
for ia = 1:numel(unique_alpha)
    for ik = 1:numel(unique_Kct)
        mask = (summary_stats.alpha_scale == unique_alpha(ia)) & ...
               (summary_stats.Kct == unique_Kct(ik));
        if any(mask)
            Z(ik, ia) = mean(summary_stats.mean_mean_pearsonR(mask), 'omitnan');
        end
    end
end

% --- Plot heatmap ---
figure('Color','w');
imagesc(unique_alpha, unique_Kct, Z);
set(gca,'YDir','normal');  % low Kct at bottom
axis tight;
colormap(parula);
cb = colorbar;
ylabel(cb, 'Mean PearsonR ($\mu$m$^2$/min$^2$)', 'Interpreter','latex');

xlabel('$\alpha$ Scale','Interpreter','latex','FontSize',16);
ylabel('$K_{ct}$ (pN/$\mu$m)','Interpreter','latex','FontSize',16);
title('Mean PearsonR Heatmap','Interpreter','latex','FontSize',16);

% Optional cosmetic tuning
set(gca, 'XTick', unique_alpha, 'YTick', unique_Kct);
box on;

% --- Optional contour overlay at specific MSV level ---
hold on;
[AA, KK] = meshgrid(unique_alpha, unique_Kct);
contour_level = 5.6;  % adjust as desired
[C, h] = contour(AA, KK, Z, [contour_level contour_level], 'k', 'LineWidth', 2);
clabel(C, h, 'FontSize', 10, 'Color', 'k');
hold off;



%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_all, 'kct_alpha_pearson_results.csv');   % CSV for Python
