%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha_Kct_Sweep_MSV.m
%
% PURPOSE
%   Run a stochastic (multi-iteration; noisy) parameter sweep of:
%     (1) AlphaScale  -> controls Alpha (velocity-dependent attachment term)
%     (2) Kct         -> centromere stiffness (pN/µm)
%   while INCLUDING inter-chromosomal mechanical coupling via a stochastic
%   spring network between chromosome centers-of-mass (CMs).
%
%   For each (AlphaScale, Kct) pair, simulate Nchromosomes sister-chromosome
%   dynamics for N_iter independent “cells” (noise realizations) and quantify
%   chromosome "activity" using chromosome_activity_measurement.
%
% IMPORTANT:
%   Avg_KE (code) ≡ MSV (paper)
%   These terms are used interchangeably to describe the same chromosome
%   activity measure derived from time-resolved chromosome velocities.
%
% HOW THIS SCRIPT IS USED in the paper: 
%   -Figure S3B: Parameter-phase heatmap of mean MSV (chromosome activityt)
%    on the (AlphaScale, Kct) grid in the presence of inter-chromosomal coupling.
%
% NOTES / IMPLEMENTATION DETAILS
%   To facilitate reproducibility without rerunning the full sweep, the
%   processed CSV output has been COMMITTED to the repository:Users can 
%   directly load this CSV file to reproduce Figure S3B without re-running
%   this script.
%
% KEY OUTPUTS (IN-MEMORY)
%  1) CSV summary (used for Python replotting / panel assembly):
%     Columns:alpha_scale, alpha, Kct, iteration, chrom_idx, MSV 
%
% DEPENDENCIES 
%   spring_connect.m
%   chromosome_activity_measurement.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC–Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% ---------------- Simulation settings ----------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
duration  = 10;             % Duration (min) used in KE analysis
% ---------------- Number of chromosome ---------------- 
Nchromosomes = 5;           % Number of chromosome pairs

% ------------------ Sweep definitions ------------------
Alpha_scales = 3:0.2:15;    % Alpha scale factors
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
    'VariableNames', {'alpha_scale','alpha','Kct','iteration','chrom_idx','MSV'});


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
            cL = zeros(Nsteps+1, 2,Nchromosomes);       % left chromosome position: (time, [x y], chromosome index)
            cR = zeros(Nsteps+1, 2,Nchromosomes);       % right chromosome position: (time, [x y], chromosome index)
            xL = zeros(Nsteps+1, Nchromosomes);         % left MT tip x-position over time
            xR = zeros(Nsteps+1, Nchromosomes);         % right MT tip x-position over time
            NL = zeros(Nsteps+1, Nchromosomes);         % left KT–MT attachment number over time 
            NR = zeros(Nsteps+1, Nchromosomes);         % right KT–MT attachment number over time 
            vL = zeros(Nsteps, Nchromosomes);           % left chromosome velocity over time 
            vR = zeros(Nsteps, Nchromosomes);           % right chromosome velocity over time

            % ---------------- Initial conditions ----------------
            % Spring connectivity state (updated each time step)
            flag_c=zeros(Nchromosomes,Nchromosomes);   % current adjacency (1 = connected)
            flag_cp=flag_c;                            % Previous state 
            flag_c_sum_accumulated = zeros(Nsteps, 1); % To accumulate sum of flag_c over time
            % Fixed y positions are assigned
            y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
            % Each Pairs have their own 
            for i = 1: Nchromosomes 
                % Initial attachment numbers (slightly perturbed to break symmetry)
                NR(1,i)=Nbar(i)*(1-epsilon);
                NL(1,i)=Nbar(i)*(1+epsilon); 
                % Initial chromosome positions (x in µm; y fixed) 
                cL(1,:,i) =-I0/2-epsilon*(2*rand(1)-1);  % left sister initial x-position
                cR(1,:,i) =I0/2+epsilon*(2*rand(1)-1);   % right sister initial x-position
                cL(:,2,i) =  y_positions(i);             % y fixed
                cR(:,2,i) =  y_positions(i);             % y fixed
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
            skipped = false;
            if any(isnan(cL), 'all') || any(isnan(cR), 'all')
                skipped = true;
            else
                % ---- ---------------- Activity metrics----------- ----
                activity_table = chromosome_activity_measurement(cL, cR, dt);
                % get numeric vector of KE per chromosome
                T = activity_table(:, {'PairID','Avg_KE'});
                T.Properties.VariableNames = {'chrom_idx','MSV'};
             
                % Add run-level columns (same value for all rows in this iteration)
                T.alpha_scale = repmat(a,               height(T), 1);
                T.alpha       = repmat(alpha_effective, height(T), 1);
                T.Kct         = repmat(Kct,             height(T), 1);
                T.iteration   = repmat(double(iter),    height(T), 1);
                
                % Put run-level columns first to match results_all order
                T = movevars(T, {'alpha_scale','alpha','Kct','iteration'}, 'Before', 'chrom_idx');
                results_all = [results_all; T];

                % --- --------------------Progress note ---------------
                counter = counter + 1;
                fprintf('Progress: %d/%d (%.1f%%) | a=%.1f, Kct=%.2f, iter=%d\n', ...
                    counter, nTotal, 100*counter/nTotal, a, Kct, iter);
            end
        end
    end
end

%% ========================================================================
%   FIGURE S3B: HEATMAP  
% =========================================================================

% --- Load sweep results from CSV (exported by MATLAB) ---
csvFile = "Output_Sweep_Bifurcation/kct_alpha_MSV_results.csv";
results_all = readtable(csvFile);

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
            mean_KE_map(iA, iK) = mean(subset.MSV);
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


%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_all, 'kct_alpha_MSV.csv');   % CSV for Python
