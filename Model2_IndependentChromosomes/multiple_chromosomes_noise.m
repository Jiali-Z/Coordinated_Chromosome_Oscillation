%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiple_chromosomes_noise.m
%
% PURPOSE
%   Simulate N independent sister-chromosome pairs in 1D (x), each with
%   chromosome-specific parameter noise and measurement-like positional jitter.
%   This script is a negative control: because chromosome pairs are independent,
%   cross-chromosome coordination metrics (Pearson r, displacement cross-corr)
%   should be ~0 on average.
%
% WHAT'S DIFFERENT FROM one_oscillating_chromosome_noise.m: 
%   - Simulates N chromosome pairs (indexed by i = 1..Nchromosomes) 
%   - Draws kinetic parameters independently for each chromosome i
%   - Assigns fixed y-positions for visualization only (dynamics remain 1D in x)
%   - Adds coordination analyses across chromosome pairs (Pearson + displacement cross-corr)
%
% OUTPUTS
%   Figures:
%     (1) CM traces over time for all chromosomes (stacked/offset view)
%     (2+) Peaks and valleys identified by oscillation_measurement for each chromosome (if plotting enabled)
%     (2) Activity metrics bar plot per chromosome (Avg_KE, Avg|v_L|, Avg|v_R|)
%     (3) Pearson correlation distribution across chromosome pairs
%     (4) Cross-correlation vs neighbor order (two-point microrheology analysis of chromosome embedded in spindle)
%     (5) Cross-correlation vs lag time (two-point microrheology analysis of chromosome embedded in spindle)
%
%   Printed tables:
%     - Per-chromosome CM/KK oscillation summary (mean ± SD)
%     - Per-chromosome activity summary
%   
%   Optional CSV:
%     – Per-chromosome cL/cR trajectories vs time for Python-based plotting
%
% DEPENDENCIES
%   oscillation_measurement.m
%   chromosome_activity_measurement.m
%   pearsons_correlation_analysis.m
%   crosscorrelation_measurement.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

% ---------------- Time information ---------------- 
dt = 2e-3;         % time step (min)
Nsteps = 5000;     % number of simulation steps 

% ---------------- Number of chromosome ---------------- 
Nchromosomes=5;    

% ---------------- Noise parameters ----------------
noise=0.05;        % magnitude of uniform parameter perturbations for n_dot and Nbar
noise_b=0.05;      % magnitude of common-mode positional jitter added each time step

% ---------------- Model parameters ----------------
Kct = 12.30;       % pN/µm  centromere spring constant  
I0 = 2;            % µm     centromere rest length
Kkt = 1;           % pN/µm  KT–MT spring constant 
Gamma = 0.1;       %        effective drag coefficient 
Beta = 0.7;        %        scaling for MT tip velocity
Nmax = 25;         % count  max number of KT–MT attachments
epsilon =0.1;      % small symmetry-breaking perturbation (dimensionless) 

% ---------------- Parameter noise per chromosome ----------------
% Each chromosome i gets its own kinetic parameters (chromosome-to-chromosome variability).
% These are drawn ONCE (at initialization) and held fixed during the simulation.
Nbar   = zeros(Nchromosomes, 1);  % steady-state attachment number (initialization reference)
n_dot  = zeros(Nchromosomes, 1);  % baseline gain term for attachments
Lambda = zeros(Nchromosomes, 1);  % detachment rate constant
Alpha  = zeros(Nchromosomes, 1);  % velocity-attachment coupling coefficient
for i = 1:Nchromosomes
    Nbar(i)  = 20 * (1 + noise * (2*rand(1) - 1));
    n_dot(i) =  1 * (1 + noise * (2*rand(1) - 1));
    Lambda(i) = n_dot(i) / Nbar(i);
    Alpha(i)  = n_dot(i) * 6.2 / (1 - Beta);
end
% Print parameter draws
fprintf("Ndot = %4.4f\n", n_dot);
fprintf("Nbar = %4.4f\n", Nbar);
fprintf("Lambda = %4.4f\n", Lambda);
fprintf("Alpha = %4.4f\n", Alpha);

% ---------------- Pre-allocate state variables ----------------
cL = zeros(Nsteps+1, 2,Nchromosomes);       % left chromosome position: (time, [x y], chromosome index)
cR = zeros(Nsteps+1, 2,Nchromosomes);       % right chromosome position: (time, [x y], chromosome index)
xL = zeros(Nsteps+1, Nchromosomes);         % left MT tip x-position over time
xR = zeros(Nsteps+1, Nchromosomes);         % right MT tip x-position over time
NL = zeros(Nsteps+1, Nchromosomes);         % left KT–MT attachment number over time 
NR = zeros(Nsteps+1, Nchromosomes);         % right KT–MT attachment number over time 
vL = zeros(Nsteps, Nchromosomes);           % left chromosome velocity over time
vR = zeros(Nsteps, Nchromosomes);           % right chromosome velocity over time

% ---------------- Initial conditions ----------------
% Fixed y positions are assigned
y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
% Each Pairs have their own 
for i = 1: Nchromosomes 
    % Initial attachment numbers (slightly perturbed to break symmetry)
    NR(1,i)=Nbar(i)*(1-epsilon);
    NL(1,i)=Nbar(i)*(1+epsilon); 
    % Initial chromosome positions (x in µm; y fixed) 
    cL(1,1,i) =-I0/2-epsilon*(2*rand(1)-1); % left sister initial x-position
    cR(1,1,i) =I0/2+epsilon*(2*rand(1)-1);  % right sister initial x-position
    cL(:,2,i) = y_positions(i);             % y fixed
    cR(:,2,i) = y_positions(i);             % y fixed
    % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
    xL(1,i) = cL(1,1,i)-0.3;
    xR(1,i) = cR(1,1,i)+0.3;
    % Initial chromosome velocities
    vL(1,i)=0;
    vR(1,i)=0;
end 

% ---------------- Time integration loop (explicit Euler) ----------------
for t = 1:Nsteps
    for i=1:Nchromosomes 
    % ----- Force calculations -----
    % Define centromere and MT forces on the given chromosome
    F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));  % force on left chromosome from KT–MT bundle
    F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0);      % Centromere forces on Left Chromosome
    F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));  % force on right chromosome from KT–MT bundle
    F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0);     % Centromere forces on Right Chromosome
    % ----- Velocities (overdamped dynamics) -----
    vL(t,i) = (F_KT_L + F_CT_L)/Gamma;
    vR(t,i) = (F_KT_R + F_CT_R)/Gamma;
    
    % ----- Update chromosome positions -----
    % Common-mode positional jitter: the same dx_diff is added to both sisters,
    dx_diff=noise_b *(randn(1))*sqrt(dt);
    cL(t+1,1,i) = cL(t,1,i) + dt * vL(t,i)+dx_diff ;
    cR(t+1,1,i) = cR(t,1,i) + dt * vR(t,i)+dx_diff ;
    % Maintain fixed Y positions
    cL(t+1,2,i) = cL(1,2,i);
    cR(t+1,2,i) = cR(1,2,i);
    % ----- Update MT tip positions -----
    xL(t+1,i) = xL(t,i) + dt * Beta * (vL(t,i));
    xR(t+1,i) = xR(t,i) + dt * Beta * (vR(t,i));
    % ----- Update attachment numbers -----
    NL(t+1,i) = NL(t,i) + dt * (...
        n_dot(i) +(-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
        - Lambda(i)* NL(t,i));
    NR(t+1,i) = NR(t,i) + dt * (...
        n_dot(i) + (Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
        - Lambda(i)* NR(t,i));  
    end
end

% ---------------- Derived observables and plotting ----------------
% Figure1: side-by-side CM positions overtime 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2+i;
    plot(time, CMx, 'LineWidth', 2);
end
xlabel('Time (min)');
ylabel('Center of Mass X Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% ---------------- Oscillation and activity analysis ----------------
% Oscillation Analysis 
[~, ~, summary] = oscillation_measurement(cL(:,:,:), cR(:,:,:), dt, 1, false);
disp('Per-Chromosome Oscillation Summary:');
disp(summary);
% Chromosome Activity Analysis
activity_all= chromosome_activity_measurement(cL(:,:,:), cR(:,:,:), dt);
disp('Per-Chromosome Activity Summary:');
disp(activity_all);
% Figure 2: activity metrics bar plot for each chromosome 
figure;
bar_data = [activity_all.Avg_KE, activity_all.vL, activity_all.vR];
bar(bar_data, 'grouped');
legend({'Avg KE', '|v_L|', '|v_R|'});
xlabel('Chromosome Index');
ylabel('Value');
title('Chromosome Activity');
grid on;

% ---------------- inter-chromosomal correlation analysis ----------------
% Pearon's correlation analysis 
duration = 10;  % minutes
pearson_table = pearsons_correlation_analysis(cL, cR, dt, duration);
% Figure 3: distribution of Pearson's coefficients
figure; 
histogram(pearson_table.pearson_r, 20);                        % 20 bins; adjust if you like
xlabel('Pearson''s r (CM vs CM)');
ylabel('Count');
title('Distribution of Pearson''s Correlations (5 s sampled)');
grid on; box on; hold on;
xline(0, '--k');                                % reference at r=0
xline(mean(pearson_table.pearson_r,'omitnan'), 'r-', 'LineWidth', 1.5);  % mean r
hold off;

% Cross-correlation Analysis 
% Create a structure array to store results
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s         
tau_groups = reshape(tau_values, 10, [])';  % 3 groups of 10 τs each (3×10 matrix)
% Calculate CM in the x axis 
CMx = (cL(:,1,:) + cR(:,1,:)) / 2;
iteration_id=1;
% Compute correlations at selected tau values
for tau = tau_values
    [C, N] = crosscorrelation_measurement(CMx, tau, Nsteps, Nchromosomes);

    for neighbor_idx = 1:(Nchromosomes - 1)
        if N(neighbor_idx) > 0
            C_norm = C(neighbor_idx) / N(neighbor_idx);
        else
            C_norm = NaN;
        end
        records(end+1) = struct( ...
            'iteration', iteration_id, ...
            'tau', tau, ...
            'neighbor', neighbor_idx, ...
            'C_value', C_norm, ...
            'N_value', N(neighbor_idx) ...
        );
    end
end
records_table = struct2table(records);

% Figure 4: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
ylim([-1,1]);
xticks(1:1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

% Figure 5: C vs Time Lag for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:Nchromosomes-1
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2);
end
ylim([-1,1]);
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for 1st, 2nd, 3rd Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);

%% Save .csv for plotting in python 
Step   = (0:Nsteps)';            % 0,1,2,...,Nsteps 
T_min = (0:Nsteps)' * dt;        % time in minutes
T_sec = T_min * 60;              % time in seconds 

cL_x = squeeze(cL(:,1,:));      % (Nsteps+1) x Nchromosomes
cR_x = squeeze(cR(:,1,:));      % (Nsteps+1) x Nchromosomes
% Build data matrix
Chromosome = repelem((1:Nchromosomes)', numel(Step), 1);   
Step_all   = repmat(Step, Nchromosomes, 1);
Tmin_all   = repmat(T_min, Nchromosomes, 1);
Tsec_all   = repmat(T_sec, Nchromosomes, 1);
cL_col     = reshape(cL_x, [], 1);
cR_col     = reshape(cR_x, [], 1);

cLcR_tbl = table(Chromosome, Step_all, Tmin_all, Tsec_all, cL_col, cR_col, ...
    'VariableNames', {'Chromosome','Step','T_min','T_sec','cL','cR'});
writetable(cLcR_tbl, 'interdependent_chromosome_trajectory_sample01.csv');