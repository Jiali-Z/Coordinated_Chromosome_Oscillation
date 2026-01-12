%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiple_connected_chromosomes_noise.m
%
% PURPOSE
%   Simulate N sister-chromosome pairs that oscillate in 1D (x) and are
%   MECHANICALLY COUPLED through stochastic inter-chromosomal springs.
%
%   Each chromosome pair follows the same 1D dynamics as the single-chromosome
%   model, but chromosomes can additionally connect to each other via springs 
%   that stochastically form/break based on distance. When connected, chromosomes
%   experience coupling forces that can generate nonzero coordination. 
%
% MODEL OVERVIEW 
%   This script extends the single-chromosome force-balance oscillation model
%   to N chromosome pairs by adding two ingredients:
%
%   (A) Chromosome-to-chromosome variability 
%       - Each pair i draws its own kinetic parameters at initialization:
%         Nbar(i), n_dot(i), which set Lambda(i) and Alpha(i)
%       - These parameters remain fixed over the entire simulation
%
%   (B) Inter-chromosomal mechanical coupling (dynamic spring network)
%       - Each time step, spring_connect updates a binary adjacency matrix
%         flag_c(t) indicating which chromosome CMs are connected
%       - Connections stochastically form/break with kon_0 and koff_0, with
%         distance dependence set by the spring rest length l0
%       - When connected, chromosome i experiences an additional coupling
%         force from neighbors j:
%             F_couple,i = -k * Σ_j flag_c(i,j) * (CM_i - CM_j)
%         which is added identically to both sisters (acts on the pair CM)
%
%   All other force terms and attachment dynamics per chromosome pair follow
%   the single-chromosome model described in one_oscillating_chromosome*.m
%
% OUTPUTS
%   Figures:
%     (1) CM x-position vs time for all chromosomes (stacked/offset view)
%     (2) CM x-position vs time for all chromosomes (unshifted; raw CM traces)
%     (3) 2D CM trajectories (CM_x vs fixed y) for visualization only
%     (4) Total number of spring connections vs time
%         (sum of flag_c upper triangle; unique connections)
%     (5) Sliding-window average of total connections (5 s window)
%     (6) Chromosome activity metrics bar plot:
%         Avg_KE, Avg |v_L|, Avg |v_R| for each chromosome
%     (7) Histogram of Pearson correlation coefficients across chromosome pairs
%         (CM vs CM, sampled per pearsons_correlation_analysis settings)
%     (8) Cross-correlation vs neighbor level (averaged within τ groups)
%     (9) Cross-correlation vs lag time τ for neighbor levels 1..(Nchromosomes-1)
%
%   Command-window tables:
%     - Per-chromosome oscillation summary (CM and KK):
%       mean ± SD of amplitude and period
%     - Per-chromosome activity summary:
%       Avg_KE, Avg |v_L|, Avg |v_R|
%
%   Saved data:
%     - CSV of per-chromosome trajectories vs time:
%       (columns: Chromosome, Step, T_min, T_sec, cL, cR)
%
%
% DEPENDENCIES
%   oscillation_measurement.m
%   chromosome_activity_measurement.m
%   pearsons_correlation_analysis.m
%   crosscorrelation_measurement.m
%   spring_connect.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear; close all;

% ---------------- Time information ---------------- 
dt = 2e-3;         % time step (min)
Nsteps = 5000;     % number of simulation steps 

% ---------------- Number of chromosome ---------------- 
Nchromosomes=5; %Number of chromosome pairs 

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
% ---------------- Inter-chromosome coupling parameters ----------------
l0=3;              % µm     inter-chromosomal spring effective length 
koff_0=40;         % 1/s    base connection rate
kon_0=20;          % 1/s    base disconnection rate
Nspring = 1;       % Number of springs per chromosome not used 
k = 0.5;           % pN/µm  coupling strength 

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
    cL(1,1,i) =-I0/2-epsilon*(2*rand(1)-1);  % left sister initial x-position
    cR(1,1,i) =I0/2+epsilon*(2*rand(1)-1);   % right sister initial x-position
    cL(:,2,i) = y_positions(i);              % y fixed
    cR(:,2,i) = y_positions(i);              % y fixed
    % MT tip position 
    xL(1,i) = cL(1,1,i)-0.3;
    xR(1,i) = cR(1,1,i)+0.3;
    vL(1,i)=0;
    vR(1,i)=0;
end 

% ODE solver loop
for t = 1:Nsteps
    %Update connection state 
    CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);  % N x 1 vector
    YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);  % N x 1 vector
    [flag_c,dist]=spring_connect(CM,YM,flag_cp,Nchromosomes,dt,l0,koff_0,kon_0);
    flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
    flag_cp = flag_c;  % Update prior state
    for i=1:Nchromosomes 
        % Force Calculations 
        %Coupling force 
        F_coupling_x=0;
        for j = 1:Nchromosomes
            if j ~= i  % Skip self
                delta_x = CM(i) - CM(j);
                F_coupling_x = F_coupling_x - k * delta_x * flag_c(i, j);
            end
        end
        % Random Force 
        F_noise = 0*noise_b *(randn(1,1)) / sqrt(dt); % Brownian Force
        % Define centromere and MT forces on the given chromosome
        F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i)); % MT forces on Left Chromosome
        F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Left Chromosome
        F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i)); % MT forces on Right Chromosome
        F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Right Chromosome
    
        % Calculate velocities 
        vL(t,i) = (F_KT_L + F_CT_L+F_noise+F_coupling_x)/Gamma;
        vR(t,i) = (F_KT_R + F_CT_R+F_noise+F_coupling_x)/Gamma;
        
        % Update Positions 
        % Calculate the current chromosome position in both directions
        dx_diff=noise_b *(randn(1))*sqrt(dt);
        cL(t+1,1,i) = cL(t,1,i) + dt * vL(t,i)+dx_diff ;
        cR(t+1,1,i) = cR(t,1,i) + dt * vR(t,i)+dx_diff ;
        % Maintain fixed Y positions
        cL(t+1,2,i) = cL(1,2,i);
        cR(t+1,2,i) = cR(1,2,i);
        % Calculate the current microtubule tip position in both directions
        xL(t+1,i) = xL(t,i) + dt * Beta * (vL(t,i));
        xR(t+1,i) = xR(t,i) + dt * Beta * (vR(t,i));
        % Update number of attachments 
        NL(t+1,i) = NL(t,i) + dt * (...
            n_dot(i) +(-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
            - Lambda(i)* NL(t,i));
        NR(t+1,i) = NR(t,i) + dt * (...
            n_dot(i) + (Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
            - Lambda(i)* NR(t,i));  
    end
end

%Draw afterwards
% side-by-side CM positions overtime 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = ((cL(:,1,i) + cR(:,1,i)) / 2)+i;
    plot(time, CMx, 'LineWidth', 2);
end
xlabel('Time (min)');
ylabel('Center of Mass X Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

%Draw afterwards
% side-by-side CM positions overtime 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    plot(time, CMx, 'LineWidth', 2);
end
xlabel('Time (min)');
ylabel('Center of Mass X Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% Movement Trajectories 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    CMy = (cL(:,2,i) + cR(:,2,i)) / 2;
    plot(CMx, CMy, 'LineWidth', 2);
end
xlabel('Center of Mass X Position (µm)');
ylabel('Center of Mass Y Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% Plot total number of connections over time
figure;
plot((1:Nsteps)*dt, flag_c_sum_accumulated, 'LineWidth', 2);
xlabel('Time (min)', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Total Connections', 'FontSize', 18, 'FontName', 'Arial');
title('Number of Chromosome Connections Over Time', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Arial');
grid on;


% Plot sliding average of total number of connections over time
figure;
window_size = round(5 / (dt * 60));  % steps per 5-second interval
flag_c_sliding_avg = movmean(flag_c_sum_accumulated, window_size);
% Plot sliding average over time
plot((1:Nsteps)*dt, flag_c_sliding_avg, 'LineWidth', 2, 'Color', 'b');
xlabel('Time (min)', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Total Connections (5s average)', 'FontSize', 18, 'FontName', 'Arial');
title('Number of Chromosome Connections (Sliding Window Average)', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Arial');
ylim([0,12]);
grid on;
hold on; 
% % Movement Trajectories 
% figure;
% hold on;
% time = (0:Nsteps) * dt;
% for i = 1:Nchromosomes
%     CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
%     CMy = (cL(:,2,i) + cR(:,2,i)) / 2;
%     plot(CMx, CMy, 'LineWidth', 2);
% end
% xlabel('Center of Mass X Position (µm)');
% ylabel('Center of Mass Y Position (µm)');
% title('Chromosome CM Oscillations');
% legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
% set(gca, 'FontSize', 14);
% grid on;

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
writetable(cLcR_tbl, 'activity_fit_strong_spring_chromosome_trajectory.csv');
%% Analysis 
% Oscillation Analysis 
[~, ~, summary] = oscillation_measurement(cL(:,:,:), cR(:,:,:), dt, 1, true);
disp('Per-Chromosome Oscillation Summary:');
disp(summary);
% Chromosome Activity Analysis
activity_all= chromosome_activity_measurement(cL(:,:,:), cR(:,:,:), dt);
disp('Per-Chromosome Activity Summary:');
disp(activity_all);

figure;
bar_data = [activity_all.Avg_KE, activity_all.vL, activity_all.vR];
bar(bar_data, 'grouped');
legend({'Avg KE', '|v_L|', '|v_R|'});
xlabel('Chromosome Index');
ylabel('Value');
title('Chromosome Activity');
grid on;

%% Pearon's correlation analysis 
duration = 10;  % minutes
pearson_table = pearsons_correlation_analysis(cL, cR, dt, duration);
% Plot the distribution of Pearson's coefficients
figure; 
histogram(pearson_table.pearson_r, 20);                        % 20 bins; adjust if you like
xlabel('Pearson''s r (CM vs CM)');
ylabel('Count');
title('Distribution of Pearson''s Correlations (5 s sampled)');
grid on; box on; hold on;
xline(0, '--k');                                % reference at r=0
xline(mean(pearson_table.pearson_r,'omitnan'), 'r-', 'LineWidth', 1.5);  % mean r
hold off;

%% Cross-correlation Analysis 
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

% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
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

% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
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

