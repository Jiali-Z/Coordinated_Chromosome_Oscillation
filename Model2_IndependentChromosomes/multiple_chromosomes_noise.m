%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multi_chromosomes_noise.m
%
% PURPOSE
%   N independent sister-chromosome pairs with small stochasticity added to:
%     (i) kinetic parameters (n_dot, Nbar) per chromosome, and
%     (ii) chromosome x-positions via TWO noise modes that drive BOTH
%         center-of-mass (COM) and inter-kinetochore (KK) jitter.
%
% What's added to the main script: 
%   • Parameter noise (cell-to-cell variability), drawn independently for
%     each chromosome i 
%   • Chromsomal positional jitter in 1D for each chromosome i 
%   • Motion is strictly 1D; y is fixed (layout uses evenly spaced y’s).
%
% OUTPUTS 
%  Figures:
%   – Side-by-side CM traces over time for all chromosomes
%   – 2D CM (x vs fixed y) trajectories
%  Command-window tables from:
%   - Oscillation summary (mean ± SD of amplitude and period for CM and KK oscillation)
%   – Chromosome activity summary: Avg_KE, Avg |v_L|, Avg |v_R|
%  Bar chart:
%   – “Chromosome Activity Metrics” (Avg_KE, Avg |v_L|, Avg |v_R|) for each
%   chromosome i 
%  Cross-correlation analysis: 
%  - Distribution of Pearson's coefficients 
%  Microrheology analysis: 
%  - Cross-correlation VS neighbor order in the first 100sec
%  - Crpss-correlation VS lag time for neighbor order of 1,2,3 
%
% DEPENDENCIES
%   oscillation_measurement.m
%   chromosome_activity_measurement.m
%   crosscorrelation_measurement.m   
%
%
% AUTHOR
%   Jiali Zhu, 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% Time Information 
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
Nchromosomes=5; %Number of chromosome pairs 

%Parameter Values for chromosome movement
noise=0.05;
noise_b=0.05;
Kct = 14.80; % pN/µm centromee spring constant, Harasymiw et al, 2019
I0 = 2; % µm Rest length of centromere spring
Kkt = 1;% pN/m kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.1; % kg/s Drag coefficient
Beta = 0.7; % Scaling factor
Nmax = 25; % Maximum number of attachments
% n_dot = 1*(1+noise*(2*rand(1)-1)); %Add noise 
% Nbar=20*(1+noise*(2*rand(1)-1)); % Steady state number of MTs when Ch is centered
% Lambda=n_dot/(Nbar);% s^-2 KMT detach rate Akioshi et al. 2010
%Alpha=n_dot*6.3/(1-Beta); 
epsilon =0.1; % Small perturbation factor

%Varying parameters for each chromoosme
Nbar = zeros(Nchromosomes,1); %The steady-state number of MTs when the Ch is centered
for i = 1:Nchromosomes
    Nbar(i) = 20*(1+noise*(2*rand(1)-1));
end
n_dot = zeros(Nchromosomes,1); %The steady-state number of MTs when the Ch is centered
for i = 1:Nchromosomes
    n_dot(i) = 1*(1+noise*(2*rand(1)-1));
end
Lambda = zeros(Nchromosomes, 1);% s^-2 KMT detach rate Akioshi et al. 2010
for i = 1:Nchromosomes
    Lambda(i)=n_dot(i)/Nbar(i);
end
Alpha = zeros(Nchromosomes, 1);
for i = 1:Nchromosomes
    Alpha(i)=n_dot(i)*6.2/(1-Beta);
end
fprintf("Ndot = %4.4f\n", n_dot);
fprintf("Nbar = %4.4f\n", Nbar);
fprintf("Lambda = %4.4f\n", Lambda);
fprintf("Alpha = %4.4f\n", Alpha);

% Set up the vectors for results
cL = zeros(Nsteps+1, 2,Nchromosomes); % Left chromosome position (time,axis,Nchromosome) 
cR = zeros(Nsteps+1, 2,Nchromosomes); % Right chromosome position (time,axis,Nchromosome) 
xL = zeros(Nsteps+1, Nchromosomes); % Left MT tip position 
xR = zeros(Nsteps+1, Nchromosomes); % Right MT tip position 
NL = zeros(Nsteps+1, Nchromosomes); % Left MT-KT attachment amount 
NR = zeros(Nsteps+1, Nchromosomes); % Right MT-KT attachment amount 
vL = zeros(Nsteps, Nchromosomes); % Left chromosome position 
vR = zeros(Nsteps, Nchromosomes); % Right chromosome position



% Initial Condition 
% Fixed y positions
y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
% Each Pairs have their own 
for i = 1: Nchromosomes 
    % Number of attachments
    NR(1,i)=Nbar(i)*(1-epsilon);
    NL(1,i)=Nbar(i)*(1+epsilon); 
    % Ch position 
    cL(1,1,i) =-I0/2-epsilon*(2*rand(1)-1);%[x,y]
    cR(1,1,i) =I0/2+epsilon*(2*rand(1)-1);%[x,y]
    cL(:,2,i) = y_positions(i);
    cR(:,2,i) = y_positions(i);
    % MT tip position 
    xL(1,i) = cL(1,1,i)-0.3;
    xR(1,i) = cR(1,1,i)+0.3;
    vL(1,i)=0;
    vR(1,i)=0;
end 

% ODE solver loop
for t = 1:Nsteps
    for i=1:Nchromosomes 
    % Force Calculations 
    % Random Force 
    F_noise = 0*noise_b *(randn(1,1)) / sqrt(dt); % Brownian Force
    % Define centromere and MT forces on the given chromosome
    F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i)); % MT forces on Left Chromosome
    F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Left Chromosome
    F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i)); % MT forces on Right Chromosome
    F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Right Chromosome

    % Calculate velocities 
    vL(t,i) = (F_KT_L + F_CT_L+F_noise)/Gamma;
    vR(t,i) = (F_KT_R + F_CT_R+F_noise)/Gamma;
    
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
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2+i;
    plot(time, CMx, 'LineWidth', 2);
end
xlabel('Time (min)');
ylabel('Center of Mass X Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

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
writetable(cLcR_tbl, 'interdependent_chromosome_trajectory_sample01.csv');

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
