clc; clear; close all;

% Simulation settings
N_iter = 100;              % Number of cells
dt = 2e-3;                % Timestep (min)
Nsteps = 5000;            % Number of time steps
results = [];             % Store summary per iteration

for iter = 1:N_iter
    % Fixed noise value for this "cell"
    noise = 0.15; %* (1 + 0.1*(2*rand - 1)); % e.g. ±40% variation

    % Parameters
    noise_b = 0.2; % Brownian noise is constant across iterations (still random within)
    Kct = 12.3; % pN/µm centromee spring constant, Harasymiw et al, 2019
    I0 = 2; % µm Rest length of centromere spring
    Kkt = 1; % pN/m kinetocore spring constant, Cojoc et al. 2016
    Gamma = 0.1; % kg/s Drag coefficient
    Beta = 0.7;% Scaling factor
    Nmax = 25;% Maximum number of attachments
    n_dot = 1*(1 + noise*(2*rand(1) - 1));
    Nbar = 20*(1 + noise*(2*rand(1) - 1));% Steady state number of MTs when Ch is centered
    Lambda = n_dot / (Nbar);% s^-2 KMT detach rate Akioshi et al. 2010
    Alpha = n_dot * 6.2 / (1 - Beta);
    epsilon = 0.1; % Small perturbation factor

    % Initialize arrays
    xL = zeros(Nsteps, 1); xR = zeros(Nsteps, 1);
    NL = zeros(Nsteps, 1); NR = zeros(Nsteps, 1);
    cL = zeros(Nsteps+1, 2, 1); cR = zeros(Nsteps+1, 2, 1);
    vL = zeros(Nsteps, 1); vR = zeros(Nsteps, 1);

    % Initial conditions
    NR(1) = Nbar*(1 - epsilon);
    NL(1) = Nbar*(1 + epsilon);
    cL(1,:,1) = [-I0/2 - epsilon*(2*rand(1) - 1), 0];
    cR(1,:,1) = [ I0/2 + epsilon*(2*rand(1) - 1), 0];
    xbar = (Kct * I0) / (Nbar * Kkt);
    xL(1) = cL(1,1,1) - 0.3;
    xR(1) = cR(1,1,1) + 0.3;
    vL(1) = 0; vR(1) = 0;

    % ODE loop
    for t = 1:Nsteps
        F_noise = 0 * noise_b * randn(1,1) / sqrt(dt);  % Set to 0 if not using
        % Define centromere and MT forces on the given chromosome
        F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));
        F_CT_L = Kct * (cR(t,1,1) - cL(t,1,1) - I0);
        F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));
        F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0);
        % Calculate velocities 
        vL(t) = (F_KT_L + F_CT_L + F_noise) / Gamma;
        vR(t) = (F_KT_R + F_CT_R + F_noise) / Gamma;
        % Calculate the current chromosome position in both directions
        dx_diff = noise_b * randn(1,1) * sqrt(dt);
        cL(t+1,:,1) = cL(t,:,1) + dt * vL(t) + dx_diff;
        cR(t+1,:,1) = cR(t,:,1) + dt * vR(t) + dx_diff;
         % Calculate the current microtubule tip position in both directions
        xL(t+1) = xL(t) + dt * Beta * vL(t);
        xR(t+1) = xR(t) + dt * Beta * vR(t);
        % Update number of attachments 
        NL(t+1) = NL(t) + dt * (...
            n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
            - Lambda * NL(t));
        NR(t+1) = NR(t) + dt * (...
            n_dot + (Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
            - Lambda * NR(t));
    end

    % Oscillation analysis
    [~, ~, table_summary] = oscillation_measurement(cL, cR, dt, iter, false);
    activity_table = chromosome_activity_measurement(cL, cR, dt);

    % Add Avg_KE, vL, vR from activity_table to each row of table_summary
    table_summary.Avg_KE = repmat(activity_table.Avg_KE, height(table_summary), 1);
    table_summary.vL     = repmat(activity_table.vL,     height(table_summary), 1);
    table_summary.vR     = repmat(activity_table.vR,     height(table_summary), 1);
    
    % Append to results
    results = [results; table_summary];
end

%% Summary Statistics and Boxplots for All 5 Measurements

% Summary Statistics and Boxplots for All 5 Measurements
fprintf('\nSummary Statistics Across Simulated Cells:\n');
fprintf('%-20s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n', ...
    'Metric', 'Mean', 'Std', 'Min', '25%', 'Median', '75%', 'Max');

% Extract unique iterations
unique_iters = unique(results.iteration);
N_iter = numel(unique_iters);

% Define metrics and values
metric_defs = {
    'CM Amplitude', results.mean_amplitude(strcmp(results.metric, 'CM'));
    'CM Period',    results.mean_period(strcmp(results.metric, 'CM'));
    'KK Amplitude', results.mean_amplitude(strcmp(results.metric, 'KK'));
    'KK Period',    results.mean_period(strcmp(results.metric, 'KK'));
    'Average KE',   results.Avg_KE(strcmp(results.metric, 'CM'));  % KE per chromosome
};

% Loop through metrics and display stats + plot boxplot
for i = 1:size(metric_defs,1)
    label = metric_defs{i,1};
    raw_values = metric_defs{i,2};

    if strcmp(label, 'Average KE')     
        KE_per_iter = results.Avg_KE(strcmp(results.metric, 'CM'));
        % Pad with 0s if any iterations missing
        if numel(KE_per_iter) < N_iter
            KE_per_iter(end+1:N_iter) = 0;
        end
        
        avg_ke_normalized = mean(KE_per_iter);
        std_ke = std(KE_per_iter);
        fprintf('%-20s %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f\n', ...
            label, mean(KE_per_iter),  std(KE_per_iter), min(KE_per_iter), ...
            prctile(KE_per_iter,25), median(KE_per_iter), ...
            prctile(KE_per_iter,75), max(KE_per_iter));

        % Boxplot for KE
        values_to_plot = KE_per_iter;
    else
        fprintf('%-20s %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f\n', ...
            label, mean(raw_values), std(raw_values), min(raw_values), ...
            prctile(raw_values,25), median(raw_values), ...
            prctile(raw_values,75), max(raw_values));

        values_to_plot = raw_values;
    end

    % Boxplot
    figure;
    boxplot(values_to_plot);
    title(label, 'FontSize', 16);
    ylabel(label, 'FontSize', 14);
    set(gca, 'FontSize', 14);
    set(gcf, 'Color', 'w');
    grid on;
end

