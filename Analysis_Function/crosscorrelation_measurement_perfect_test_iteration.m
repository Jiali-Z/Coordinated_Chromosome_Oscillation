%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscorrelation_measurement_perfect_control.M
%
% PURPOSE:
%   This script generates a **perfect correlation control** where all
%   chromosomes follow the *exact same trajectory*, and thus their 
%   normalized cross-correlation should approach 1 across all τ and 
%   neighbor levels.
%
%   The output allows direct comparison with noisy stochastic simulations.
%
% MODEL:
%   - All chromosomes follow an identical Brownian-harmonic trajectory.
%   - No inter-chromosomal interaction is needed — the motion is cloned.
%
% OUTPUTS:
%   - Same as positive control:
%       (1) C vs Neighbor Level (grouped by τ)
%       (2) C vs τ (for Neighbor Levels 1–3)
%   - `records_table`: cross-correlation results
%
% AUTHOR:
%   Jiali Zhu, August,2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Parameters
Nchromosomes = 5;
dt = 2e-3;
Nsteps = 5000;
K = 1;
Gamma = 0.01;
noise_strength = 0.1;
tau_step = round(5 / (dt * 60));  % steps per 5 sec
tau_values = tau_step * (1:30);   % From 5s to 150s       
tau_groups = reshape(tau_values, 10, [])';  % 3x10 groups

% Generate shared trajectory
x_shared = zeros(Nsteps, 1);
x_shared(1) = randn();  % initial condition

for t = 1:Nsteps-1
    F_det = -K * x_shared(t);  % harmonic force
    F_noise = noise_strength * randn() / sqrt(dt);  % Brownian
    dx = dt * (F_det + F_noise) / Gamma;
    x_shared(t+1) = x_shared(t) + dx;
end

% Duplicate shared trajectory across chromosomes
x = repmat(x_shared, 1, Nchromosomes);

% === Cross-correlation measurement ===
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
iteration_id = 1;

for tau = tau_values
    [C, N] = crosscorrelation_measurement(x, tau, Nsteps, Nchromosomes);

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

%% Plot 1: Avg Cross-Correlation vs Neighbor Level (Grouped by τ)
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value, grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, ...
                                             (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
xticks(1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('Perfect Correlation: C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

%% Plot 2: C vs τ for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:min(Nchromosomes-1, 3)
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");

    errorbar((grouped.tau / tau_step) * 5, grouped.mean_C_value, grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2);
end
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('Perfect Correlation: C vs \tau (Neighbors 1–3)', 'FontSize', 16);
legend('Location', 'best'); grid on; box on;
set(gca, 'FontSize', 12);

%% Plot 3: Trajectories of All Chromosomes (should perfectly overlap)
figure;
plot((1:Nsteps) * dt * 60, x, 'LineWidth', 1.5);  % time in seconds
xlabel('Time (s)');
ylabel('Position (x)');
title('Perfect Correlation: Identical Chromosome Trajectories');
legend(arrayfun(@(i) sprintf('Ch %d', i), 1:Nchromosomes, 'UniformOutput', false));
grid on; box on;
set(gca, 'FontSize', 12);

%% Plot 4: Overlay One Trajectory + Differences (Should be Zero)
figure;
plot((1:Nsteps) * dt * 60, x(:,1), 'k-', 'LineWidth', 2); hold on;
for i = 2:Nchromosomes
    plot((1:Nsteps) * dt * 60, x(:,i) - x(:,1), '--', 'LineWidth', 1.2);
end
xlabel('Time (s)');
ylabel('Deviation from x_1');
title('Deviation from Shared Trajectory (Should be Zero)');
legend(["x_1", arrayfun(@(i) sprintf('x_%d - x_1', i), 2:Nchromosomes, 'UniformOutput', false)]);
grid on; box on;
set(gca, 'FontSize', 12);