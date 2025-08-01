%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscorrelation_measurement_null_test_iteration.M
%
% PURPOSE:
%   This script runs the **null model** for chromosome oscillation for 
%   **multiple independent iterations**, to assess the averaged behavior of 
%   the cross-correlation metric 
%
%   Compared to the single-run version, this script repeats the simulation
%   for multiple valid rounds (`target_iterations`) and accumulates statistics
%   to better quantify expected variability in a non-interacting system.
%
% Model Detail:
%   - Each particle experiences a deterministic restoring force (-Kx)
%     and stochastic Brownian noise in the x-direction.
%   - The y-position is fixed. 
%   - No inter-particle interactions 
%   - The simulation is repeated for N iterations, excluding any run with NaNs.
%   - Cross-correlations are computed for each run, at multiple lag times (τ),
%     and averaged across valid iterations.
%
% OUTPUTS:
%   - Two plots:
%       (1) Averaged normalized correlation vs. neighbor level (grouped by τ)
%       (2) Averaged normalized correlation vs. τ (for neighbor levels 1–4)
%   - Summary table `records_table` containing:
%       iteration, tau, neighbor level, accumulated correlation (C_value), and count (N_value)
%
% DEPENDENCIES:
%   - Requires `crosscorrelation_measurement.m` on the MATLAB path
%
% PARAMETERS:
%   - Nchromosomes      : Number of particles (default = 5)
%   - dt                : Time step in minutes (default = 2e-3)
%   - Nsteps            : Number of time steps per iteration (default = 5000)
%   - K                 : Harmonic spring constant
%   - Gamma             : Drag coefficient
%   - noise_strength    : Brownian noise magnitude
%   - target_iterations : Number of valid simulation runs to perform
%
% NOTE:
%   - τ values range from 5 to 150 seconds in 5-second steps
%   - τ values are grouped into 3 bins of 10 τ values each for plotting
%   - Correlation values are normalized by the number of contributing terms
%
% AUTHOR:
%   Jiali Zhu, 2025/08/01, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;
% Brownian oscillator simulation
Nchromosomes = 5;         % Number of particles
dt = 2e-3;                % Time step
Nsteps = 5000;            % Number of steps
K = 1;                    % Spring constant
Gamma = 0.01;                % Drag coefficient
noise_strength = 0.1;     % Noise magnitude
target_iterations = 100; % Number of valid rounds of data needed
valid_iterations = 0; % Counter for valid rounds

% Create a structure array to store results
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s       
tau_groups = reshape(tau_values, 10, [])';  % 3 groups of 10 τs each (3×10 matrix)

while valid_iterations < target_iterations
   % Simulate trajectories
    x = zeros(Nsteps, Nchromosomes);
    y = zeros(Nsteps, Nchromosomes);
    x(1,:) = randn(1, Nchromosomes);  % initial condition
    y(1,:) = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
    for t = 1:Nsteps-1
        % x-direction: harmonic potential + noise
        F_det_x = -K * x(t, :);
        F_noise_x = noise_strength * randn(1, Nchromosomes) / sqrt(dt);
        dx = dt * (F_det_x + F_noise_x) / Gamma;
        x(t+1, :) = x(t, :) + dx;
        % y-direction: Brownian noise only 
        F_noise_y = 0*noise_strength * randn(1, Nchromosomes) / sqrt(dt);
        %F_noise_y = 0; 
        % Apply Brownian step
        dy = dt * F_noise_y / Gamma;
        y(t+1, :)= y(t, :) + dy;
    end
    
    % Check for NaNs in cL_x, cL_y, cR_x, cR_y
    if any(isnan(x(:)))
        fprintf('NaN detected, skipping this round and starting a new one\n');
        continue; % Skip this round and start a new iteration
    end

    % Cross-correlation Calculation
    % Cross-correlation versus Tau_group 
    CM = x;
    iteration_id = valid_iterations + 1;
    % Compute correlations at selected tau values
    for tau = tau_values
        [C, N] = crosscorrelation_measurement(CM, tau, Nsteps, Nchromosomes);
    
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
    
    % If no NaNs, proceed with correlation calculationcompu
    fprintf('Iteration %d complete\n', valid_iterations + 1);
    valid_iterations = valid_iterations + 1; % Increment valid rounds counter
end

records_table = struct2table(records);

%% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value, grouped.std_C_value,'-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 3, 'Color', colors(g, :));
end
ylim([-0.1,0.1]);
xticks(1:1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

%% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3, 4
figure; hold on;
for n = 1:Nchromosomes-1
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 3);
end
ylim([-0.1,0.1])
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for 1st, 2nd, 3rd, 4th Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);





