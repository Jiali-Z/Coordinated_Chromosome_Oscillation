%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscorrelation_measurement_null.M
%
% PURPOSE:
%   This scipr impements a **null model** for chromosome osciilation, 
%   designed to test and validate the 'crosscorrelation_measurement'
%   function
% 
%   Simulate the dynamics of multiple Brownian particles under a harmonic 
%   potential. There are no inter-particle particles - any obsered
%   correlation should arise purely from random coincidence, providing a
%   baseline for interpreting cross-correlation values in the actual
%   simulation. 
%   
%
% Model Detail:
%   - Each particle experiences a deterministic restoring force (-Kx)
%     and stochastic Brownian noise in the x-direction.
%   - The y-position is fixed 
%   - Cross-correlation is calculated using a custom function 
%     `crosscorrelation_measurement.m` for multiple time lags.
%
% OUTPUTS:
%   - Two plots:
%       (1) Normalized cross-correlation vs. neighbor level, grouped by τ
%       (2) Normalized cross-correlation vs. τ, for each neighbor level
%   - Summary table `records_table` containing:
%       iteration, tau, neighbor level, accumulated correlation (C_value), and count (N_value)
%
% DEPENDENCIES:
%   - Requires `crosscorrelation_measurement.m` in the same directory
%
% PARAMETERS (can be tuned):
%   - Nchromosomes  : Number of particles (default = 5)
%   - dt            : Time step in minutes (default = 2e-3)
%   - Nsteps        : Number of time steps (default = 5000)
%   - K             : Spring constant (harmonic trap strength)
%   - noise_strength: Magnitude of Brownian noise
%   - tau_step      : Step size in τ (computed from 5s intervals)
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

figure; 
plot(x, y, 'LineWidth', 1.5);
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
title('Brownian Particles in Harmonic Potential', 'FontSize', 16);
grid on;
set(gca, 'LineWidth', 2, 'FontSize', 14);

% Time vector matches the number of steps
time = (0:Nsteps-1) * dt;

% Plot CMx (oscillation) over time for each chromosome
figure;
hold on;
for i = 1:Nchromosomes
    plot(time, x(:, i), 'LineWidth', 2);  % Plot one chromosome's x-position
end
xlabel('Time (min)');
ylabel('X Position (µm)');
title('Brownian Particle X Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;


%% Cross-correlation Calculation
% Create a structure array to store results
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s       
tau_groups = reshape(tau_values, 10, [])';  % 3 groups of 10 τs each (3×10 matrix)
% Cross-correlation versus Tau_group 
CM = x;
iteration_id = 1; 
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
            'C_value', C_norm, ... % normalied crosscorrelation value 
            'N_value', N(neighbor_idx) ... % number of N used to calculate C_norm 
        );
    end
end

records_table = struct2table(records);
fprintf('Final table size: %d rows\n', height(records_table));

%Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");% Average from different tau groups

    errorbar(grouped.neighbor, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
ylim([-0.1,0.1]);
xticks(1:1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

%Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:Nchromosomes-1
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");  

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2);
end
ylim([-0.1,0.1])
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for 1st, 2nd, 3rd Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);




