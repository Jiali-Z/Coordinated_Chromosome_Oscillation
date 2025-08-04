%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscorrelation_measurement_positive_Test_iteration.M
%
% PURPOSE:
%   This script simulates a **positive control model** for chromosome 
%   oscillations for **multiple independent iterations**, to assess the 
%   averaged behavior of the cross-correlation metric 
%
%   Compared to the single-run version, this script repeats the simulation
%   for multiple valid rounds (`target_iterations`) and aggregate
%   correlation data across iterations to assess average treads 
%
% Model Detail:
%   - Each chromosome experiences:
%       (1) Harmonic potential: F= -K * x
%       (2) Inter-chromosomal spring force: to adjacent neighbors
%       (3) Brownian noise
%   - The y-position is fixed or unused.
%   - Forces are tracked over time for analysis and visualization.
%   - Cross-correlations are computed over a range of τ values (5–150s),
%   grouped by neighbor level. 
%
% OUTPUTS:
%   - Two plots:
%       (1) Avg cross-correlation vs. neighbor level (grouped by τ)
%       (2) Avg cross-correlation vs. τ (for each neighbor level)
%- `records_table`: correlation values for all τ and neighbor levels
%
% DEPENDENCIES:
%   - Requires `crosscorrelation_measurement.m` on the MATLAB path
%
% PARAMETERS:
%   - Nchromosomes   : Number of particles (default = 5)
%   - dt             : Time step (minutes)
%   - Nsteps         : Number of simulation steps
%   - K              : Harmonic spring constant
%   - K_inter        : Inter-chromosome spring constant (controls coupling strength)
%   - Gamma          : Drag coefficient
%   - noise_strength : Noise amplitude
%   - target_iterations : Number of successful simulation runs to accumulate
%
% AUTHOR:
%   Jiali Zhu, 2025/08/01, UNC-CHapel Hill 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
% Brownian oscillator simulation
Nchromosomes = 5;         % Number of particles
dt = 2e-3;                % Time step
Nsteps = 5000;            % Number of steps
K = 1;                    % Spring constant
Gamma = 0.01;                % Drag coefficient
noise_strength = 0.1;     % Noise magnitude
% Add inter-particle springs 
K_inter = 1;  % Inter-particle spring constant 
target_iterations = 100; % Number of valid rounds of data needed
valid_iterations = 0; % Counter for valid rounds

% Create a structure array to store results
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s       
tau_groups = reshape(tau_values, 10, [])';  % 3 groups of 10 τs each (3×10 matrix)

while valid_iterations < target_iterations
    % Intialization 
    x = zeros(Nsteps, Nchromosomes);
    y = zeros(Nsteps, Nchromosomes);
    x(1,:) = randn(1, Nchromosomes);  % initial condition
    y(1,:) = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
    % Preallocate force tracking ===
    F_harmonic_all = zeros(Nsteps, Nchromosomes);
    F_inter_all    = zeros(Nsteps, Nchromosomes);
    F_noise_all    = zeros(Nsteps, Nchromosomes);
    for t = 1:Nsteps-1
        % x-direction: harmonic potential + noise
        F_harmonic = -K * x(t, :); %harmonic force from potential 
        % Inter-chromosome spring coupling (nearest neighbors)
        F_inter = zeros(1, Nchromosomes);
        for i = 1:Nchromosomes
            if i > 1
                F_inter(i) = F_inter(i) + K_inter * (x(t, i-1) - x(t, i));
            end
            if i < Nchromosomes
                F_inter(i) = F_inter(i) + K_inter * (x(t, i+1) - x(t, i));
            end
        end
        % Brownian Force
        F_noise= noise_strength * randn(1, Nchromosomes) / sqrt(dt);
        dx = dt * (F_harmonic + F_inter+F_noise) / Gamma;
        x(t+1, :) = x(t, :) + dx;
        % y-direction: Brownian noise only 
        dy = dt * 0 * randn(1, Nchromosomes) / Gamma;  % no movement in y
        y(t+1, :)= y(t, :) + dy;
        %Store force
        F_harmonic_all(t, :) = F_harmonic;
        F_inter_all(t, :)    = F_inter;
        F_noise_all(t, :)    = F_noise;
    end
    % Check for NaNs in cL_x, cL_y, cR_x, cR_y
    if any(isnan(x(:)))
        fprintf('NaN detected, skipping this round and starting a new one\n');
        continue; % Skip this round and start a new iteration
    end
    
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
                'C_value', C_norm, ... % normalied crosscorrelation value 
                'N_value', N(neighbor_idx) ... % number of N used to calculate C_norm 
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
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");% Average from different tau groups

    errorbar(grouped.neighbor, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
ylim([-0.3,0.3]);
xticks(1:1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

%% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:Nchromosomes-1
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");  

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2);
end
ylim([-0.3,0.3])
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for 1st, 2nd, 3rd Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);




