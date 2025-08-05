%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscorrelation_measurement_positive_test_iteration_kinter.M
%
% PURPOSE:
%   This script extends the multiple-iteration **positive control model** 
%   by sweeping across a range of inter-chromosomal spring constants (K_inter)
%   to assess how the strength of inter-chromosome interactions influences 
%   the normalized cross-correlation between chromosome movements.
%
%   It builds on the previous model where chromosomes were coupled via
%   spring-like forces and simulated with Brownian noise. The goal here is 
%   to systematically test whether increasing K_inter leads to higher 
%   pairwise correlation, as measured by `crosscorrelation_measurement`.
%
%
% MODEL DETAILS:
%   - Each chromosome experiences:
%       (1) A harmonic restoring force:        F = -K * x
%       (2) Inter-chromosomal spring coupling: to adjacent neighbors (via K_inter)
%       (3) Brownian noise
%   - Motion is simulated in 1D (x-axis); y is fixed and unused
%   - No movement is simulated in y-direction (dy = 0)
%
% OUTPUTS:
%   - `all_results(k_idx).records_table`: table of normalized cross-correlation
%     for each value of K_inter across all τ and neighbor levels
%   - Final figure: Overlay plot showing C vs. neighbor level for different K_inter
%     using the middle τ group (55–100 s)
%
% DEPENDENCIES:
%   - Requires `crosscorrelation_measurement.m` on MATLAB path
%
% PARAMETERS:
%   - Nchromosomes      : Number of chromosomes (default = 5)
%   - dt                : Time step (minutes)
%   - Nsteps            : Number of simulation steps per run
%   - K                 : Harmonic potential spring constant
%   - K_inter_values    : Array of K_inter values to sweep (e.g., 0:0.1:1)
%   - Gamma             : Drag coefficient
%   - noise_strength    : Brownian noise amplitude
%   - target_iterations : Number of valid simulation runs per K_inter value
%
% NOTES:
%   - This script provides a **quantitative test** of whether cross-correlation
%     increases with coupling strength, and how strongly neighbor level influences
%     the effect of K_inter.
%   - Sweeping restuls showed that increasing k_inter **does NOT** always
%     lead to higher correlation. In fact, beyon a certain point (1.0),
%     stronger inter-chromosomal springs can "suppress" cross-correlation. 
%   - There is a **nonlinear interplay** between k_inter(inter-particle coupling 
%     strength) and k(harmonic potential strength anchoring each particle)
%   - When K_inter is too large relative to K, the system becomes overly
%     stiff and loses flexibility for correlated fluctuations - resulting 
%     in reduced cross-correlation 
%   - Those observation motivate the next step: a **2D parameter sweep**
%     over both K and k_inter to identify where corr-correlation peaks 
%
% AUTHOR:
%   Jiali Zhu, 2025/08/01, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% === Parameters ===
Nchromosomes = 5; % Number of particles
dt = 2e-3;        % Time step
Nsteps = 5000;    % Number of steps
K = 1;            % Spring constant
Gamma = 0.01;     % Drag coefficient
noise_strength = 0.1; % Noise magnitude 
target_iterations = 100; % Number of valid rounds of data needed

% Create a structure array to store results
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s       
tau_groups = reshape(tau_values, 10, [])';  % 3×10

% === Sweep over K_inter ===
K_inter_values = 0:0.5:2;
all_results = struct();

for k_idx = 1:length(K_inter_values)
    K_inter = K_inter_values(k_idx);
    fprintf('\nRunning K_inter = %.1f...\n', K_inter);

    valid_iterations = 0;
    records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});

    while valid_iterations < target_iterations
        % === Initialization ===
        x = zeros(Nsteps, Nchromosomes);
        y = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
        x(1,:) = randn(1, Nchromosomes);

        for t = 1:Nsteps-1
            F_harmonic = -K * x(t, :);
            F_inter = zeros(1, Nchromosomes);
            for i = 1:Nchromosomes
                if i > 1
                    F_inter(i) = F_inter(i) + K_inter * (x(t, i-1) - x(t, i));
                end
                if i < Nchromosomes
                    F_inter(i) = F_inter(i) + K_inter * (x(t, i+1) - x(t, i));
                end
            end
            F_noise = noise_strength * randn(1, Nchromosomes) / sqrt(dt);
            dx = dt * (F_harmonic + F_inter + F_noise) / Gamma;
            x(t+1, :) = x(t, :) + dx;
        end

        if any(isnan(x(:)))
            continue;
        end

        % === Cross-correlation ===
        CM = x;
        iteration_id = valid_iterations + 1;
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
                    'N_value', N(neighbor_idx));
            end
        end
        valid_iterations = valid_iterations + 1;
    end

    % === Save Table ===
    all_results(k_idx).K_inter = K_inter;
    all_results(k_idx).records_table = struct2table(records);
end
% === Overlay Plot: Avg C vs Neighbor Level (for middle τ group) ===
figure; hold on;
colors = lines(length(K_inter_values));
g = 2;  % middle τ group

for k_idx = 1:length(K_inter_values)
    K_inter = K_inter_values(k_idx);
    records_table = all_results(k_idx).records_table;
    current_group = tau_groups(g, :);

    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean", "std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value, grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('K_{inter} = %.1f', K_inter), ...
        'LineWidth', 2, 'Color', colors(k_idx, :));
end

xlabel('Neighbor Level');
ylabel('Avg Norm Cross-Corr');
title('Overlay: C vs Neighbor Level at τ = 55–100 s');
legend('Location', 'best'); grid on; box on;
set(gca, 'FontSize', 14);
xticks(1:Nchromosomes-1);
ylim([-0.3, 0.3]);
