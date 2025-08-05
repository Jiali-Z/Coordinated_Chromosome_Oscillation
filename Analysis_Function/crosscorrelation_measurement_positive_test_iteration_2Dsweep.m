%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosscorrelation_measurement_positive_test_iteration_2Dsweep.m
%
% PURPOSE:
%   Perform a 2D parameter sweep across:
%       (1) K: harmonic potential strength
%       (2) K_inter: inter-particle spring coupling
%   The goal is to determine how these two parameters jointly influence
%   the 1st-neighbor normalized cross-correlation in a system of 
%   Brownian-coupled oscillators.
%
%   This is an extension of the null positive control model of chromosome
%   oscillation. Previous simulations showed that increasing `K_inter`
%   enhances coordination to a point, after which too strong coupling
%   can suppress correlation due to overly rigid connections.
%
%   This sweep helps identify the balance between confinement (K) and 
%   interconnection (K_inter) that maximizes correlated movement.
%
% MODEL DETAILS:
%   - Nchromosomes move in 1D under:
%       • Harmonic potential:     F = -K*x
%       • Inter-particle springs: nearest-neighbor coupling (K_inter)
%       • Brownian force:         random noise
%   - For each parameter pair, run multiple iterations and average
%     the normalized cross-correlation of the 1st neighbors only.
%
% OUTPUT:
%   - A smoothed 3D surface plot:
%         Z = mean normalized C (1st neighbor)
%         vs K and K_inter
%   - Peak (maximum) correlation value and its location
% Note: 
%   - Max correlation = 2.2187 at K = 0.10, K_inter = 0.50
%
% DEPENDENCY:
%   - Requires: crosscorrelation_measurement.m on MATLAB path
%
% AUTHOR:
%   Jiali Zhu, August 2025, UNC-Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc; clear; close all;

% === Parameter Grid ===
K_values       = [0.1,0.5,1, 2];
K_inter_values = 0:0.5:2;
target_iterations = 50;  % reduced for faster sweep

% === System Settings ===
Nchromosomes = 5;        % Number of particles        
dt = 2e-3;               % Time step
Nsteps = 5000;           % reduce to 2000 for speed
Gamma = 0.01;            % Drag coefficient
noise_strength = 0.1;    % Noise magnitude 
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s

% === Output Storage ===
mean_corr_surface = NaN(length(K_values), length(K_inter_values));

% === Main Loop ===
for ki = 1:length(K_values)
    K = K_values(ki);
    for kj = 1:length(K_inter_values)
        K_inter = K_inter_values(kj);
        fprintf('K = %.2f, K_inter = %.2f\n', K, K_inter);

        C_vals_first_neighbor = [];
        valid_iterations = 0;

        while valid_iterations < target_iterations
            % Initialize position
            x = zeros(Nsteps, Nchromosomes);
            x(1,:) = randn(1, Nchromosomes);

            for t = 1:Nsteps-1
                F_harmonic = -K * x(t,:);
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
                x(t+1,:) = x(t,:) + dx;
            end

            if any(isnan(x(:))), continue; end

            % Cross-correlation
            C_sum_first_neighbor = 0; count=0;
            for tau = tau_values
                [C, N] = crosscorrelation_measurement(x, tau, Nsteps, Nchromosomes);
                if N(1) >0
                    C1=C(1)./N(1);
                    C_sum_first_neighbor=C_sum_first_neighbor+C1;
                    count = count + 1;
                end
            end

            if count > 0
                C_vals_first_neighbor(end+1) = C_sum_first_neighbor / count;
                valid_iterations = valid_iterations + 1;
            end
        end

        mean_corr_surface(ki, kj) = mean(C_vals_first_neighbor, 'omitnan');
    end
end

% === Plot 3D Surface ===
figure;
[KINTER, KGRID] = meshgrid(K_inter_values, K_values);
surf(KINTER, KGRID, mean_corr_surface, 'EdgeColor', 'none');
xlabel('K_{inter}'); ylabel('K'); zlabel('Mean Norm Cross-Corr');
title('1st neighbor Cross-Correlation vs K and K_{inter}');
colorbar;
view(45, 30);
set(gca, 'FontSize', 12);

% === Peak Detection ===
[max_val, max_idx] = max(mean_corr_surface(:));
[row, col] = ind2sub(size(mean_corr_surface), max_idx);
fprintf('\n Max correlation = %.4f at K = %.2f, K_inter = %.2f\n', ...
        max_val, K_values(row), K_inter_values(col));
