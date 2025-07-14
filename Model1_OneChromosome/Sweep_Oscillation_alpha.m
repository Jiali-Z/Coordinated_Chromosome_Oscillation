% Sweep_Oscillation_Kct.m
% Sweep Kct values and compute oscillation amplitude/period summary stats

clc; clear; close all;

% Sweep Settings
Alpha_scale = [4,5,6,7,8,9,10,11,12,13,14,15];  % Values to test
num_iterations = 1;  % Number of simulations per Kct

% Time and system parameters
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
n_dot = 1;
Kct=15.4; % pN/µm centromere spring constant
I0 = 2; % µm Rest length of centromere spring
Kkt = 1; % pN/µm kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.1; % kg/s Drag coefficient
Beta = 0.7; % Scaling factor
Nmax = 25; % Maximum number of attachments
Nbar = 20;  % Steady state number of MTs when Ch is centered 
Lambda = n_dot / Nbar; % s^-2 KMT detach rate Akioshi et al. 2010
epsilon = 0.1; % Small perturbation factor

% Storage for output
summary_results = table();

for k_idx = 1:length(Alpha_scale)
    Alpha = n_dot * Alpha_scale(k_idx) / (1 - Beta);
    
    for iter = 1:num_iterations
        % Allocate arrays
        cL = zeros(Nsteps, 2, 1);
        cR = zeros(Nsteps, 2, 1);
        xL = zeros(Nsteps, 1);
        xR = zeros(Nsteps, 1);
        vL = zeros(Nsteps, 1);
        vR = zeros(Nsteps, 1);
        NL = zeros(Nsteps, 1);
        NR = zeros(Nsteps, 1);

        % Initial conditions
        NL(1) = Nbar * (1 + epsilon);
        NR(1) = Nbar * (1 - epsilon);
        cL(1,:,1) = [-I0/2 - epsilon, 0];
        cR(1,:,1) = [ I0/2 + epsilon, 0];
        xL(1) = cL(1,1,1) - 0.3;
        xR(1) = cR(1,1,1) + 0.3;

        for t = 1:Nsteps
            F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));
            F_CT_L =  Kct * (cR(t,1,1) - cL(t,1,1) - I0);
            F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));
            F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0);

            vL(t) = (F_KT_L + F_CT_L ) / Gamma;
            vR(t) = (F_KT_R + F_CT_R ) / Gamma;

            cL(t+1,:,1) = cL(t,:,1) + dt * vL(t);
            cR(t+1,:,1) = cR(t,:,1) + dt * vR(t);
            xL(t+1) = xL(t) + dt * Beta * vL(t);
            xR(t+1) = xR(t) + dt * Beta * vR(t);

            NL(t+1) = NL(t) + dt * (n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) - Lambda * NL(t));
            NR(t+1) = NR(t) + dt * (n_dot + ( Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) - Lambda * NR(t));
        end

        % Run oscillation analysis and tag result with Kct
        [~, ~, table_summary] = oscillation_measurement(cL, cR, dt, iter, true);
        table_summary.Kct = repmat(Kct, height(table_summary), 1);
        table_summary.Alpha = repmat(Alpha, height(table_summary), 1);
        summary_results = [summary_results; table_summary];
    end
end

% Extract CM and KK data
cm_data = summary_results(strcmp(summary_results.metric, 'CM'), :);
kk_data = summary_results(strcmp(summary_results.metric, 'KK'), :);

% Unique Alpha values
alpha_vals = unique(cm_data.Alpha);

% Initialize storage
cm_amp = zeros(size(alpha_vals));
cm_per = zeros(size(alpha_vals));
kk_amp = zeros(size(alpha_vals));
kk_per = zeros(size(alpha_vals));

for i = 1:length(alpha_vals)
    a = alpha_vals(i);
    cm_amp(i) = mean(cm_data.mean_amplitude(cm_data.Alpha == a));
    cm_per(i) = mean(cm_data.mean_period(cm_data.Alpha == a));
    kk_amp(i) = mean(kk_data.mean_amplitude(kk_data.Alpha == a));
    kk_per(i) = mean(kk_data.mean_period(kk_data.Alpha == a));
end

% Plotting
figure;
subplot(2,2,1)
plot(alpha_vals, cm_amp, '-o', 'LineWidth', 1.8)
xlabel('\alpha'); ylabel('CM Amplitude'); title('CM Amplitude vs Alpha'); grid on

subplot(2,2,2)
plot(alpha_vals, cm_per, '-o', 'LineWidth', 1.8)
xlabel('\alpha'); ylabel('CM Period'); title('CM Period vs Alpha'); grid on

subplot(2,2,3)
plot(alpha_vals, kk_amp, '-o', 'LineWidth', 1.8)
xlabel('\alpha'); ylabel('KK Amplitude'); title('KK Amplitude vs Alpha'); grid on

subplot(2,2,4)
plot(alpha_vals, kk_per, '-o', 'LineWidth', 1.8)
xlabel('\alpha'); ylabel('KK Period'); title('KK Period vs Alpha'); grid on
