% Sweep_Oscillation_Kct.m
% Sweep Kct values and compute oscillation amplitude/period summary stats

clc; clear; close all;

% Sweep Settings
Kct_values = [12.5,13,14,15,16,17,18,19,20, 22, 25, 30, 35, 40];  % Values to test
num_iterations = 1;  % Number of simulations per Kct

% Time and system parameters
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
n_dot = 1; 
I0 = 2; % µm Rest length of centromere spring
Kkt = 1; % pN/µm kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.1; % kg/s Drag coefficient
Beta = 0.7; % Scaling factor
Nmax = 25; % Maximum number of attachments
Nbar = 20;  % Steady state number of MTs when Ch is centered 
Lambda = n_dot / Nbar; % s^-2 KMT detach rate Akioshi et al. 2010
Alpha = n_dot * 5.6 / (1 - Beta);
epsilon = 0.1; % Small perturbation factor


% Storage for output
summary_results = table();

for k_idx = 1:length(Kct_values)
    Kct = Kct_values(k_idx);
    
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
            F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t)); % MT forces on Left Chromosome
            F_CT_L =  Kct * (cR(t,1,1) - cL(t,1,1) - I0); % Centromere forces on Left Chromosome
            F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t)); % MT forces on Right Chromosome
            F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0); % Centromere forces on Right Chromosome

            vL(t) = (F_KT_L + F_CT_L) / Gamma;
            vR(t) = (F_KT_R + F_CT_R ) / Gamma;
            % Calculate the current chromosome position in both directions
            cL(t+1,:,1) = cL(t,:,1) + dt * vL(t);
            cR(t+1,:,1) = cR(t,:,1) + dt * vR(t);
            % Calculate the current microtubule tip position in both directions
            xL(t+1) = xL(t) + dt * Beta * vL(t);
            xR(t+1) = xR(t) + dt * Beta * vR(t);
            % Update number of attachments 
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

% Plotting
% Make sure summary_results is already loaded

% Extract CM and KK separately
cm_data = summary_results(strcmp(summary_results.metric, 'CM'), :);
kk_data = summary_results(strcmp(summary_results.metric, 'KK'), :);

% Use unique Kct values for x-axis
Kct_vals = unique(summary_results.Kct);

% Preallocate arrays
cm_amp = zeros(size(Kct_vals));
cm_per = zeros(size(Kct_vals));
kk_amp = zeros(size(Kct_vals));
kk_per = zeros(size(Kct_vals));

% Fill arrays
for i = 1:length(Kct_vals)
    k = Kct_vals(i);
    cm_amp(i) = mean(cm_data.mean_amplitude(cm_data.Kct == k));
    cm_per(i) = mean(cm_data.mean_period(cm_data.Kct == k));
    kk_amp(i) = mean(kk_data.mean_amplitude(kk_data.Kct == k));
    kk_per(i) = mean(kk_data.mean_period(kk_data.Kct == k));
end

% Plotting
figure;
subplot(2,2,1)
plot(Kct_vals, cm_amp, '-o', 'LineWidth', 1.8)
xlabel('Kct'); ylabel('CM Amplitude'); title('CM Amplitude vs Kct'); grid on

subplot(2,2,2)
plot(Kct_vals, cm_per, '-o', 'LineWidth', 1.8)
xlabel('Kct'); ylabel('CM Period'); title('CM Period vs Kct'); grid on

subplot(2,2,3)
plot(Kct_vals, kk_amp, '-o', 'LineWidth', 1.8)
xlabel('Kct'); ylabel('KK Amplitude'); title('KK Amplitude vs Kct'); grid on

subplot(2,2,4)
plot(Kct_vals, kk_per, '-o', 'LineWidth', 1.8)
xlabel('Kct'); ylabel('KK Period'); title('KK Period vs Kct'); grid on
