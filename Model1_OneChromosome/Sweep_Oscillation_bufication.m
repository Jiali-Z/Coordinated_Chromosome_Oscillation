clc; clear; close all;

% ------------------ Simulation Parameters ------------------
dt = 2e-3;              % Timestep (min)
Nsteps = 5000;          % Number of steps
n_dot = 1;
Beta = 0.7;
duration = 10;          % Duration (min) used in KE analysis
Alpha_scales = 3:15;    % Alpha scale factors
Kct_values = 10:0.1:20; % Kct values to sweep

% ------------------ Storage for results ------------------
results = [];

% ------------------ Parameter Sweep ------------------
for a = Alpha_scales
    Alpha = n_dot * a / (1 - Beta);  % Recalculate Alpha

    for Kct = Kct_values
        % Fixed parameters
        I0 = 2;
        Kkt = 1;
        Gamma = 0.1;
        Nmax = 25;
        Nbar = 20;
        Lambda = n_dot / Nbar;
        epsilon = 0.1;

        % Initialization
        xL = zeros(Nsteps, 1); xR = zeros(Nsteps, 1);
        NL = zeros(Nsteps, 1); NR = zeros(Nsteps, 1);
        cL = zeros(Nsteps, 2, 1); cR = zeros(Nsteps, 2, 1);
        vL = zeros(Nsteps, 1); vR = zeros(Nsteps, 1);

        NR(1) = Nbar * (1 - epsilon);
        NL(1) = Nbar * (1 + epsilon);
        cL(1,:,1) = [-I0/2 - epsilon, 0];
        cR(1,:,1) = [ I0/2 + epsilon, 0];
        xL(1) = cL(1,1,1) - 0.3;
        xR(1) = cR(1,1,1) + 0.3;

        % ------------------ Run Simulation ------------------
        for t = 1:Nsteps
            F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));
            F_CT_L = Kct * (cR(t,1,1) - cL(t,1,1) - I0);
            F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));
            F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0);

            vL(t) = (F_KT_L + F_CT_L) / Gamma;
            vR(t) = (F_KT_R + F_CT_R) / Gamma;

            cL(t+1,:,1) = cL(t,:,1) + dt * vL(t);
            cR(t+1,:,1) = cR(t,:,1) + dt * vR(t);
            xL(t+1) = xL(t) + dt * Beta * vL(t);
            xR(t+1) = xR(t) + dt * Beta * vR(t);

            NL(t+1) = NL(t) + dt * (n_dot - Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax) - Lambda * NL(t));
            NR(t+1) = NR(t) + dt * (n_dot + Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax) - Lambda * NR(t));
        end

        % ------------------ Activity Measurement ------------------
        sampling_interval = round(5 / 60 / dt);
        total_steps = min(round(duration / dt), Nsteps);

        signalL = squeeze(cL(1:sampling_interval:total_steps, 1, 1));
        signalR = squeeze(cR(1:sampling_interval:total_steps, 1, 1));
        time_sampled = (0:length(signalL)-1)' * sampling_interval * dt;

        vL_list = []; vR_list = [];

        for j = 1:(length(signalL) - 2)
            t_chunk = time_sampled(j:j+2);
            sL_chunk = signalL(j:j+2);
            sR_chunk = signalR(j:j+2);

            pL = polyfit(t_chunk, sL_chunk, 1);
            pR = polyfit(t_chunk, sR_chunk, 1);

            vL_list(end+1) = pL(1);
            vR_list(end+1) = pR(1);
        end

        KE = vL_list.^2 + vR_list.^2;
        avg_KE = mean(KE);

        % ------------------ Store result ------------------
        results = [results; struct( ...
            'AlphaScale', a, ...
            'Kct', Kct, ...
            'Avg_KE', avg_KE ...
        )];
    end
end

results_table = struct2table(results);

% ------------------ Plot: Bifurcation in Avg KE ------------------
figure;
hold on;
colors = lines(length(Alpha_scales));

for i = 1:length(Alpha_scales)
    a = Alpha_scales(i);
    subset = results_table(results_table.AlphaScale == a, :);
    
    % Assign marker by Alpha category
    if a == 3
        marker_style = 'o';
    elseif a >= 4 && a <= 10
        marker_style = 's';
    else
        marker_style = '^';
    end

    plot(subset.Kct, subset.Avg_KE, ...
        'LineStyle', '-', ...
        'Marker', marker_style, ...
        'Color', colors(i,:), ...
        'LineWidth', 1, ...
        'MarkerSize', 5, ...
        'DisplayName', ['$\alpha = $' num2str(a)]);
end

xlabel('$K_{ct}$ (pN/$\mu$m)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Avg. Kinetic Energy ($\mu$m$^2$/min$^2$)', 'Interpreter', 'latex', 'FontSize', 18);
title('Bifurcation Plot: Avg KE vs $K_{ct}$', 'Interpreter', 'latex', 'FontSize', 20);
legend('Location', 'eastoutside', 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 14);
