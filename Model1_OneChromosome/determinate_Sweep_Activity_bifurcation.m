clc; clear; close all;

% ------------------ Simulation Parameters ------------------
dt = 2e-3;              % Timestep (min)
Nsteps = 5000;          % Number of steps
n_dot = 1;
Beta = 0.7;
duration = 10;          % Duration (min) used in KE analysis
Alpha_scales = 3:0.2:15;    % Alpha scale factors
Kct_values = 10:0.1:20; % Kct values to sweep
nTotal=numel(Alpha_scales)*numel(Kct_values);
% ------------------ Storage for results ------------------
results = [];
counter = 0;
% ------------------ Parameter Sweep ------------------
for a = Alpha_scales
    Alpha = n_dot * a / (1 - Beta);  % Recalculate Alpha

    for Kct = Kct_values
        counter = counter + 1;
        % Fixed parameters
        I0 = 2; % µm
        Kkt = 1; % pN/µm
        Gamma = 0.1; % kg/s
        Nmax = 25;
        Nbar = 20;
        Lambda = n_dot / Nbar;
        epsilon = 0.1;

        % Initialization
        xL = zeros(Nsteps+1, 1); xR = zeros(Nsteps+1, 1);
        NL = zeros(Nsteps+1, 1); NR = zeros(Nsteps+1, 1);
        cL = zeros(Nsteps+1, 2, 1); cR = zeros(Nsteps+1, 2, 1);
        vL = zeros(Nsteps, 1); vR = zeros(Nsteps, 1);

        NR(1) = Nbar * (1 - epsilon);
        NL(1) = Nbar * (1 + epsilon);
        cL(1,:,1) = [-I0/2 - epsilon, 0];
        cR(1,:,1) = [ I0/2 + epsilon, 0];
        xL(1) = cL(1,1,1) - 0.3;
        xR(1) = cR(1,1,1) + 0.3;

        % ------------------ Run Simulation ------------------
        for t = 1:Nsteps
            % Forces
            F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));
            F_CT_L = Kct * (cR(t,1,1) - cL(t,1,1) - I0);
            F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));
            F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0);
            % Velocities
            vL(t) = (F_KT_L + F_CT_L) / Gamma;
            vR(t) = (F_KT_R + F_CT_R) / Gamma;
            % Update chromosome x positions
            cL(t+1,1,1) = cL(t,1,1) + dt * vL(t);
            cR(t+1,1,1) = cR(t,1,1) + dt * vR(t);
            % Update MT tip positions
            xL(t+1) = xL(t) + dt * Beta * vL(t);
            xR(t+1) = xR(t) + dt * Beta * vR(t);
            % Update attachments 
            NL(t+1) = NL(t) + dt * ( ...
                n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
                - Lambda * NL(t) );

            NR(t+1) = NR(t) + dt * ( ...
                n_dot + ( Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
                - Lambda * NR(t) );
        end

        % ------------------ Activity Measurement ------------------
        activity_tbl = chromosome_activity_measurement(cL, cR, dt);
        avg_KE = activity_tbl.Avg_KE;   % assumes single-row output

        % ------------------ Store result ------------------
        results = [results; struct( ...
            'AlphaScale', a, ...
            'Alpha',      Alpha, ...
            'Kct',        Kct, ...
            'Avg_KE',     avg_KE )];
        % ------------------ Progress Print ------------------
        fprintf('Progress: %d / %d (%.1f%%) | Alpha = %.1f, Kct = %.2f\n', ...
            counter, nTotal, 100*counter/nTotal, a, Kct);
    end
end

results_table = struct2table(results);

%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_table, 'MSV_determinate_Alpha_3_0.2_15.csv');      % per-point CSV
%%
% ------------------ Plot: Avg KE vs Kct (grouped by AlphaScale) ------------------
figure; hold on;
alpha_values = unique(results_table.AlphaScale);
colors = turbo(length(alpha_values));

for i = 1:length(alpha_values)
    a_val = alpha_values(i);
    subset = results_table(results_table.AlphaScale == a_val, :);
    plot(subset.Kct, subset.Avg_KE, ...
        'LineStyle','-', 'Marker','o', 'Color', colors(i,:), ...
        'LineWidth', 2, 'MarkerSize', 5, ...
        'DisplayName', sprintf('$\\alpha = %.1f$', a_val));
end

xlabel('\bf $K_{ct}$ (pN/$\mu$m)','Interpreter','latex', 'FontName','Arial','FontSize',18);
ylabel('\bf Avg. MSV ($\mu$m$^2$/min$^2$)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
legend('Location','best', 'Interpreter','latex');
grid on; set(gca,'FontSize',14);

%%
% ------------------ Plot: Avg KE vs alpha scale (grouped by Kct) ------------------
figure; hold on;
kct_values = unique(results_table.Kct);
colors = turbo(length(kct_values));

for i = 1:length(kct_values)
    kct = kct_values(i);
    subset = results_table(results_table.Kct == kct, :);
    plot(subset.AlphaScale, subset.Avg_KE, ...
        'LineStyle','-', 'Marker','o', 'Color', colors(i,:), ...
        'LineWidth', 2, 'MarkerSize', 5, ...
        'DisplayName', sprintf('$K_{ct} = %.2f$', kct));
end

xlabel('\bf$\alpha$ (scale factor)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
ylabel('\bf Avg. MSV ($\mu$m$^2$/min$^2$)', 'Interpreter','latex', 'FontName','Arial','FontSize',18);
legend('Location','best', 'Interpreter','latex');
grid on; set(gca,'FontSize',14);

%% 
% --- 1) Make a regular grid (rows = AlphaScale, cols = Kct) ---
aVals = unique(results_table.AlphaScale, 'stable');    % row coordinates
kVals = unique(results_table.Kct,        'stable');    % column coordinates

% Robust indexing in case there are repeated (AlphaScale,Kct) pairs:
[~, ia] = ismember(results_table.AlphaScale, aVals);
[~, ik] = ismember(results_table.Kct,   kVals);

% Average duplicates into a matrix (NaN where missing)
Z = accumarray([ia ik], results_table.Avg_KE, [numel(aVals) numel(kVals)], @mean, NaN);

% --- 2) Normalize KE to [0,1] for sharp contrast ---
% minZ = min(Z(:), [], 'omitnan');
% maxZ = max(Z(:), [], 'omitnan');
% Znorm = (Z - minZ) ./ (maxZ - minZ);

% --- 3) Heatmap (2D) ---
figure('Color','w'); 
imagesc(aVals, kVals, Z);
set(gca, 'YDir','normal');                 % so low AlphaScale is at bottom
axis tight
colormap(parula);                           % or 'turbo' if you prefer
c = colorbar; c.Label.String = 'Average MSV (0–1)';
ylabel('K_{ct} (pN/\mu m)');
xlabel('\alpha_{scale}');
title('Avg KE Phase Map (normalized)');




