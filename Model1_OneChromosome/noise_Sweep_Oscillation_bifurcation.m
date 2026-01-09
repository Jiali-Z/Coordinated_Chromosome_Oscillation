clc; clear; close all;

% ------------------ Global Settings ------------------
dt        = 2e-3;           % Timestep (min)
Nsteps    = 5000;           % Number of steps
Beta      = 0.7;
duration  = 10;             % Duration (min) used in KE analysis
Alpha_scales = 3:1:10;    % Alpha scale factors
Kct_values   = 10:0.1:20;   % Kct values to sweep

% Fixed (per-run) constants
I0   = 2;       % µm, rest length centromere spring
Kkt  = 1;       % pN/µm, kinetochore spring constant
Gamma= 0.1;     % kg/s, drag
Nmax = 25;      % Max attachments


% Noise / cell-to-cell variability
N_iter  = 250;              % Number of cells per (alpha, Kct)
noise   = 0.05;             % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (position + force)
results_osc = table();   % one row per iteration with CM/KK stats
nTotal  = numel(Alpha_scales) * numel(Kct_values) * N_iter;
counter = 0;

% ------------------ Parameter Sweep with Noise ------------------
for a = Alpha_scales
    for Kct = Kct_values
        for iter = 1:N_iter
            % ---- cell-specific parameters with noise ----
            n_dot  = 1 * (1 + noise*(2*rand - 1));  % production rate
            Nbar   = 20   * (1 + noise*(2*rand - 1));  % steady-state MTs
            Lambda = n_dot / Nbar;                             % detach rate
            Alpha  = n_dot * a / (1 - Beta);                   % scaled Alpha
            epsilon = 0.1;                                     % small perturb.

            % ---- allocate & initialize ----
            % Note: Nsteps+1 because we use t and t+1
            xL = zeros(Nsteps+1,1); xR = zeros(Nsteps+1,1);
            NL = zeros(Nsteps+1,1); NR = zeros(Nsteps+1,1);
            cL = zeros(Nsteps+1,2,1); cR = zeros(Nsteps+1,2,1);
            vL = zeros(Nsteps+1,1);  vR = zeros(Nsteps+1,1);

            NR(1)     = Nbar*(1 - epsilon);
            NL(1)     = Nbar*(1 + epsilon);
            cL(1,1,1) =-I0/2-epsilon;
            cR(1,1,1) =I0/2+epsilon;
            cL(:,2,1) = 0;
            cR(:,2,1) = 0;

            % MT tips start slightly offset from centromeres
            xL(1) = cL(1,1,1) - 0.3;
            xR(1) = cR(1,1,1) + 0.3;
            vL(1) = 0; vR(1) = 0;

            % ---- time stepping ----
            for t = 1:Nsteps
                % Forces on L and R centromeres
                F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t));
                F_CT_L =  Kct   * (cR(t,1,1) - cL(t,1,1) - I0);
                F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t));
                F_CT_R = -Kct   * (cR(t,1,1) - cL(t,1,1) - I0);

                % Velocities
                vL(t) = (F_KT_L + F_CT_L ) / Gamma;
                vR(t) = (F_KT_R + F_CT_R ) / Gamma;

                % Brownian step (shared to keep the pair coherent; change if desired)
                dx_diff = noise_b * randn(1,1) * sqrt(dt);
                % Update centromeres
                cL(t+1,1,1) = cL(t,1,1) + dt*vL(t) + dx_diff;
                cR(t+1,1,1) = cR(t,1,1) + dt*vR(t) + dx_diff;
                % Update MT tips
                xL(t+1) = xL(t) + dt * Beta * vL(t);
                xR(t+1) = xR(t) + dt * Beta * vR(t);
                % Attachments
                NL(t+1) = NL(t) + dt * ( ...
                    n_dot + (-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
                    - Lambda * NL(t));
                NR(t+1) = NR(t) + dt * ( ...
                    n_dot + ( Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
                    - Lambda * NR(t));
            end

            % ---- Oscillation Analysis + Safe Append ----
             [table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, false);
            
            if ~isempty(table_summary)
                table_summary.alpha_scale = repmat(a,     height(table_summary), 1);
                table_summary.alpha       = repmat(Alpha, height(table_summary), 1);
                table_summary.Kct         = repmat(Kct,   height(table_summary), 1);
                table_summary = movevars(table_summary, {'alpha_scale','alpha','Kct'}, 'Before', 1);
            
                if isempty(results_osc)
                    results_osc = table_summary;
                else
                    results_osc = [results_osc; table_summary(:, results_osc.Properties.VariableNames)];
                end
            else
                fprintf('Iter %d: no oscillation rows — skipping append.\n', iter);
            end

            % --- Progress note ---
            counter = counter + 1;
            fprintf('Progress: %d/%d (%.1f%%) | a=%.1f, Kct=%.2f, iter=%d\n', ...
                counter, nTotal, 100*counter/nTotal, a, Kct, iter);

        end
    end 
end

%% Export the data for python plots 
% results_all has columns: alpha_scale, alpha, Kct, iteration, mean_KE
writetable(results_osc, 'alpha3-10_noise_results_oscillation.csv');   % CSV for Python


%%

% Group by alpha_scale, Kct, and metric (CM/KK)
[G, key_alpha, key_Kct, key_metric] = findgroups(results_osc.alpha_scale, results_osc.Kct, results_osc.metric);

% Helper to compute mean while omitting NaNs
mean_omit = @(x) mean(x, 'omitnan');

% Aggregate
avg_amp1   = splitapply(mean_omit, results_osc.mean_amplitude1, G);
avg_amp2   = splitapply(mean_omit, results_osc.mean_amplitude2, G);
avg_period = splitapply(mean_omit, results_osc.mean_period,    G);
N          = splitapply(@(x) sum(~isnan(x)), results_osc.mean_amplitude1, G);  % count used

% Build the summary table (long format: one row per (alpha_scale, Kct, metric))
summary_tbl = table(key_alpha, key_Kct, key_metric, N, ...
    avg_amp1, avg_amp2, avg_period, ...
    'VariableNames', {'alpha_scale','Kct','metric','N', ...
                      'avg_amplitude1','avg_amplitude2','avg_period'});

summary_tbl = sortrows(summary_tbl,"metric","ascend");


%% ---------- INPUT ----------
T = summary_tbl;  % alpha_scale, Kct, metric, N, avg_amplitude1, avg_amplitude2, avg_period

% Axes (sorted)
a = unique(T.alpha_scale(:));     % x-axis (alpha_scale)
k = unique(T.Kct(:));             % y-axis (Kct)

% Map each row to (row=Kct, col=alpha_scale)
[~, ia] = ismember(T.alpha_scale, a);   % columns
[~, ik] = ismember(T.Kct,        k);    % rows

% Helper: put a field into a (Kct × alpha_scale) grid for a given metric
gridify = @(field, met) accumarray( ...
    [ik(T.metric==met), ia(T.metric==met)], ...
    T.(field)(T.metric==met), ...
    [numel(k) numel(a)], @mean, NaN);

% ---------- BUILD 6 MATRICES ----------
CM_amp1   = gridify('avg_amplitude1','CM');
CM_amp2   = gridify('avg_amplitude2','CM');
CM_period = gridify('avg_period',    'CM');

KK_amp1   = gridify('avg_amplitude1','KK');
KK_amp2   = gridify('avg_amplitude2','KK');
KK_period = gridify('avg_period',    'KK');

% ---------- PLOT ----------
figure('Color','w','Position',[100 100 1200 650]);
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

plotspec = { ...
    CM_amp1,   'CM amplitude1', 0.7868;   % (1)
    CM_amp2,   'CM amplitude2', 0.7868;   % (1)
    CM_period, 'CM period',     2.035;    % (2)
    KK_amp1,   'KK amplitude1', 0.32;     % (3)
    KK_amp2,   'KK amplitude2', 0.32;     % (3)
    KK_period, 'KK period',     1.076     % (4)
};

% For contour coordinates
[AA, KKgrid] = meshgrid(a, k);

for i = 1:size(plotspec,1)
    nexttile
    Z = plotspec{i,1};
    level = plotspec{i,3};

    hImg = imagesc(a, k, Z);
    set(gca,'YDir','normal'); colormap(parula); colorbar
    xlabel('\alpha scale'); ylabel('K_{ct} (pN/\mu m)')
    title(plotspec{i,2})
    set(gca,'Color',[0.95 0.95 0.95]);
    set(hImg,'AlphaData',~isnan(Z));

    hold on
    Zc = Z; Zc(isnan(Zc)) = -Inf;

    % dashed line for period plots, solid for amplitudes
    isPeriod = contains(lower(plotspec{i,2}), 'period');
    ls = '--'; if ~isPeriod, ls = '-'; end

    [C, hC] = contour(AA, KKgrid, Zc, [level level], ...
                      'LineColor','k','LineWidth',2,'LineStyle',ls);

    % IMPORTANT: don't capture outputs from clabel when supplying hC
    clabel(C, hC, 'LabelSpacing', 500);  % optional label, no output

    hold off
end
