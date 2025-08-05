function [table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, visualize)
% Extract CM and KK oscillation amplitude and period per iteration

% Time vector
Nsteps = size(cL, 1);
t = (0:Nsteps-1) * dt;
Nchromosomes = size(cL, 3);

table_raw_all = [];
table_cycles_all = [];

% Parameters
buffer_time =1; % min to exclude peaks near edges
buffer_idx = round(buffer_time / dt);

for chr = 1:Nchromosomes
    % Get CM and KK
    CM = (cL(:,1,chr) + cR(:,1,chr)) / 2;
    KK = cR(:,1,chr) - cL(:,1,chr);

    % --- CM peak/valley detection ---
    [CM_peaks, CM_peak_locs] = findpeaks(CM, 'MinPeakProminence', 1,'MinPeakDistance',round(1/dt));
    [CM_valleys, CM_valley_locs] = findpeaks(-CM, 'MinPeakProminence', 1,'MinPeakDistance',round(1/dt));
    CM_valleys = -CM_valleys;

    % Exclude near-boundary peaks/valleys
    valid_CM_peaks = CM_peak_locs > buffer_idx & CM_peak_locs < Nsteps - buffer_idx;
    valid_CM_valleys = CM_valley_locs > buffer_idx & CM_valley_locs < Nsteps - buffer_idx;
    CM_peaks = CM_peaks(valid_CM_peaks);
    CM_peak_locs = CM_peak_locs(valid_CM_peaks);
    CM_valleys = CM_valleys(valid_CM_valleys);
    CM_valley_locs = CM_valley_locs(valid_CM_valleys);

    % --- KK peak/valley detection ---
    [KK_peaks, KK_peak_locs] = findpeaks(KK, 'MinPeakProminence', 0.35,'MinPeakDistance',round(0.5/dt));
    [KK_valleys, KK_valley_locs] = findpeaks(-KK, 'MinPeakProminence', 0.35,'MinPeakDistance',round(0.5/dt));
    KK_valleys = -KK_valleys;

    valid_KK_peaks = KK_peak_locs > buffer_idx & KK_peak_locs < Nsteps - buffer_idx;
    valid_KK_valleys = KK_valley_locs > buffer_idx & KK_valley_locs < Nsteps - buffer_idx;
    KK_peaks = KK_peaks(valid_KK_peaks);
    KK_peak_locs = KK_peak_locs(valid_KK_peaks);
    KK_valleys = KK_valleys(valid_KK_valleys);
    KK_valley_locs = KK_valley_locs(valid_KK_valleys);

    % --- Construct raw peak/valley table ---
    T_Raw_CM = table([t(CM_peak_locs)'; t(CM_valley_locs)'], ...
        [CM_peaks; CM_valleys], ...
        [repmat("peak", length(CM_peaks), 1); repmat("valley", length(CM_valleys), 1)], ...
        repmat("CM", length(CM_peaks)+length(CM_valleys), 1), ...
        repmat(chr, length(CM_peaks)+length(CM_valleys), 1), ...
        repmat(iter, length(CM_peaks)+length(CM_valleys), 1), ...
        'VariableNames', {'time', 'value', 'type', 'metric', 'chromosome', 'iteration'});
    T_Raw_CM = sortrows(T_Raw_CM, 'time');

    T_Raw_KK = table([t(KK_peak_locs)'; t(KK_valley_locs)'], ...
        [KK_peaks; KK_valleys], ...
        [repmat("peak", length(KK_peaks), 1); repmat("valley", length(KK_valleys), 1)], ...
        repmat("KK", length(KK_peaks)+length(KK_valleys), 1), ...
        repmat(chr, length(KK_peaks)+length(KK_valleys), 1), ...
        repmat(iter, length(KK_peaks)+length(KK_valleys), 1), ...
        'VariableNames', {'time', 'value', 'type', 'metric', 'chromosome', 'iteration'});
    T_Raw_KK = sortrows(T_Raw_KK, 'time');

    table_raw_all = [table_raw_all; T_Raw_CM; T_Raw_KK];

    % --- Pairing for CM cycles ---
    CM_events = [CM_valley_locs(:), zeros(length(CM_valley_locs),1); ...
                 CM_peak_locs(:), ones(length(CM_peak_locs),1)];
    CM_events = sortrows(CM_events, 1);

    paired_amp_CM = [];
    paired_period_CM = [];
    for i = 1:(size(CM_events,1)-1)
        if CM_events(i,2) == 0 && CM_events(i+1,2) == 1 % valley → peak
            vi = CM_events(i,1); pi = CM_events(i+1,1);
            dt_cycle = t(pi) - t(vi);
            if dt_cycle <= 10.0 % CM: max valley-to-peak time = 2 min
                paired_amp_CM(end+1,1) = (CM(pi) - CM(vi)) / 2;
                paired_period_CM(end+1,1) = dt_cycle * 2; % full period
            end
        end
    end

    % --- Pairing for KK cycles ---
    KK_events = [KK_valley_locs(:), zeros(length(KK_valley_locs),1); ...
                 KK_peak_locs(:), ones(length(KK_peak_locs),1)];
    KK_events = sortrows(KK_events, 1);

    paired_amp_KK = [];
    paired_period_KK = [];
    for i = 1:(size(KK_events,1)-1)
        if KK_events(i,2) == 0 && KK_events(i+1,2) == 1 % valley → peak
            vi = KK_events(i,1); pi = KK_events(i+1,1);
            dt_cycle = t(pi) - t(vi);
            if dt_cycle <= 10.0 % KK: max valley-to-peak time = 1 min
                paired_amp_KK(end+1,1) = (KK(pi) - KK(vi)) / 2;
                paired_period_KK(end+1,1) = dt_cycle * 2; % full period
            end
        end
    end

    % Cycle table: CM
    T_cycles_CM = table((1:length(paired_amp_CM))', ...
        repmat(chr, length(paired_amp_CM), 1), ...
        repmat(iter, length(paired_amp_CM), 1), ...
        repmat("CM", length(paired_amp_CM), 1), ...
        paired_amp_CM, paired_period_CM, ...
        'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'amplitude', 'period'});

    % Cycle table: KK
    T_cycles_KK = table((1:length(paired_amp_KK))', ...
        repmat(chr, length(paired_amp_KK), 1), ...
        repmat(iter, length(paired_amp_KK), 1), ...
        repmat("KK", length(paired_amp_KK), 1), ...
        paired_amp_KK, paired_period_KK, ...
        'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'amplitude', 'period'});

    table_cycles_all = [table_cycles_all; T_cycles_CM; T_cycles_KK];

    % --- Visualization ---
    if visualize
        figure;
        plot(t, CM, 'k', 'LineWidth', 2); hold on;
        plot(t, KK, 'm', 'LineWidth', 2);
        plot(t(CM_peak_locs), CM(CM_peak_locs), 'ko', 'MarkerFaceColor', 'g');
        plot(t(CM_valley_locs), CM(CM_valley_locs), 'ks', 'MarkerFaceColor', 'y');
        plot(t(KK_peak_locs), KK(KK_peak_locs), 'mo', 'MarkerFaceColor', 'c');
        plot(t(KK_valley_locs), KK(KK_valley_locs), 'ms', 'MarkerFaceColor', 'w');
        legend('CM', 'KK', 'CM peak', 'CM valley', 'KK peak', 'KK valley');
        xlabel('Time (min)');
        ylabel('Position (µm)');
        xticks(0:1:20);
        title(sprintf('Chromosome %d (Iter %d)', chr, iter));
        grid on;
    end
end

% Summary
table_summary = groupsummary(table_cycles_all, ...
    {'iteration', 'chromosome', 'metric'}, ...
    ["mean", "std"], ["amplitude", "period"]);

% Outputs
table_raw = table_raw_all;
table_cycle = table_cycles_all;
end



















% function [table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, visualize)
% % ANALYZE_CM_KK_OSCILLATIONS 
% % Extracts CM and KK oscillation features (amplitude, period) from one iteration of cL and cR.
% % 
% % INPUTS:
% %   cL, cR     - [Nsteps × 2 × Nchromosomes] left and right chromatid positions
% %   dt         - time step in minutes
% %   iter       - iteration ID (integer)
% %   visualize  - true/false flag to plot traces and peak/valley locations
% %
% % OUTPUTS:
% %   table_raw     - table of all CM and KK peaks/valleys with time and value
% %   table_cycle   - per-cycle amplitude and period data
% %   table_summary - mean/std of amplitude and period grouped by chromosome and metric
% 
% % Time vector
% Nsteps = size(cL, 1);
% t = (0:Nsteps-1) * dt;
% Nchromosomes = size(cL, 3);
% 
% table_raw_all = [];
% table_cycles_all = [];
% 
% for chr = 1:Nchromosomes
%     % Get CM and KK
%     CM = (cL(:,1,chr) + cR(:,1,chr)) / 2;
%     KK = cR(:,1,chr) - cL(:,1,chr);
% 
%     % Find peaks and valleys
%     [CM_peaks, CM_peak_locs] = findpeaks(CM, 'MinPeakProminence', 1,'MinPeakDistance',1);
%     [CM_valleys, CM_valley_locs] = findpeaks(-CM, 'MinPeakProminence', 1,'MinPeakDistance',1); 
%     CM_valleys = -CM_valleys;
% 
%     [KK_peaks, KK_peak_locs] = findpeaks(KK, 'MinPeakProminence', 0.4,'MinPeakDistance',0.5);
%     [KK_valleys, KK_valley_locs] = findpeaks(-KK, 'MinPeakProminence',0.4,'MinPeakDistance',0.5); 
%     KK_valleys = -KK_valleys;
% 
%     % --- Raw table ---
%     T_Raw_CM = table(...
%         [t(CM_peak_locs)'; t(CM_valley_locs)'], ...
%         [CM_peaks; CM_valleys], ...
%         [repmat("peak", length(CM_peaks), 1); repmat("valley", length(CM_valleys), 1)], ...
%         repmat("CM", length(CM_peaks)+length(CM_valleys), 1), ...
%         repmat(chr, length(CM_peaks)+length(CM_valleys), 1), ...
%         repmat(iter, length(CM_peaks)+length(CM_valleys), 1), ...
%         'VariableNames', {'time', 'value', 'type', 'metric', 'chromosome', 'iteration'});
%     T_Raw_CM = sortrows(T_Raw_CM, 'time');
% 
%     T_Raw_KK = table(...
%         [t(KK_peak_locs)'; t(KK_valley_locs)'], ...
%         [KK_peaks; KK_valleys], ...
%         [repmat("peak", length(KK_peaks), 1); repmat("valley", length(KK_valleys), 1)], ...
%         repmat("KK", length(KK_peaks)+length(KK_valleys), 1), ...
%         repmat(chr, length(KK_peaks)+length(KK_valleys), 1), ...
%         repmat(iter, length(KK_peaks)+length(KK_valleys), 1), ...
%         'VariableNames', {'time', 'value', 'type', 'metric', 'chromosome', 'iteration'});
%     T_Raw_KK = sortrows(T_Raw_KK, 'time');
% 
%     table_raw_all = [table_raw_all; T_Raw_CM; T_Raw_KK];
% 
%     % --- Amplitude table ---
%     N_CM_amp = min(length(CM_peaks), length(CM_valleys));
%     T_amplitude_CM = table(...
%         (1:N_CM_amp)', ...
%         repmat(chr, N_CM_amp, 1), ...
%         repmat(iter, N_CM_amp, 1), ...
%         repmat("CM", N_CM_amp, 1), ...
%         (CM_peaks(1:N_CM_amp) - CM_valleys(1:N_CM_amp)) / 2, ...
%         'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'amplitude'});
% 
%     N_KK_amp = min(length(KK_peaks), length(KK_valleys));
%     T_amplitude_KK = table(...
%         (1:N_KK_amp)', ...
%         repmat(chr, N_KK_amp, 1), ...
%         repmat(iter, N_KK_amp, 1), ...
%         repmat("KK", N_KK_amp, 1), ...
%         (KK_peaks(1:N_KK_amp) - KK_valleys(1:N_KK_amp)) / 2, ...
%         'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'amplitude'});
% 
%     % --- Period table ---
%     N_CM_period = length(CM_peak_locs) - 1;
%     T_period_CM = table(...
%         (1:N_CM_period)', ...
%         repmat(chr, N_CM_period, 1), ...
%         repmat(iter, N_CM_period, 1), ...
%         repmat("CM", N_CM_period, 1), ...
%         diff(t(CM_peak_locs)'), ...
%         'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'period'});
% 
%     N_KK_period = length(KK_peak_locs) - 1;
%     T_period_KK = table(...
%         (1:N_KK_period)', ...
%         repmat(chr, N_KK_period, 1), ...
%         repmat(iter, N_KK_period, 1), ...
%         repmat("KK", N_KK_period, 1), ...
%         diff(t(KK_peak_locs)'), ...
%         'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'period'});
% 
%     % --- Merge amplitude and period ---
%     T_cycles_CM = innerjoin(T_amplitude_CM, T_period_CM, ...
%         'Keys', {'cycle_id', 'chromosome', 'iteration', 'metric'});
%     T_cycles_KK = innerjoin(T_amplitude_KK, T_period_KK, ...
%         'Keys', {'cycle_id', 'chromosome', 'iteration', 'metric'});
% 
%     table_cycles_all = [table_cycles_all; T_cycles_CM; T_cycles_KK];
% 
%     % --- Visualization ---
%     if visualize
%         figure;
%         plot(t, CM, 'k', 'LineWidth', 2); hold on;
%         plot(t, KK, 'm', 'LineWidth', 2);
%         plot(t(CM_peak_locs), CM_peaks, 'ko', 'MarkerFaceColor', 'g');
%         plot(t(CM_valley_locs), CM_valleys, 'ks', 'MarkerFaceColor', 'y');
%         plot(t(KK_peak_locs), KK_peaks, 'mo', 'MarkerFaceColor', 'c');
%         plot(t(KK_valley_locs), KK_valleys, 'ms', 'MarkerFaceColor', 'w');
%         legend('CM', 'KK', 'CM peak', 'CM valley', 'KK peak', 'KK valley');
%         xlabel('Time (min)');
%         ylabel('Position (µm)');
%         xticks(0:1:20);  % Adjust the range as needed
%         title(sprintf('Chromosome %d (Iter %d)', chr, iter));
%         grid on;
%     end
% end
% 
% % Summary table
% table_summary = groupsummary(table_cycles_all, ...
%     {'iteration', 'chromosome', 'metric'}, ...
%     ["mean", "std"], ["amplitude", "period"]);
% 
% % Final outputs
% table_raw = table_raw_all;
% table_cycle = table_cycles_all;
% 
% end
