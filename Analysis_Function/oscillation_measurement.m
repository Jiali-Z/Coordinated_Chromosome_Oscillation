function [table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, visualize)
%   Extract oscillatory features of chromosome motion from a single simulation
%   iteration by analyzing:
%     (1) Center-of-mass (CM) oscillations
%     (2) Inter-kinetochore distance (KK) oscillations
%
%   Oscillation cycles are identified using a three-extrema (peak–valley–peak
%   or valley–peak–valley) detector with metric-specific amplitude and period
%   constraints. Only VALID cycles are retained for summary statistics.
%
%   The function outputs raw extrema, per-cycle measurements, and per-chromosome
%   summary statistics suitable for population-level pooling across iterations.
%
% INPUTS
%   cL, cR     : [Nsteps × 2 × Nchromosomes] arrays
%                Left and right chromatid positions.
%                Dimension 2 corresponds to [x, y]; analysis is performed on x.

%   dt         : Scalar (minutes)
%                Simulation time step.
%
%   iter       : Scalar integer
%                Iteration index (used for bookkeeping in output tables).
%
%   visualize  : (optional) logical, default = false
%                If true, plots CM and KK trajectories for each chromosome and
%                overlays only VALID oscillation cycles used in the analysis.
%
% OUTPUTS
%   table_raw
%       Table of all detected CM and KK extrema (peaks and valleys), including:
%         - time        : time of extrema (minutes)
%         - value       : signal value at extrema
%         - type        : "peak" or "valley"
%         - metric      : "CM" or "KK"
%         - chromosome  : chromosome index
%         - iteration   : iteration ID
%
%   table_cycle
%       Table of per-cycle oscillation measurements extracted from consecutive
%       extrema triplets, including:
%         - amplitude1, amplitude2 : half peak-to-valley distances
%         - period                 : time between first and third extrema
%         - is_cycle               : "Y" (valid) or "N" (rejected)
%         - metric, chromosome, iteration
%
%   table_summary
%       Summary table computed ONLY from valid cycles (is_cycle == "Y"),
%       reporting mean oscillation amplitude and period for each:
%         (iteration × chromosome × metric) group.
%       Returns an empty table if no valid cycles are detected.
%
% NOTES
%   - CM and KK oscillations use different detection thresholds reflecting
%     their characteristic amplitudes and timescales.
%   - This function operates on a SINGLE simulation iteration and is intended
%     to be called inside multi-iteration simulation loops.

if nargin < 5, visualize = false; end

% --- Time & sizes ---
Nsteps       = size(cL, 1);
t            = (0:Nsteps-1) * dt;
Nchromosomes = size(cL, 3);

% --- Collectors ---
table_raw_all    = table();
table_cycles_all = table();

for chr = 1:Nchromosomes
    % ---- Signals ----
    CM = (cL(:,1,chr) + cR(:,1,chr)) / 2;
    KK =  cR(:,1,chr) - cL(:,1,chr);

    % ---- CM peak/valley detection (distances in samples) ----
    [CM_peaks,   CM_peak_locs]   = findpeaks(CM,  'MinPeakProminence', 1.00, 'MinPeakDistance', round(1/dt));
    [CM_valleys, CM_valley_locs] = findpeaks(-CM, 'MinPeakProminence', 1.00, 'MinPeakDistance', round(1/dt));
    CM_valleys = -CM_valleys;

    % ---- KK peak/valley detection ----
    [KK_peaks,   KK_peak_locs]   = findpeaks(KK,  'MinPeakProminence', 0.20, 'MinPeakDistance', round(0.5/dt));
    [KK_valleys, KK_valley_locs] = findpeaks(-KK, 'MinPeakProminence', 0.20, 'MinPeakDistance', round(0.5/dt));
    KK_valleys = -KK_valleys;

    % ---- Raw tables (with metric) ----
    T_Raw_CM = table( ...
        [t(CM_peak_locs)'; t(CM_valley_locs)'], ...
        [CM_peaks; CM_valleys], ...
        [repmat("peak",   numel(CM_peaks),   1); repmat("valley", numel(CM_valleys), 1)], ...
        repmat("CM", numel(CM_peaks)+numel(CM_valleys), 1), ...
        repmat(chr,  numel(CM_peaks)+numel(CM_valleys), 1), ...
        repmat(iter, numel(CM_peaks)+numel(CM_valleys), 1), ...
        'VariableNames', {'time','value','type','metric','chromosome','iteration'});
    T_Raw_CM = sortrows(T_Raw_CM,'time');

    T_Raw_KK = table( ...
        [t(KK_peak_locs)'; t(KK_valley_locs)'], ...
        [KK_peaks; KK_valleys], ...
        [repmat("peak",   numel(KK_peaks),   1); repmat("valley", numel(KK_valleys), 1)], ...
        repmat("KK", numel(KK_peaks)+numel(KK_valleys), 1), ...
        repmat(chr,  numel(KK_peaks)+numel(KK_valleys), 1), ...
        repmat(iter, numel(KK_peaks)+numel(KK_valleys), 1), ...
        'VariableNames', {'time','value','type','metric','chromosome','iteration'});
    T_Raw_KK = sortrows(T_Raw_KK,'time');

    table_raw_all = [table_raw_all; T_Raw_CM; T_Raw_KK]; %#ok<AGROW>

    % ---- Cycle tables (minJump per metric) ----
    cycle_CM = extract_cycles(T_Raw_CM, 0.76, 3.18, 1.33,3.58, "CM");
    cycle_KK = extract_cycles(T_Raw_KK, 0.20,1.74, 0.50,1.7 ,"KK");
    table_cycles_all = [table_cycles_all; cycle_CM; cycle_KK]; %#ok<AGROW>

    % ---- Optional viz (ONLY valid cycles) ----
    if visualize
        valid_CM = cycle_CM(cycle_CM.is_cycle=="Y", :);
        valid_KK = cycle_KK(cycle_KK.is_cycle=="Y", :);

        figure('Name',sprintf('Chr %d (Iter %d)', chr, iter));
        plot(t, CM, 'k', 'LineWidth', 1.5); hold on;
        plot(t, KK, 'm', 'LineWidth', 1.5);

        % mark the 3 points for each valid cycle
        for k = 1:height(valid_CM)
            pts = valid_CM.times_3pt{k};
            plot(pts, interp1(t, CM, pts), 'o', 'MarkerFaceColor','g');
        end
        for k = 1:height(valid_KK)
            pts = valid_KK.times_3pt{k};
            plot(pts, interp1(t, KK, pts), 's', 'MarkerFaceColor','c');
        end

        xlabel('Time (min)');
        ylabel('Position (\mum)');
        legend('CM','KK','Valid CM cycle','Valid KK cycle','Location','best');
        grid on;
        title(sprintf('Chromosome %d (Iter %d)', chr, iter));
    end
end

% ---- Summaries (only valid cycles) ----
valid_cycles = table_cycles_all(table_cycles_all.is_cycle=="Y", :);
if isempty(valid_cycles)
    table_summary = table();  % clean empty summary
else
    table_summary = groupsummary(valid_cycles, ...
        {'iteration','chromosome','metric'}, 'mean', ...
        {'amplitude1','amplitude2','period'});
end

% ---- Outputs ----
table_raw   = table_raw_all;
table_cycle = table_cycles_all;

% ================== helper ==================
function cycle_table = extract_cycles(T_raw, minJump,maxJump,minTime,maxTime,metricName)
% Identify cycles from alternating peak/valley triplets.
T_raw = sortrows(T_raw,'time');
n = height(T_raw);
% --- Early return if fewer than 3 extrema (cannot form a cycle) ---
if n < 3
    cycle_table = table( ...
        cell(0,1), strings(0,1), [], [], [], strings(0,1), [], [], strings(0,1), ...
        'VariableNames', {'times_3pt','is_cycle','amplitude1','amplitude2','period', ...
                          'start_type','chromosome','iteration','metric'});
    return;
end
times_cell  = {};
is_cycle    = strings(0,1);
start_type  = strings(0,1);
amp1        = [];
amp2        = [];
period      = [];
chromosomes = [];
iterations  = [];
metrics     = strings(0,1);

i = 1;
while i <= n-2
    t1 = T_raw.time(i);   y1 = T_raw.value(i);   s1 = T_raw.type(i);
    t2 = T_raw.time(i+1); y2 = T_raw.value(i+1); s2 = T_raw.type(i+1);
    t3 = T_raw.time(i+2); y3 = T_raw.value(i+2); s3 = T_raw.type(i+2);

    isAlt = (s1 ~= s2) && (s1 == s3); % P-V-P or V-P-V
    if isAlt
        d1 = abs(y2 - y1);
        d2 = abs(y3 - y2);
        p = t3-t1;
        pass = (d1 >= minJump && d1 <= maxJump) && ...
               (d2 >= minJump && d2 <= maxJump) && ...
               (p  >= minTime && p  <= maxTime);

        times_cell{end+1,1}  = [t1, t2, t3];
        start_type(end+1,1)  = s1;
        chromosomes(end+1,1) = T_raw.chromosome(i);
        iterations(end+1,1)  = T_raw.iteration(i);
        metrics(end+1,1)     = metricName;

        if pass
            is_cycle(end+1,1) = "Y";
            amp1(end+1,1)     = d1/2;
            amp2(end+1,1)     = d2/2;
            period(end+1,1)   = t3 - t1;
        else
            is_cycle(end+1,1) = "N";
            amp1(end+1,1)     = d1/2;
            amp2(end+1,1)     = d2/2;
            period(end+1,1)   = t3 - t1;
        end
    end
    i = i + 1; % overlapping cycles allowed
end

cycle_table = table( ...
    times_cell, is_cycle, amp1, amp2, period, start_type, ...
    chromosomes, iterations, metrics, ...
    'VariableNames', {'times_3pt','is_cycle','amplitude1','amplitude2','period', ...
                      'start_type','chromosome','iteration','metric'});
end

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
