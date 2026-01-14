clc; clear; close all;

% ----------------------- Generate Synthetic Data -----------------------
dt = 2e-3;
Nsteps = 5000;
t = (0:Nsteps-1) * dt;

f_CM = 1/2;
f_KK = 1/1;

amp_CM = 0.8;
amp_KK = 0.3;

CM = amp_CM * sin(2*pi*f_CM*t) ;
KK = amp_KK * sin(2*pi*f_KK*t) ;

% Back-calculate chromatid positions
cL(:,1,1) = CM - KK/2;
cL(:,2,1) = 0;
cR(:,1,1) = CM + KK/2;
cR(:,2,1) = 0;

% ----------------------- Calculate Trajectories -----------------------
CM_x = (cL(:,1,1) + cR(:,1,1)) / 2;
KK_x = cR(:,1,1) - cL(:,1,1);

% ----------------------- Plot Oscillations -----------------------
figure;
plot(t, CM_x, 'k', 'LineWidth', 2); hold on;
plot(t, KK_x, 'm', 'LineWidth', 2);
legend('CM', 'KK');
xlabel('Time (min)');
ylabel('Position (µm)');
title('Synthetic CM and KK Oscillations');
grid on;

% ----------------------- Peak and Valley Detection -----------------------
[CM_peaks, CM_peak_locs] = findpeaks(CM_x, 'MinPeakProminence', 1,'MinPeakDistance',round(1/dt));
[CM_valleys, CM_valley_locs] = findpeaks(-CM_x, 'MinPeakProminence', 1,'MinPeakDistance',round(1/dt)); 
CM_valleys = -CM_valleys;

[KK_peaks, KK_peak_locs] = findpeaks(KK_x, 'MinPeakProminence', 0.3,'MinPeakDistance',round(1/dt));
[KK_valleys, KK_valley_locs] = findpeaks(-KK_x, 'MinPeakProminence', 0.3,'MinPeakDistance',round(1/dt));
KK_valleys = -KK_valleys;

% ----------------------- Visualize Peaks and Valleys -----------------------
plot(t(CM_peak_locs), CM_peaks, 'ko', 'MarkerFaceColor', 'g');
plot(t(CM_valley_locs), CM_valleys, 'ks', 'MarkerFaceColor', 'y');
plot(t(KK_peak_locs), KK_peaks, 'mo', 'MarkerFaceColor', 'c');
plot(t(KK_valley_locs), KK_valleys, 'ms', 'MarkerFaceColor', 'w');

% ----------------------- Set Meta Info -----------------------
iter = 1;
chr = 1;

% ----------------------- Raw Table -----------------------
T_Raw_CM = table(...
    [t(CM_peak_locs)'; t(CM_valley_locs)'], ...
    [CM_peaks; CM_valleys], ...
    [repmat("peak", length(CM_peaks), 1); repmat("valley", length(CM_valleys), 1)], ...
    repmat("CM", length(CM_peaks)+length(CM_valleys), 1), ...
    repmat(chr, length(CM_peaks)+length(CM_valleys), 1), ...
    repmat(iter, length(CM_peaks)+length(CM_valleys), 1), ...
    'VariableNames', {'time', 'value', 'type', 'metric', 'chromosome', 'iteration'});
T_Raw_CM = sortrows(T_Raw_CM , 'time');

T_Raw_KK = table(...
    [t(KK_peak_locs)'; t(KK_valley_locs)'], ...
    [KK_peaks; KK_valleys], ...
    [repmat("peak", length(KK_peaks), 1); repmat("valley", length(KK_valleys), 1)], ...
    repmat("KK", length(KK_peaks)+length(KK_valleys), 1), ...
    repmat(chr, length(KK_peaks)+length(KK_valleys), 1), ...
    repmat(iter, length(KK_peaks)+length(KK_valleys), 1), ...
    'VariableNames', {'time', 'value', 'type', 'metric', 'chromosome', 'iteration'});
T_Raw_KK = sortrows(T_Raw_KK , 'time');
table_raw_all = [T_Raw_CM; T_Raw_KK];

% Cycle_table 
cycle_CM = extract_cycles(T_Raw_CM, 0.76);
cycle_KK = extract_cycles(T_Raw_KK, 0.76);
% Summarty_Table
summary_CM = summarize_cycle_table(cycle_CM);           
summary_KK = summarize_cycle_table(cycle_KK);      


function cycle_table = extract_cycles(T_raw, minJump)
%EXTRACT_CYCLES  Identify oscillation cycles from peak/valley raw table
    % Ensure sorted by time
    T_raw = sortrows(T_raw, 'time');
    n = height(T_raw);
    i = 1;
    % Input for the table 
    times_cell = {};
    is_cycle   = strings(0,1);
    start_type = strings(0,1);
    amp1 = [];amp2 = [];period = [];
    chromosomes = []; iterations = [];
    while i <= n-2
        t1 = T_raw.time(i);   y1 = T_raw.value(i);   s1 = T_raw.type(i);
        t2 = T_raw.time(i+1); y2 = T_raw.value(i+1); s2 = T_raw.type(i+1);
        t3 = T_raw.time(i+2); y3 = T_raw.value(i+2); s3 = T_raw.type(i+2);
        % Must alternate and return to same type
        isAltCycle = (s1 ~= s2) && (s1 == s3);
        % The difference in the jump needs to pass minJump to be considered
        % as an actual cycle 
        if isAltCycle
            d1 = abs(y2 - y1);
            d2 = abs(y3 - y2);
            pass = (d1 >= minJump) && (d2 >= minJump);
            times_cell{end+1,1} = [t1, t2, t3];
            start_type(end+1,1) = s1;
            chromosomes(end+1,1) = T_raw.chromosome(i);
            iterations(end+1,1)  = T_raw.iteration(i);
            if pass
                is_cycle(end+1,1) = "Y";
                amp1(end+1,1)   = d1/2;   % amplitude half jump
                amp2(end+1,1)   = d2/2;
                period(end+1,1) = t3 - t1;
            else
                is_cycle(end+1,1) = "N";
                amp1(end+1,1)   = NaN;
                amp2(end+1,1)   = NaN;
                period(end+1,1) = NaN;
            end
        end
        i = i + 1;  % overlapping cycles
    end
        cycle_table = table( ...
            times_cell, is_cycle, amp1, amp2, period, start_type,chromosomes, iterations, ...
            'VariableNames', {'times_3pt','is_cycle','amplitude1','amplitude2','period', ...
                              'start_type','chromosome','iteration'});
end

function summary_table = summarize_cycle_table(cycle_tbl)
    % Keep the valid cycle 
    valid = cycle_tbl(cycle_tbl.is_cycle=="Y", :);
    % ---- Use groupsummary ----
    summary_table = groupsummary(valid, {'iteration','chromosome'}, {'mean'}, {'amplitude1','amplitude2','period'});
end

%% Scripts 
% ----------------------- Cycle Table -----------------------
% CM osciilation 
i=1; n=height(T_Raw_CM);minJump=0.76;
% Input for the table 
times_cell = {};
is_cycle = strings(0,1);start_type = strings(0,1);    % "peak" or "valley"
amp1= [];amp2= [];period= [];
while i <= n-2
    t1 = T_Raw_CM.time(i);   y1 = T_Raw_CM.value(i);   s1 = T_Raw_CM.type(i);
    t2 = T_Raw_CM.time(i+1); y2 = T_Raw_CM.value(i+1); s2 = T_Raw_CM.type(i+1);
    t3 = T_Raw_CM.time(i+2); y3 = T_Raw_CM.value(i+2); s3 = T_Raw_CM.type(i+2);
    % Must be alternating and start/end on the same type
    isAltCycle = (s1 ~= s2) && (s1 == s3) && (s2 ~= s3);
    if isAltCycle
        % Compute side jump to determine if it's an actual cycle 
        d1 = abs(y2 - y1);    % first half-cycle jump
        d2 = abs(y3 - y2);    % second half-cycle jump
        pass = (d1 >= minJump) && (d2 >= minJump);
        times_cell{end+1,1} = [t1, t2, t3]; 
        start_type(end+1,1) = s1;
        if pass 
            is_cycle(end+1,1)='Y';
            amp1(end+1,1)=d1/2;
            amp2(end+1,1)=d2/2;
            period(end+1,1)=t3-t1;
        else
            is_cycle(end+1,1)='N';
            amp1(end+1,1)=NaN;
            amp2(end+1,1)=NaN;
            period(end+1,1)=NaN;
        end
    end
    % OVERLAPPING windows to capture offset cycles
    i = i + 1;
end

cycle_table = table( ...
    times_cell, is_cycle, amp1, amp2, period, start_type, ...
    'VariableNames', {'times_3pt','is_cycle','amplitude1','amplitude2','period','start_type'});

%%
% ----------------------- Cycle Table -----------------------
%  Amplitude Tables 
N_CM_amp = min(length(CM_peaks), length(CM_valleys));
T_amplitude_CM = table(...
    (1:N_CM_amp)', ...
    repmat(chr, N_CM_amp, 1), ...
    repmat(iter, N_CM_amp, 1), ...
    repmat("CM", N_CM_amp, 1), ...
    (CM_peaks(1:N_CM_amp) - CM_valleys(1:N_CM_amp)) / 2, ...
    'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'amplitude'});

N_KK_amp = min(length(KK_peaks), length(KK_valleys));
T_amplitude_KK = table(...
    (1:N_KK_amp)', ...
    repmat(chr, N_KK_amp, 1), ...
    repmat(iter, N_KK_amp, 1), ...
    repmat("KK", N_KK_amp, 1), ...
    (KK_peaks(1:N_KK_amp) - KK_valleys(1:N_KK_amp)) / 2, ...
    'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'amplitude'});

%  Period Tables 
N_CM_period = length(CM_peak_locs) - 1;
T_period_CM = table(...
    (1:N_CM_period)', ...
    repmat(chr, N_CM_period, 1), ...
    repmat(iter, N_CM_period, 1), ...
    repmat("CM", N_CM_period, 1), ...
    diff(t(CM_peak_locs)'), ...
    'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'period'});

N_KK_period = length(KK_peak_locs) - 1;
T_period_KK = table(...
    (1:N_KK_period)', ...
    repmat(chr, N_KK_period, 1), ...
    repmat(iter, N_KK_period, 1), ...
    repmat("KK", N_KK_period, 1), ...
    diff(t(KK_peak_locs)'), ...
    'VariableNames', {'cycle_id', 'chromosome', 'iteration', 'metric', 'period'});

% Merge CM amplitude and period tables
T_cycles_CM = innerjoin(T_amplitude_CM, T_period_CM, ...
    'Keys', {'cycle_id', 'chromosome', 'iteration', 'metric'});

% Merge KK amplitude and period tables
T_cycles_KK = innerjoin(T_amplitude_KK, T_period_KK, ...
    'Keys', {'cycle_id', 'chromosome', 'iteration', 'metric'});

% Combine into one cycles table
table_cycles = [T_cycles_CM; T_cycles_KK];


% ----------------------- Summary Table -----------------------
table_summary = groupsummary(table_cycles, ...
    {'iteration', 'chromosome', 'metric'}, ...
    ["mean", "std"], ["amplitude", "period"]);
