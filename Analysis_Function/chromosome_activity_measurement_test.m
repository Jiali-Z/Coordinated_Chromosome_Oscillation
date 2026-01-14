clc; clear; close all;

% ----------------------- Generate Synthetic Data -----------------------
dt = 2e-3;
Nsteps = 5000;
t_min = (0:Nsteps-1) * dt; % time in minuts; 
t_sec = t_min * 60;        % time in seconds;

f_CM = 1/5; % cycles/mins (this is frequency not period!)
f_KK = 1/1; % cycles/mins 

amp_CM = 1.5; %amplitue in um
CM = amp_CM * sin(2*pi*f_CM*t_min);

KK_mean = 2.0;        % target separation (µm)
KK_amp  = 0.2;        % ± variation (keep < KK_mean to stay positive)
KK = KK_mean + KK_amp * sin(2*pi*f_KK*t_min);  % always positive

% Back-calculate chromatid positions
cL(:,1,1) = CM - KK/2;
cL(:,2,1) = 0;
cR(:,1,1) = CM + KK/2;
cR(:,2,1) = 0;

%-------------------------Plot chromosome Position Over Time-------------
figure;
plot(t_sec, squeeze(cL(:,1,1)), 'co-', 'LineWidth', 1, 'Markersize',2,'MarkerFaceColor','c');hold on; 
plot(t_sec, squeeze(cR(:,1,1)), 'mo-', 'LineWidth', 1, 'Markersize',2,'MarkerFaceColor','m');
xlabel('Time (sec)');
ylabel('X Position (\mum)');
legend('Left Chromatid', 'Right Chromatid');
title('Chromatid X-Positions from Synthetic Data');
grid on;

duration=10; %default 10mins
[Nsteps, ~, Nchromosomes] = size(cL);
sampling_interval=round((5/60)/dt);  % every 5second, ~21
total_steps = min(round(duration / dt), Nsteps); % limit to Nsteps

for idx = 1:Nchromosomes
    t_samp  = t_sec(1:sampling_interval:total_steps);
    xL_samp = squeeze(cL(1:sampling_interval:total_steps, 1, idx));;
    xR_samp = squeeze(cR(1:sampling_interval:total_steps, 1, idx));;

    % Plot
    figure('Color','w');
    % --- Plot all points (neutral colors), no lines ---
    plot(t_samp, xL_samp, 'co-', 'Markersize',5,'MarkerFaceColor','c');hold on;
    plot(t_samp, xR_samp, 'mo-', 'Markersize',5,'MarkerFaceColor','m');               
    
    xlabel('Time (min)');
    ylabel('X Position (\mum)');
    legend('L (all points)','R (all points)','L (sampled)','R (sampled)','Location','best');
    title(sprintf('Sampling check (Pair %d): stride = %d frames (~%.1f s)', ...
          idx, sampling_interval, sampling_interval*dt*60));
    grid on;
end

% ----------------------- Calculate chromosome activity -----------------------
results = [];
duration=10; %default 10mins

[Nsteps, ~, Nchromosomes] = size(cL);
sampling_interval=round(5/60/dt);  % every 5second, ~21
total_steps = min(round(duration / dt), Nsteps); % limit to Nsteps

for idx = 1:Nchromosomes
    signalL = squeeze(cL(1:sampling_interval:total_steps, 1, idx)); % x only
    signalR = squeeze(cR(1:sampling_interval:total_steps, 1, idx)); 
    time_sampled_min = (0:numel(signalL)-1)' * (sampling_interval * dt);  % minutes
    time_sampled_sec = time_sampled_min * 60;  % seconds 
    vL_list = [];
    vR_list = [];
    % Sliding window 3-point slope calculation
    for j = 1:(length(signalL) - 2)
        t_chunk = time_sampled_min(j:j+2);
        sL_chunk = signalL(j:j+2);
        sR_chunk = signalR(j:j+2);

        % Fit line: first-order polynomial (slope only)
        pL = polyfit(t_chunk, sL_chunk, 1);
        pR = polyfit(t_chunk, sR_chunk, 1);

        vL_list(end+1) = pL(1);  % slope = velocity
        vR_list(end+1) = pR(1);
    end

    % Kinetic energy and velocity stats
    KE = vL_list.^2 + vR_list.^2;
    avg_KE = mean(KE, 'omitnan');
    avg_vL = mean((vL_list), 'omitnan');
    avg_vR = mean((vR_list), 'omitnan');

    % Store result
    results = [results; struct( ...
        'PairID', idx, ...
        'Avg_KE', avg_KE, ...
        'vL', avg_vL, ...
        'vR', avg_vR)];
end
results_table = struct2table(results);
disp(results_table);

figure;
bar_data = [results_table.Avg_KE; results_table.vL; results_table.vR];
bar_labels = {'Avg KE', 'Avg |v_L|', 'Avg |v_R|'};
bar(bar_data);
set(gca, 'xticklabel', bar_labels);
ylabel('Value');
title('Chromosome Activity Metrics');
grid on;