clc; clear; close all;

% ----------------------- Generate Synthetic Data -----------------------
dt = 4e-3;
Nsteps = 2500;
t = (0:Nsteps-1) * dt;

f_CM = 1;
f_KK = 2;

amp_CM = 0.8;
amp_KK = 0.3;

CM = amp_CM * sin(2*pi*f_CM*t);
KK = amp_KK * sin(2*pi*f_KK*t);

% Back-calculate chromatid positions
cL(:,1,1) = CM - KK/2;
cL(:,2,1) = 0;
cR(:,1,1) = CM + KK/2;
cR(:,2,1) = 0;

%-------------------------Plot chromosome Position Over Time-------------
figure;
plot(t, squeeze(cL(:,1,1)), 'c', 'LineWidth', 2); hold on;
plot(t, squeeze(cR(:,1,1)), 'm', 'LineWidth', 2);
xlabel('Time (min)');
ylabel('X Position (\mum)');
legend('Left Chromatid', 'Right Chromatid');
title('Chromatid X-Positions from Synthetic Data');
grid on;

% ----------------------- Calculate chromosome activity -----------------------
results = [];
duration=10; %default 10mins

[Nsteps, ~, Nchromosomes] = size(cL);
sampling_interval=round(5/60/dt);  % every 5second, ~21
total_steps = min(round(duration / dt), Nsteps); % limit to Nsteps

for idx = 1:Nchromosomes
    signalL = squeeze(cL(1:sampling_interval:total_steps, 1, idx)); % x only
    signalR = squeeze(cR(1:sampling_interval:total_steps, 1, idx)); 
    time_sampled = (0:length(signalL)-1)' * sampling_interval * dt;
    vL_list = [];
    vR_list = [];
    % Sliding window 3-point slope calculation
    for j = 1:(length(signalL) - 2)
        t_chunk = time_sampled(j:j+2);
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
    avg_KE = mean(KE);
    avg_vL = mean(abs(vL_list));
    avg_vR = mean(abs(vR_list));

    % Store result
    results = [results; struct( ...
        'PairID', idx, ...
        'Avg_KE', avg_KE, ...
        'vL', avg_vL, ...
        'vR', avg_vR)];
end
results_table = struct2table(results);

figure;
plot(time_sampled(2:end-1), vL_list, '-k*', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Velocity (\mum/min)');
title('Instantaneous Left Chromosome Velocity');
grid on;

figure;
plot(time_sampled(2:end-1), vR_list, '-k*', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Velocity (\mum/min)');
title('Instantaneous Right Chromosome Velocity');
grid on;

figure;
plot(time_sampled(2:end-1), vL_list.^2 + vR_list.^2, '-k*', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Instantaneous KE');
title('Instantaneous Kinetic Energy Over Time');
grid on;


figure;
bar_data = [results_table.Avg_KE; results_table.vL; results_table.vR];
bar_labels = {'Avg KE', 'Avg |v_L|', 'Avg |v_R|'};
bar(bar_data);
set(gca, 'xticklabel', bar_labels);
ylabel('Value');
title('Chromosome Activity Metrics');
grid on;