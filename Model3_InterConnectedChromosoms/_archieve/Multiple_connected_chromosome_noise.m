% Inter-connected chromosomes 

clc;clear;close all;
% Time Information 
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
Nchromosomes=5; %Number of chromosome pairs 

%Parameter Values for chromosome movement
noise=0.05;
noise_b=0.01;
Kct = 11.2; % pN/µm centromee spring constant, Harasymiw et al, 2019
I0 = 2; % µm Rest length of centromere spring
Kkt = 1;% pN/m kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.1; % kg/s Drag coefficient
Beta = 0.7; % Scaling factor
Nmax = 25; % Maximum number of attachments
% n_dot = 1*(1+noise*(2*rand(1)-1)); %Add noise 
% Nbar=20*(1+noise*(2*rand(1)-1)); % Steady state number of MTs when Ch is centered
% Lambda=n_dot/(Nbar);% s^-2 KMT detach rate Akioshi et al. 2010
%Alpha=n_dot*6.3/(1-Beta); 
epsilon =0.1; % Small perturbation factor

%Parameter Values for coupling sping dynamics
l0=2; % natural length of the spring 
koff_0=40; % base disconnection rate
kon_0=80; % base connection rate 
Nspring = 1; % Number of springs per chromosome not used 
k = 2; % interchromosome coupling value, make weak?

%Varying parameters for each chromoosme
Nbar = zeros(Nchromosomes,1); %The steady-state number of MTs when the Ch is centered
for i = 1:Nchromosomes
    Nbar(i) = 20*(1+noise*(2*rand(1)-1));
end
n_dot = zeros(Nchromosomes,1); %The steady-state number of MTs when the Ch is centered
for i = 1:Nchromosomes
    n_dot(i) = 1*(1+noise*(2*rand(1)-1));
end
Lambda = zeros(Nchromosomes, 1);% s^-2 KMT detach rate Akioshi et al. 2010
for i = 1:Nchromosomes
    Lambda(i)=n_dot(i)/Nbar(i);
end
Alpha = zeros(Nchromosomes, 1);
for i = 1:Nchromosomes
    Alpha(i)=n_dot(i)*3.6/(1-Beta);
end
fprintf("Ndot = %4.4f\n", n_dot);
fprintf("Nbar = %4.4f\n", Nbar);
fprintf("Lambda = %4.4f\n", Lambda);
fprintf("Alpha = %4.4f\n", Alpha);

% Set up the vectors for results
cL = zeros(Nsteps+1, 2,Nchromosomes); % Left chromosome position (time,axis,Nchromosome) 
cR = zeros(Nsteps+1, 2,Nchromosomes); % Right chromosome position (time,axis,Nchromosome) 
xL = zeros(Nsteps+1, Nchromosomes); % Left MT tip position 
xR = zeros(Nsteps+1, Nchromosomes); % Right MT tip position 
NL = zeros(Nsteps+1, Nchromosomes); % Left MT-KT attachment amount 
NR = zeros(Nsteps+1, Nchromosomes); % Right MT-KT attachment amount 
vL = zeros(Nsteps, Nchromosomes); % Left chromosome position 
vR = zeros(Nsteps, Nchromosomes); % Right chromosome position



% Initial Condition 
% Fixed y positions
y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
% Coupling state
flag_c=zeros(Nchromosomes,Nchromosomes);  % current state, update at each timestep 
flag_cp=flag_c; % Previous state 
flag_c_sum_accumulated = zeros(Nsteps, 1); % To accumulate sum of flag_c over time
% Each Pairs have their own 
for i = 1: Nchromosomes 
    % Number of attachments
    NR(1,i)=Nbar(i)*(1-epsilon);
    NL(1,i)=Nbar(i)*(1+epsilon); 
    % Ch position 
    cL(1,1,i) =-I0/2-epsilon*(2*rand(1)-1);%[x,y]
    cR(1,1,i) =I0/2+epsilon*(2*rand(1)-1);%[x,y]
    cL(:,2,i) = y_positions(i);
    cR(:,2,i) = y_positions(i);
    % MT tip position 
    xL(1,i) = cL(1,1,i)-0.3;
    xR(1,i) = cR(1,1,i)+0.3;
    vL(1,i)=0;
    vR(1,i)=0;
end 

% Save video 
% video_filename = 'Connected_chromosome_oscillation_simulation.mp4';
% v = VideoWriter(video_filename, 'MPEG-4'); % You can also use 'Motion JPEG AVI'
% v.FrameRate = 30;  % Adjust frame rate as needed
% open(v);
% ODE solver loop
for t = 1:Nsteps
    %Update connection state 
    CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);  % N x 1 vector
    YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);  % N x 1 vector
    [flag_c,dist]=spring_connect(CM,YM,flag_cp,Nchromosomes,dt,l0,koff_0,kon_0);
    flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
    flag_cp = flag_c;  % Update prior state
    for i=1:Nchromosomes 
        % Force Calculations
        %Coupling force 
        F_coupling_x=0;
        for j = 1:Nchromosomes
            if j ~= i  % Skip self
                delta_x = CM(i) - CM(j);
                F_coupling_x = F_coupling_x - k * delta_x * flag_c(i, j);
            end
        end
        % Define centromere and MT forces on the given chromosome
        F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i)); % MT forces on Left Chromosome
        F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Left Chromosome
        F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i)); % MT forces on Right Chromosome
        F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Right Chromosome
    
        % Calculate velocities 
        vL(t,i) = (F_KT_L + F_CT_L+F_coupling_x)/Gamma;
        vR(t,i) = (F_KT_R + F_CT_R+F_coupling_x)/Gamma;
        
        % Update Positions 
        % Calculate the current chromosome position in both directions
        dx_diff=noise_b *(randn(1))*sqrt(dt);
        cL(t+1,1,i) = cL(t,1,i) + dt * vL(t,i)+dx_diff ;
        cR(t+1,1,i) = cR(t,1,i) + dt * vR(t,i)+dx_diff ;
        % Maintain fixed Y positions
        cL(t+1,2,i) = cL(1,2,i);
        cR(t+1,2,i) = cR(1,2,i);
        % Calculate the current microtubule tip position in both directions
        xL(t+1,i) = xL(t,i) + dt * Beta * (vL(t,i));
        xR(t+1,i) = xR(t,i) + dt * Beta * (vR(t,i));
        % Update number of attachments 
        NL(t+1,i) = NL(t,i) + dt * (...
            n_dot(i) +(-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
            - Lambda(i)* NL(t,i));
        NR(t+1,i) = NR(t,i) + dt * (...
            n_dot(i) + (Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
            - Lambda(i)* NR(t,i));  
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Video %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if (mod(t,20)==0)
    %     cla; % Clear the plot content but keep the axis limits
    %     hold on; % Hold on to allow multiple plots
    %     for idx=1:Nchromosomes
    %         %Plot the chromosome location 
    %         plot(cL(t,1,idx), cL(t,2,idx), 'ro'); hold on; % Left Chromosome
    %         plot(cR(t,1,idx), cR(t,2,idx), 'bo'); hold on; % Right chromosome
    %         %Plot the microtubule tips (xL and xR)
    %         plot(xL(t,idx), cL(t,2,idx), 'rx', 'MarkerSize', 10); % Left microtubule tip
    %         plot(xR(t,idx), cR(t,2,idx), 'bx', 'MarkerSize', 10); % Right microtubule tip
    %         %Draw a single line and label NL if NL is not zero
    %         if NL(t, idx) > 0
    %             line([cL(t,1,idx), xL(t, idx)], [cL(t,2,idx), cL(t,2,idx)], 'Color', 'r', 'LineWidth', 1.5);
    %             text((cL(t,1,idx) + xL(t, idx)) / 2, cL(t,2,idx), sprintf('%.1f', NL(t, idx)), ...
    %                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'r', 'FontSize', 10);
    %         end
    %         %Draw a single line and label NR if NR is not zero
    %         if NR(t, idx) > 0
    %             line([cR(t,1,idx), xR(t, idx)], [cR(t,2,idx), cR(t,2,idx)], 'Color', 'b', 'LineWidth', 1.5);
    %             text((cR(t,1,idx) + xR(t, idx)) / 2, cR(t,2,idx), sprintf('%.1f', NR(t, idx)), ...
    %                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Color', 'b', 'FontSize', 10);
    %         end
    %         %Draw coupling springs between connected chromosomes 
    %         for jdx = 1:Nchromosomes
    %             if flag_c(idx, jdx) == 1 && idx < jdx  % avoid drawing each connection twice
    %                 % Compute center-of-mass coordinates
    %                 CM_i_x = (cL(t,1,idx) + cR(t,1,idx)) / 2;
    %                 CM_i_y = (cL(t,2,idx) + cR(t,2,idx)) / 2;
    %                 CM_j_x = (cL(t,1,jdx) + cR(t,1,jdx)) / 2;
    %                 CM_j_y = (cL(t,2,jdx) + cR(t,2,jdx)) / 2;
    % 
    %                 % Draw the connection spring
    %                 plot([CM_i_x, CM_j_x], [CM_i_y, CM_j_y], 'k-', 'LineWidth', 2);
    %             end
    %         end
    % 
    %     end
    %     set(gca, 'XLim', [-5, 5], 'YLim', [-5, 5]);
    %     drawnow;
    %     %Capture frame and write to video
    %     frame=getframe(gcf);
    %     writeVideo(v,frame);
    %     hold off;
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%Draw afterwards
% side-by-side CM positions overtime 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    plot(time, CMx, 'LineWidth', 2);
end
xlabel('Time (min)');
ylabel('Center of Mass X Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% Movement Trajectories 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    CMy = (cL(:,2,i) + cR(:,2,i)) / 2;
    plot(CMx, CMy, 'LineWidth', 2);
end
xlabel('Center of Mass X Position (µm)');
ylabel('Center of Mass Y Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% close(v);
% fprintf('Video saved to %s\n', video_filename);

%%
%Draw afterwards
% side-by-side CM positions overtime 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    plot(time, CMx, 'LineWidth', 2);
end
xlabel('Time (min)');
ylabel('Center of Mass X Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% Movement Trajectories 
figure;
hold on;
time = (0:Nsteps) * dt;
for i = 1:Nchromosomes
    CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    CMy = (cL(:,2,i) + cR(:,2,i)) / 2;
    plot(CMx, CMy, 'LineWidth', 2);
end
xlabel('Center of Mass X Position (µm)');
ylabel('Center of Mass Y Position (µm)');
title('Chromosome CM Oscillations');
legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
set(gca, 'FontSize', 14);
grid on;

% Plot total number of connections over time
figure;
plot((1:Nsteps)*dt, flag_c_sum_accumulated, 'LineWidth', 2);
xlabel('Time (min)', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Total Connections', 'FontSize', 18, 'FontName', 'Arial');
title('Number of Chromosome Connections Over Time', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Arial');
grid on;


% Plot sliding average of total number of connections over time
figure;
window_size = round(5 / (dt * 60));  % steps per 5-second interval
flag_c_sliding_avg = movmean(flag_c_sum_accumulated, window_size);
% Plot sliding average over time
plot((1:Nsteps)*dt, flag_c_sliding_avg, 'LineWidth', 2, 'Color', 'b');
xlabel('Time (min)', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Total Connections (5s average)', 'FontSize', 18, 'FontName', 'Arial');
title('Number of Chromosome Connections (Sliding Window Average)', 'FontSize', 20);
set(gca, 'FontSize', 16, 'FontName', 'Arial');
ylim([0,12]);
grid on;
hold on; 

%% Output data for plotting in python 
window_size = round(5 / (dt * 60));  % steps per 5-second interval
flag_c_sliding_avg = movmean(flag_c_sum_accumulated, window_size);
time_min = (1:Nsteps)' * dt;  % time vector aligns with flag_c_* which runs 1..Nsteps
T = table( ...
    time_min, ...
    flag_c_sum_accumulated(:), ...
    flag_c_sliding_avg(:), ...
    'VariableNames', {'Time_min','TotalConnections','TotalConnections_5sAvg'} ...
);

out_csv = 'kon80_koff_40_connections_over_time.csv';
writetable(T, out_csv);

%% Output data for plotting side-by-side CM postions 
time_min = (0:Nsteps)' * dt;     % includes t=0
nT = numel(time_min);

% Build long-form arrays: one row per (time, chromosome)
Time_col  = repmat(time_min, Nchromosomes, 1);
PairID_col = repelem((1:Nchromosomes)', nT, 1);

CM_X_all = [];
for i = 1:Nchromosomes
    CMx_i = (cL(:,1,i) + cR(:,1,i)) / 2;     % (Nsteps+1) x 1
    CM_X_all = [CM_X_all; CMx_i];
end

T_cm = table(Time_col, PairID_col, CM_X_all, ...
    'VariableNames', {'Time_min','PairID','CM_X'});

out_csv_cm = 'sim_CM_over_time.csv';
writetable(T_cm, out_csv_cm);
fprintf('Wrote %s with columns: Time_min, PairID, CM_X\n', out_csv_cm);


%% Analysis 
   
% Oscillation Analysis 
all_summary = [];
for i = 1:Nchromosomes
    [~, ~, summary_i] = oscillation_measurement(cL(:,:,i), cR(:,:,i), dt, i, true);
    all_summary = [all_summary; summary_i];
end
disp('Per-Chromosome Oscillation Summary:');
disp(all_summary);

% Chromosome Activity Analysis
activity_all = [];
for i = 1:Nchromosomes
    table_i = chromosome_activity_measurement(cL(:,:,i), cR(:,:,i), dt);
    activity_all = [activity_all; table_i];
end
disp('Per-Chromosome Activity Summary:');
disp(activity_all);

figure;
bar_data = [activity_all.Avg_KE, activity_all.vL, activity_all.vR];
bar(bar_data, 'grouped');
legend({'Avg KE', '|v_L|', '|v_R|'});
xlabel('Chromosome Index');
ylabel('Value');
title('Chromosome Activity');
grid on;

% Cross-correlation Analysis 
% Create a structure array to store results
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s         
tau_groups = reshape(tau_values, 10, [])';  % 3 groups of 10 τs each (3×10 matrix)
% Calculate CM in the x axis 
CMx = (cL(:,1,:) + cR(:,1,:)) / 2;
iteration_id=1;
% Compute correlations at selected tau values
for tau = tau_values
    [C, N] = crosscorrelation_measurement(CMx, tau, Nsteps, Nchromosomes);

    for neighbor_idx = 1:(Nchromosomes - 1)
        if N(neighbor_idx) > 0
            C_norm = C(neighbor_idx) / N(neighbor_idx);
        else
            C_norm = NaN;
        end

        records(end+1) = struct( ...
            'iteration', iteration_id, ...
            'tau', tau, ...
            'neighbor', neighbor_idx, ...
            'C_value', C_norm, ...
            'N_value', N(neighbor_idx) ...
        );
    end
end
records_table = struct2table(records);

% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', colors(g, :));
end
ylim([-1,1]);
xticks(1:1:Nchromosomes-1);
xlabel('Neighbor Level'); ylabel('Avg Norm Cross-Corr');
title('C vs Neighbor Level (Grouped by \tau)');
legend; grid on; box on;

% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:Nchromosomes-1
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2);
end
ylim([-1,1]);
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
title('C vs \tau for 1st, 2nd, 3rd Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
set(gca, 'FontSize', 12);
