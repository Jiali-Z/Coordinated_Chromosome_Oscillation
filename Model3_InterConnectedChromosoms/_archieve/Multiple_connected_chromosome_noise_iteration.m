% Inter-connected chromosomes 

clc;clear;close all;
% Time Information 
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
Nchromosomes=5; %Number of chromosome pairs 

%Parameter Values for chromosome movement
noise=0.05;
noise_b=0.05;
Kct = 12.3; % pN/µm centromee spring constant, Harasymiw et al, 2019
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
l0=0.5; % natural length of the spring 
koff_0=40; % base disconnection rate
kon_0=20; % base connection rate 
Nspring = 1; % Number of springs per chromosome not used 
k = 0.1; % interchromosome coupling value, make weak?

%Iteration setup 
valid_iterations = 0; % Counter for valid rounds
target_iterations = 100; % Number of valid rounds of data needed

% Create a structure array to store results
records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});
tau_step = round(5 / (dt * 60));  % Number of steps per 5-second interval
tau_values = tau_step * (1:30);   % From 5s to 150s         
tau_groups = reshape(tau_values, 10, [])';  % 3 groups of 10 τs each (3×10 matrix)

while valid_iterations < target_iterations
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
        Alpha(i)=n_dot(i)*6.2/(1-Beta);
    end
    
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
    % Number of attachments
    for i = 1: Nchromosomes 
        NR(1,i)=Nbar(i)*(1-epsilon);
        NL(1,i)=Nbar(i)*(1+epsilon); 
        % Ch position 
        cL(1,:,i) =[-I0/2-epsilon*(2*rand(1)-1),y_positions(i)];%[x,y]
        cR(1,:,i) =[ I0/2+epsilon*(2*rand(1)-1),y_positions(i)];%[x,y]
        % MT tip position 
        xL(1,i) = cL(1,1,i)-0.3;
        xR(1,i) = cR(1,1,i)+0.3;
        vL(1,i)=0;
        vR(1,i)=0;
    end 
    
    
    % ODE solver loop
    for t = 1:Nsteps
        %Update connection state 
        CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);  % N x 1 vector
        YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);  % N x 1 vector
        [flag_c,dist]=spring_connect(CM,YM,flag_cp,Nchromosomes,dt,l0,koff_0,kon_0);
        flag_c_sum_accumulated(t) = sum(sum(flag_c));
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
            % Random Force 
            dx_diff=noise_b *(randn(1))*sqrt(dt);
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
    end
     % Check for NaNs in cL and cR
    if any(isnan(cL),'all') || any(isnan(cR),'all')
        fprintf('NaN detected, skipping this round and starting a new one\n');
        continue; % Skip this round and start a new iteration
    end
    % if (mod(valid_iterations,100)==0)
    %     %Draw afterwards
    %     % side-by-side CM positions overtime 
    %     figure;
    %     hold on;
    %     time = (0:Nsteps) * dt;
    %     for i = 1:Nchromosomes
    %         CMx = (cL(:,1,i) + cR(:,1,i)) / 2;
    %         plot(time, CMx, 'LineWidth', 2);
    %     end
    %     xlabel('Time (min)');
    %     ylabel('Center of Mass X Position (µm)');
    %     title('Chromosome CM Oscillations');
    %     legend(arrayfun(@(i) sprintf('Chr %d', i), 1:Nchromosomes, 'UniformOutput', false));
    %     set(gca, 'FontSize', 14);
    %     grid on;
    % 
    %     % Plot total number of connections over time
    %     figure;
    %     plot((1:Nsteps)*dt, flag_c_sum_accumulated, 'LineWidth', 2);
    %     xlabel('Time (min)', 'FontSize', 18, 'FontName', 'Arial');
    %     ylabel('Total Connections', 'FontSize', 18, 'FontName', 'Arial');
    %     title('Number of Chromosome Connections Over Time', 'FontSize', 20);
    %     set(gca, 'FontSize', 16, 'FontName', 'Arial');
    %     grid on;
    % end
    % Cross-correlation Analysis 
    % Calculate CM in the x axis 
    CMx = (cL(:,1,:) + cR(:,1,:)) / 2;
   iteration_id=valid_iterations + 1;
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
                'C_value', C_norm, ... #Already normalized 
                'N_value', N(neighbor_idx) ...
            );
        end
    end
    % If no NaNs, proceed with correlation calculation
    fprintf('Iteration %d complete\n', valid_iterations + 1);
    valid_iterations = valid_iterations + 1; % Increment valid rounds counte
end
records_table = struct2table(records);

%%
% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
figure; hold on;
colors = lines(size(tau_groups, 1));
for g = 1:size(tau_groups, 1)
    if g ~= 1, continue; end
    current_group = tau_groups(g, :);
    subset = records_table(ismember(records_table.tau, current_group), :);
    grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");

    errorbar(grouped.neighbor, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        "MarkerSize",10,...
        'DisplayName', sprintf('Avg. over:\\tau = %d–%d', (current_group(1)/tau_step)*5, (current_group(end)/tau_step)*5), ...
        'LineWidth', 2, 'Color', "#2171b5");
end
ylim([-0.2, 0.30]);
yticks(-0.2:0.05:0.30);    % y-axis ticks every 0.05
xlim([0.5,4.5]);
xticks(1:1:Nchromosomes-1);
xticklabels({'1st','2nd','3rd','4th'});   
xlabel('Neighbor Level'); ylabel('Avg Cross-Correlation');
%title('C vs Neighbor Level (Grouped by \tau)');
grid on; box on;
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;  
lg = legend('Location','northeast');

% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
figure; hold on;
for n = 1:Nchromosomes-1
    if n ~= 1, continue; end
    subset = records_table(records_table.neighbor == n, :);
    grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, '-o', ...
        'DisplayName', sprintf('%d^{st} neighbor', n), 'LineWidth', 2,...
        'Color',"#238b45");
end
ylim([-0.2, 0.30]);
yticks(-0.2:0.05:0.30);   % y-axis ticks every 0.05
xlim([0,100]);
xticks(0:10:120);
xlabel('\tau (seconds)', 'FontSize', 14);
ylabel('Avg Norm Cross-Corr', 'FontSize', 14);
%title('C vs \tau for 1st, 2nd, 3rd Neighbors', 'FontSize', 16);
legend('Location', 'best');
grid on; box on;
ax = gca;
ax.FontSize = 12;
ax.LineWidth = 1.2;
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.GridLineStyle = '--';
ax.GridAlpha = 0.5;  
lg = legend('Location','northeast');

%% Save .csv to plot in pytyhon 
% ----- Identify this simulation run -----
condition_name = 'TSA';  % change per run: 'WT', 'TSA', 'Noc'

% Convert tau to seconds for easy Python aggregation
records_table.tau_sec = (records_table.tau ./ tau_step) * 5;  % 5,10,...,150 s
records_table.condition = repmat(string(condition_name), height(records_table), 1);

% (optional but recommended) keep tau_index too
records_table.tau_index = records_table.tau;

% Save a single tidy file for this condition
writetable(records_table, sprintf('sim_%s_records.csv', condition_name));
