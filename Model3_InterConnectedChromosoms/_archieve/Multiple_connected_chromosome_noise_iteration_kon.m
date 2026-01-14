clc; clear; close all;

% Simulation Parameters
dt = 2e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 
Nchromosomes = 5; %Number of chromosome pairs 
noise = 0.05;
noise_b = 0.05;
Kct = 12.3;  % pN/µm centromee spring constant 
I0 = 2; % µm Rest length of centromere spring
Kkt = 1; % pN/m kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.1;  % kg/s Drag coefficient
Beta = 0.7;  % Scaling factor
Nmax = 25; % Maximum number of attachments
epsilon = 0.1; % Small perturbation factor

% l0 values to scan
r_values = [0.1 0.2 0.5 1 2 5 10];   % ratio = kon/koff
Ksum = 20;                           % keep kon+koff constant (same timescale)
l0 = 2; % natural length of the spring 
k = 1;  % interchromosome coupling value 
Nspring = 1; % Number of springs per chromosome not used  

% Cross-correlation parameters
tau_step = round(5 / (dt * 60));  % 5s steps
tau_values = tau_step * (1:30);   % From 5s to 150s
tau_groups = reshape(tau_values, 10, [])';

% Store results for all l0
all_records = struct();           % store each table here
labels = strings(numel(r_values),1);              % legend labels
for ri = 1:numel(r_values)
    r = r_values(ri);
    kon_0  = Ksum * (r / (1 + r));
    koff_0 = Ksum * (1 / (1 + r));
    labels(ri) = sprintf('r=%.3g (kon=%.3g, koff=%.3g)', r, kon_0, koff_0);
    fprintf('Running %s...\n', labels(ri));
    valid_iterations = 0;
    target_iterations = 50;
    records = struct('iteration', {}, 'tau', {}, 'neighbor', {}, 'C_value', {}, 'N_value', {});

    while valid_iterations < target_iterations
        % Initialize chromosome parameters
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
        cL = zeros(Nsteps+1, 2, Nchromosomes);
        cR = zeros(Nsteps+1, 2, Nchromosomes);
        xL = zeros(Nsteps+1, Nchromosomes);
        xR = zeros(Nsteps+1, Nchromosomes);
        NL = zeros(Nsteps+1, Nchromosomes);
        NR = zeros(Nsteps+1, Nchromosomes);
        vL = zeros(Nsteps, Nchromosomes);
        vR = zeros(Nsteps, Nchromosomes);
        % Coupling state
        flag_c_sum_accumulated = zeros(Nsteps, 1);
        flag_c = zeros(Nchromosomes);
        flag_cp = flag_c;

        % Initial Condition 
        y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
        for i = 1:Nchromosomes
            % Number of attachments
            NR(1,i) = Nbar(i)*(1 - epsilon);
            NL(1,i) = Nbar(i)*(1 + epsilon);
            % Ch position 
            cL(1,:,i) =-I0/2-epsilon*(2*rand(1)-1);%[x,y]
            cR(1,:,i) =I0/2+epsilon*(2*rand(1)-1);%[x,y]
            cL(:,2,i) =  y_positions(i);
            cR(:,2,i) =  y_positions(i);
            % MT tip position
            xL(1,i) = cL(1,1,i) - 0.3;
            xR(1,i) = cR(1,1,i) + 0.3;
            vL(1,i)=0;
            vR(1,i)=0;
        end
        % ODE solver loop
        for t = 1:Nsteps
            CM = squeeze((cL(t,1,:) + cR(t,1,:)) / 2);
            YM = squeeze((cL(t,2,:) + cR(t,2,:)) / 2);
            [flag_c, ~] = spring_connect(CM, YM, flag_cp, Nchromosomes, dt, l0, koff_0, kon_0);
            flag_c_sum_accumulated(t) = sum(triu(flag_c,1),'all');
            flag_cp = flag_c;
            for i = 1:Nchromosomes
                % Force Calculations 
                %Coupling force 
                F_coupling_x = 0;
                for j = 1:Nchromosomes
                    if j ~= i
                        delta_x = CM(i) - CM(j);
                        F_coupling_x = F_coupling_x - k * delta_x * flag_c(i,j);
                    end
                end
                 % Define centromere and MT forces on the given chromosome
                F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i)); % MT forces on Left Chromosome
                F_CT_L = Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Left Chromosome
                F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i)); % MT forces on Right Chromosome
                F_CT_R = -Kct * (cR(t,1,i) - cL(t,1,i) - I0); % Centromere forces on Right Chromosome
                % Calculate velocities 
                vL(t,i) = (F_KT_L + F_CT_L + F_coupling_x) / Gamma;
                vR(t,i) = (F_KT_R + F_CT_R + F_coupling_x) / Gamma;
                % Update Positions 
                % Calculate the current chromosome position 
                dx_diff=noise_b *(randn(1))*sqrt(dt);
                cL(t+1,1,i) = cL(t,1,i) + dt * vL(t,i) + dx_diff;
                cR(t+1,1,i) = cR(t,1,i) + dt * vR(t,i) + dx_diff;
                cL(t+1,2,i) = cL(1,2,i);  % fixed y
                cR(t+1,2,i) = cR(1,2,i);  % fixed y
                % Calculate the current microtubule tip position in both directions
                xL(t+1,i) = xL(t,i) + dt * Beta * vL(t,i);
                xR(t+1,i) = xR(t,i) + dt * Beta * vR(t,i);
                % Update number of attachments 
                NL(t+1,i) = NL(t,i) + dt * (n_dot(i) + ...
                    (-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) - Lambda(i)* NL(t,i));
                NR(t+1,i) = NR(t,i) + dt * (n_dot(i) + ...
                    (Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) - Lambda(i)* NR(t,i));
            end
        end
        %%%%%%%%%%%%%%%%%% Analysis at the end of each iteration %%%%%%%%%%%%%%
        % Check for NaNs in cL and cR
        if any(isnan(cL), 'all') || any(isnan(cR), 'all')
            fprintf('NaN detected for l0 = %d, skipping\n', l0);
            continue;
        end
        % Crosscorrelation Analysis
        CMx = (cL(:,1,:) + cR(:,1,:)) / 2;
        iteration_id = valid_iterations + 1;
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
                    'N_value', N(neighbor_idx));
            end
        end
        valid_iterations = valid_iterations + 1;
        fprintf('Iteration %d complete for %s\n', valid_iterations, labels(ri));
    end
    fieldname = matlab.lang.makeValidName(sprintf('r_%g', r));
    all_records.(fieldname) = struct2table(records);
end

%% Final Plot — One Plot per τ Group, Overlaying l0 Values
cmap = turbo(numel(r_values));
group_labels = ["τ = 5–50 s", "τ = 55–100 s", "τ = 105–150 s"];
% Plot 1: Avg Cross-Correlation vs Neighbor Level for Each τ Group
for g = 1:size(tau_groups, 1)
    current_group = tau_groups(g, :);
    figure; hold on;
    for li = 1:length(r_values)
        r_val = r_values(li);
        r_key = matlab.lang.makeValidName(sprintf('r_%g', r));;
        
        records_table = all_records.(r_key);
        subset = records_table(ismember(records_table.tau, current_group), :);
        grouped = groupsummary(subset, "neighbor", ["mean","std"], "C_value");
        if isempty(grouped), continue; end

        errorbar(grouped.neighbor, grouped.mean_C_value, grouped.std_C_value, ...
            '-o', 'Color', cmap(li, :), 'MarkerSize',7,...
            'DisplayName', sprintf('r = %d', r_val), ...
            'LineWidth', 2);
    end
    xlabel('Neighbor Level');
    ylabel('Avg Norm Cross-Corr');
    title(sprintf('C vs Neighbor Level (%s)', group_labels(g)));
    legend('Location', 'eastoutside');
    ylim([-0.2, 0.3]);
    grid on; box on;
    set(gca, 'FontSize', 14);
end

% Plot 2: C vs Time Lag for Neighbor Level = 1, 2, 3
for n = 1:Nchromosomes-1
    figure; hold on;
    for li = 1:length(r_values)
        r_val = r_values(li);
        r_key = matlab.lang.makeValidName(sprintf('r_%g', r));;
        
        records_table = all_records.(r_key);
        subset = records_table(records_table.neighbor == n, :);
        grouped = groupsummary(subset, "tau", ["mean","std"], "C_value");
        if isempty(grouped), continue; end

    errorbar((grouped.tau/tau_step)*5, grouped.mean_C_value,grouped.std_C_value, ...
        '-o', 'Color', cmap(li, :),'MarkerSize',7,...
        'DisplayName', sprintf('r = %d', r_val), 'LineWidth', 2);
    end
    xlabel('Time lag (s)');
    ylabel('Avg. Norm Cross-Corr');
    title(sprintf('C vs time lag (neighbor %d)', n));
    legend('Location', 'eastoutside');
    ylim([-0.2, 0.3]);
    grid on; box on;
    set(gca, 'FontSize', 14);
end

%% ===== EXPORT RAW all_records TO CSV =====
% all_records has fields like 'l0_1','l0_2',... each a table with:
% iteration, tau, neighbor, C_value, N_value

raw = table();
fn = fieldnames(all_records);

for f = 1:numel(fn)
    key = fn{f};                 % e.g., 'l0_3'
    T = all_records.(key);
    l0_val = sscanf(key, 'l0_%d');  % parse l0 number from field name
    T.l0 = repmat(l0_val, height(T), 1);
    
    % Optional: write per-l0 files too
    writetable(T, sprintf('raw_records_%s.csv', key));
    
    raw = [raw; T]; %#ok<AGROW>
end

% One combined long CSV (recommended for Python)
writetable(raw, 'effect_of_l0_records_all.csv');
disp('Wrote raw_records_all.csv (and raw_records_l0_*.csv)');

