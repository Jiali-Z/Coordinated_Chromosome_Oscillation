clc; clear; close all;

% --------------------- USER SETTINGS ---------------------
% Kon and Koff grid 
koff0_list  = 0:5:100;     % base disconnection rate(s)
kon0_list   = 0:5:100; % base connection rate(s)

avg_win_sec = 5;            % sliding window (seconds) for time-averaged connectivity
ylim_max    = 10;           % 5 chromosome have maximum of 10 connected pair 

% ---------- Fixed model parameters (your values) ----------
dt = 2e-3;         % min
Nsteps = 5000;     % steps
Nchromosomes = 5;

noise   = 0.05;
noise_b = 0.05;
Kct  = 12.30;      % pN/µm
I0   = 2;          % µm
Kkt  = 1;          % pN/µm
Gamma= 0.1;        % kg/s
Beta = 0.7;
Nmax  = 25;
epsilon = 0.1;

% Inter-chromosomal spring parameters that vary per combo
l0 = 2;    % µm
k  = 1;    % coupling strength

% Time vector (minutes)
tvec = (1:Nsteps)*dt;

% Precompute sliding window size (# of steps in avg_win_sec seconds)
win_steps = round(avg_win_sec / (dt*60));  % dt in min 

% Storage for heatmap (rows=koff_0, cols=kon_0) 
avg_conn = nan(numel(koff0_list), numel(kon0_list));
%progress counter 
nTotal  = numel(kon0_list) * numel(koff0_list);
counter = 0;

% --------- Sweep across (kon_0, koff_0) combinations ---------
for iKon = 1:numel(kon0_list)
    for iKoff = 1:numel(koff0_list)
        kon_0  = kon0_list(iKon);
        koff_0 = koff0_list(iKoff);

        % Storage 
        cL = zeros(Nsteps+1, 2, Nchromosomes); % [time, (x,y), chr]
        cR = zeros(Nsteps+1, 2, Nchromosomes);
        xL = zeros(Nsteps+1, Nchromosomes);
        xR = zeros(Nsteps+1, Nchromosomes);
        NL = zeros(Nsteps+1, Nchromosomes);
        NR = zeros(Nsteps+1, Nchromosomes);
        vL = zeros(Nsteps,   Nchromosomes);
        vR = zeros(Nsteps,   Nchromosomes);

        % Per-chromosome parameter variability
        Nbar   = zeros(Nchromosomes,1);
        n_dot  = zeros(Nchromosomes,1);
        Lambda = zeros(Nchromosomes,1);
        Alpha  = zeros(Nchromosomes,1);

        for i = 1:Nchromosomes
            Nbar(i)   = 20 * (1 + noise*(2*rand-1));
            n_dot(i)  =  1 * (1 + noise*(2*rand-1));
            Lambda(i) = n_dot(i)/Nbar(i);
            Alpha(i)  = n_dot(i) * 6.2 / (1 - Beta);
        end

        % Fixed y layout
        y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);

        % Spring connection flags
        flag_c  = zeros(Nchromosomes,Nchromosomes);
        flag_cp = flag_c;
        flag_c_sum_accumulated = zeros(Nsteps,1);

        % Initial conditions
        for i = 1:Nchromosomes
            NR(1,i) = Nbar(i)*(1 - epsilon);
            NL(1,i) = Nbar(i)*(1 + epsilon);

            cL(1,1,i) = -I0/2 - epsilon*(2*rand-1);  % x
            cR(1,1,i) =  I0/2 + epsilon*(2*rand-1);  % x
            cL(:,2,i) =  y_positions(i);            % y fixed
            cR(:,2,i) =  y_positions(i);

            xL(1,i) = cL(1,1,i) - 0.3;
            xR(1,i) = cR(1,1,i) + 0.3;

            vL(1,i) = 0;
            vR(1,i) = 0;
        end

        % ODE solver
        for t = 1:Nsteps
            % Update connections
            CM = squeeze((cL(t,1,:) + cR(t,1,:))/2);
            YM = squeeze((cL(t,2,:) + cR(t,2,:))/2);

            % Requires your existing spring_connect(...)
            [flag_c, dist] = spring_connect(CM, YM, flag_cp, Nchromosomes, dt, l0, koff_0, kon_0);

            % Track total active connections (upper-triangular to count pairs once)
            flag_c_sum_accumulated(t) = sum(triu(flag_c, 1), 'all');

            % Carry to next step
            flag_cp = flag_c;

            % Per-chromosome forces & updates
            for i = 1:Nchromosomes
                % Inter-chromosomal coupling along x
                F_coupling_x = 0;
                for j = 1:Nchromosomes
                    if j ~= i
                        delta_x = CM(i) - CM(j);
                        F_coupling_x = F_coupling_x - k * delta_x * flag_c(i, j);
                    end
                end

                % Brownian force (kept as in your script; multiplier 0 means disabled)
                F_noise = 0 * noise_b * (randn(1)) / sqrt(dt);

                % Forces
                F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));
                F_CT_L =  Kct      * (cR(t,1,i) - cL(t,1,i) - I0);

                F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));
                F_CT_R = -Kct      * (cR(t,1,i) - cL(t,1,i) - I0);

                % Velocities
                vL(t,i) = (F_KT_L + F_CT_L + F_noise + F_coupling_x) / Gamma;
                vR(t,i) = (F_KT_R + F_CT_R + F_noise + F_coupling_x) / Gamma;

                % Positions (shared jitter in x for both sisters)
                dx_diff = noise_b * (randn(1)) * sqrt(dt);
                cL(t+1,1,i) = cL(t,1,i) + dt*vL(t,i) + dx_diff;
                cR(t+1,1,i) = cR(t,1,i) + dt*vR(t,i) + dx_diff;

                % Keep y fixed
                cL(t+1,2,i) = cL(1,2,i);
                cR(t+1,2,i) = cR(1,2,i);

                % MT tip positions
                xL(t+1,i) = xL(t,i) + dt*Beta*vL(t,i);
                xR(t+1,i) = xR(t,i) + dt*Beta*vR(t,i);

                % Attachments
                NL(t+1,i) = NL(t,i) + dt*( n_dot(i) ...
                    + (-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
                    - Lambda(i)*NL(t,i) );

                NR(t+1,i) = NR(t,i) + dt*( n_dot(i) ...
                    + ( +Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
                    - Lambda(i)*NR(t,i) );
            end
        end

        % Time-averaged connectivity (sliding mean) 
        flag_c_sliding_avg = movmean(flag_c_sum_accumulated, win_steps);
        % Store mean over entire series for heatmap
        avg_conn(iKoff, iKon) = mean(flag_c_sum_accumulated);
        counter = counter + 1;
        fprintf('Progress: %d/%d (%.1f%%) | kon_0=%g, koff_0=%g\n', ...
            counter, nTotal, 100*counter/nTotal, kon_0, koff_0);
        
    end
end

%%
% Heatmap of mean connections over kon_0 × koff_0 
figure('Color','w');
imagesc(kon0_list, koff0_list, avg_conn);   % X = kon_0, Y = koff_0
set(gca,'YDir','normal');                   % low koff_0 at bottom
colormap(turbo); colorbar;
xlabel('kon_0','FontSize',16,'FontName','Arial');
ylabel('koff_0','FontSize',16,'FontName','Arial');
title('Average connections over entire time series','FontSize',18);
set(gca,'FontSize',14,'FontName','Arial',...
    'XTick',kon0_list,'YTick',koff0_list,'TickDir','out','Box','on');
%% Save .csv for python plotting 
[KON, KOFF] = meshgrid(kon0_list, koff0_list);   % size = [numel(koff0_list) x numel(kon0_list)]
T_avg = table( ...
    KON(:), ...
    KOFF(:), ...
    avg_conn(:), ...
    'VariableNames', {'kon_0','koff_0','avg_connections'} ...
);

writetable(T_avg, 'avg_conn_long.csv');   % tidy CSV
save('Kon_Koff_avg_conn_results.mat');  % MATLAB file