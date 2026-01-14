%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kon_Koff_Sweep_Connectivity.m
%
% PURPOSE
%   Run a stochastic (multi-iteration; noisy) parameter sweep over the
%   inter-chromosomal spring network *kinetics*:
%       (1) kon_0  -> base connection rate (formation) for inter-chromosomal links
%       (2) koff_0 -> base disconnection rate (breakage) for inter-chromosomal links
%   while simulating N sister-chromosome pairs that oscillate in 1D (x) using
%   the same force-balance oscillation dynamics as the single-chromosome model.
%
%   For each (kon_0, koff_0) pair, the script run simulation once and  updates
%   the stochastic spring connectivity matrix over time via spring_connect(...), 
%   tracks the total number of active chromosome–chromosome connections 
%   at each time step, and summarizes connectivity as a time-averaged metric 
%   for heatmap visualization.
%
% NOTE / PAPER NOTATION
%   Same notations are used in paper 
%
% HOW THIS SCRIPT IS USED IN THE PAPER
%   - Figure S4A: Parameter-phase heatmap of  number of active 
%     chromosome–chromosome connection on the (kon_0, koff_0) grid
%
% NOTES / IMPLEMENTATION DETAILS
%   To facilitate reproducibility without rerunning the full sweep, the
%   processed CSV output has been COMMITTED to the repository:Users can 
%   directly load this CSV file to reproduce Figure 4B without re-running
%   this script.
%
%
% OUTPUTS
%   1) avg_conn heatmap matrix 
%        avg_conn(iKoff, iKon) = mean # connected chromosome pairs over time
%   2) CSV export (tidy / long format) for downstream plotting:
%        avg_conn_long.csv with columns: kon_0, koff_0, avg_connections
%
%
% DEPENDENCIES
%   spring_connect.m
%
% AUTHOR
%   Jiali Zhu, 2025, UNC–Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;


% ------------------ Sweep definitions ------------------
koff0_list  = 0:5:100;         % base disconnection rate(s)
kon0_list   = 0:5:100;         % base connection rate(s)
nTotal  = numel(kon0_list) * numel(koff0_list);
counter = 0;

% ---------------- Simulation settings ----------------
dt = 2e-3;                     % Timestep (min)
Nsteps = 5000;                 % Number of steps
tvec = (1:Nsteps)*dt;          % Time vector (minutes)

% ---------------- Number of chromosome ---------------- 
Nchromosomes = 5;              % Number of chromosome pairs

% ------------------ Stochasticity parameters ------------------
noise   = 0.05;            % variability for n_dot and Nbar
noise_b = 0.05;             % Brownian noise strength (diffusive displacement)

% ------------------ Fixed model constants ------------------
Kct  = 12.30;                % pN/µm: centromere spring constant
I0   = 2;                    % µm: centromere rest length
Kkt  = 1;                    % pN/µm: KT–MT spring constant
Gamma= 0.1;                  % effective drag coefficient
Beta = 0.7;                  % MT tip motion coupling factor (x tip follows chromosome)
Nmax  = 25;                  % max attachments per kinetochore
epsilon = 0.1;               % small symmetry-breaking perturbation (dimensionless)   

% ---------------- Inter-chromosome coupling parameters ----------------
% koff_0=40;                  % 1/s    base disconnection rate 
% kon_0=20;                   % 1/s    base disconnection rate
l0 = 2;                       % µm     inter-chromosomal spring effective length 
k  = 1;                       % pN/µm: inter-chromosomal spring constant 


% ------------------ analysis parameter ------------------
avg_win_sec = 5;                             % sliding window (seconds) for time-averaged connectivity
win_steps = round(avg_win_sec / (dt*60));    % # of steps in avg_win_sec seconds
ylim_max    = 10;                            % 5 chromosome have maximum of 10 connected pair 


% ------------------ Storage for results ------------------
avg_conn = nan(numel(koff0_list), numel(kon0_list));


% ========================================================================
%                         PARAMETER SWEEP LOOP
% =========================================================================
for iKon = 1:numel(kon0_list)
    for iKoff = 1:numel(koff0_list)
        kon_0  = kon0_list(iKon);
        koff_0 = koff0_list(iKoff);

        % ---- cell-specific parameters with noise ----
        Nbar   = zeros(Nchromosomes, 1);  % steady-state attachment number (initialization reference)
        n_dot  = zeros(Nchromosomes, 1);  % baseline gain term for attachments
        Lambda = zeros(Nchromosomes, 1);  % detachment rate constant
        Alpha  = zeros(Nchromosomes, 1);  % velocity-attachment coupling coefficient
        for i = 1:Nchromosomes
            Nbar(i)  = 20 * (1 + noise * (2*rand(1) - 1));
            n_dot(i) =  1 * (1 + noise * (2*rand(1) - 1));
            Lambda(i) = n_dot(i) / Nbar(i);
            Alpha(i)  = n_dot(i) * 6.2 / (1 - Beta); alpha_effective = mean(Alpha); 
        end 

        % ------------------ Allocate state variables ------------------
        cL = zeros(Nsteps+1, 2, Nchromosomes);       % left chromosome position: (time, [x y], chromosome index)
        cR = zeros(Nsteps+1, 2, Nchromosomes);       % right chromosome position: (time, [x y], chromosome index)
        xL = zeros(Nsteps+1, Nchromosomes);          % left MT tip x-position over time
        xR = zeros(Nsteps+1, Nchromosomes);          % right MT tip x-position over time
        NL = zeros(Nsteps+1, Nchromosomes);          % left KT–MT attachment number over time 
        NR = zeros(Nsteps+1, Nchromosomes);          % right KT–MT attachment number over time 
        vL = zeros(Nsteps, Nchromosomes);            % left chromosome velocity over time 
        vR = zeros(Nsteps, Nchromosomes);            % right chromosome velocity over time

        % ---------------- Initial conditions ----------------
        % Spring connectivity state (updated each time step)
        flag_c = zeros(Nchromosomes);                   % current adjacency (1 = connected)
        flag_cp = flag_c;                               % Previous state 
        flag_c_sum_accumulated = zeros(Nsteps, 1);      % To accumulate sum of flag_c over time
        % Fixed y positions are assigned 
        y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
        % Each Pairs have their own 
        for i = 1:Nchromosomes
            % Initial attachment numbers (slightly perturbed to break symmetry)
            NR(1,i) = Nbar(i)*(1 - epsilon);
            NL(1,i) = Nbar(i)*(1 + epsilon);
            % Initial chromosome positions (x in µm; y fixed) 
            cL(1,:,i) = -I0/2 - epsilon*(2*rand(1)-1);   % left sister initial x-position
            cR(1,:,i) =  I0/2 + epsilon*(2*rand(1)-1);   % right sister initial x-position
            cL(:,2,i) =  y_positions(i);                 % y fixed
            cR(:,2,i) =  y_positions(i);                 % y fixed
            % Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
            xL(1,i) = cL(1,1,i) - 0.3;
            xR(1,i) = cR(1,1,i) + 0.3;
            % Initial chromosome velocities
            vL(1,i)=0; vR(1,i)=0;
        end

       % ---------------- Time integration loop (explicit Euler) ----------------
        for t = 1:Nsteps
            % ----- Update stochastic spring network  -----
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
                % ----- Force calculations -----
                % Inter-chromosomal Coupling force 
                F_coupling_x = 0;
                for j = 1:Nchromosomes
                    if j ~= i
                        delta_x = CM(i) - CM(j);
                        F_coupling_x = F_coupling_x - k * delta_x * flag_c(i, j);
                    end
                end
                % Define centromere and MT forces on the given chromosome
                F_KT_L = -NL(t,i) * Kkt * (cL(t,1,i) - xL(t,i));     % force on left chromosome from KT–MT bundle
                F_CT_L =  Kct     * (cR(t,1,i) - cL(t,1,i) - I0);    % Centromere forces on Left Chromosome
                F_KT_R = -NR(t,i) * Kkt * (cR(t,1,i) - xR(t,i));     % force on right chromosome from KT–MT bundle
                F_CT_R = -Kct     * (cR(t,1,i) - cL(t,1,i) - I0);    % Centromere forces on Right Chromosome
                 % ----- Velocities (overdamped dynamics) -----
                vL(t,i) = (F_KT_L + F_CT_L + F_coupling_x) / Gamma;
                vR(t,i) = (F_KT_R + F_CT_R + F_coupling_x) / Gamma;
                % ----- Update chromosome positions -----
                % Common-mode positional jitter: the same dx_diff is added to both sisters,
                dx_diff = noise_b * (randn(1)) * sqrt(dt);
                cL(t+1,1,i) = cL(t,1,i) + dt*vL(t,i) + dx_diff;
                cR(t+1,1,i) = cR(t,1,i) + dt*vR(t,i) + dx_diff;
                % Keep y fixed
                cL(t+1,2,i) = cL(1,2,i);
                cR(t+1,2,i) = cR(1,2,i);
                % ----- Update MT tip positions -----
                xL(t+1,i) = xL(t,i) + dt*Beta*vL(t,i);
                xR(t+1,i) = xR(t,i) + dt*Beta*vR(t,i);
                % ----- Update attachment numbers -----
                NL(t+1,i) = NL(t,i) + dt*( n_dot(i) ...
                    + (-Alpha(i) * NL(t,i) * (1 - Beta) * vL(t,i) * (1 - NL(t,i)/Nmax)) ...
                    - Lambda(i)*NL(t,i) );

                NR(t+1,i) = NR(t,i) + dt*( n_dot(i) ...
                    + ( +Alpha(i) * NR(t,i) * (1 - Beta) * vR(t,i) * (1 - NR(t,i)/Nmax)) ...
                    - Lambda(i)*NR(t,i) );
            end
        end
        %%%%%%%%%%%%%%%%%% Analysis at the end of each iteration %%%%%%%%%%%%%%
        % ---------------Time-averaged connectivity (sliding mean)-------- 
        flag_c_sliding_avg = movmean(flag_c_sum_accumulated, win_steps);
         % Append the results 
        avg_conn(iKoff, iKon) = mean(flag_c_sum_accumulated);
        % --- --------------------Progress note ---------------
        counter = counter + 1;
        fprintf('Progress: %d/%d (%.1f%%) | kon_0=%g, koff_0=%g\n', ...
            counter, nTotal, 100*counter/nTotal, kon_0, koff_0);    
    end
end

%% ========================================================================
%   FIGURE S4A: HEATMAP  
% =========================================================================
% Heatmap of mean connections over kon_0 × koff_0 

% --- Load sweep results from CSV (exported earlier) ---
csvFile = "Output_Sweep_Bifurcation/Kon_Koff_avg_conn_sweep.csv";
results_all = readtable(csvFile);
% Extract unique parameter grids
kon0_list  = unique(results_all.kon_0,  'sorted');
koff0_list = unique(results_all.koff_0, 'sorted');
% Preallocate heatmap matrix
avg_conn = nan(numel(koff0_list), numel(kon0_list));
% Fill matrix: rows = koff_0, cols = kon_0
for i = 1:height(results_all)
    iKon  = find(kon0_list  == results_all.kon_0(i));
    iKoff = find(koff0_list == results_all.koff_0(i));
    avg_conn(iKoff, iKon) = results_all.avg_connections(i);
end
% --- Plot heatmap ---
figure('Color','w');
imagesc(kon0_list, koff0_list, avg_conn);   % X = kon_0, Y = koff_0
set(gca,'YDir','normal');                   % low koff_0 at bottom
colormap(parula); colorbar;
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
