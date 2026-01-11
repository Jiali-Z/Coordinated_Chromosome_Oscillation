%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one_oscillating_chromosome_main.m
%
% PURPOSE
%   Simulate oscillatory dynamics of a single sister-chromosome pair in 1D (x).
%   The sisters are connected by a centromere spring and each sister is
%   attached to a kinetochore–microtubule (KT–MT) bundle modeled as a spring
%   with a dynamic number of attachments N(t). 
%
%   This script can be used to locate bifurcation/transition points: the
%   critical centromere stiffness (Kct) at which sustained oscillations
%   emerge or vanish as a function of the velocity–attachment coupling
%   parameter (Alpha).
%
% BIFURCATION/TRANSITION POINTS (for reference)
% NOTE: The values below refer to the "scale factor" used in Alpha = n_dot*scale/(1-Beta).
%   Alpha scale = 3 : No Kct recovers oscillation
%   Alpha scale = 4 : Kct = 15.3 → 15.4
%   Alpha scale = 5 : Kct = 12.9 → 13.0
%   Alpha scale = 6 : Kct = 12.3 → 12.4
%   Alpha scale = 7 : Kct = 12.2 → 12.3
%   Alpha scale = 8 : Kct = 12.4 → 12.5
%   Alpha scale = 9 : Kct = 12.8 → 12.9
%   Alpha scale = 10: Kct = 13.4 → 13.2
%   Alpha scale = 11: transition becomes less sharp
%
% MODEL OVERVIEW 
%   Force balance (overdamped) on each sister kinetochore:
%     (1) Centromere spring force (F_CT)
%     (2) KT–MT spring force (F_KT) scaled by attachment number N(t)
%     (3) Overdamped drag (v = F/Gamma) 
%
%   Attachment dynamics: 
%     N(t) evolves with a baseline growth term (n_dot), a velocity-coupled term
%     controlled by Alpha, and a detachment term controlled by Lambda.
%  
%   Assumptions: 
%     - Microtubule tip positions (xL, xR) move with scaled chromosome velocity (factor Beta)
%     - Initial conditions include a small perturbation (epsilon) to break symmetry
% 
% OUTPUTS
%   Figures:
%     (1) Inter-kinetochore distance (KK = xR_ch − xL_ch) and center of mass
%         (COM = 0.5*(xL_ch + xR_ch)) vs time
%     (2) KT–MT attachments: NL, NR, and NL+NR vs time
%     (3) Left/right chromosome x-position trajectories: xL_ch and xR_ch vs time
%     (4) Peaks and valleys identified by the chromosome oscillation analysis function
%     (5) Bar plot of chromosome activity metrics: Avg_KE, Avg |v_L|, Avg |v_R|
%    
%   Printed tables: 
%     - Oscillation amplitude/period summary (mean ± STD) via oscillation_measurement
%     - Chromosome activity metrics via chromosome_activity_measurement 
%
%   Optional export:
%      - CSV with columns: Step, T (sec), X_Left, X_Right
%     
% DEPENDENCIES
%   - `oscillation_measurement.m`
%   - `chromosome_activity_measurement.m`
%
% NOTES
%   - Motion is 1D along x; y is fixed at 0.
%   - Alpha is calculated from n_dot and Beta unless manually overridden
%
% AUTHOR: 
%   Jiali Zhu, August, 2025, UNC-Chapel Hill 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close all;

% ---------------- Time information ----------------
dt = 2e-3;  % time step (min)
Nsteps = 5000; % number of simulation steps 

% ---------------- Model parameters ----------------
Kct = 12.4;       % pN/µm  centromere spring constant (Transition near Kct ≈ 15.3–15.4 for alpha scale of 4)  
I0 = 2;           % µm     centromere rest length
Kkt = 1;          % pN/µm  KT–MT spring constant 
Gamma = 0.1;      %        effective drag coefficient 
Beta = 0.7;       % dimensionless scaling for MT tip velocity 
Nmax = 25;        % count  max number of KT–MT attachments 
n_dot = 1;        %        baseline attachment gain term
Nbar=20;          % count  attachment number at steady state (used for initialization)
Lambda=n_dot/(Nbar);  %            detachment rate constant
Alpha=n_dot*6/(1-Beta); %          velocity-attachment coupling coefficient (scale factor = 6)
epsilon =0.1;     %         small symmetry-breaking perturbation (dimensionless)

% ---------------- Pre-allocate state variables ----------------
xL = zeros(Nsteps+1, 1);     % left MT tip x-position over time
xR = zeros(Nsteps+1, 1);     % right MT tip x-position over time
NL = zeros(Nsteps+1, 1);     % left KT–MT attachment number over time
NR = zeros(Nsteps+1, 1);     % right KT–MT attachment number over time
cL = zeros(Nsteps+1, 2,1);   % left chromosome position: (time, [x y], chromosome index)
cR = zeros(Nsteps+1, 2,1);   % right chromosome position: (time, [x y], chromosome index)
vL = zeros(Nsteps, 1);       % left chromosome velocity over time
vR = zeros(Nsteps, 1);       % right chromosome velocity over time

% ---------------- Initial conditions ----------------
% Initial attachment numbers (slightly perturbed to break symmetry)
NR(1)=Nbar*(1-epsilon);
NL(1)=Nbar*(1+epsilon); 
% Initial chromosome positions (x in µm; y fixed at 0) 
cL(1,1,1) =-I0/2-epsilon;    % left sister initial x-position
cR(1,1,1) =I0/2+epsilon;     % right sister initial x-position
cL(:,2,1) = 0;               % y fixed
cR(:,2,1) = 0;               % y fixed
% Initial MT tip positions (offset from chromosomes; chosen to start dynamics)
xbar=(Kct*I0)/(Nbar*Kkt);    % predicted steady-state tip offset (for reference; not used later)
xL(1) = cL(1,1,1)-0.3;       % left MT tip initial x-position
xR(1) = cR(1,1,1)+0.3;       % right MT tip initial x-position
% Initial chromosome velocities
vL(1)=0;
vR(1)=0;

% ---------------- Time integration loop (explicit Euler) ----------------
for t = 1:Nsteps
    % ----- Force calculations -----
    % KT–MT spring forces: proportional to extension (chromosome - tip), scaled by N(t)
    F_KT_L = -NL(t) * Kkt * (cL(t,1,1) - xL(t)); % force on left chromosome from KT–MT bundle
    F_KT_R = -NR(t) * Kkt * (cR(t,1,1) - xR(t)); % force on right chromosome from KT–MT bundle
    % Centromere spring forces: proportional to deviation from rest length I0
    F_CT_R = -Kct * (cR(t,1,1) - cL(t,1,1) - I0); % Centromere forces on Right Chromosome
    F_CT_L = Kct * (cR(t,1,1) - cL(t,1,1) - I0);  % Centromere forces on Left Chromosome
    % ----- Velocities (overdamped dynamics) -----
    vL(t) = (F_KT_L + F_CT_L)/Gamma;
    vR(t) = (F_KT_R + F_CT_R)/Gamma;
    
    % ----- Update chromosome positions -----
    cL(t+1,1,1) = cL(t,1,1) + dt * vL(t) ;
    cR(t+1,1,1) = cR(t,1,1) + dt * vR(t) ;
    % ----- Update MT tip positions -----
    % Tip motion follows scaled chromosome motion: dx_tip/dt = Beta*v
    xL(t+1) = xL(t) + dt * Beta * (vL(t));
    xR(t+1) = xR(t) + dt * Beta * (vR(t));
    % ----- Update attachment numbers -----
    % N(t) evolves via baseline growth + velocity-dependent term (with saturation) - detachment
    NL(t+1) = NL(t) + dt * (...
        n_dot +(-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
        - Lambda* NL(t));
    NR(t+1) = NR(t) + dt * (...
        n_dot + (Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
        - Lambda* NR(t));  
    % Optional safety bounds (uncomment to enforce 0 <= N <= Nmax)
    % NL(t+1) = min(max(NL(t+1), 0), Nmax);
    % NR(t+1) = min(max(NR(t+1), 0), Nmax);
end

% ---------------- Derived observables and plotting ----------------
tmin = (0:Nsteps-1)*dt;              % time vector (min)
xL_ch = squeeze(cL(1:Nsteps,1,1));   % left chromosome x(t)
xR_ch = squeeze(cR(1:Nsteps,1,1));   % right chromosome x(t)
KK    = xR_ch - xL_ch;               % inter-kinetochore distance(µm)
COM   = 0.5*(xL_ch + xR_ch);         % center of mass position (µm)
% Figure 1: KK & COM vs time
figure(1); clf; hold on;
plot(tmin, KK, 'c', 'LineWidth', 2);
plot(tmin, COM, 'm', 'LineWidth', 2);
hold off;
xlabel('Time (min)', 'Interpreter','latex');
ylabel('$\mathrm{KK},\,\mathrm{COM}\;(\mu\mathrm{m})$', 'Interpreter','latex');
legend({'KK','COM'}, 'Location','best', 'Box','off');
set(gca,'FontSize',18); grid on; set(gcf,'color','w');

% Figure 2: NL, NR, and total attachments vs time
figure(2); clf; hold on;
plot(tmin, NL(1:Nsteps), 'c', 'LineWidth', 2);
plot(tmin, NR(1:Nsteps), 'm', 'LineWidth', 2);
plot(tmin, NL(1:Nsteps)+NR(1:Nsteps), 'y', 'LineWidth', 2);
hold off;
xlabel('Time (min)', 'Interpreter','latex');
ylabel('Attachments (count)', 'Interpreter','latex');
legend({'NL','NR','NL+NR'}, 'Location','best', 'Box','off');
set(gca,'FontSize',18); grid on; set(gcf,'color','w');


% Figure 3: left/right chromosome x trajectories vs time
figure(3); clf; hold on;
plot(tmin, xL_ch, 'o-', 'LineWidth', 2, 'MarkerSize', 1);
plot(tmin, xR_ch, 'o-', 'LineWidth', 2, 'MarkerSize', 1);
hold off;
xlabel('Time (min)', 'Interpreter','latex');
ylabel('$x_{\mathrm{ch}}\;(\mu\mathrm{m})$', 'Interpreter','latex');
legend({'Left','Right'}, 'Location','best', 'Box','off');
set(gca,'FontSize',18); grid on; set(gcf,'color','w');


% ---------------- Oscillation and activity analysis ----------------
% Chromosome Oscillation Analysis
iter = 1; % since this is one simulation run
[table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, true);
disp('Summary of Oscillation Amplitude and Period (Mean ± STD):');% Print Summary
disp(table_summary);

% Chromosome Activity Analysis
activty_table = chromosome_activity_measurement(cL, cR, dt);
disp('Chromosome Activity Summary:');
disp(activty_table);

% Figure 4: activity metrics bar plot
figure;
bar_data = [activty_table.Avg_KE; activty_table.vL; activty_table.vR];
bar_labels = {'Avg KE', 'Avg |v_L|', 'Avg |v_R|'};

bar(bar_data, 'FaceColor', [0.3 0.3 0.3]);
set(gca, 'XTickLabel', bar_labels, 'FontSize', 14);
ylabel('Value', 'FontSize', 16);
ylim([0,40])
title('Chromosome Activity Metrics', 'FontSize', 18);
grid on;
%% Export to CSV for Python plotting (optional) 
% ---- Save xL and xR to CSV for Python plotting ----
Step   = (0:Nsteps)';            % 0,1,2,...,Nsteps 
t_min = (0:Nsteps)' * dt;        % time in minutes
T_sec = t_min * 60;              % time in seconds 

X_Left  = cL(:,1,1);             % left chromosome x
X_Right = cR(:,1,1);             % right chromosome x

export_tbl = table(Step,T_sec, X_Left, X_Right, ...
    'VariableNames', {'Step','T','X_Left','X_Right'});

writetable(export_tbl, 'chromsome_position_alpha_6.csv');   % write CSV for Python