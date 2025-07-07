% One Oscillating Chromosome 
% See overleaf document "Chromosome Oscillation Models(cleanup)" for model
% details.



clc;clear;close all;

% Time Information 
dt = 4e-3; % Timestep
Nsteps = 5000; % Number of steps Unit:min 

%Parameter Values for chromosome movement
n_dot = 1;
Kct = 30; % pN/µm centromere spring constant, Harasymiw et al, 2019
I0 = 2; % µm Rest length of centromere spring
Kkt = 1; % pN/µm kinetocore spring constant, Cojoc et al. 2016
Gamma = 0.25; % kg/s Drag coefficient
Beta = 0.7; % Scaling factor
Nmax = 25; % Maximum number of attachments
Nbar=20; % Steady state number of MTs when Ch is centered
Lambda=n_dot/(Nbar); % s^-2 KMT detach rate Akioshi et al. 2010
Alpha=n_dot*5/(1-Beta); 
epsilon =0.1; % Small perturbation factor
noise=0;

% Set up the vectors for results
xL = zeros(Nsteps, 1); % Left MT tip position 
xR = zeros(Nsteps, 1); % Right MT tip position 
NL = zeros(Nsteps, 1); % Left MT-KT attachment amount 
NR = zeros(Nsteps, 1); % Right MT-KT attachment amount 
cL = zeros(Nsteps, 2,1); % Left chromosome position (time,axis,Nchromosome) 
cR = zeros(Nsteps, 2,1); % Right chromosome position (time,axis,Nchromosome) 
vL = zeros(Nsteps, 1); % Left chromosome position 
vR = zeros(Nsteps, 1); % Right chromosome position



% Initial Condition 
% Number of attachments
NR(1)=Nbar*(1-epsilon);
NL(1)=Nbar*(1+epsilon); 
% Ch position 
cL(1,:,1) =[-I0/2-epsilon,0];%[x,y]
cR(1,:,1) =[ I0/2+epsilon,0];%[x,y]
% MT tip position 
xbar=(Kct*I0)/(Nbar*Kkt);%Steady state MT tip position 
xL(1) = cL(1)-0.3;
xR(1) = cR(1)+0.3;
vL(1)=0;
vR(1)=0;


% ODE solver loop
for t = 1:Nsteps
    % Force Calculations 
    % Random Force 
    F_noise = noise * randn(1,1) / sqrt(dt); % Brownian Force
    % Define centromere and MT forces on the given chromosome
    F_KT_L = -NL(t) * Kkt * (cL(t) - xL(t)); % MT forces on Left Chromosome
    F_CT_L = Kct * (cR(t) - cL(t) - I0); % Centromere forces on Left Chromosome
    F_KT_R = -NR(t) * Kkt * (cR(t) - xR(t)); % MT forces on Right Chromosome
    F_CT_R = -Kct * (cR(t) - cL(t) - I0); % Centromere forces on Right Chromosome

    % Calculate velocities 
    vL(t) = (F_KT_L + F_CT_L+F_noise)/Gamma;
    vR(t) = (F_KT_R + F_CT_R+F_noise)/Gamma;
    
    % Update Positions 
    % Calculate the current chromosome position in both directions
    cL(t+1,:,1) = cL(t,:,1) + dt * vL(t) ;
    cR(t+1,:,1) = cR(t,:,1) + dt * vR(t) ;
    % Calculate the current microtubule tip position in both directions
    xL(t+1) = xL(t) + dt * Beta * (vL(t));
    xR(t+1) = xR(t) + dt * Beta * (vR(t));
    % Update number of attachments 
    NL(t+1) = NL(t) + dt * (...
        n_dot +(-Alpha * NL(t) * (1 - Beta) * vL(t) * (1 - NL(t)/Nmax)) ...
        - Lambda* NL(t));
    NR(t+1) = NR(t) + dt * (...
        n_dot + (Alpha * NR(t) * (1 - Beta) * vR(t) * (1 - NR(t)/Nmax)) ...
        - Lambda* NR(t));   
end

% % Draw afterwards
% figure(1);
% hold on;
% plot((0:1:Nsteps-1)*dt,cR(1:Nsteps)-cL(1:Nsteps),'c','LineWidth',2);
% plot((0:1:Nsteps-1)*dt,(cL(1:Nsteps)+cR(1:Nsteps))/2,'m','LineWidth',2);
% hold off;
% xlabel({'$\tau$'}, 'Interpreter', 'latex');
% ylabel('kk distance versus Center of Mass');
% set(gca, 'FontSize', 18);
% set(gcf, 'color', 'w');
% grid on; 
% xticks(0:1:20);  % Adjust the range as needed
% 
% figure(2);
% hold on;
% plot((0:Nsteps-1)*dt,NL(1:Nsteps),'c','LineWidth',2);
% plot((0:Nsteps-1)*dt,NR(1:Nsteps),'m','LineWidth',2);
% plot((0:Nsteps-1)*dt,NL(1:Nsteps)+NR(1:Nsteps),'y','LineWidth',2);
% hold off;
% xlabel({'$\tau$'}, 'Interpreter', 'latex');
% ylabel({'$\sigma{L}$'}, 'Interpreter', 'latex');
% set(gca, 'FontSize', 18);
% set(gcf, 'color', 'w');
% grid on; 
% 
% figure(3);
% hold on;
% plot((0:Nsteps-1)*dt,xL(1:Nsteps),'b','LineWidth',2);
% plot((0:Nsteps-1)*dt,xR(1:Nsteps),'r','LineWidth',2);
% hold off;
% xlabel({'$\tau$'}, 'Interpreter', 'latex');
% ylabel({'$\delta{L,R}$'}, 'Interpreter', 'latex');
% set(gca, 'FontSize', 18);
% set(gcf, 'color', 'w');
% 
% figure(4);
% hold on;
% plot(vR(1:Nsteps),vL(1:Nsteps),'b','LineWidth',2);
% hold off;
% xlabel({'$\tau$'}, 'Interpreter', 'latex');
% ylabel({'$\delta{L,R}$'}, 'Interpreter', 'latex');
% set(gca, 'FontSize', 18);
% set(gcf, 'color', 'w');

   
% Oscillation Analysis 
iter = 1; % since this is one simulation run
[table_raw, table_cycle, table_summary] = oscillation_measurement(cL, cR, dt, iter, true);
disp('Summary of Oscillation Amplitude and Period (Mean ± STD):');% Print Summary
disp(table_summary);
