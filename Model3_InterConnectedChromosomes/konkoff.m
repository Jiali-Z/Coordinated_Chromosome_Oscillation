
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE (Figure S4B)
%   For a given choice of base inter-chromosomal coupling kinetics
%   (kon_0, koff_0), this script visualizes how the connection (kon) 
%   and disconnection (koff) rates vary as a function of inter-chromosomal
%   distance. 
%
% USAGE
%   Input any desired kon_0 and koff_0 values to generate the corresponding
%   kon(distance) and koff(distance) curves shown in Figure S4B.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the parameters
l0 = 2;
koff_0=20; % base disconnection rate
kon_0=40; % base connection rate 
dt=2e-3;
% Create a range of distances
distances = linspace(0, 5, 100);  % From 0 to 2 with 100 points

% Calculate kon and koff based on the distance
kon = kon_0 * exp(-(distances ./ l0).^2)*dt;
koff = koff_0 * exp((distances ./ l0).^2)*dt;

% Plotting the results
figure;
plot(distances, kon, 'b-', 'LineWidth', 2);
xlabel('Distance');
ylabel('kon');
title('kon vs. Distance');
grid on;

figure;
plot(distances, koff, 'r-', 'LineWidth', 2);
xlabel('Distance');
ylabel('koff');
title('koff vs. Distance');
grid on;
