function [C, N] = crosscorrelation_measurement(CM, tau, Nsteps, Nchromosomes)
% CROSSCORRELATION_MEASUREMENT Compute lagged cross-correlation between chromosome pairs 
%
% PURPOSE:
%   Computes the time-lagged cross-correlation between the x-positions of a pair of the centroids  
%   of chromosome pair at a given lag (tau), grouped by their neighbor level 
%   (i.e., distance in chromosome index: 1st neighbor, 2nd neighbor, etc.).
%   This is used to quantify how the movement of chromosomes is coordinated
%   as a function of their relative spatial arrangement.
%
% INPUTS:
%   CM              double matrix of size [Nsteps × 2 x Nchromosomes]
%                   X-position of the center of mass of each chromosome pair over time.               
%
%   tau             scalar integer
%                   Time lag (in steps) to compute the correlation at.
%
%   Nsteps          scalar integer
%                   Total number of time steps 
%
%   Nchromosomes    scalar integer
%                   Total number of chromosomes.
%
% OUTPUTS:
%   C               1 × (Nchromosomes - 1) vector
%                   Accumulated cross-correlation values for each neighbor level.
%                  
%   N               1 × (Nchromosomes - 1) vector
%                   Number of contributions used to compute and accumulate C at each neighbor level. 

% NOTE:
%   - Output C is **not normalized**; average cross-correlation:  C ./ N 
%   - Only pairs where j > i are considered (upper triangle).
%   - At each time step t, the function computes:
%         (CM(t+tau, i) - CM(t, i)) * (CM(t+tau, j) - CM(t, j))
%     which captures whether two chromosomes move in the same or opposite direction
%     over the time interval [t → t+tau].
    
    CM = squeeze(CM);
    % Initialize outputs
    max_neighbor = Nchromosomes - 1;
    C = zeros(1, max_neighbor);
    N = zeros(1, max_neighbor);

    % Loop over time steps
    for t = 1:(Nsteps - tau)
        % Loop over chromosome pairs
        for i = 1:(Nchromosomes - 1)
            for j = (i + 1):Nchromosomes
                neighbor_level = j - i;  % 1st, 2nd, ..., N-1
                idx = neighbor_level;    % MATLAB is 1-based
                corr = (CM(t + tau, i) - CM(t, i)) * (CM(t + tau, j) - CM(t, j));
                C(idx) = C(idx) + corr;
                N(idx) = N(idx) + 1;
            end
        end
    end
end
