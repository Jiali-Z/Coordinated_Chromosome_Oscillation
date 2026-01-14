function pearson_table = pearsons_correlation_analysis(cL, cR, dt,duration)
% PEARSONS_CORRELATION_ANALYSIS
%   Downsample to every 5 seconds, compute CM for each chromosome, and
%   return Pearson's r for all unique chromosome pairs as a table.
%
% INPUTS
%   cL        : [T x 2 x Nchromosomes] left sister positions
%   cR        : [T x 2 x Nchromosomes] right sister positions
%   dt        : timestep in minutes
%   duration  : duration (minutes) to analyze (cap if longer than data)
%
% OUTPUT
%   pearson_table : table with columns
%                   - pair       : [i j] chromosome index pairs (i<j)
%                   - pearson_r  : Pearson correlation coefficient in [-1, 1]

    [Nsteps, ~, Nchromosomes] = size(cL);
    sampling_interval=round((5/60)/dt);  % every 5second, ~21
    total_steps = min(round(duration / dt), Nsteps); % limit to Nsteps
    CMx_sampled=[];
    for idx = 1:Nchromosomes
        signalL = cL(1:sampling_interval:total_steps, 1, idx); % x only
        signalR = cR(1:sampling_interval:total_steps, 1, idx); 
        CMx_sampled(:, idx)=(signalR+signalL)/2;
    end
    % Build the results table 
    nPairs = Nchromosomes*(Nchromosomes-1)/2;
    pearsons = zeros(nPairs, 1);            
    pairs = zeros(nPairs, 2);
    k = 0;
    for i = 1:Nchromosomes-1
        for j = i+1:Nchromosomes
            k = k + 1;
            pairs(k, :) = [i, j];
            R = corrcoef(CMx_sampled(:, i), CMx_sampled(:, j));
            r = R(1, 2);
            pearsons(k)=r;
        end
    end
    pearson_table = table(pairs, pearsons, ...
        'VariableNames', {'pair', 'pearson_r'});
end
