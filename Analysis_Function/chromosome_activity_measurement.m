function results_table=chromosome_activity_measurement(cL,cR,dt,duration)
% Chromosome_activity_measurement measures instaneous velovities for each chromosome
% at the experimental sampling rate(every 5seconds) and their kinetic
% activity
%Inout: 
%   cL,cR: Nsteps*2*Nchromosomes arrays of left and right chromosome pos
%   dt: time steps(mins)
%   duration: duration of the trajectoy to be analyzed (default=10min)
%Output:
%   results_taboe: table with Avg_KE,cL,cR for each chromosome


% Set duration default to 10 if not provided 
if nargin < 4
    duration = 10; % default 10 minutes
end

results = [];
[Nsteps, ~, Nchromosomes] = size(cL);
sampling_interval=round((5/60)/dt);  % every 5second, ~21
total_steps = min(round(duration / dt), Nsteps); % limit to Nsteps

for idx = 1:Nchromosomes
    signalL = squeeze(cL(1:sampling_interval:total_steps, 1, idx)); % x only
    signalR = squeeze(cR(1:sampling_interval:total_steps, 1, idx)); 
    time_sampled_min = (0:length(signalL)-1)' * sampling_interval * dt;  % minutes
    time_sampled_sec = time_sampled_min * 60;  % seconds 
    vL_list = [];
    vR_list = [];
    % Sliding window 3-point slope calculation
    for j = 1:(length(signalL) - 2)
        t_chunk = time_sampled_min(j:j+2);
        sL_chunk = signalL(j:j+2);
        sR_chunk = signalR(j:j+2);

        % Fit line: first-order polynomial (slope only)
        pL = polyfit(t_chunk, sL_chunk, 1);
        pR = polyfit(t_chunk, sR_chunk, 1);

        vL_list(end+1) = pL(1);  % slope = velocity
        vR_list(end+1) = pR(1);
    end

    % Kinetic energy and velocity stats
    avg_KE = mean(vL_list.^2 + vR_list.^2,'omitnan');
    
    avg_vL = mean((vL_list), 'omitnan');
    avg_vR = mean((vR_list), 'omitnan');

    % Store result
    results = [results; struct( ...
        'PairID', idx, ...
        'Avg_KE', avg_KE, ...
        'vL', avg_vL, ...
        'vR', avg_vR)];
end
results_table = struct2table(results);
end




