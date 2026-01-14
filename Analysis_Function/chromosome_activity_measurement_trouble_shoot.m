%% KE / activity from experimental CSV
clc; clear; close all;

% ----------------------- User settings -----------------------
% Open a dialog to pick a CSV
[fn, fp] = uigetfile({'*.csv','CSV files (*.csv)'; '*.*','All files'}, ...
                     'Select the CSV file');
if isequal(fn,0)
    error('No file selected.');
end
csv_file = fullfile(fp, fn);   % full path to the chosen file

%%
% ----------------------- Load & standardize ------------------
opts = detectImportOptions(csv_file);
opts = setvartype(opts,'Pole','string');   % keep 'L'/'R' as strings
TBL  = readtable(csv_file, opts);
TBL.Properties.VariableNames = matlab.lang.makeValidName(TBL.Properties.VariableNames);


% Units: convert time to minutes (so velocities are µm/min)
TBL.t_min = TBL.T / 60;
% 
% % Compute SIGNED, pole-corrected positions (don’t use abs!)
% TBL.x_rel = TBL.X_k_cor - TBL.X_p_cor;
% TBL.y_rel = TBL.Y_k_cor - TBL.Y_p_cor;

% ----------------------- Per-Pair workflow -------------------
results = [];
pairIDs = unique(TBL.PairID);

for idx = 1:numel(pairIDs)
    pid = pairIDs(idx);

    % Split L / R for this pair, sort by time
    L = sortrows(TBL(TBL.PairID==pid & TBL.Pole=="L", :), 't_min');
    R = sortrows(TBL(TBL.PairID==pid & TBL.Pole=="R", :), 't_min');
    
    % ---------- keep only overlapping times ----------
    % 1) exact intersection
    [t_common, iL, iR] = intersect(L.t_min, R.t_min);
    Lc = L(iL, :);
    Rc = R(iR, :);
    % Now Lc and Rc have the same times & same length
    t_sync = Lc.t_min;       % == Rc.t_min
    xL = Lc.X_Pole;
    xR = Rc.X_Pole;

    % Your sliding 3-point slope calculation
    vL_list = [];
    vR_list = [];
    for j = 1:(length(t_sync) - 2)
        t_chunk = t_sync(j:j+2);
        sL_chunk = xL(j:j+2);
        sR_chunk = xR(j:j+2);

        pL = polyfit(t_chunk, sL_chunk, 1);
        pR = polyfit(t_chunk, sR_chunk, 1);

        vL_list(end+1) = pL(1);  % slope = velocity (µm/min)
        vR_list(end+1) = pR(1);
    end

    % Kinetic energy and velocity stats (your definitions)
    KE     = vL_list.^2 + vR_list.^2;
    avg_KE = mean(KE, 'omitnan');
    avg_vL = mean((vL_list), 'omitnan');
    avg_vR = mean((vR_list), 'omitnan');

    % Store per-pair result
    results = [results; struct( ...
        'PairID', pid, ...
        'Avg_KE', avg_KE, ...
        'vL',     avg_vL, ...
        'vR',     avg_vR)];
end

% Results table across all pairs + summary bar
results_table = struct2table(results);
disp(results_table);

figure;
bar(results_table.PairID, results_table.Avg_KE);
xlabel('PairID'); ylabel('Avg KE (v_L^2 + v_R^2)');
title('Average KE per sister pair'); grid on;