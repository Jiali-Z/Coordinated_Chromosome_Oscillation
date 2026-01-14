clc; clear; 
%close all;

% --------------------- USER SETTINGS ---------------------
koff0_list  = [40];
kon0_list   = [20 ];

avg_win_sec = 5;
ylim_max    = 12;

line_width = 2;
marker_every = 400;

% ---------- Fixed model parameters ----------
dt = 2e-3;  Nsteps = 5000;  Nchromosomes = 5;
noise = 0.05; noise_b = 0.05; Kct = 12.30; I0 = 2; Kkt = 1;
Gamma = 0.1; Beta = 0.7; Nmax = 25; epsilon = 0.1;

% Inter-chromosomal spring
l0 = 3;    % µm
k  = 1;    % coupling strength

tvec = (1:Nsteps)*dt;
win_steps = round(avg_win_sec / (dt*60));

figure('Color','w'); hold on; grid on;
xlabel('Time (min)','FontSize',18,'FontName','Arial');
ylabel(sprintf('Total Connections (%ds average)', avg_win_sec),'FontSize',18,'FontName','Arial');
title('Time-Averaged Chromosome Connectivity vs. Time','FontSize',20);
set(gca,'FontSize',16,'FontName','Arial');
ylim([0, ylim_max]);

colors = lines(numel(kon0_list)*numel(koff0_list));

% === Accumulator for tidy CSV (initialize ONCE, before loops) ===
all_rows = table();

% --------- Sweep across (kon_0, koff_0) combinations ---------
for iKon = 1:numel(kon0_list)
    kon_0 = kon0_list(iKon);
    colors = lines(numel(koff0_list));

    for iKoff = 1:numel(koff0_list)
        koff_0 = koff0_list(iKoff);

        % ====== Initialize per-run state ======
        cL = zeros(Nsteps+1, 2, Nchromosomes);
        cR = zeros(Nsteps+1, 2, Nchromosomes);
        xL = zeros(Nsteps+1, Nchromosomes);
        xR = zeros(Nsteps+1, Nchromosomes);
        NL = zeros(Nsteps+1, Nchromosomes);
        NR = zeros(Nsteps+1, Nchromosomes);
        vL = zeros(Nsteps,   Nchromosomes);
        vR = zeros(Nsteps,   Nchromosomes);

        Nbar   = zeros(Nchromosomes,1);
        n_dot  = zeros(Nchromosomes,1);
        Lambda = zeros(Nchromosomes,1);
        Alpha  = zeros(Nchromosomes,1);

        for ii = 1:Nchromosomes
            Nbar(ii)   = 20 * (1 + noise*(2*rand-1));
            n_dot(ii)  =  1 * (1 + noise*(2*rand-1));
            Lambda(ii) = n_dot(ii)/Nbar(ii);
            Alpha(ii)  = n_dot(ii) * 6.2 / (1 - Beta);
        end

        y_positions = linspace(-Nchromosomes/2, Nchromosomes/2, Nchromosomes);
        flag_c  = zeros(Nchromosomes,Nchromosomes);
        flag_cp = flag_c;
        flag_c_sum_accumulated = zeros(Nsteps,1);

        for ii = 1:Nchromosomes
            NR(1,ii) = Nbar(ii)*(1 - epsilon);
            NL(1,ii) = Nbar(ii)*(1 + epsilon);
            cL(1,1,ii) = -I0/2 - epsilon*(2*rand-1);
            cR(1,1,ii) =  I0/2 + epsilon*(2*rand-1);
            cL(:,2,ii) =  y_positions(ii);
            cR(:,2,ii) =  y_positions(ii);
            xL(1,ii) = cL(1,1,ii) - 0.3;
            xR(1,ii) = cR(1,1,ii) + 0.3;
            vL(1,ii) = 0; vR(1,ii) = 0;
        end

        % ====== ODE loop ======
        for t = 1:Nsteps
            CM = squeeze((cL(t,1,:) + cR(t,1,:))/2);
            YM = squeeze((cL(t,2,:) + cR(t,2,:))/2);

            [flag_c, ~] = spring_connect(CM, YM, flag_cp, Nchromosomes, dt, l0, koff_0, kon_0);
            flag_c_sum_accumulated(t) = sum(triu(flag_c, 1), 'all');
            flag_cp = flag_c;

            for ii = 1:Nchromosomes
                F_coupling_x = 0;
                for jj = 1:Nchromosomes
                    if jj ~= ii
                        delta_x = CM(ii) - CM(jj);
                        F_coupling_x = F_coupling_x - k * delta_x * flag_c(ii, jj);
                    end
                end

                F_KT_L = -NL(t,ii) * Kkt * (cL(t,1,ii) - xL(t,ii));
                F_CT_L =  Kct      * (cR(t,1,ii) - cL(t,1,ii) - I0);
                F_KT_R = -NR(t,ii) * Kkt * (cR(t,1,ii) - xR(t,ii));
                F_CT_R = -Kct      * (cR(t,1,ii) - cL(t,1,ii) - I0);

                vL(t,ii) = (F_KT_L + F_CT_L + F_coupling_x) / Gamma;
                vR(t,ii) = (F_KT_R + F_CT_R + F_coupling_x) / Gamma;

                dx_diff = noise_b * (randn(1)) * sqrt(dt);
                cL(t+1,1,ii) = cL(t,1,ii) + dt*vL(t,ii) + dx_diff;
                cR(t+1,1,ii) = cR(t,1,ii) + dt*vR(t,ii) + dx_diff;

                cL(t+1,2,ii) = cL(1,2,ii);
                cR(t+1,2,ii) = cR(1,2,ii);

                xL(t+1,ii) = xL(t,ii) + dt*Beta*vL(t,ii);
                xR(t+1,ii) = xR(t,ii) + dt*Beta*vR(t,ii);

                NL(t+1,ii) = NL(t,ii) + dt*( n_dot(ii) ...
                    + (-Alpha(ii) * NL(t,ii) * (1 - Beta) * vL(t,ii) * (1 - NL(t,ii)/Nmax)) ...
                    - Lambda(ii)*NL(t,ii) );
                NR(t+1,ii) = NR(t,ii) + dt*( n_dot(ii) ...
                    + ( +Alpha(ii) * NR(t,ii) * (1 - Beta) * vR(t,ii) * (1 - NR(t,ii)/Nmax)) ...
                    - Lambda(ii)*NR(t,ii) );
            end
        end

        % ====== Time-averaged connectivity and plot ======
        flag_c_sliding_avg = movmean(flag_c_sum_accumulated, win_steps);

        lbl = sprintf('koff_0=%g', koff_0);
        plot(tvec, flag_c_sliding_avg, 'LineWidth', line_width, 'Color', colors(iKoff,:),'DisplayName',lbl);

        if marker_every > 0
            mk_idx = 1:marker_every:numel(tvec);
            plot(tvec(mk_idx), flag_c_sliding_avg(mk_idx), 'o', 'MarkerSize', 3, ...
                 'Color', colors(iKoff,:), 'HandleVisibility','off');
        end

        % === append this run's series to the tidy table (INSIDE inner loop) ===
        T = table( ...
            repmat(kon_0, numel(tvec), 1), ...
            repmat(koff_0, numel(tvec), 1), ...
            tvec(:), ...
            flag_c_sliding_avg(:), ...
            'VariableNames', {'kon_0','koff_0','time_min','avg_conn_win'});
        all_rows = [all_rows; T];   %#ok<AGROW>

    end

    legend('Location','bestoutside');
    ylim([0, ylim_max]);
end

% --- write once, AFTER the loops ---
writetable(all_rows, 'connectivity_timeseries.csv');
fprintf('Wrote connectivity_timeseries.csv with %d rows\n', height(all_rows));
