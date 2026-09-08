%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  RPE uniquement : G1 vs G2 (Borg CR-10)                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% PATHS
PathCycle = 'J:\Piano_Fatigue\Data_Exported\';
PathInfo  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
PathSave  = fileparts(mfilename('fullpath'));

N_bins   = 10;
T        = (1:N_bins) * (100/N_bins);

G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];
G2_Liszt = [3,8,10,11,12,13,14,16,18,19,24,25,26,27,33,34,36,40,41,43,47,48,50];

cycles_info = load(fullfile(PathCycle, 'cycles_Li.mat'), 'participants');
tmp_info    = load(fullfile(PathInfo, 'info_participants_corrected.mat'), 'Info_participants');

%% CHARGEMENT RPE
function [Mat_RPE, valid] = loadRPE(GroupList, cycles_info, tmp_info, N_bins)
    n       = length(GroupList);
    Mat_RPE = NaN(n, N_bins);
    valid   = false(n, 1);

    for iG = 1:n
        iP     = GroupList(iG);
        t_cyc  = cycles_info.participants(iP).t_do(:);
        if isempty(t_cyc), continue; end

        % Nombre de timepoints
        tp_s = max(1, floor(t_cyc(1)));
        tp_e = floor(t_cyc(end));
        if tp_e <= tp_s, continue; end
        N_tp = tp_e - tp_s + 1;

        % RPE interpolée
        rpe_raw = tmp_info.Info_participants(iP).Liszt(:);
        rpe_raw(isnan(rpe_raw)) = [];
        if iP == 22, rpe_raw = rpe_raw(1:end-1); end
        t_rpe = (1:length(rpe_raw))' * 30;
        t_tp  = (0:N_tp-1)';
        t_rpe = min(t_rpe, t_tp(end));
        rpe_tp = interp1(t_rpe, double(rpe_raw), t_tp, 'linear', 'extrap');
        rpe_tp = max(0, min(10, rpe_tp));

        % Bins
        be = zeros(1, N_bins+1); be(1) = 1;
        for i = 1:N_bins, be(i+1) = round(N_tp * i / N_bins); end
        for k = 1:N_bins
            Mat_RPE(iG, k) = mean(rpe_tp(be(k):be(k+1)));
        end
        valid(iG) = true;
    end
end

[MRPE_G1, val_G1] = loadRPE(G1_Liszt, cycles_info, tmp_info, N_bins);
[MRPE_G2, val_G2] = loadRPE(G2_Liszt, cycles_info, tmp_info, N_bins);

Vi1 = find(val_G1);
Vi2 = find(val_G2);

%% STATS
col_G1 = [0.85 0.15 0.15];
col_G2 = [0.15 0.45 0.85];

Mu_rpe1 = nanmean(MRPE_G1(Vi1,:), 1);
SD_rpe1 = nanstd(MRPE_G1(Vi1,:), 0, 1);
Mu_rpe2 = nanmean(MRPE_G2(Vi2,:), 1);
SD_rpe2 = nanstd(MRPE_G2(Vi2,:), 0, 1);

%% FIGURE
figure('Name','RPE_G1vsG2','NumberTitle','off','Position',[100 100 700 420]);
hold on;

fill([T fliplr(T)], [Mu_rpe1+SD_rpe1 fliplr(Mu_rpe1-SD_rpe1)], col_G1, ...
    'FaceAlpha', 0.20, 'EdgeColor', 'none');
fill([T fliplr(T)], [Mu_rpe2+SD_rpe2 fliplr(Mu_rpe2-SD_rpe2)], col_G2, ...
    'FaceAlpha', 0.20, 'EdgeColor', 'none');

h1 = plot(T, Mu_rpe1, '-o', 'Color', col_G1, 'LineWidth', 2.5, 'MarkerSize', 5);
h2 = plot(T, Mu_rpe2, '-o', 'Color', col_G2, 'LineWidth', 2.5, 'MarkerSize', 5);

legend([h1 h2], {['G1 (n=' num2str(length(Vi1)) ')'], ...
                  ['G2 (n=' num2str(length(Vi2)) ')']}, ...
    'Location', 'northwest', 'FontSize', 10);

xlabel('Temps normalisé (%)', 'FontSize', 11);
ylabel('RPE (Borg CR-10)', 'FontSize', 11);
title('Évolution du RPE - G1 vs G2', 'FontSize', 13, 'FontWeight', 'bold');
xlim([5 105]); ylim([0 10]);
xticks(T); yticks(0:2:10);
grid on; box on; hold off;

saveas(gcf, fullfile(PathSave, 'Fig_RPE_G1vsG2.png'));
fprintf('Figure sauvegardée.\n');