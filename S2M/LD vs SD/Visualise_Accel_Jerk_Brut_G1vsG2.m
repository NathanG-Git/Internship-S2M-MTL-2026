%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Acceleration & Jerk BRUTS par segment - G1 vs G2                %%%%
%%%%  Source : Features3_XSENS_S1XX_Li.mat                             %%%%
%%%%  Variables : Mean.Module_Acceleration, Mean.Module_Jerk           %%%%
%%%%  G1 : 26 participants | G2 : 23 participants                      %%%%
%%%%  Valeurs BRUTES sans normalisation (m/s^2 et m/s^3)               %%%%
%%%%  1 figure par signal : 2 lignes (Mean/Median) x 7 cols (segments) %%%%
%%%%  G1 = rouge | G2 = bleu | RPE G1 = tirets rouge | RPE G2 = tirets bleu %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS & PARAMETRES
%  -----------------------------------------------------------------------
PathFeat  = 'J:\Piano_Fatigue\Data_Exported\Features_XSENS\Li\';
PathCycle = 'J:\Piano_Fatigue\Data_Exported\';
PathInfo  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
PathSave  = fileparts(mfilename('fullpath'));

N_bins   = 10;
T        = (1:N_bins) * (100/N_bins);

G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];
G2_Liszt = [3,8,10,11,12,13,14,16,18,19,24,25,26,27,33,34,36,40,41,43,47,48,50];

Segments = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
N_seg    = length(Segments);

cycles_info = load(fullfile(PathCycle, 'cycles_Li.mat'), 'participants');
tmp_info    = load(fullfile(PathInfo, 'info_participants_corrected.mat'), 'Info_participants');

%% -----------------------------------------------------------------------
%  FONCTION DE CHARGEMENT (commune G1 et G2)
%  -----------------------------------------------------------------------
function [Mat_Accel, Mat_Jerk, Mat_RPE, valid] = loadGroup(GroupList, PathFeat, cycles_info, tmp_info, N_bins, N_seg)
    n = length(GroupList);
    Mat_Accel = NaN(n, N_bins, N_seg);
    Mat_Jerk  = NaN(n, N_bins, N_seg);
    Mat_RPE   = NaN(n, N_bins);
    valid     = false(n, 1);

    for iG = 1:n
        iP     = GroupList(iG);
        SubjID = sprintf('S1%02d', iP);
        fname  = fullfile(PathFeat, sprintf('Features3_XSENS_%s_Li.mat', SubjID));
        if ~exist(fname, 'file'), fprintf('  MISSING: %s\n', fname); continue; end

        tmp  = load(fname, 'Features_XSENS');
        feat = tmp.Features_XSENS;
        N_f  = size(feat.Mean.Module_Acceleration, 1);

        t_cyc = cycles_info.participants(iP).t_do(:);
        if isempty(t_cyc), fprintf('  %s: t_do vide\n', SubjID); continue; end

        tp_s = max(1, floor(t_cyc(1)));
        tp_e = min(N_f, floor(t_cyc(end)));
        if tp_e <= tp_s, continue; end
        N_tp = tp_e - tp_s + 1;

        % RPE
        rpe_raw = tmp_info.Info_participants(iP).Liszt(:);
        rpe_raw(isnan(rpe_raw)) = [];
        if iP == 22, rpe_raw = rpe_raw(1:end-1); end
        t_rpe = (1:length(rpe_raw))' * 30;
        t_tp  = (0:N_tp-1)';
        t_rpe = min(t_rpe, t_tp(end));
        rpe_tp = interp1(t_rpe, double(rpe_raw), t_tp, 'linear', 'extrap');
        rpe_tp = max(0, min(10, rpe_tp));

        % Bins
        be    = zeros(1, N_bins+1); be(1) = 1;
        for i = 1:N_bins, be(i+1) = round(N_tp * i / N_bins); end

        mod_a = feat.Mean.Module_Acceleration;
        mod_j = feat.Mean.Module_Jerk;

        for k = 1:N_bins
            ks = be(k); ke = be(k+1);
            if ke < ks, continue; end
            Mat_RPE(iG, k) = mean(rpe_tp(ks:ke));
            gs = tp_s+ks-1; ge = min(tp_s+ke-1, N_f);
            for iSeg = 1:N_seg
                Mat_Accel(iG, k, iSeg) = mean(mod_a(gs:ge, iSeg));
                Mat_Jerk(iG,  k, iSeg) = mean(mod_j(gs:ge, iSeg));
            end
        end
        valid(iG) = true;
        fprintf('  %s: OK\n', SubjID);
    end
end

%% -----------------------------------------------------------------------
%  CHARGEMENT G1 et G2
%  -----------------------------------------------------------------------
fprintf('=== Chargement G1 (%d participants) ===\n', length(G1_Liszt));
[MA_G1, MJ_G1, MRPE_G1, val_G1] = loadGroup(G1_Liszt, PathFeat, cycles_info, tmp_info, N_bins, N_seg);

fprintf('\n=== Chargement G2 (%d participants) ===\n', length(G2_Liszt));
[MA_G2, MJ_G2, MRPE_G2, val_G2] = loadGroup(G2_Liszt, PathFeat, cycles_info, tmp_info, N_bins, N_seg);

Vi1 = find(val_G1);
Vi2 = find(val_G2);
fprintf('\nG1 valides: %d | G2 valides: %d\n', length(Vi1), length(Vi2));

[Mu_rpe1, SD_rpe1] = gStats(MRPE_G1(Vi1,:));
[Mu_rpe2, SD_rpe2] = gStats(MRPE_G2(Vi2,:));

%% -----------------------------------------------------------------------
%  FIGURES : 1 par signal (Accel / Jerk)
%  Chaque figure : 2 lignes (Mean, Median) x 7 colonnes (segments)
%  Sur chaque subplot : G1 rouge + G2 bleu + RPE G1 tirets rouge + RPE G2 tirets bleu
%  -----------------------------------------------------------------------
col_G1     = [0.85 0.15 0.15];  % rouge G1
col_G2     = [0.15 0.45 0.85];  % bleu  G2
col_rpe_G1 = [0.70 0.10 0.10];  % rouge fonce RPE G1
col_rpe_G2 = [0.10 0.25 0.70];  % bleu fonce  RPE G2

sig_names  = {'Acceleration', 'Jerk'};
lbl_units  = {'Acceleration (m/s^2)', 'Jerk (m/s^3)'};
MatG1_all  = {MA_G1, MJ_G1};
MatG2_all  = {MA_G2, MJ_G2};

for iSig = 1:2
    MatG1 = MatG1_all{iSig};
    MatG2 = MatG2_all{iSig};
    lbl_y = lbl_units{iSig};

    fig_name = ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_Brut_G1vsG2'];
    figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [30 30 1800 700]);

    for iScore = 1:2  % 1=Mean, 2=Median
        for iSeg = 1:N_seg
            subplot(2, N_seg, (iScore-1)*N_seg + iSeg);
            hold on;

            Mv1 = squeeze(MatG1(Vi1, :, iSeg));
            Mv2 = squeeze(MatG2(Vi2, :, iSeg));

            if iScore == 1
                [Mu1, SD1] = gStats(Mv1);
                [Mu2, SD2] = gStats(Mv2);
                score_str  = 'Mean';
            else
                [Mu1, SD1] = gStatsMedian(Mv1);
                [Mu2, SD2] = gStatsMedian(Mv2);
                score_str  = 'Median';
            end

            yyaxis left

            % Bande SD G1
            fill([T fliplr(T)], [Mu1+SD1 fliplr(Mu1-SD1)], col_G1, ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none');
            % Bande SD G2
            fill([T fliplr(T)], [Mu2+SD2 fliplr(Mu2-SD2)], col_G2, ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none');
            % Courbes
            h1 = plot(T, Mu1, '-',  'Color', col_G1, 'LineWidth', 2.5);
            h2 = plot(T, Mu2, '-',  'Color', col_G2, 'LineWidth', 2.5);

            ylabel(lbl_y, 'FontSize', 8);
            ax = gca; ax.YColor = [0.2 0.2 0.2];

            yyaxis right

            % RPE G1 et G2 superposes
            fill([T fliplr(T)], [Mu_rpe1+SD_rpe1 fliplr(Mu_rpe1-SD_rpe1)], col_rpe_G1, ...
                'FaceAlpha', 0.08, 'EdgeColor', 'none');
            fill([T fliplr(T)], [Mu_rpe2+SD_rpe2 fliplr(Mu_rpe2-SD_rpe2)], col_rpe_G2, ...
                'FaceAlpha', 0.08, 'EdgeColor', 'none');
            plot(T, Mu_rpe1, ':', 'Color', col_rpe_G1, 'LineWidth', 1.5);
            plot(T, Mu_rpe2, ':', 'Color', col_rpe_G2, 'LineWidth', 1.5);
            ylabel('RPE (Borg CR-10)', 'FontSize', 8, 'Color', [0.4 0.1 0.6]);
            ax.YColor = [0.4 0.1 0.6];
            ylim([0 10]);

            % Legende (1er subplot seulement)
            if iSeg == 1 && iScore == 1
                legend([h1 h2], ...
                    {['G1 n=' num2str(length(Vi1))], ['G2 n=' num2str(length(Vi2))]}, ...
                    'Location', 'northwest', 'FontSize', 7);
            end

            xlabel('Temps normalise (%)', 'FontSize', 8);
            title([score_str ' - ' Segments{iSeg}], 'FontSize', 9, ...
                'FontWeight', 'bold', 'Interpreter', 'none');
            xlim([5 105]); xticks(T); grid on; box on; hold off;
        end
    end

    sgtitle([sig_names{iSig} ' BRUT par segment - G1 (rouge) vs G2 (bleu) - Temps normalise (%)'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, [PathSave '\' fig_name '.png']);
    fprintf('Figure : %s.png\n', fig_name);
end

fprintf('\nTermine.\n');


%% -----------------------------------------------------------------------
%  FIGURES DELTA : evolution relative par rapport au bin 1
%  Pour chaque participant : delta(k) = signal(k) - signal(bin1)
%  Moyenne de groupe sur ces deltas -> tout le monde part de 0
%  Pas d'ecart-type RPE affiche
%  -----------------------------------------------------------------------

sig_names_d = {'Acceleration', 'Jerk'};
lbl_units_d = {'Delta Acceleration (m/s^2)', 'Delta Jerk (m/s^3)'};
MatG1_all_d = {MA_G1, MJ_G1};
MatG2_all_d = {MA_G2, MJ_G2};

% Centrage RPE individuel (bin1 soustrait)
RPE1_delta = MRPE_G1(Vi1,:) - MRPE_G1(Vi1,1);
RPE2_delta = MRPE_G2(Vi2,:) - MRPE_G2(Vi2,1);
[Mu_rpe1_d, ~] = gStats(RPE1_delta);
[Mu_rpe2_d, ~] = gStats(RPE2_delta);

for iSig = 1:2
    MatG1_raw = MatG1_all_d{iSig};
    MatG2_raw = MatG2_all_d{iSig};
    lbl_y     = lbl_units_d{iSig};

    % Centrage individuel : signal(k) - signal(bin1)
    MatG1_d = MatG1_raw(Vi1,:,:) - MatG1_raw(Vi1,1,:);
    MatG2_d = MatG2_raw(Vi2,:,:) - MatG2_raw(Vi2,1,:);

    fig_name = ['Fig_TimeNorm_BySegment_' sig_names_d{iSig} '_Delta_G1vsG2'];
    figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [30 30 1800 700]);

    for iScore = 1:2
        for iSeg = 1:N_seg
            subplot(2, N_seg, (iScore-1)*N_seg + iSeg);
            hold on;

            Mv1 = squeeze(MatG1_d(:, :, iSeg));
            Mv2 = squeeze(MatG2_d(:, :, iSeg));

            if iScore == 1
                [Mu1, SD1] = gStats(Mv1);
                [Mu2, SD2] = gStats(Mv2);
                score_str  = 'Mean';
            else
                [Mu1, SD1] = gStatsMedian(Mv1);
                [Mu2, SD2] = gStatsMedian(Mv2);
                score_str  = 'Median';
            end

            yyaxis left
            fill([T fliplr(T)], [Mu1+SD1 fliplr(Mu1-SD1)], col_G1, 'FaceAlpha',0.15,'EdgeColor','none');
            fill([T fliplr(T)], [Mu2+SD2 fliplr(Mu2-SD2)], col_G2, 'FaceAlpha',0.15,'EdgeColor','none');
            h1 = plot(T, Mu1, '-', 'Color', col_G1, 'LineWidth', 2.5);
            h2 = plot(T, Mu2, '-', 'Color', col_G2, 'LineWidth', 2.5);
            plot(T, zeros(1,length(T)), '--k', 'LineWidth', 0.8);
            ylabel(lbl_y, 'FontSize', 8);
            ax = gca; ax.YColor = [0.2 0.2 0.2];

            yyaxis right
            plot(T, Mu_rpe1_d, ':', 'Color', col_rpe_G1, 'LineWidth', 1.5);
            plot(T, Mu_rpe2_d, ':', 'Color', col_rpe_G2, 'LineWidth', 1.5);
            ylabel('Delta RPE (Borg CR-10)', 'FontSize', 8, 'Color', [0.4 0.1 0.6]);
            ax.YColor = [0.4 0.1 0.6];

            if iSeg == 1 && iScore == 1
                legend([h1 h2], {['G1 n=' num2str(length(Vi1))], ...
                    ['G2 n=' num2str(length(Vi2))]}, ...
                    'Location','northwest','FontSize',7);
            end

            xlabel('Temps normalise (%)', 'FontSize', 8);
            title([score_str ' - ' Segments{iSeg}], 'FontSize', 9, ...
                'FontWeight','bold','Interpreter','none');
            xlim([5 105]); xticks(T); grid on; box on; hold off;
        end
    end

    sgtitle(['Delta ' sig_names_d{iSig} ' par segment - G1 (rouge) vs G2 (bleu) - Evolution depuis bin1'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, [PathSave '\\' fig_name '.png']);
    fprintf('Figure : %s.png\n', fig_name);
end

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------

function [Mu, SD] = gStats(Mat)
    N  = sum(~isnan(Mat), 1);
    Mu = nanmean(Mat, 1);
    SD = nanstd(Mat, 0, 1);
    Mu(N < 3) = NaN;
    SD(N < 3) = NaN;
end

function [Md, SD] = gStatsMedian(Mat)
    N  = sum(~isnan(Mat), 1);
    Md = nanmedian(Mat, 1);
    SD = nanstd(Mat, 0, 1);
    Md(N < 3) = NaN;
    SD(N < 3) = NaN;
end