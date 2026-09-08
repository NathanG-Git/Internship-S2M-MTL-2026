%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Variables Goubault BRUTES - G1 vs G2 superposes                 %%%%
%%%%  Source : Features3_XSENS_S1XX_Li.mat                             %%%%
%%%%  10 variables Goubault 2023 (Chord task, Table III)               %%%%
%%%%  Valeurs BRUTES sans normalisation                                 %%%%
%%%%  1 figure : 2 lignes x 5 colonnes (10 variables)                  %%%%
%%%%  G1 = rouge | G2 = bleu | RPE G1/G2 = pointilles                 %%%%
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

%% -----------------------------------------------------------------------
%  DEFINITION DES 10 VARIABLES GOUBAULT
%  {nom_affichage, feature_field, signal_field, col_1based, unite}
%  -----------------------------------------------------------------------
Vars = {
    'Hand MedianFreq Accel Y',           'MedianFreq',             'Acceleration',           20, 'Hz';
    'Hand Prc90 Accel X',                'Percentile90',           'Acceleration',           19, 'm/s^2';
    'Hand PeakPower Accel Mod',          'PeakPower',              'Module_Acceleration',     7, 'W';
    'Hand SpectralEntropy Accel Y',      'SpectralEntropy',        'Acceleration',           20, 'u.a.';
    'Head PeakPower AngVel X',           'PeakPower',              'Angular_Velocity',        7, 'W';
    'Forearm Mean Accel Y',              'Mean',                   'Acceleration',           17, 'm/s^2';
    'Forearm SpectralEntropy AngVel Mod','SpectralEntropy',        'Module_Angular_Velocity', 6, 'u.a.';
    'Forearm PeakPowerFreq AngVel Mod',  'PeakPower_Freq',         'Module_Angular_Velocity', 6, 'Hz';
    'Shoulder Mean AngVel X',            'Mean',                   'Angular_Velocity',       10, 'rad/s';
    'Shoulder Mean AngVel Y',            'Mean',                   'Angular_Velocity',       11, 'rad/s';
};
N_vars = size(Vars, 1);

cycles_info = load(fullfile(PathCycle, 'cycles_Li.mat'), 'participants');
tmp_info    = load(fullfile(PathInfo, 'info_participants_corrected.mat'), 'Info_participants');

%% -----------------------------------------------------------------------
%  FONCTION DE CHARGEMENT
%  -----------------------------------------------------------------------
function [Mat_Vars, Mat_RPE, valid] = loadGroup(GroupList, PathFeat, cycles_info, tmp_info, N_bins, Vars)
    N_vars = size(Vars, 1);
    n      = length(GroupList);
    Mat_Vars = NaN(n, N_bins, N_vars);
    Mat_RPE  = NaN(n, N_bins);
    valid    = false(n, 1);

    for iG = 1:n
        iP     = GroupList(iG);
        SubjID = sprintf('S1%02d', iP);
        fname  = fullfile(PathFeat, sprintf('Features3_XSENS_%s_Li.mat', SubjID));
        if ~exist(fname, 'file'), fprintf('  MISSING: %s\n', fname); continue; end

        tmp  = load(fname, 'Features_XSENS');
        feat = tmp.Features_XSENS;
        N_f  = size(feat.Mean.Acceleration, 1);

        t_cyc = cycles_info.participants(iP).t_do(:);
        if isempty(t_cyc), continue; end
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
        be = zeros(1, N_bins+1); be(1) = 1;
        for i = 1:N_bins, be(i+1) = round(N_tp * i / N_bins); end

        for k = 1:N_bins
            ks = be(k); ke = be(k+1);
            if ke < ks, continue; end
            Mat_RPE(iG, k) = mean(rpe_tp(ks:ke));
            gs = tp_s+ks-1; ge = min(tp_s+ke-1, N_f);

            for iV = 1:N_vars
                feat_name = Vars{iV, 2};
                sig_name  = Vars{iV, 3};
                col       = Vars{iV, 4};
                try
                    blk = feat.(feat_name).(sig_name)(gs:ge, col);
                    if strcmp(feat_name, 'Mean')
                        blk = abs(blk);
                    end
                    Mat_Vars(iG, k, iV) = mean(blk);
                catch
                end
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
[MV_G1, MRPE_G1, val_G1] = loadGroup(G1_Liszt, PathFeat, cycles_info, tmp_info, N_bins, Vars);

fprintf('\n=== Chargement G2 (%d participants) ===\n', length(G2_Liszt));
[MV_G2, MRPE_G2, val_G2] = loadGroup(G2_Liszt, PathFeat, cycles_info, tmp_info, N_bins, Vars);

Vi1 = find(val_G1);
Vi2 = find(val_G2);
fprintf('\nG1 valides: %d | G2 valides: %d\n', length(Vi1), length(Vi2));

[Mu_rpe1, SD_rpe1] = gStats(MRPE_G1(Vi1,:));
[Mu_rpe2, SD_rpe2] = gStats(MRPE_G2(Vi2,:));

%% -----------------------------------------------------------------------
%  COULEURS
%  -----------------------------------------------------------------------
col_G1     = [0.85 0.15 0.15];   % rouge G1
col_G2     = [0.15 0.45 0.85];   % bleu  G2
col_rpe_G1 = [0.70 0.10 0.10];
col_rpe_G2 = [0.10 0.25 0.70];

%% -----------------------------------------------------------------------
%  FIGURE : 2 lignes x 5 colonnes (Mean ligne 1, Median ligne 2)
%  -----------------------------------------------------------------------
for iScore = 1:2  % fenetre 1 = Mean, fenetre 2 = Median

    if iScore == 1
        fig_name  = 'Fig_Goubault_Brut_G1vsG2_Mean';
        score_str = 'Mean';
    else
        fig_name  = 'Fig_Goubault_Brut_G1vsG2_Median';
        score_str = 'Median';
    end

    figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [30 30 1900 700]);

    for iV = 1:N_vars
        subplot(2, 5, iV);
        hold on;

        Mv1 = squeeze(MV_G1(Vi1, :, iV));
        Mv2 = squeeze(MV_G2(Vi2, :, iV));

        if iScore == 1
            [Mu1, SD1] = gStats(Mv1);
            [Mu2, SD2] = gStats(Mv2);
        else
            [Mu1, SD1] = gStatsMedian(Mv1);
            [Mu2, SD2] = gStatsMedian(Mv2);
        end

        yyaxis left
        fill([T fliplr(T)], [Mu1+SD1 fliplr(Mu1-SD1)], col_G1, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        fill([T fliplr(T)], [Mu2+SD2 fliplr(Mu2-SD2)], col_G2, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        h1 = plot(T, Mu1, '-', 'Color', col_G1, 'LineWidth', 2.5);
        h2 = plot(T, Mu2, '-', 'Color', col_G2, 'LineWidth', 2.5);
        ylabel(Vars{iV, 5}, 'FontSize', 7);
        ax = gca; ax.YColor = [0.2 0.2 0.2];

        yyaxis right
        fill([T fliplr(T)], [Mu_rpe1+SD_rpe1 fliplr(Mu_rpe1-SD_rpe1)], col_rpe_G1, 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        fill([T fliplr(T)], [Mu_rpe2+SD_rpe2 fliplr(Mu_rpe2-SD_rpe2)], col_rpe_G2, 'FaceAlpha', 0.07, 'EdgeColor', 'none');
        plot(T, Mu_rpe1, ':', 'Color', col_rpe_G1, 'LineWidth', 1.5);
        plot(T, Mu_rpe2, ':', 'Color', col_rpe_G2, 'LineWidth', 1.5);
        ylabel('RPE (Borg CR-10)', 'FontSize', 7, 'Color', [0.4 0.1 0.6]);
        ax.YColor = [0.4 0.1 0.6];
        ylim([0 10]);

        if iV == 1
            legend([h1 h2], {['G1 n=' num2str(length(Vi1))], ['G2 n=' num2str(length(Vi2))]}, ...
                'Location', 'northwest', 'FontSize', 7);
        end

        xlabel('Temps normalise (%)', 'FontSize', 7);
        title(Vars{iV,1}, 'FontSize', 7, 'FontWeight', 'bold', 'Interpreter', 'none');
        xlim([5 105]); xticks(T); grid on; box on; hold off;
    end

    sgtitle(['Variables Goubault BRUTES - ' score_str ' - G1 (rouge) vs G2 (bleu) - Temps normalise (%)'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, [PathSave '\' fig_name '.png']);
    fprintf('Figure : %s.png\n', fig_name);
end

fprintf('Termine.\n');


%% -----------------------------------------------------------------------
%  FIGURES DELTA : evolution relative par rapport au bin 1
%  Pour chaque participant : delta(k) = signal(k) - signal(bin1)
%  Moyenne de groupe sur ces deltas -> tout le monde part de 0
%  2 fenetres : Mean et Median | Pas d'ecart-type RPE affiche
%  -----------------------------------------------------------------------

% Centrage RPE individuel
RPE1_delta = MRPE_G1(Vi1,:) - MRPE_G1(Vi1,1);
RPE2_delta = MRPE_G2(Vi2,:) - MRPE_G2(Vi2,1);
[Mu_rpe1_d, ~] = gStats(RPE1_delta);
[Mu_rpe2_d, ~] = gStats(RPE2_delta);

% Centrage individuel des variables Goubault
MV_G1_delta = MV_G1(Vi1,:,:) - MV_G1(Vi1,1,:);
MV_G2_delta = MV_G2(Vi2,:,:) - MV_G2(Vi2,1,:);

for iScore = 1:2  % fenetre 1 = Mean, fenetre 2 = Median

    if iScore == 1
        fig_name  = 'Fig_Goubault_Delta_G1vsG2_Mean';
        score_str = 'Mean';
    else
        fig_name  = 'Fig_Goubault_Delta_G1vsG2_Median';
        score_str = 'Median';
    end

    figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [30 30 1900 700]);

    for iV = 1:N_vars
        subplot(2, 5, iV);
        hold on;

        Mv1 = squeeze(MV_G1_delta(:, :, iV));
        Mv2 = squeeze(MV_G2_delta(:, :, iV));

        if iScore == 1
            [Mu1, SD1] = gStats(Mv1);
            [Mu2, SD2] = gStats(Mv2);
        else
            [Mu1, SD1] = gStatsMedian(Mv1);
            [Mu2, SD2] = gStatsMedian(Mv2);
        end

        yyaxis left
        fill([T fliplr(T)], [Mu1+SD1 fliplr(Mu1-SD1)], col_G1, 'FaceAlpha',0.15,'EdgeColor','none');
        fill([T fliplr(T)], [Mu2+SD2 fliplr(Mu2-SD2)], col_G2, 'FaceAlpha',0.15,'EdgeColor','none');
        h1 = plot(T, Mu1, '-', 'Color', col_G1, 'LineWidth', 2.5);
        h2 = plot(T, Mu2, '-', 'Color', col_G2, 'LineWidth', 2.5);
        plot(T, zeros(1,length(T)), '--k', 'LineWidth', 0.8);
        ylabel(['Delta ' Vars{iV,5}], 'FontSize', 7);
        ax = gca; ax.YColor = [0.2 0.2 0.2];

        yyaxis right
        plot(T, Mu_rpe1_d, ':', 'Color', col_rpe_G1, 'LineWidth', 1.5);
        plot(T, Mu_rpe2_d, ':', 'Color', col_rpe_G2, 'LineWidth', 1.5);
        ylabel('Delta RPE (Borg CR-10)', 'FontSize', 7, 'Color', [0.4 0.1 0.6]);
        ax.YColor = [0.4 0.1 0.6];

        if iV == 1
            legend([h1 h2], {['G1 n=' num2str(length(Vi1))], ...
                ['G2 n=' num2str(length(Vi2))]}, ...
                'Location','northwest','FontSize',7);
        end

        xlabel('Temps normalise (%)', 'FontSize', 7);
        title(Vars{iV,1}, 'FontSize', 7, 'FontWeight','bold','Interpreter','none');
        xlim([5 105]); xticks(T); grid on; box on; hold off;
    end

    sgtitle(['Variables Goubault DELTA - ' score_str ' - G1 (rouge) vs G2 (bleu) - Evolution depuis bin1'], ...
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