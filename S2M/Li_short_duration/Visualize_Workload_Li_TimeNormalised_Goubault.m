%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Visualization — Variables Goubault vs RPE — Temps normalisé    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pipeline identique à Visualize_Workload_Li_TimeNormalised.m :
%   1. Charger Features3_XSENS_S10X_Li.mat (1 timepoint/seconde)
%   2. Extraire les 10 variables Goubault 2023 (Chord task, Table III)
%   3. Découper en 10 bins temporels normalisés
%   4. Normaliser par la moyenne des bins 1-3 (baseline)
%   5. Moyenne ± SD de groupe + droite de tendance linéaire rouge
%
% 10 variables (Goubault 2023, Chord task — Table III) :
%   Colonnes dans Features3 (21 cols = 7 seg × 3 axes) :
%     Segments : L5(1), T8(2), Head(3), Shoulder(4), Arm(5), Forearm(6), Hand(7)
%     Axe X = seg*3-2, Y = seg*3-1, Z = seg*3  (1-based)
%   Colonnes Module (7 cols) : même ordre segments
%
%   #  Segment   Feature          Signal                 Col
%   1  Hand      MedianFreq       Acceleration Y         20
%   2  Hand      Percentile90     Acceleration X         19
%   3  Hand      PeakPower        Module_Acceleration     7
%   4  Hand      SpectralEntropy  Acceleration Y         20
%   5  Head      PeakPower        Angular_Velocity X      7
%   6  Forearm   Mean             Acceleration Y         17
%   7  Forearm   SpectralEntropy  Module_AngVel           6
%   8  Forearm   PeakPower_Freq   Module_AngVel           6
%   9  Shoulder  Mean             Angular_Velocity X     10
%  10  Shoulder  Mean             Angular_Velocity Y     11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS & PARAMETERS
%  -----------------------------------------------------------------------
PathData  = 'J:\Piano_Fatigue\Data_Exported\Features_XSENS\Li\';
PathInfo  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
PathCycle = 'J:\Piano_Fatigue\Data_Exported\';
PathSave  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';

N_bins   = 10;

G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];

%% -----------------------------------------------------------------------
%  DÉFINITION DES 10 VARIABLES GOUBAULT
%  {nom_affichage, feature_field, signal_field, col_1based}
%  -----------------------------------------------------------------------
Vars = {
    'Hand — MedianFreq — Accel Y',            'MedianFreq',      'Acceleration',              20;
    'Hand — Prc90 — Accel X',                 'Percentile90',    'Acceleration',              19;
    'Hand — PeakPower — Accel Module',         'PeakPower',       'Module_Acceleration',        7;
    'Hand — SpectralEntropy — Accel Y',        'SpectralEntropy', 'Acceleration',              20;
    'Head — PeakPower — AngVel X',             'PeakPower',       'Angular_Velocity',           7;
    'Forearm — Mean — Accel Y',                'Mean',            'Acceleration',              17;
    'Forearm — SpectralEntropy — AngVel Mod',  'SpectralEntropy', 'Module_Angular_Velocity',    6;
    'Forearm — PeakPowerFreq — AngVel Mod',    'PeakPower_Freq',  'Module_Angular_Velocity',    6;
    'Shoulder — Mean — AngVel X',              'Mean',            'Angular_Velocity',          10;
    'Shoulder — Mean — AngVel Y',              'Mean',            'Angular_Velocity',          11;
};

N_vars = size(Vars, 1);

%% -----------------------------------------------------------------------
%  MATRICES DE STOCKAGE
%  -----------------------------------------------------------------------
Mat_Vars   = NaN(length(G1_Liszt), N_bins, N_vars);
Mat_RPE    = NaN(length(G1_Liszt), N_bins);
valid_subj = false(length(G1_Liszt), 1);

cycles_info = load(fullfile(PathCycle, 'cycles_Li.mat'), 'participants');
tmp_info    = load(fullfile(PathInfo,  'info_participants_corrected.mat'), 'Info_participants');

%% -----------------------------------------------------------------------
%  MAIN LOOP
%  -----------------------------------------------------------------------
for iG = 1:length(G1_Liszt)
    iP     = G1_Liszt(iG);
    SubjID = sprintf('S1%02d', iP);
    fname  = fullfile(PathData, sprintf('Features3_XSENS_%s_Li.mat', SubjID));

    if ~exist(fname, 'file')
        fprintf('  MISSING: %s\n', fname);
        continue;
    end
    fprintf('Processing %s ...\n', SubjID);

    %% Charger Features3
    tmp  = load(fname, 'Features_XSENS');
    feat = tmp.Features_XSENS;
    N_feat_total = size(feat.Mean.Acceleration, 1);  % durée totale en secondes

    %% Délimiter la portion Liszt
    t_cycles_p = cycles_info.participants(iP).t_do(:);
    t_start_li = t_cycles_p(1);
    t_end_li   = t_cycles_p(end);

    % Indices dans Features3 (1 timepoint = 1 seconde)
    tp_start = max(1, floor(t_start_li));
    tp_end   = min(N_feat_total, floor(t_end_li));  % pas de marge post-jeu

    if tp_end <= tp_start
        fprintf('  %s: indices invalides, skip\n', SubjID);
        continue;
    end
    N_tp = tp_end - tp_start + 1;
    fprintf('  %s: %d timepoints (%.0f-%.0fs)\n', SubjID, N_tp, t_start_li, t_end_li);

    %% RPE discret → interpoler sur les N_tp timepoints
    rpe_raw = tmp_info.Info_participants(iP).Liszt(:);
    rpe_raw(isnan(rpe_raw)) = [];

    if iP == 22
        rpe_raw = rpe_raw(1:end-1);
    end

    N_rpe    = length(rpe_raw);
    t_rpe    = (1:N_rpe)' * 30;
    t_tp     = (0:N_tp-1)';
    t_rpe    = min(t_rpe, t_tp(end));

    rpe_tp   = interp1(t_rpe, double(rpe_raw), t_tp, 'linear', 'extrap');
    rpe_tp   = max(0, min(10, rpe_tp));

    %% Découper en N_bins intervalles
    bin_edges    = zeros(1, N_bins+1);
    bin_edges(1) = 1;
    for i = 1:N_bins
        bin_edges(i+1) = round(N_tp * i / N_bins);
    end

    rpe_bins = NaN(N_bins, 1);
    var_bins = NaN(N_bins, N_vars);

    for k = 1:N_bins
        k_s = bin_edges(k);
        k_e = bin_edges(k+1);
        if k_e < k_s, continue; end

        rpe_bins(k) = mean(rpe_tp(k_s:k_e));

        % Indices globaux dans Features3
        g_s = tp_start + k_s - 1;
        g_e = tp_start + k_e - 1;
        g_e = min(g_e, N_feat_total);

        for iV = 1:N_vars
            feat_name = Vars{iV, 2};
            sig_name  = Vars{iV, 3};
            col       = Vars{iV, 4};

            try
                sig_data = feat.(feat_name).(sig_name);  % (N_feat × N_cols)
                blk      = sig_data(g_s:g_e, col);

                % Pour les variables Mean (signal signé) : valeur absolue
                % afin d'éviter l'annulation de signe entre participants
                if strcmp(feat_name, 'Mean')
                    blk = abs(blk);
                end

                var_bins(k, iV) = mean(blk);
            catch
                % champ manquant
            end
        end
    end

    %% Stocker
    Mat_RPE(iG, :)     = rpe_bins';
    Mat_Vars(iG, :, :) = var_bins;
    valid_subj(iG)     = true;

    fprintf('  %s: OK — RPE [%.1f-%.1f]\n', SubjID, ...
        nanmin(rpe_bins), nanmax(rpe_bins));
end

fprintf('\n%d / %d participants traités\n', sum(valid_subj), length(G1_Liszt));

%% -----------------------------------------------------------------------
%  SAVE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'Workload_Li_TimeNormalised_Goubault.mat'), ...
    'Mat_Vars', 'Mat_RPE', 'valid_subj', 'Vars', 'N_bins', 'G1_Liszt');
fprintf('Sauvegardé : Workload_Li_TimeNormalised_Goubault.mat\n');

%% -----------------------------------------------------------------------
%  NORMALISATION PAR LA MOYENNE DES BINS 1-3
%  Uniquement pour les variables dont le signal n'est pas centré sur 0.
%  Les variables Mean (signées, oscillent autour de 0) ne sont PAS normalisées
%  car diviser par une valeur proche de 0 produit des ratios aberrants.
%  -----------------------------------------------------------------------
Vi = find(valid_subj);

normByBin1_3 = @(M) M ./ mean(M(:, 1:3), 2);

% Flag : 1 = normaliser par bins 1-3, 0 = brut
% Toutes les variables sont maintenant normalisables :
% - Mean a été passée en valeur absolue → plus de division par ~0
DoNorm = ones(N_vars, 1);

Norm_Vars = NaN(length(Vi), N_bins, N_vars);
for iV = 1:N_vars
    M = squeeze(Mat_Vars(Vi,:,iV));
    if DoNorm(iV)
        Norm_Vars(:,:,iV) = normByBin1_3(M);
    else
        Norm_Vars(:,:,iV) = M;  % unités brutes
    end
end

%% -----------------------------------------------------------------------
%  FIGURES — 1 figure, 2 lignes × 5 colonnes
%  -----------------------------------------------------------------------
T = (1:N_bins) * (100/N_bins);  % 10, 20, ..., 100

[Mu_rpe, SD_rpe] = gStats(Mat_RPE(Vi, :));

% Couleurs par segment
SegColors = [
    0.75 0.10 0.10;   % 1  Hand      rouge
    0.75 0.10 0.10;   % 2  Hand
    0.75 0.10 0.10;   % 3  Hand
    0.75 0.10 0.10;   % 4  Hand
    0.50 0.10 0.80;   % 5  Head      violet
    0.10 0.50 0.80;   % 6  Forearm   bleu
    0.10 0.50 0.80;   % 7  Forearm
    0.10 0.50 0.80;   % 8  Forearm
    0.10 0.70 0.30;   % 9  Shoulder  vert
    0.10 0.70 0.30;   % 10 Shoulder
];

figure('Name','Fig_TimeNorm_Goubault','NumberTitle','off', ...
    'Position',[30 30 1900 820]);

for iV = 1:N_vars
    subplot(2, 5, iV);
    Mat_v = squeeze(Norm_Vars(:,:,iV));
    [Mu_v, SD_v] = gStats(Mat_v);
    if DoNorm(iV)
        if strcmp(Vars{iV,2}, 'Mean')
            lbl = '|Mean| normalisé (ratio / bins 1-3)';
        else
            lbl = 'Amplitude normalisée (ratio / bins 1-3)';
        end
    else
        lbl = 'Amplitude brute (unités originales)';
    end
    plotDual(T, Mu_v, SD_v, Mu_rpe, SD_rpe, SegColors(iV,:), lbl, Vars{iV,1});
end

sgtitle('Goubault 2023 — Top 10 variables Chord task — Temps normalisé (%)', ...
    'FontSize', 12, 'FontWeight', 'bold');

saveas(gcf, fullfile(PathSave, 'Fig_TimeNorm_Goubault.png'));
fprintf('Figure sauvegardée : Fig_TimeNorm_Goubault.png\n');

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

function plotDual(T, Mu_v, SD_v, Mu_r, SD_r, col_v, lbl_v, ttl)
    T    = T(:)';    Mu_v = Mu_v(:)';  SD_v = SD_v(:)';
    Mu_r = Mu_r(:)'; SD_r = SD_r(:)';
    col_r     = [0.5 0.1 0.5];
    col_trend = [0.85 0.1 0.1];

    valid_idx = find(~isnan(Mu_v));
    if length(valid_idx) >= 2
        p     = polyfit(T(valid_idx), Mu_v(valid_idx), 1);
        trend = polyval(p, T(valid_idx));
    else
        trend = [];
    end

    yyaxis left
    hold on
    fill([T fliplr(T)], [Mu_v+SD_v fliplr(Mu_v-SD_v)], col_v, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(T, Mu_v, '-', 'Color', col_v, 'LineWidth', 2.5);
    if ~isempty(trend)
        plot(T(valid_idx), trend, '-', 'Color', col_trend, 'LineWidth', 2);
    end
    ylabel(lbl_v, 'FontSize', 7, 'Color', col_v);
    ax = gca; ax.YColor = col_v;

    yyaxis right
    fill([T fliplr(T)], [Mu_r+SD_r fliplr(Mu_r-SD_r)], col_r, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(T, Mu_r, '--', 'Color', col_r, 'LineWidth', 1.5);
    ylabel('RPE (Borg CR-10)', 'FontSize', 7, 'Color', col_r);
    ax.YColor = col_r;
    ylim([0 10]);

    xlabel('Temps normalisé (%)', 'FontSize', 7);
    title(ttl, 'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlim([5 105]); xticks(T); grid on; box on; hold off
end