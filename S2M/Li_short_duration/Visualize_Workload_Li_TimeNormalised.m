%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Visualization — Workload & RPE vs Temps normalisé (%) — BRUT   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calcul depuis les données brutes frame par frame :
%   1. Charger S10X_Li_Accel.mat → signal brut (N_frames × 69)
%   2. Filtrer (Butterworth 0.5-29 Hz, comme Goubault)
%   3. Calculer |accel| et jerk = diff(accel)/dt par frame
%   4. Agréger par segment/axe/total (SOMMES, cohérent avec Workload_Li_G1)
%   5. Rééchantillonner sur 10 intervalles normalisés (mean par intervalle)
%   6. Interpoler RPE sur les 10 intervalles
%   7. Normalisation par la moyenne individuelle (ratio / baseline)
%   8. Moyenne ± SD de groupe + droite de tendance linéaire rouge
%
% Coupure RPE : chaque participant est tronqué dès que son RPE discret
%   atteint 7 deux fois consécutives (critère d'arrêt G1).
%   Les intervalles au-delà de ce point sont mis à NaN.
%
% Convention d'agrégation (identique à Workload_Li_Cycles.m) :
%   ByAxis   : SOMME des |accel| sur les 7 segments pour chaque axe
%   BySegMod : norme √(X²+Y²+Z²) pour chaque segment
%   Total    : SOMME sur les 21 colonnes = sum(ByAxis)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS & PARAMETERS
%  -----------------------------------------------------------------------
PathAccel = 'J:\Piano_Fatigue\Data_Exported\Xsens\';
PathInfo  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
PathSave  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';
PathCycle = 'J:\Piano_Fatigue\Data_Exported\';
PathData26 = 'J:\Piano_Fatigue\Data_Exported\Features_XSENS\Li\';  % pour BySegMod 26 participants

fs     = 60;    % Hz
dt     = 1/fs;
N_bins = 10;    % 10 intervalles (méthode identique au script de référence)

G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];

Segments = {'L5', 'T8', 'Head', 'Shoulder', 'Arm', 'Forearm', 'Hand'};
SegCols  = [4,5,6; 13,14,15; 19,20,21; 22,23,24; 25,26,27; 28,29,30; 31,32,33];

N_seg        = length(Segments);
Axes         = {'X', 'Y', 'Z'};
AxCols_X     = SegCols(:, 1);
AxCols_Y     = SegCols(:, 2);
AxCols_Z     = SegCols(:, 3);
AxColsByAxis = {AxCols_X, AxCols_Y, AxCols_Z};
AllSegCols   = SegCols(:);

[b_filt, a_filt] = butter(2, [0.5 29] / (fs/2), 'bandpass');

%% -----------------------------------------------------------------------
%  MATRICES DE STOCKAGE
%  -----------------------------------------------------------------------
Mat_Total_Accel = NaN(length(G1_Liszt), N_bins);
Mat_Total_Jerk  = NaN(length(G1_Liszt), N_bins);
Mat_Axis_Accel  = NaN(length(G1_Liszt), N_bins, 3);
Mat_Axis_Jerk   = NaN(length(G1_Liszt), N_bins, 3);
Mat_Seg_Accel   = NaN(length(G1_Liszt), N_bins, N_seg);
Mat_Seg_Jerk    = NaN(length(G1_Liszt), N_bins, N_seg);
Mat_SegAxis_Accel = NaN(length(G1_Liszt), N_bins, N_seg, 3);  % 7 seg × 3 axes
Mat_SegAxis_Jerk  = NaN(length(G1_Liszt), N_bins, N_seg, 3);
Mat_RPE         = NaN(length(G1_Liszt), N_bins);

cycles_info = load(fullfile(PathCycle, 'cycles_Li.mat'), 'participants');
valid_subj  = false(length(G1_Liszt), 1);

% Précharger RPE pour les 26 participants (utilisé aussi pour BySegMod_26)
tmp_info26 = load(fullfile(PathInfo, 'info_participants_corrected.mat'), 'Info_participants');

%% -----------------------------------------------------------------------
%  MAIN LOOP
%  -----------------------------------------------------------------------
for iG = 1:length(G1_Liszt)
    iP     = G1_Liszt(iG);
    SubjID = sprintf('S1%02d', iP);
    fname  = fullfile(PathAccel, sprintf('%s_Li_Accel.mat', SubjID));

    if ~exist(fname, 'file')
        fprintf('  MISSING: %s\n', fname);
        continue;
    end
    fprintf('Processing %s ...\n', SubjID);

    %% Charger et découper sur la portion Liszt
    tmp      = load(fname, 'Accel');
    data_raw = tmp.Accel.data;
    fs_file  = tmp.Accel.freq;

    t_cycles_p = cycles_info.participants(iP).t_do(:);
    t_start_li = t_cycles_p(1);
    t_end_li   = t_cycles_p(end);

    idx_s    = max(1, round(t_start_li * fs_file));
    idx_e    = min(size(data_raw,1), round(t_end_li * fs_file));
    data     = data_raw(idx_s:idx_e, :);
    N_frames = size(data, 1);

    fprintf('  %s: total=%d frames, Liszt=%d frames (%.1f-%.1fs)\n', ...
        SubjID, size(data_raw,1), N_frames, t_start_li, t_end_li);

    if N_frames < 7, fprintf('  Signal too short, skipping\n'); continue; end

    %% Filtrage et jerk
    data_filt = filtfilt(b_filt, a_filt, data);

    % Jerk — différence centrée à 3 points (plus précis que la différence forward)
    % jerk(t) = (accel(t+1) - accel(t-1)) / (2*dt)
    % Bords : différence forward (t=1) et backward (t=end)
    jerk_filt = zeros(size(data_filt));
    jerk_filt(2:end-1, :) = (data_filt(3:end, :) - data_filt(1:end-2, :)) / (2*dt);
    jerk_filt(1, :)       = (data_filt(2, :) - data_filt(1, :)) / dt;
    jerk_filt(end, :)     = (data_filt(end, :) - data_filt(end-1, :)) / dt;

    %% Charger RPE discret
    try
        rpe_raw  = tmp_info26.Info_participants(iP).Liszt(:);
        rpe_raw(isnan(rpe_raw)) = [];
    catch
        fprintf('  Cannot load RPE for %s, skipping\n', SubjID);
        continue;
    end

    % Correction S122 : supprimer la dernière valeur RPE
    if iP == 22
        rpe_raw = rpe_raw(1:end-1);
    end

    %% Interpoler RPE sur les frames
    t_frames   = (0:N_frames-1)' / fs;
    N_rpe      = length(rpe_raw);
    t_rpe      = (1:N_rpe)' * 30;
    t_rpe      = min(t_rpe, t_frames(end));
    rpe_frames = interp1(t_rpe, double(rpe_raw), t_frames, 'linear', 'extrap');
    rpe_frames = max(0, min(10, rpe_frames));

    %% Découper en N_bins intervalles (méthode du script de référence)
    % Borne i = round(N_frames * i / N_bins)
    bin_edges    = zeros(1, N_bins+1);
    bin_edges(1) = 1;
    for i = 1:N_bins
        bin_edges(i+1) = round(N_frames * i / N_bins);
    end

    rpe_bins    = NaN(N_bins, 1);
    total_accel = NaN(N_bins, 1);
    total_jerk  = NaN(N_bins, 1);
    axis_accel  = NaN(N_bins, 3);
    axis_jerk   = NaN(N_bins, 3);
    seg_accel   = NaN(N_bins, N_seg);
    seg_jerk    = NaN(N_bins, N_seg);
    segaxis_accel = NaN(N_bins, N_seg, 3);
    segaxis_jerk  = NaN(N_bins, N_seg, 3);

    for k = 1:N_bins
        k_s = bin_edges(k);
        k_e = bin_edges(k+1);
        if k_e < k_s, continue; end

        rpe_bins(k) = mean(rpe_frames(k_s:k_e));

        blk_a = data_filt(k_s:k_e, :);
        blk_j = jerk_filt(k_s:k_e, :);

        total_accel(k) = mean(sum(abs(blk_a(:, AllSegCols)), 2));
        total_jerk(k)  = mean(sum(abs(blk_j(:, AllSegCols)), 2));

        for iAx = 1:3
            ax_cols = AxColsByAxis{iAx};
            axis_accel(k, iAx) = mean(sum(abs(blk_a(:, ax_cols)), 2));
            axis_jerk(k, iAx)  = mean(sum(abs(blk_j(:, ax_cols)), 2));
        end

        for iSeg = 1:N_seg
            cols = SegCols(iSeg, :);
            seg_accel(k, iSeg) = mean(sqrt(sum(blk_a(:, cols).^2, 2)));
            seg_jerk(k, iSeg)  = mean(sqrt(sum(blk_j(:, cols).^2, 2)));
        end

        % By Segment × Axis : mean(|accel|) pour chaque segment et chaque axe
        % → toujours positif car abs() appliqué sur frames brutes
        for iSeg = 1:N_seg
            for iAx = 1:3
                col = SegCols(iSeg, iAx);
                segaxis_accel(k, iSeg, iAx) = mean(abs(blk_a(:, col)));
                segaxis_jerk(k, iSeg, iAx)  = mean(abs(blk_j(:, col)));
            end
        end
    end

    %% Stocker
    Mat_RPE(iG, :)            = rpe_bins';
    Mat_Total_Accel(iG, :)    = total_accel';
    Mat_Total_Jerk(iG, :)     = total_jerk';
    Mat_Axis_Accel(iG, :, :)  = axis_accel;
    Mat_Axis_Jerk(iG, :, :)   = axis_jerk;
    Mat_Seg_Accel(iG, :, :)   = seg_accel;
    Mat_Seg_Jerk(iG, :, :)    = seg_jerk;
    Mat_SegAxis_Accel(iG, :, :, :) = segaxis_accel;
    Mat_SegAxis_Jerk(iG, :, :, :)  = segaxis_jerk;
    valid_subj(iG) = true;

    fprintf('  %s: %d bins valides, RPE [%.1f-%.1f]\n', ...
        SubjID, sum(~isnan(rpe_bins)), nanmin(rpe_bins), nanmax(rpe_bins));
end

fprintf('\n%d / %d participants traités\n', sum(valid_subj), length(G1_Liszt));

%% -----------------------------------------------------------------------
%  SAVE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'Workload_Li_TimeNormalised.mat'), ...
    'Mat_RPE', 'Mat_Total_Accel', 'Mat_Total_Jerk', ...
    'Mat_Axis_Accel', 'Mat_Axis_Jerk', ...
    'Mat_Seg_Accel',  'Mat_Seg_Jerk', ...
    'Mat_SegAxis_Accel', 'Mat_SegAxis_Jerk', ...
    'valid_subj', 'Segments', 'Axes', 'N_bins', 'G1_Liszt');
fprintf('Données sauvegardées dans %s\n', PathSave);

%% -----------------------------------------------------------------------
%  NORMALISATION PAR LA MOYENNE DES BINS 1-3 (premier tiers = baseline)
%  signal_norm(t) = signal(t) / mean(signal(bins 1-3))
%  → 1.0 = niveau moyen en début de session (fatigue minimale)
%  → réduit la variance du dénominateur vs normalisation par bin 1 seul
%  Appliquée uniquement pour la visualisation.
%  -----------------------------------------------------------------------
Vi = find(valid_subj);

normByBin1_3 = @(M) M ./ mean(M(:, 1:3), 2);

Norm_Total_Accel = normByBin1_3(Mat_Total_Accel(Vi, :));
Norm_Total_Jerk  = normByBin1_3(Mat_Total_Jerk(Vi, :));

Norm_Axis_Accel = NaN(size(Mat_Axis_Accel(Vi,:,:)));
Norm_Axis_Jerk  = NaN(size(Mat_Axis_Jerk(Vi,:,:)));
for iAx = 1:3
    Norm_Axis_Accel(:,:,iAx) = normByBin1_3(squeeze(Mat_Axis_Accel(Vi,:,iAx)));
    Norm_Axis_Jerk(:,:,iAx)  = normByBin1_3(squeeze(Mat_Axis_Jerk(Vi,:,iAx)));
end

Norm_Seg_Accel = NaN(size(Mat_Seg_Accel(Vi,:,:)));
Norm_Seg_Jerk  = NaN(size(Mat_Seg_Jerk(Vi,:,:)));
for iSeg = 1:N_seg
    Norm_Seg_Accel(:,:,iSeg) = normByBin1_3(squeeze(Mat_Seg_Accel(Vi,:,iSeg)));
    Norm_Seg_Jerk(:,:,iSeg)  = normByBin1_3(squeeze(Mat_Seg_Jerk(Vi,:,iSeg)));
end

% Normalisation BySegAxis (7 seg × 3 axes)
Norm_SegAxis_Accel = NaN(length(Vi), N_bins, N_seg, 3);
Norm_SegAxis_Jerk  = NaN(length(Vi), N_bins, N_seg, 3);
for iSeg = 1:N_seg
    for iAx = 1:3
        Norm_SegAxis_Accel(:,:,iSeg,iAx) = normByBin1_3(squeeze(Mat_SegAxis_Accel(Vi,:,iSeg,iAx)));
        Norm_SegAxis_Jerk(:,:,iSeg,iAx)  = normByBin1_3(squeeze(Mat_SegAxis_Jerk(Vi,:,iSeg,iAx)));
    end
end

%% -----------------------------------------------------------------------
%  FIGURES
%  -----------------------------------------------------------------------
% Abscisse en % : 10, 20, ..., 100
T = (1:N_bins) * (100/N_bins);

[Mu_rpe, SD_rpe] = gStats(Mat_RPE(Vi, :));

SegColors = lines(7);
AxColors  = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.1 0.7 0.3]};
sig_names = {'Acceleration', 'Jerk'};
lbl_norm  = 'Amplitude normalisée (ratio / bins 1-3)';

%% --- FIG 1 : TOTAL ---
figure('Name','Fig_TimeNorm_Total','NumberTitle','off','Position',[50 50 1200 700]);
for iSig = 1:2
    if iSig == 1, MatN = Norm_Total_Accel;
    else,         MatN = Norm_Total_Jerk; end
    [Mu_v, SD_v]  = gStats(MatN);
    [Md_v, SD_md] = gStatsMedian(MatN);
    for iScore = 1:2
        subplot(2, 2, (iScore-1)*2 + iSig);
        if iScore == 1
            plotDual(T, Mu_v, SD_v, Mu_rpe, SD_rpe, [0 0 0], lbl_norm, ...
                ['Mean groupe — Total ' sig_names{iSig}]);
        else
            plotDual(T, Md_v, SD_md, Mu_rpe, SD_rpe, [0 0 0], lbl_norm, ...
                ['Median groupe — Total ' sig_names{iSig}]);
        end
    end
end
sgtitle('Total Workload normalisé vs RPE — Temps normalisé (%)', ...
    'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf, fullfile(PathSave, 'Fig_TimeNorm_Total.png'));

%% --- FIG 2-3 : BY AXIS ---
for iSig = 1:2
    if iSig==1, MatN3 = Norm_Axis_Accel; else, MatN3 = Norm_Axis_Jerk; end
    figure('Name', ['Fig_TimeNorm_ByAxis_' sig_names{iSig}], ...
        'NumberTitle', 'off', 'Position', [50 50 1400 700]);
    for iScore = 1:2
        for iAx = 1:3
            subplot(2, 3, (iScore-1)*3 + iAx);
            Mat_v = squeeze(MatN3(:, :, iAx));
            if iScore == 1
                [Mu_v, SD_v] = gStats(Mat_v);
                plotDual(T, Mu_v, SD_v, Mu_rpe, SD_rpe, AxColors{iAx}, lbl_norm, ...
                    ['Mean groupe — Axe ' Axes{iAx}]);
            else
                [Md_v, SD_md] = gStatsMedian(Mat_v);
                plotDual(T, Md_v, SD_md, Mu_rpe, SD_rpe, AxColors{iAx}, lbl_norm, ...
                    ['Median groupe — Axe ' Axes{iAx}]);
            end
        end
    end
    sgtitle([sig_names{iSig} ' By Axis normalisé vs RPE — Temps normalisé (%)'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig_TimeNorm_ByAxis_' sig_names{iSig} '.png']));
end

%% --- FIG 4-5 : BY SEGMENT MODULE — 26 participants (Features3_XSENS) ---
% On charge Mean.Module_Acceleration et Mean.Module_Jerk depuis Features3
% pour avoir les 26 participants G1 au lieu de 14.
fprintf('\nChargement des 26 participants pour BySegMod...\n');

Mat_Seg_Accel_26 = NaN(length(G1_Liszt), N_bins, N_seg);
Mat_Seg_Jerk_26  = NaN(length(G1_Liszt), N_bins, N_seg);
Mat_RPE_26       = NaN(length(G1_Liszt), N_bins);
valid_26         = false(length(G1_Liszt), 1);

for iG26 = 1:length(G1_Liszt)
    iP26   = G1_Liszt(iG26);
    SubjID26 = sprintf('S1%02d', iP26);
    fname26  = fullfile(PathData26, sprintf('Features3_XSENS_%s_Li.mat', SubjID26));

    if ~exist(fname26, 'file'), continue; end

    tmp26  = load(fname26, 'Features_XSENS');
    feat26 = tmp26.Features_XSENS;
    N_f26  = size(feat26.Mean.Module_Acceleration, 1);

    t_cyc26    = cycles_info.participants(iP26).t_do(:);
    tp_s26 = max(1, floor(t_cyc26(1)));
    tp_e26 = min(N_f26, floor(t_cyc26(end)));
    if tp_e26 <= tp_s26, continue; end
    N_tp26 = tp_e26 - tp_s26 + 1;

    rpe26 = tmp_info26.Info_participants(iP26).Liszt(:);
    rpe26(isnan(rpe26)) = [];
    if iP26 == 22, rpe26 = rpe26(1:end-1); end

    t_rpe26 = (1:length(rpe26))' * 30;
    t_tp26  = (0:N_tp26-1)';
    t_rpe26 = min(t_rpe26, t_tp26(end));
    rpe_tp26 = interp1(t_rpe26, double(rpe26), t_tp26, 'linear', 'extrap');
    rpe_tp26 = max(0, min(10, rpe_tp26));

    be26    = zeros(1, N_bins+1); be26(1) = 1;
    for i = 1:N_bins, be26(i+1) = round(N_tp26 * i / N_bins); end

    mod_a26 = feat26.Mean.Module_Acceleration;
    mod_j26 = feat26.Mean.Module_Jerk;

    for k = 1:N_bins
        ks = be26(k); ke = be26(k+1);
        if ke < ks, continue; end
        Mat_RPE_26(iG26, k) = mean(rpe_tp26(ks:ke));
        gs = tp_s26+ks-1; ge = min(tp_s26+ke-1, N_f26);
        for iSeg = 1:N_seg
            Mat_Seg_Accel_26(iG26, k, iSeg) = mean(mod_a26(gs:ge, iSeg));
            Mat_Seg_Jerk_26(iG26, k, iSeg)  = mean(mod_j26(gs:ge, iSeg));
        end
    end
    valid_26(iG26) = true;
end
% Calculer Total_26 = somme des normes de segments (Features3)
Mat_Total_Accel_26 = sum(Mat_Seg_Accel_26, 3);  % N_part x N_bins
Mat_Total_Jerk_26  = sum(Mat_Seg_Jerk_26,  3);  % N_part x N_bins
fprintf('%d / 26 participants charges pour BySegMod_26\n', sum(valid_26));

%% Sauvegarder les matrices _26 dans le meme fichier .mat
save(fullfile(PathSave, 'Workload_Li_TimeNormalised.mat'), ...
    'Mat_Seg_Accel_26', 'Mat_Seg_Jerk_26', ...
    'Mat_Total_Accel_26', 'Mat_Total_Jerk_26', ...
    'Mat_RPE_26', 'valid_26', '-append');
fprintf('Matrices _26 ajoutees au fichier .mat\n');

Vi26 = find(valid_26);
[Mu_rpe26, SD_rpe26] = gStats(Mat_RPE_26(Vi26, :));

Norm_Seg_Accel_26 = NaN(length(Vi26), N_bins, N_seg);
Norm_Seg_Jerk_26  = NaN(length(Vi26), N_bins, N_seg);
for iSeg = 1:N_seg
    Norm_Seg_Accel_26(:,:,iSeg) = normByBin1_3(squeeze(Mat_Seg_Accel_26(Vi26,:,iSeg)));
    Norm_Seg_Jerk_26(:,:,iSeg)  = normByBin1_3(squeeze(Mat_Seg_Jerk_26(Vi26,:,iSeg)));
end

for iSig = 1:2
    if iSig==1, MatN3 = Norm_Seg_Accel_26; Mu_r = Mu_rpe26; SD_r = SD_rpe26;
    else,        MatN3 = Norm_Seg_Jerk_26;  Mu_r = Mu_rpe26; SD_r = SD_rpe26; end
    figure('Name', ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_26'], ...
        'NumberTitle', 'off', 'Position', [50 50 1800 700]);
    for iScore = 1:2
        for iSeg = 1:N_seg
            subplot(2, N_seg, (iScore-1)*N_seg + iSeg);
            Mat_v = squeeze(MatN3(:, :, iSeg));
            if iScore == 1
                [Mu_v, SD_v] = gStats(Mat_v);
                plotDual(T, Mu_v, SD_v, Mu_r, SD_r, SegColors(iSeg,:), lbl_norm, ...
                    ['Mean groupe — ' Segments{iSeg}]);
            else
                [Md_v, SD_md] = gStatsMedian(Mat_v);
                plotDual(T, Md_v, SD_md, Mu_r, SD_r, SegColors(iSeg,:), lbl_norm, ...
                    ['Median groupe — ' Segments{iSeg}]);
            end
        end
    end
    sgtitle([sig_names{iSig} ' By Segment normalisé vs RPE — Temps normalisé (%) — n=26'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_26.png']));
end

fprintf('\nToutes les figures sauvegardées dans %s\n', PathSave);

%% --- FIG 4b : BY SEGMENT MODULE — 26 participants — Normalisation Bin 1 ---
normBin1 = @(M) M ./ M(:, 1);

Norm_B1_Seg_Accel_26 = NaN(length(Vi26), N_bins, N_seg);
Norm_B1_Seg_Jerk_26  = NaN(length(Vi26), N_bins, N_seg);
for iSeg = 1:N_seg
    Norm_B1_Seg_Accel_26(:,:,iSeg) = normBin1(squeeze(Mat_Seg_Accel_26(Vi26,:,iSeg)));
    Norm_B1_Seg_Jerk_26(:,:,iSeg)  = normBin1(squeeze(Mat_Seg_Jerk_26(Vi26,:,iSeg)));
end

for iSig = 1:2
    if iSig==1, MatN3 = Norm_B1_Seg_Accel_26; Mu_r = Mu_rpe26; SD_r = SD_rpe26;
    else,        MatN3 = Norm_B1_Seg_Jerk_26;  Mu_r = Mu_rpe26; SD_r = SD_rpe26; end
    figure('Name', ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_26_Bin1'], ...
        'NumberTitle', 'off', 'Position', [50 50 1800 700]);
    for iScore = 1:2
        for iSeg = 1:N_seg
            subplot(2, N_seg, (iScore-1)*N_seg + iSeg);
            Mat_v = squeeze(MatN3(:, :, iSeg));
            if iScore == 1
                [Mu_v, SD_v] = gStats(Mat_v);
                plotDual(T, Mu_v, SD_v, Mu_r, SD_r, SegColors(iSeg,:), ...
                    'Amplitude normalisée (ratio / bin 1)', ...
                    ['Mean groupe — ' Segments{iSeg}]);
            else
                [Md_v, SD_md] = gStatsMedian(Mat_v);
                plotDual(T, Md_v, SD_md, Mu_r, SD_r, SegColors(iSeg,:), ...
                    'Amplitude normalisée (ratio / bin 1)', ...
                    ['Median groupe — ' Segments{iSeg}]);
            end
        end
    end
    sgtitle([sig_names{iSig} ' By Segment — Normalisation Bin 1 — n=26'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_26_Bin1.png']));
end

%% --- FIG 4c : BY SEGMENT MODULE — 26 participants — Normalisation Max ---
normMaxInd = @(M) M ./ max(M, [], 2);

Norm_Mx_Seg_Accel_26 = NaN(length(Vi26), N_bins, N_seg);
Norm_Mx_Seg_Jerk_26  = NaN(length(Vi26), N_bins, N_seg);
for iSeg = 1:N_seg
    Norm_Mx_Seg_Accel_26(:,:,iSeg) = normMaxInd(squeeze(Mat_Seg_Accel_26(Vi26,:,iSeg)));
    Norm_Mx_Seg_Jerk_26(:,:,iSeg)  = normMaxInd(squeeze(Mat_Seg_Jerk_26(Vi26,:,iSeg)));
end

for iSig = 1:2
    if iSig==1, MatN3 = Norm_Mx_Seg_Accel_26; Mu_r = Mu_rpe26; SD_r = SD_rpe26;
    else,        MatN3 = Norm_Mx_Seg_Jerk_26;  Mu_r = Mu_rpe26; SD_r = SD_rpe26; end
    figure('Name', ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_26_Max'], ...
        'NumberTitle', 'off', 'Position', [50 50 1800 700]);
    for iScore = 1:2
        for iSeg = 1:N_seg
            subplot(2, N_seg, (iScore-1)*N_seg + iSeg);
            Mat_v = squeeze(MatN3(:, :, iSeg));
            if iScore == 1
                [Mu_v, SD_v] = gStats(Mat_v);
                plotDual(T, Mu_v, SD_v, Mu_r, SD_r, SegColors(iSeg,:), ...
                    'Amplitude normalisée (ratio / max)', ...
                    ['Mean groupe — ' Segments{iSeg}]);
            else
                [Md_v, SD_md] = gStatsMedian(Mat_v);
                plotDual(T, Md_v, SD_md, Mu_r, SD_r, SegColors(iSeg,:), ...
                    'Amplitude normalisée (ratio / max)', ...
                    ['Median groupe — ' Segments{iSeg}]);
            end
        end
    end
    sgtitle([sig_names{iSig} ' By Segment — Normalisation Max — n=26'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig_TimeNorm_BySegment_' sig_names{iSig} '_26_Max.png']));
end

%% --- FIG 6-11 : BY SEGMENT × AXIS (3 fenêtres par signal) ---
% Une fenêtre par axe (X, Y, Z), chaque fenêtre = 2 lignes × 7 colonnes
% (Mean groupe ligne 1, Median groupe ligne 2)
% → 3 axes × 2 signaux = 6 figures au total

AxColors_full = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.1 0.7 0.3]};

for iSig = 1:2
    if iSig == 1
        MatSA = Norm_SegAxis_Accel;
        sig_label = 'Acceleration';
    else
        MatSA = Norm_SegAxis_Jerk;
        sig_label = 'Jerk';
    end

    for iAx = 1:3
        fig_name = sprintf('Fig_TimeNorm_BySegAxis_%s_%s', sig_label, Axes{iAx});
        figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [30 30 1800 700]);

        for iScore = 1:2
            for iSeg = 1:N_seg
                subplot(2, N_seg, (iScore-1)*N_seg + iSeg);
                Mat_v = squeeze(MatSA(:, :, iSeg, iAx));
                if iScore == 1
                    [Mu_v, SD_v] = gStats(Mat_v);
                    ttl = sprintf('Mean — %s', Segments{iSeg});
                else
                    [Mu_v, SD_v] = gStatsMedian(Mat_v);
                    ttl = sprintf('Median — %s', Segments{iSeg});
                end
                plotDual(T, Mu_v, SD_v, Mu_rpe, SD_rpe, ...
                    AxColors_full{iAx}, lbl_norm, ttl);
            end
        end

        sgtitle(sprintf('%s By Segment — Axe %s normalisé vs RPE — Temps normalisé (%%)', ...
            sig_label, Axes{iAx}), 'FontSize', 12, 'FontWeight', 'bold');
        saveas(gcf, fullfile(PathSave, [fig_name '.png']));
    end
end

fprintf('Figures BySegAxis sauvegardées.\n');

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

function plotDual(T, Mu_v, SD_v, Mu_r, SD_r, col_v, lbl_v, ttl)
% Double axe Y : variable cinématique normalisée (gauche) + RPE (droite).
% Pas de smoothing (10 points = déjà lisible).
% Droite de tendance linéaire en rouge (polyfit degré 1).
    T    = T(:)';    Mu_v = Mu_v(:)';  SD_v = SD_v(:)';
    Mu_r = Mu_r(:)'; SD_r = SD_r(:)';
    col_r     = [0.5 0.1 0.5];
    col_trend = [0.85 0.1 0.1];

    % Droite de régression sur les points valides
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
        plot(T(valid_idx), trend, '-', 'Color', col_trend, 'LineWidth', 2.5);
    end
    ylabel(lbl_v, 'FontSize', 8, 'Color', col_v);
    ax = gca; ax.YColor = col_v;

    yyaxis right
    fill([T fliplr(T)], [Mu_r+SD_r fliplr(Mu_r-SD_r)], col_r, ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(T, Mu_r, '--', 'Color', col_r, 'LineWidth', 2);
    ylabel('RPE (Borg CR-10)', 'FontSize', 8, 'Color', col_r);
    ax.YColor = col_r;
    ylim([0 10]);

    xlabel('Temps normalisé (%)', 'FontSize', 8);
    title(ttl, 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlim([5 105]); xticks(T); grid on; box on; hold off
end