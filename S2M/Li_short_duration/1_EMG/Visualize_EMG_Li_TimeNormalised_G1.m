%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Visualization - EMG vs RPE vs Temps normalise (%) - G1         %%%%
%%%%   v2 : Normalisation %MVC + debug Activity + outliers SE          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGEMENTS v2 :
%   1. Amplitude normalisee en %MVC
%      Fichier : EMG_Max_S1XX_MVC.mat  champ : EMG_max (minuscule)
%      Vecteur 47x1 - indices 43:47 = Biceps Triceps DeltAnt DeltMed SupTrap
%   2. Activity : diagnostic unites dans la console (1er participant)
%      + normalisation par mediane bin1 pour rendre l'echelle lisible
%   3. SampleEntropy : nettoyage outliers +-3SD sur la figure Median
%      + commandes console en Section 9 pour identifier les participants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS & PARAMETERS
%  -----------------------------------------------------------------------
PathBase  = 'J:\Piano_Fatigue\Data_Exported\';
PathAmp   = fullfile(PathBase, 'Amplitude_Li\');
PathAct   = fullfile(PathBase, 'Activity_Li\');
PathMob   = fullfile(PathBase, 'Mobility_Li\');
PathSE    = fullfile(PathBase, 'SampleEntropy_Li\');
PathTFRmf = fullfile(PathBase, 'TFR_MedianFreq_Li\');
PathTFRse = fullfile(PathBase, 'TFR_SpectralEntropy_Li\');
PathMVC   = fullfile(PathBase, 'EMG_Max_MVC\');
PathCycle = 'J:\Piano_Fatigue\Data_Exported\';
PathInfo  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
PathSave  = fileparts(mfilename('fullpath'));  % dossier du script

N_bins   = 10;
G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];

Channels = {'Biceps', 'Triceps', 'DeltAnt', 'DeltMed', 'SupTrap'};
N_ch     = length(Channels);

% Indices MATLAB (1-based) des 5 canaux surface dans le vecteur EMG_max (47 elements)
% Ordre : 42 HD-EMG puis Biceps(43) Triceps(44) DeltAnt(45) DeltMed(46) SupTrap(47)
MVC_idx = [43, 44, 45, 46, 47];

VarsCycle = {
    'Amplitude',     PathAmp,   'Amplitude_%s_Li_EMG.mat',     'Amplitude';
    'Activity',      PathAct,   'Activity%s_Li_EMG.mat',        'Activity';
    'Mobility',      PathMob,   'Mobility%s_Li_EMG.mat',        'Mobility';
    'SampleEntropy', PathSE,    'SampleEntropy_%s_Li_EMG.mat',  'SampleEntropy';
};
VarsTFR = {
    'TFR_MedianFreq',       PathTFRmf, 'TFR_MedianFreq_%s_Li_EMG.mat',       'TFR_MedianFreq',       'MedianFreq';
    'TFR_SpectralEntropy',  PathTFRse, 'TFR_SpectralEntropy_%s_Li_EMG.mat',  'TFR_SpectralEntropy',  'SpectralEntropy';
};

N_varCycle = size(VarsCycle, 1);
N_varTFR   = size(VarsTFR, 1);
N_vars     = N_varCycle + N_varTFR;

%% -----------------------------------------------------------------------
%  MATRICES DE STOCKAGE
%  -----------------------------------------------------------------------
Mat_EMG    = NaN(length(G1_Liszt), N_bins, N_ch, N_vars);
Mat_RPE    = NaN(length(G1_Liszt), N_bins);
valid_subj = false(length(G1_Liszt), 1);
Mat_MVC    = NaN(length(G1_Liszt), N_ch);   % MVC par participant x canal

cycles_info = load(fullfile(PathCycle, 'cycles_Li.mat'), 'participants');
tmp_info    = load(fullfile(PathInfo,  'info_participants_corrected.mat'), 'Info_participants');

Activity_debug_done = false;

%% -----------------------------------------------------------------------
%  MAIN LOOP
%  -----------------------------------------------------------------------
for iG = 1:length(G1_Liszt)
    iP     = G1_Liszt(iG);
    SubjID = sprintf('S1%02d', iP);
    fprintf('Processing %s ...\n', SubjID);

    %% Cycles Liszt
    t_do = cycles_info.participants(iP).t_do(:);
    if isempty(t_do)
        fprintf('  %s: pas de cycles Liszt, skip\n', SubjID); continue;
    end
    N_cycles = length(t_do) - 1;
    t_start  = t_do(1);
    t_end    = t_do(end);
    duree    = t_end - t_start;

    %% RPE -> interpolation sur N_cycles
    rpe_raw = tmp_info.Info_participants(iP).Liszt(:);
    rpe_raw(isnan(rpe_raw)) = [];
    if iP == 22, rpe_raw = rpe_raw(1:end-1); end
    N_rpe   = length(rpe_raw);
    t_rpe   = (1:N_rpe)' * 30;
    t_cyc   = linspace(0, duree, N_cycles)';
    t_rpe   = min(t_rpe, t_cyc(end));
    rpe_cyc = interp1(t_rpe, double(rpe_raw), t_cyc, 'linear', 'extrap');
    rpe_cyc = max(0, min(10, rpe_cyc));

    %% Bins temporels
    bin_edges    = zeros(1, N_bins+1);
    bin_edges(1) = 1;
    for i = 1:N_bins
        bin_edges(i+1) = round(N_cycles * i / N_bins);
    end
    rpe_bins = NaN(N_bins, 1);
    for k = 1:N_bins
        ks = bin_edges(k); ke = bin_edges(k+1);
        if ke < ks, continue; end
        rpe_bins(k) = mean(rpe_cyc(ks:ke));
    end

    %% ---------------------------------------------------------------
    %  CHARGEMENT MVC
    %  Fichier : EMG_Max_S1XX_MVC.mat  |  champ : EMG_max (minuscule)
    %  Vecteur 47x1, indices 43-47 = 5 muscles de surface
    %% ---------------------------------------------------------------
    fname_mvc = fullfile(PathMVC, sprintf('EMG_Max_%s_MVC.mat', SubjID));
    mvc_ok    = false;
    if exist(fname_mvc, 'file')
        tmp_mvc = load(fname_mvc, 'EMG_max');
        if isfield(tmp_mvc, 'EMG_max')
            vec_mvc = tmp_mvc.EMG_max(:);
            if length(vec_mvc) >= 47
                for iC = 1:N_ch
                    Mat_MVC(iG, iC) = vec_mvc(MVC_idx(iC));
                end
                mvc_ok = true;
                fprintf('  MVC OK : [%.1f %.1f %.1f %.1f %.1f] uV\n', ...
                    Mat_MVC(iG,1), Mat_MVC(iG,2), Mat_MVC(iG,3), ...
                    Mat_MVC(iG,4), Mat_MVC(iG,5));
            else
                fprintf('  MVC: vecteur trop court (%d) pour %s\n', length(vec_mvc), SubjID);
            end
        else
            fprintf('  MVC: champ EMG_max absent dans %s\n', fname_mvc);
        end
    else
        fprintf('  MVC MISSING: %s\n', fname_mvc);
    end

    %% --- Variables par cycle ---
    any_ok = false;
    for iV = 1:N_varCycle
        var_name  = VarsCycle{iV, 1};
        var_path  = VarsCycle{iV, 2};
        fname_pat = VarsCycle{iV, 3};
        field     = VarsCycle{iV, 4};
        fname = fullfile(var_path, sprintf(fname_pat, SubjID));
        if ~exist(fname, 'file')
            fprintf('  MISSING: %s\n', fname); continue;
        end

        tmp      = load(fname, field);
        data_raw = tmp.(field);
        data     = data_raw(1,1);

        for iC = 1:N_ch
            ch  = Channels{iC};
            sig = data.(ch)(:);  % 136x1

            %% DIAGNOSTIC ACTIVITY - affiche les unites (1 seule fois)
            if strcmp(var_name, 'Activity') && ~Activity_debug_done
                fprintf('\n  === DIAGNOSTIC ACTIVITY (%s - %s) ===\n', SubjID, ch);
                fprintf('  min=%.2f  max=%.2f  mean=%.2f  median=%.2f\n', ...
                    min(sig), max(sig), mean(sig), median(sig));
                if strcmp(ch, Channels{end})
                    Activity_debug_done = true;
                    fprintf('  => Valeurs brutes : pas un ratio [0,1].\n');
                    fprintf('  => Normalisation par mediane bin1 appliquee sur la figure.\n\n');
                end
            end

            if length(sig) ~= N_cycles
                fprintf('  %s %s %s: %d valeurs != %d cycles\n', ...
                    SubjID, var_name, ch, length(sig), N_cycles);
                continue;
            end

            for k = 1:N_bins
                ks = bin_edges(k); ke = bin_edges(k+1);
                if ke < ks, continue; end
                raw_val = mean(sig(ks:ke));

                %% NORMALISATION %MVC pour Amplitude uniquement
                if strcmp(var_name, 'Amplitude') && mvc_ok && ...
                        ~isnan(Mat_MVC(iG, iC)) && Mat_MVC(iG, iC) > 0
                    Mat_EMG(iG, k, iC, iV) = (raw_val / Mat_MVC(iG, iC)) * 100;
                else
                    Mat_EMG(iG, k, iC, iV) = raw_val;
                end
            end
        end
        any_ok = true;
    end

    %% --- Variables TFR ---
    for iV = 1:N_varTFR
        iVar      = N_varCycle + iV;
        var_name  = VarsTFR{iV, 1};
        var_path  = VarsTFR{iV, 2};
        fname_pat = VarsTFR{iV, 3};
        field_top = VarsTFR{iV, 4};
        field_in  = VarsTFR{iV, 5};
        fname = fullfile(var_path, sprintf(fname_pat, SubjID));
        if ~exist(fname, 'file')
            fprintf('  MISSING: %s\n', fname); continue;
        end

        tmp      = load(fname, field_top);
        data_top = tmp.(field_top);
        data_mid = data_top(1,1).(field_in);
        data_tfr = data_mid(1,1);
        N_tfr    = size(data_tfr.(Channels{1}), 1);
        fs_tfr   = N_tfr / duree;

        for iC = 1:N_ch
            ch      = Channels{iC};
            sig     = data_tfr.(ch)(:);
            sig_cyc = NaN(N_cycles, 1);
            for cyc = 1:N_cycles
                t_s = t_do(cyc)   - t_start;
                t_e = t_do(cyc+1) - t_start;
                f_s = max(1, round(t_s * fs_tfr) + 1);
                f_e = min(N_tfr, round(t_e * fs_tfr));
                if f_e >= f_s
                    sig_cyc(cyc) = mean(sig(f_s:f_e));
                end
            end
            for k = 1:N_bins
                ks = bin_edges(k); ke = bin_edges(k+1);
                if ke < ks, continue; end
                Mat_EMG(iG, k, iC, iVar) = mean(sig_cyc(ks:ke), 'omitnan');
            end
        end
        any_ok = true;
    end

    if any_ok
        Mat_RPE(iG, :) = rpe_bins';
        valid_subj(iG) = true;
        fprintf('  %s: OK - RPE [%.1f-%.1f]\n', SubjID, ...
            nanmin(rpe_bins), nanmax(rpe_bins));
    end
end

fprintf('\n%d / %d participants traites\n', sum(valid_subj), length(G1_Liszt));

%% -----------------------------------------------------------------------
%  SAVE
%  -----------------------------------------------------------------------
VarNames = [VarsCycle(:,1); VarsTFR(:,1)];
save([PathSave '\EMG_Li_TimeNormalised_G1.mat'], ...
    'Mat_EMG', 'Mat_RPE', 'valid_subj', 'Channels', 'VarNames', ...
    'N_bins', 'G1_Liszt', 'Mat_MVC');
fprintf('Sauvegarde : EMG_Li_TimeNormalised_G1.mat\n');

%% -----------------------------------------------------------------------
%  FIGURES
%  -----------------------------------------------------------------------
T = (1:N_bins) * (100/N_bins);
Vi = find(valid_subj);
[Mu_rpe, SD_rpe] = gStats(Mat_RPE(Vi, :));

MuscleColors = {
    [0.85 0.15 0.15],   % Biceps  rouge
    [0.15 0.45 0.85],   % Triceps bleu
    [0.10 0.70 0.30],   % DeltAnt vert
    [0.80 0.45 0.05],   % DeltMed orange
    [0.55 0.10 0.75],   % SupTrap violet
};

VarLabels = {
    'Amplitude RMS (%MVC)',       % normalise par MVC
    'Activity (u.a., norm. bin1)',% valeurs brutes normalisees par bin1
    'Mobility (ratio)',
    'Sample Entropy',
    'Median Frequency (Hz)',
    'Spectral Entropy',
};

for iV = 1:N_vars
    fig_name = ['Fig_EMG_' VarNames{iV}];
    figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [30 30 1600 700]);

    for iScore = 1:2
        for iC = 1:N_ch
            subplot(2, N_ch, (iScore-1)*N_ch + iC);

            Mat_v = squeeze(Mat_EMG(Vi, :, iC, iV));  % N_valid x N_bins

            %% NORMALISATION ACTIVITY par mediane du bin1 (sur la figure seulement)
            %  Les valeurs brutes restent dans Mat_EMG
            if strcmp(VarNames{iV}, 'Activity')
                ref_bin1 = nanmedian(Mat_v(:, 1));
                if ~isnan(ref_bin1) && ref_bin1 ~= 0
                    Mat_v = Mat_v / ref_bin1;
                end
            end

            %% NETTOYAGE SampleEntropy : Inf->NaN (Mean+Median) + filtre +-3SD (Median)
            if strcmp(VarNames{iV}, 'SampleEntropy')
                Mat_v(isinf(Mat_v)) = NaN;
                if iScore == 2
                    mu_tmp = nanmean(Mat_v(:));
                    sd_tmp = nanstd(Mat_v(:));
                    Mat_v(abs(Mat_v - mu_tmp) > 3 * sd_tmp) = NaN;
                end
            end

            if iScore == 1
                [Mu_v, SD_v] = gStats(Mat_v);
                score_str = 'Mean';
            else
                [Mu_v, SD_v] = gStatsMedian(Mat_v);
                score_str = 'Median';
            end

            plotDual(T, Mu_v, SD_v, Mu_rpe, SD_rpe, ...
                MuscleColors{iC}, VarLabels{iV}, ...
                [score_str ' - ' Channels{iC}]);
        end
    end

    sgtitle([strrep(VarNames{iV}, '_', ' ') ' vs RPE - Temps normalise (%) - G1 n=' num2str(length(Vi))], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, [PathSave '\' fig_name '.png']);
    fprintf('Figure : %s.png\n', fig_name);
end

fprintf('\nTermine.\n');

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

    yyaxis left; hold on
    fill([T fliplr(T)], [Mu_v+SD_v fliplr(Mu_v-SD_v)], col_v, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(T, Mu_v, '-', 'Color', col_v, 'LineWidth', 2.5);
    if ~isempty(trend)
        plot(T(valid_idx), trend, '-', 'Color', col_trend, 'LineWidth', 2);
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

    xlabel('Temps normalise (%)', 'FontSize', 8);
    title(ttl, 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
    xlim([5 105]); xticks(T); grid on; box on; hold off
end