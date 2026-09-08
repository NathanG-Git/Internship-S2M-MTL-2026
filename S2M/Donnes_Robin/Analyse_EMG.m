%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES 3 VARIABLES EMG SUR NOUVEAU JEU DE DONNÉES            %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                         %%%%
%%%%                                                                    %%%%
%%%%  Variables calculées par intervalle (10 x 10% de session) :       %%%%
%%%%    1. TFR_MedianFreq_Triceps   (col 5) — CWT Morlet               %%%%
%%%%    2. TFR_SpectralEntropy_Deltoid (col 6) — CWT Morlet            %%%%
%%%%    3. SampleEntropy_Biceps     (col 4) — entropie non-linéaire    %%%%
%%%%                                                                    %%%%
%%%%  Pipeline :                                                        %%%%
%%%%    1. Filtre Butterworth BP 10-400 Hz (Fs=2000)                   %%%%
%%%%    2. Rééchantillonnage à 1000 Hz                                 %%%%
%%%%    3. Découpage en 10 intervalles proportionnels                  %%%%
%%%%    4. CWT Morlet + SampEn par intervalle                          %%%%
%%%%    5. Analyse temporelle et ΔRPE                                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_emg  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_save = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

addpath(path_emg);

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
Fs_acq   = 2000;   % fréquence d'acquisition EMG
Fs_emg   = 1000;   % fréquence après rééchantillonnage (conforme pipeline Robin)
N_bins   = 10;     % intervalles temporels normalisés

% Mapping muscles (ordre des colonnes dans EMG_cut)
% 1=Brachioradialis, 2=Flexor digitorum, 3=Wrist flexor,
% 4=Bicep, 5=Tricep, 6=Deltoid, 7=Trap
iC_Bicep   = 4;
iC_Tricep  = 5;
iC_Deltoid = 6;
Muscle_names = {'Brachioradialis','Flexor digitorum','Wrist flexor',...
                'Bicep','Tricep','Deltoid','Trap'};

% Design expérimental
Subject_Norm_only = {'P01','P06','P09','P11','P13','P16','P19','P22','P23'};
Subject_Both      = {'P02','P03','P04','P05','P07','P08','P10','P12',...
                     'P15','P17','P18','P20','P21','P24'};
Subject_All       = [Subject_Norm_only, Subject_Both];
N_subj = length(Subject_All);

% RPE pré/post — à remplir selon les données disponibles
% Format : RPE_data.(Subject).(Keyboard).pre et .post
% Exemple : RPE_data.P02.Norm.pre = 1; RPE_data.P02.Norm.post = 5;
% Si non disponible laisser NaN
RPE_data = struct();
for iS = 1:N_subj
    RPE_data.(Subject_All{iS}).Norm.pre  = NaN;
    RPE_data.(Subject_All{iS}).Norm.post = NaN;
    if ismember(Subject_All{iS}, Subject_Both)
        RPE_data.(Subject_All{iS}).CS60.pre  = NaN;
        RPE_data.(Subject_All{iS}).CS60.post = NaN;
    end
end
% *** REMPLIR ICI LES VALEURS RPE QUAND DISPONIBLES ***
% Exemple : RPE_data.P02.Norm.pre = 2; RPE_data.P02.Norm.post = 6;

% Paramètres CWT (pour TFR_MedianFreq et TFR_SpectralEntropy)
ond       = 'cmor8-1';
fc_morlet = 1.0;
f_cwt     = 1:1:400;                          % grille 1-400 Hz à 1 Hz (conforme Robin)
scale_cwt = Fs_emg * fc_morlet ./ f_cwt;
Freq_cwt  = f_cwt;

% Bandes fréquentielles pour la fréquence médiane EMG
% Robin utilise 1-400 Hz, on garde cette plage
idx_freq_lo = 1;   % 1 Hz
idx_freq_hi = length(f_cwt);  % 400 Hz

% Paramètres SampleEntropy
SampEn_m = 2;      % dimension d'embedding
SampEn_r = 0.25;   % tolérance (proportion de l'écart-type)

%% -----------------------------------------------------------------------
%  FILTRE BUTTERWORTH 10-400 Hz
%  -----------------------------------------------------------------------
fprintf('Paramètres du filtre Butterworth 10-400 Hz...\n');
[b_filt, a_filt] = butter(2, [10 400]/(Fs_acq/2), 'bandpass');

%% -----------------------------------------------------------------------
%  CHARGEMENT DES DONNÉES EMG
%  -----------------------------------------------------------------------
fprintf('Chargement EMG_cut...\n');
load(fullfile(path_emg, 'EMG_cut'));
fprintf('  OK\n');

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE : CALCUL DES FEATURES PAR SUJET ET CONDITION
%  -----------------------------------------------------------------------
fprintf('\n=== Calcul des features EMG (%d sujets) ===\n', N_subj);

% Structure de résultats : [N_subj x N_bins] pour chaque variable et condition
Results = struct();

for iS = 1:N_subj
    subj = Subject_All{iS};
    fprintf('\n  Sujet %s (%d/%d)...\n', subj, iS, N_subj);

    % Conditions disponibles pour ce sujet
    if ismember(subj, Subject_Both)
        Keyboards = {'Norm','CS60'};
    else
        Keyboards = {'Norm'};
    end

    for iK = 1:length(Keyboards)
        kbd = Keyboards{iK};
        fprintf('    Condition : %s\n', kbd);

        %% --- Étape 1 : Filtre + Rééchantillonnage ---
        emg_raw = EMG_cut.(subj).(kbd);   % [N_frames x 7]
        N_raw   = size(emg_raw, 1);

        % Filtre passe-bande 10-400 Hz
        emg_bp = filtfilt(b_filt, a_filt, emg_raw);

        % Rééchantillonnage 2000 → 1000 Hz par interpolation spline
        t_orig = 1:N_raw;
        t_new  = linspace(1, N_raw, floor(N_raw/2));
        emg_filt = interp1(t_orig, emg_bp, t_new, 'spline');  % [N_new x 7]
        N_new = size(emg_filt, 1);

        fprintf('      Signal : %d frames → %d frames (%.1f s)\n', ...
            N_raw, N_new, N_new/Fs_emg);

        %% --- Étape 2 : Découpage en 10 intervalles proportionnels ---
        bin_size = floor(N_new / N_bins);

        % Initialisation des résultats pour ce sujet/condition
        mf_tricep  = NaN(1, N_bins);   % TFR_MedianFreq Triceps
        se_deltoid = NaN(1, N_bins);   % TFR_SpectralEntropy Deltoid
        sampen_bic = NaN(1, N_bins);   % SampleEntropy Biceps

        for iBin = 1:N_bins
            idx_start = (iBin-1)*bin_size + 1;
            idx_end   = min(iBin*bin_size, N_new);
            seg_len   = idx_end - idx_start + 1;

            if seg_len < Fs_emg * 0.5
                fprintf('      Intervalle %d trop court (%d frames), skip\n', iBin, seg_len);
                continue;
            end

            %% --- Variable 1 : TFR_MedianFreq Triceps (col 5) ---
            x_tri = emg_filt(idx_start:idx_end, iC_Tricep);
            x_tri_z = zscore(x_tri);

            try
                % CWT sur signal z-scoré
                coef = cwt(x_tri_z, scale_cwt, ond, 'ExtendSignal', 1);
                TFR  = abs(coef);  % [N_freq x N_frames]

                if size(TFR,1) >= idx_freq_hi
                    TFR_band = TFR(idx_freq_lo:idx_freq_hi, :);
                    mf = Compute_Median_Frequency(TFR_band, Freq_cwt(idx_freq_lo:idx_freq_hi));
                    mf_tricep(iBin) = mean(mf(~isnan(mf)));
                end
            catch
                % CWT échouée — NaN conservé
            end

            %% --- Variable 2 : TFR_SpectralEntropy Deltoid (col 6) ---
            x_del = emg_filt(idx_start:idx_end, iC_Deltoid);
            x_del_z = zscore(x_del);

            try
                coef = cwt(x_del_z, scale_cwt, ond, 'ExtendSignal', 1);
                TFR  = abs(coef);

                if size(TFR,1) >= idx_freq_hi
                    TFR_band = TFR(idx_freq_lo:idx_freq_hi, :);
                    se = Compute_Spectral_Entropy(TFR_band, Freq_cwt(idx_freq_lo:idx_freq_hi));
                    se_deltoid(iBin) = mean(se(~isnan(se)));
                end
            catch
                % CWT échouée — NaN conservé
            end

            %% --- Variable 3 : SampleEntropy Biceps (col 4) ---
            x_bic = emg_filt(idx_start:idx_end, iC_Bicep);
            
            try
                % Sous-échantillonnage à 10 Hz avant SampEn
                x_sub = resample(x_bic, 1, 100);
                
                % Normalisation z-score APRÈS sous-échantillonnage
                % Standard en biomécanique/EMG : r = 0.2 * std(signal_normalisé)
                % Comme le signal est z-scoré, std=1 → r = 0.2
                x_sub_z = (x_sub - mean(x_sub)) / std(x_sub);
                r_tol = 0.2 * std(x_sub_z);  % = 0.2 par construction
                
                sampen_bic(iBin) = computeSampEn(x_sub_z, SampEn_m, r_tol);
            catch
                % SampEn échouée — NaN conservé
            end
        end

        %% --- Stocker les résultats ---
        Results.(subj).(kbd).TFR_MedianFreq_Triceps    = mf_tricep;
        Results.(subj).(kbd).TFR_SpectralEntropy_Deltoid = se_deltoid;
        Results.(subj).(kbd).SampleEntropy_Biceps      = sampen_bic;
        Results.(subj).(kbd).N_frames_raw  = N_raw;
        Results.(subj).(kbd).N_frames_filt = N_new;
        Results.(subj).(kbd).duration_s    = N_new / Fs_emg;
        Results.(subj).(kbd).RPE_pre       = RPE_data.(subj).(kbd).pre;
        Results.(subj).(kbd).RPE_post      = RPE_data.(subj).(kbd).post;
        Results.(subj).(kbd).DELTA_RPE     = RPE_data.(subj).(kbd).post - RPE_data.(subj).(kbd).pre;

        fprintf('      MedianFreq Triceps   : %.2f → %.2f Hz\n', mf_tricep(1), mf_tricep(end));
        fprintf('      SpectralEnt Deltoid  : %.4f → %.4f\n',    se_deltoid(1), se_deltoid(end));
        fprintf('      SampleEntropy Biceps : %.4f → %.4f\n',    sampen_bic(1), sampen_bic(end));
    end
end

fprintf('\n=== Calcul terminé ===\n');

%% -----------------------------------------------------------------------
%  ANALYSE TEMPORELLE : variable ~ temps * condition + (1|participant)
%  Uniquement sur les 14 sujets avec les deux conditions (crossover)
%  -----------------------------------------------------------------------
fprintf('\n=== Analyse temporelle (LMM) ===\n');

var_names = {'TFR_MedianFreq_Triceps','TFR_SpectralEntropy_Deltoid','SampleEntropy_Biceps'};
var_labels= {'TFR MedianFreq Triceps (Hz)','TFR SpectralEntropy Deltoid','SampleEntropy Biceps'};

% Construction du tableau long pour LMM [N_obs x 5]
% Colonnes : Valeur, Temps (1-10), Condition (0=Norm, 1=CS60), Participant, Variable
t_bins = (1:N_bins)';

for iV = 1:length(var_names)
    vn = var_names{iV};
    fprintf('\n  Variable : %s\n', vn);

    % Construire tableau long (uniquement sujets crossover)
    Y_all = []; T_all = []; C_all = []; S_all = [];

    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        for iK = 1:2
            kbd = {'Norm','CS60'};
            cond_val = iK - 1;  % 0=Norm, 1=CS60
            y = Results.(subj).(kbd{iK}).(vn)(:);
            ok = ~isnan(y);
            Y_all = [Y_all; y(ok)];
            T_all = [T_all; t_bins(ok)];
            C_all = [C_all; cond_val * ones(sum(ok),1)];
            S_all = [S_all; iS * ones(sum(ok),1)];
        end
    end

    % LMM : variable ~ temps + condition + temps*condition + (1|participant)
    T_lmm = table(Y_all, T_all, categorical(C_all), categorical(S_all), ...
        'VariableNames', {vn, 'Temps', 'Condition', 'Participant'});

    try
        formula = [vn ' ~ Temps + Condition + Temps:Condition + (1|Participant)'];
        lme = fitlme(T_lmm, formula);
        [~,~,FE] = fixedEffects(lme, 'DFMethod','satterthwaite');

        fprintf('    %-30s  beta    p\n', 'Effet');
        fprintf('    %s\n', repmat('-',1,50));
        for ir = 1:height(FE)
            fprintf('    %-30s  %+6.4f  %s\n', ...
                FE.Name{ir}, FE.Estimate(ir), p2star(FE.pValue(ir)));
        end

        % Stocker les résultats LMM
        LMM_Results.(vn).lme  = lme;
        LMM_Results.(vn).FE   = FE;
        LMM_Results.(vn).R2m  = lme.Rsquared.Adjusted;

    catch ME
        fprintf('    LMM échoué : %s\n', ME.message);
    end
end

%% -----------------------------------------------------------------------
%  ANALYSE ΔRPE : pente temporelle ~ ΔRPE
%  Si les RPE pré/post sont disponibles
%  -----------------------------------------------------------------------
fprintf('\n=== Analyse ΔRPE ===\n');

has_rpe = false;
for iS = 1:length(Subject_Both)
    subj = Subject_Both{iS};
    if ~isnan(Results.(subj).Norm.DELTA_RPE)
        has_rpe = true; break;
    end
end

if has_rpe
    fprintf('  RPE disponibles — calcul des corrélations pente vs ΔRPE\n');

    for iV = 1:length(var_names)
        vn = var_names{iV};
        slopes = NaN(length(Subject_Both), 2);  % col1=Norm, col2=CS60
        delta_rpe = NaN(length(Subject_Both), 2);

        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            for iK = 1:2
                kbd = {'Norm','CS60'};
                y = Results.(subj).(kbd{iK}).(vn)(:);
                ok = ~isnan(y);
                if sum(ok) >= 3
                    p_fit = polyfit(t_bins(ok), y(ok), 1);
                    slopes(iS,iK) = p_fit(1);
                end
                delta_rpe(iS,iK) = Results.(subj).(kbd{iK}).DELTA_RPE;
            end
        end

        % Corrélation pente vs ΔRPE pour chaque condition
        for iK = 1:2
            kbd_name = {'Norm','CS60'};
            ok = ~isnan(slopes(:,iK)) & ~isnan(delta_rpe(:,iK));
            if sum(ok) >= 5
                [r_val, p_val] = corr(slopes(ok,iK), delta_rpe(ok,iK));
                fprintf('  %s — %s : r=%.3f, p=%.3f %s\n', ...
                    vn, kbd_name{iK}, r_val, p_val, p2star(p_val));
            end
        end
    end
else
    fprintf('  RPE non disponibles — remplir RPE_data dans le script\n');
    fprintf('  Emplacement : ligne ~60 "*** REMPLIR ICI ***"\n');
end

%% -----------------------------------------------------------------------
%  FIGURES
%  -----------------------------------------------------------------------
fprintf('\n=== Génération des figures ===\n');

col_norm = [0.18 0.37 0.64];  % bleu — piano normal
col_ergo = [0.15 0.55 0.35];  % vert — piano CS60
t_pct    = 5:10:100;          % axe x en % de session (milieu de chaque intervalle)

% Figure 1 : Évolution temporelle des 3 variables — Normal vs CS60
fig1 = figure('Name','EMG_Evolution_Temporelle','NumberTitle','off',...
    'Position',[50 50 1100 800],'Color','white');

for iV = 1:3
    vn = var_names{iV};
    subplot(3,1,iV); hold on;

    % Collecter les matrices [N_subj x N_bins] pour chaque condition
    mat_norm = NaN(length(Subject_Both), N_bins);
    mat_cs60 = NaN(length(Subject_Both), N_bins);

    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        mat_norm(iS,:) = Results.(subj).Norm.(vn);
        mat_cs60(iS,:) = Results.(subj).CS60.(vn);
    end

    % Moyennes et IC 95%
    m_n = mean(mat_norm, 1, 'omitnan');
    s_n = std(mat_norm,  0, 1, 'omitnan') / sqrt(sum(~isnan(mat_norm(:,1))));
    m_e = mean(mat_cs60, 1, 'omitnan');
    s_e = std(mat_cs60,  0, 1, 'omitnan') / sqrt(sum(~isnan(mat_cs60(:,1))));

    % Enveloppes IC 95%
    fill([t_pct, fliplr(t_pct)], [m_n+1.96*s_n, fliplr(m_n-1.96*s_n)], ...
        col_norm, 'FaceAlpha',0.15, 'EdgeColor','none');
    fill([t_pct, fliplr(t_pct)], [m_e+1.96*s_e, fliplr(m_e-1.96*s_e)], ...
        col_ergo, 'FaceAlpha',0.15, 'EdgeColor','none');

    % Courbes moyennes
    plot(t_pct, m_n, 'o-', 'Color',col_norm, 'LineWidth',2, 'MarkerSize',5);
    plot(t_pct, m_e, 's-', 'Color',col_ergo, 'LineWidth',2, 'MarkerSize',5);

    % Annoter les intervalles significatifs (interaction temps*condition)
    % si LMM disponible
    if exist('LMM_Results','var') && isfield(LMM_Results, vn)
        FE = LMM_Results.(vn).FE;
        idx_inter = find(contains(FE.Name, 'Temps:Condition'));
        if ~isempty(idx_inter) && FE.pValue(idx_inter) < 0.05
            text(50, max([m_n m_e])*0.95, ...
                sprintf('Interaction temps×condition : p=%.3f %s', ...
                FE.pValue(idx_inter), p2star(FE.pValue(idx_inter))), ...
                'FontSize',8, 'Color',[0.4 0.4 0.4]);
        end
    end

    ylabel(var_labels{iV}, 'FontSize',10);
    if iV == 1
        legend({'Normal ±IC95%','CS60 ±IC95%'}, 'Location','best','FontSize',9);
        title('Évolution temporelle des variables EMG (sujets crossover, N=14)', ...
            'FontWeight','bold','FontSize',11);
    end
    if iV == 3
        xlabel('Pourcentage de session (%)', 'FontSize',10);
    end
    xlim([0 105]); grid on; box on;
end

saveas(fig1, fullfile(path_save,'Fig_EMG_Evolution_Temporelle.png'));
fprintf('  Fig_EMG_Evolution_Temporelle.png sauvegardée\n');

% Figure 2 : Trajectoires individuelles pour chaque variable
fig2 = figure('Name','EMG_Trajectoires_Individuelles','NumberTitle','off',...
    'Position',[50 900 1100 700],'Color','white');

for iV = 1:3
    vn = var_names{iV};

    subplot(2,3,iV); hold on;
    title(['Normal — ' strrep(vn,'_',' ')], 'FontSize',9, 'FontWeight','bold');
    mat_n = NaN(length(Subject_Both), N_bins);
    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        y = Results.(subj).Norm.(vn);
        plot(t_pct, y, '-', 'Color',[col_norm 0.35], 'LineWidth',0.8);
        mat_n(iS,:) = y;
    end
    plot(t_pct, mean(mat_n,1,'omitnan'), '-', 'Color',col_norm, 'LineWidth',2.5);
    xlabel('%session'); ylabel(var_labels{iV},'FontSize',8); grid on; box on;

    subplot(2,3,iV+3); hold on;
    title(['CS60 — ' strrep(vn,'_',' ')], 'FontSize',9, 'FontWeight','bold');
    mat_e = NaN(length(Subject_Both), N_bins);
    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        y = Results.(subj).CS60.(vn);
        plot(t_pct, y, '-', 'Color',[col_ergo 0.35], 'LineWidth',0.8);
        mat_e(iS,:) = y;
    end
    plot(t_pct, mean(mat_e,1,'omitnan'), '-', 'Color',col_ergo, 'LineWidth',2.5);
    xlabel('%session'); ylabel(var_labels{iV},'FontSize',8); grid on; box on;
end

sgtitle('Trajectoires individuelles — Normal vs CS60','FontWeight','bold');
saveas(fig2, fullfile(path_save,'Fig_EMG_Trajectoires_Individuelles.png'));
fprintf('  Fig_EMG_Trajectoires_Individuelles.png sauvegardée\n');

%% -----------------------------------------------------------------------
%  INDEX DE BLOC EMG — Z-SCORE MOYEN
%  Combine les 3 variables en un score unique
%  Signe : MedianFreq ↓ = fatigue → inverser ; SampEn ↑ = fatigue ici → garder
%  -----------------------------------------------------------------------
fprintf('\n=== Calcul de l''Index de Bloc EMG (Z-score moyen) ===\n');

% Collecter toutes les valeurs pour standardisation globale
MF_all  = []; SE_all = []; SP_all = [];
for iS = 1:N_subj
    subj = Subject_All{iS};
    if ismember(subj, Subject_Both)
        Kbds = {'Norm','CS60'};
    else
        Kbds = {'Norm'};
    end
    for iK = 1:length(Kbds)
        MF_all = [MF_all; Results.(subj).(Kbds{iK}).TFR_MedianFreq_Triceps(:)];
        SE_all = [SE_all; Results.(subj).(Kbds{iK}).SampleEntropy_Biceps(:)];
        SP_all = [SP_all; Results.(subj).(Kbds{iK}).TFR_SpectralEntropy_Deltoid(:)];
    end
end

% Paramètres de standardisation globaux
mu_mf  = mean(MF_all,'omitnan'); sd_mf  = std(MF_all,'omitnan');
mu_se  = mean(SE_all,'omitnan'); sd_se  = std(SE_all,'omitnan');
mu_sp  = mean(SP_all,'omitnan'); sd_sp  = std(SP_all,'omitnan');

fprintf('  Standardisation globale :\n');
fprintf('    MedianFreq  Triceps  : mu=%.2f sd=%.2f\n', mu_mf, sd_mf);
fprintf('    SampleEnt   Biceps   : mu=%.4f sd=%.4f\n', mu_se, sd_se);
fprintf('    SpectralEnt Deltoid  : mu=%.4f sd=%.4f\n', mu_sp, sd_sp);

% Calcul de l'index par sujet et condition
% Convention de signe :
%   MedianFreq   ↓ avec fatigue → z-score inversé  (* -1)
%   SpectralEnt  pas de tendance claire → z-score inversé (* -1) par cohérence
%   SampleEntropy ↑ avec fatigue ici → z-score gardé (* +1)
% Index_EMG > 0 = plus de fatigue que la moyenne

for iS = 1:N_subj
    subj = Subject_All{iS};
    if ismember(subj, Subject_Both)
        Kbds = {'Norm','CS60'};
    else
        Kbds = {'Norm'};
    end
    for iK = 1:length(Kbds)
        kbd = Kbds{iK};
        mf = Results.(subj).(kbd).TFR_MedianFreq_Triceps;
        se = Results.(subj).(kbd).SampleEntropy_Biceps;
        sp = Results.(subj).(kbd).TFR_SpectralEntropy_Deltoid;

        z_mf = -(mf - mu_mf) / sd_mf;   % inversé : chute MF = fatigue ↑
        z_se = +(se - mu_se) / sd_se;    % gardé   : hausse SE = fatigue ↑
        z_sp = -(sp - mu_sp) / sd_sp;    % inversé : pas de tendance claire

        Results.(subj).(kbd).Index_EMG = (z_mf + z_se + z_sp) / 3;
    end
end

fprintf('  Index_EMG calculé pour tous les sujets.\n');
fprintf('  Convention : Index_EMG > 0 = fatigue supérieure à la moyenne\n');

%% --- LMM sur l''index de bloc ---
fprintf('\n  LMM : Index_EMG ~ Temps + Condition + Temps:Condition + (1|Participant)\n');

Y_idx = []; T_idx = []; C_idx = []; S_idx = [];
for iS = 1:length(Subject_Both)
    subj = Subject_Both{iS};
    % Exclure P07 Norm (session trop courte ~113s)
    if strcmp(subj,'P07')
        y_n = NaN(N_bins,1);
    else
        y_n = Results.(subj).Norm.Index_EMG(:);
    end
    y_e = Results.(subj).CS60.Index_EMG(:);

    ok_n = ~isnan(y_n); ok_e = ~isnan(y_e);
    Y_idx = [Y_idx; y_n(ok_n); y_e(ok_e)];
    T_idx = [T_idx; t_bins(ok_n); t_bins(ok_e)];
    C_idx = [C_idx; zeros(sum(ok_n),1); ones(sum(ok_e),1)];
    S_idx = [S_idx; iS*ones(sum(ok_n),1); iS*ones(sum(ok_e),1)];
end

T_idx_lmm = table(Y_idx, T_idx, categorical(C_idx), categorical(S_idx), ...
    'VariableNames', {'Index_EMG','Temps','Condition','Participant'});

try
    lme_idx = fitlme(T_idx_lmm, ...
        'Index_EMG ~ Temps + Condition + Temps:Condition + (1|Participant)');
    [~,~,FE_idx] = fixedEffects(lme_idx,'DFMethod','satterthwaite');
    [R2m_idx, R2c_idx] = computeNakagawaR2(lme_idx, T_idx_lmm);

    fprintf('\n  %-30s  beta      p\n', 'Effet');
    fprintf('  %s\n', repmat('-',1,55));
    for ir = 1:height(FE_idx)
        fprintf('  %-30s  %+7.4f   %s\n', ...
            FE_idx.Name{ir}, FE_idx.Estimate(ir), p2star(FE_idx.pValue(ir)));
    end
    fprintf('\n  R²m = %.3f | R²c = %.3f\n', R2m_idx, R2c_idx);

    LMM_Results.Index_EMG.lme  = lme_idx;
    LMM_Results.Index_EMG.FE   = FE_idx;
    LMM_Results.Index_EMG.R2m  = R2m_idx;
    LMM_Results.Index_EMG.R2c  = R2c_idx;

catch ME
    fprintf('  LMM Index échoué : %s\n', ME.message);
end

%% --- Figure Index EMG ---
fig3 = figure('Name','Index_EMG','NumberTitle','off',...
    'Position',[50 50 700 420],'Color','white');
hold on;

mat_idx_norm = NaN(length(Subject_Both), N_bins);
mat_idx_cs60 = NaN(length(Subject_Both), N_bins);
for iS = 1:length(Subject_Both)
    subj = Subject_Both{iS};
    if ~strcmp(subj,'P07')
        mat_idx_norm(iS,:) = Results.(subj).Norm.Index_EMG;
    end
    mat_idx_cs60(iS,:) = Results.(subj).CS60.Index_EMG;
end

m_n = mean(mat_idx_norm,1,'omitnan');
s_n = std(mat_idx_norm,0,1,'omitnan') / sqrt(sum(~isnan(mat_idx_norm(:,1))));
m_e = mean(mat_idx_cs60,1,'omitnan');
s_e = std(mat_idx_cs60,0,1,'omitnan') / sqrt(sum(~isnan(mat_idx_cs60(:,1))));

fill([t_pct fliplr(t_pct)],[m_n+1.96*s_n fliplr(m_n-1.96*s_n)],...
    col_norm,'FaceAlpha',0.15,'EdgeColor','none');
fill([t_pct fliplr(t_pct)],[m_e+1.96*s_e fliplr(m_e-1.96*s_e)],...
    col_ergo,'FaceAlpha',0.15,'EdgeColor','none');
plot(t_pct, m_n,'o-','Color',col_norm,'LineWidth',2,'MarkerSize',5);
plot(t_pct, m_e,'s-','Color',col_ergo,'LineWidth',2,'MarkerSize',5);
yline(0,'--k','LineWidth',0.8);

% Annoter avec résultats LMM
if exist('LMM_Results','var') && isfield(LMM_Results,'Index_EMG')
    FE = LMM_Results.Index_EMG.FE;
    p_temps = FE.pValue(strcmp(FE.Name,'Temps'));
    p_inter = FE.pValue(contains(FE.Name,'Temps:Condition'));
    text(10, max([m_n m_e])*0.85, ...
        sprintf('Temps : p=%.3f%s | Interaction : p=%.3f%s',...
        p_temps,p2star(p_temps),p_inter,p2star(p_inter)),...
        'FontSize',9,'Color',[0.4 0.4 0.4]);
end

xlabel('Pourcentage de session (%)','FontSize',11);
ylabel('Index de Bloc EMG (u.a.)','FontSize',11);
title('Index de Bloc EMG — Piano Normal vs CS60','FontWeight','bold');
legend({'Normal ±IC95%','CS60 ±IC95%'},'Location','northwest','FontSize',10);
grid on; box on;
saveas(fig3, fullfile(path_save,'Fig_Index_EMG.png'));
fprintf('\n  Fig_Index_EMG.png sauvegardée\n');
save(fullfile(path_save,'EMG_Features_NewDataset.mat'), 'Results', 'Subject_All',...
    'Subject_Both', 'Subject_Norm_only', 'var_names', 'var_labels');
fprintf('\nRésultats sauvegardés : EMG_Features_NewDataset.mat\n');

if exist('LMM_Results','var')
    save(fullfile(path_save,'LMM_Results_NewDataset.mat'), 'LMM_Results');
    fprintf('LMM sauvegardés : LMM_Results_NewDataset.mat\n');
end

fprintf('\n========================================================\n');
fprintf('  SCRIPT TERMINÉ\n');
fprintf('========================================================\n');

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------
function se = computeSampEn(x, m, r_tol)
% Sample Entropy
% x     : signal déjà normalisé (z-score)
% m     : dimension d'embedding (typiquement 2)
% r_tol : tolérance absolue (typiquement 0.2 sur signal z-scoré)
    x = x(:); N = length(x);
    A = 0; B = 0;
    for i = 1:N-m
        for j = i+1:N-m
            if max(abs(x(i:i+m-1) - x(j:j+m-1))) < r_tol
                B = B + 1;
                if abs(x(i+m) - x(j+m)) < r_tol
                    A = A + 1;
                end
            end
        end
    end
    if B == 0
        se = Inf;
    else
        se = -log(A / B);
    end
end

function sig = p2star(p)
    if isnan(p),      sig = '';
    elseif p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,             sig = '';
    end
end

function [R2m, R2c] = computeNakagawaR2(lme, T)
    fe   = fixedEffects(lme);
    X    = designMatrix(lme, 'Fixed');
    varF = var(X * fe, 0, 'omitnan');
    try
        varR = 0;
        for ir = 1:length(lme.GroupingVariableNames)
            vc = lme.covarianceParameters{ir};
            varR = varR + trace(vc);
        end
    catch
        varR = NaN;
    end
    sigma2 = lme.MSE;
    if isnan(varR)
        R2m = varF/(varF+sigma2); R2c = NaN;
    else
        denom = varF + varR + sigma2;
        if denom <= 0 || isnan(denom)
            R2m = NaN; R2c = NaN;
        else
            R2m = varF / denom;
            R2c = (varF + varR) / denom;
        end
    end
end