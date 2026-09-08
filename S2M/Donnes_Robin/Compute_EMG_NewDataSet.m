%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES VARIABLES EMG SUR NOUVEAU JEU DE DONNÉES               %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                         %%%%
%%%%                                                                    %%%%
%%%%  Source : EMG_cut.(Sujet).(Condition) = [N_frames x 7]            %%%%
%%%%  Fs_EMG ≈ 2000 Hz (confirmé : 615903 frames / ~300s)              %%%%
%%%%                                                                    %%%%
%%%%  Mapping colonnes (confirmé séance précédente) :                  %%%%
%%%%    1=Brachioradialis 2=Flexor digitorum 3=Wrist flexor            %%%%
%%%%    4=Bicep 5=Tricep 6=Deltoid 7=Trap                              %%%%
%%%%                                                                    %%%%
%%%%  Mapping avec canaux jeu 1 :                                      %%%%
%%%%    Biceps  -> col 4 (Bicep)                                       %%%%
%%%%    Triceps -> col 5 (Tricep)                                      %%%%
%%%%    DeltAnt -> col 6 (Deltoid)                                     %%%%
%%%%    DeltMed -> col 6 (Deltoid) — dupliqué (un seul capteur deltoïde)%%%%
%%%%    SupTrap -> col 7 (Trap)                                        %%%%
%%%%                                                                    %%%%
%%%%  5 FEATURES (identiques au jeu 1, hors Amplitude exclue) :         %%%%
%%%%    Activity         : variance de l'enveloppe EMG (par bin)       %%%%
%%%%    Mobility         : var(dérivée enveloppe) / var(enveloppe)     %%%%
%%%%    SampleEntropy    : resample 100Hz -> zscore -> SampEn(m=2,r=0.2)%%%%
%%%%    TFR_MedianFreq   : CWT Morlet, fréquence médiane du spectre    %%%%
%%%%    TFR_SpectralEntropy : CWT Morlet, entropie spectrale           %%%%
%%%%                                                                    %%%%
%%%%  PARAMÈTRES CWT — :                                               %%%%
%%%%    Bande EMG typique 20-450 Hz (cf. filtre passe-bande jeu 1      %%%%
%%%%    10-400 Hz, notch 60/120/180 Hz). Fs_cwt = 1000 Hz proposé      %%%%
%%%%    (sous-échantillonnage depuis ~2000 Hz, cohérent avec le        %%%%
%%%%    pipeline jeu 1 normalisé à 1024 Hz selon CheckSegmentation).   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_save = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
Fs_emg_raw = 2000;   % Fs brut EMG_cut (confirmé : 615903/300s ≈ 2053 Hz, retenu 2000)
N_bins     = 10;

% --- Mapping canaux EMG_cut (colonnes 1-based) ---
col_brachio  = 1; col_flexdig = 2; col_wristflex = 3;
col_bicep    = 4; col_tricep  = 5; col_delt      = 6; col_trap = 7;

% --- Mapping vers les 5 canaux "classiques" du jeu 1 ---
% DeltAnt et DeltMed pointent tous les deux vers le canal Deltoid unique
Channels      = {'Biceps','Triceps','DeltAnt','DeltMed','SupTrap'};
Channels_cols = [col_bicep, col_tricep, col_delt, col_delt, col_trap];
N_ch          = length(Channels);

% --- Filtre passe-bande EMG (cohérent jeu 1 : 10-400 Hz Butterworth 2nd ordre) ---
bp_low  = 10; bp_high = 400;
[b_bp, a_bp] = butter(2, [bp_low bp_high]/(Fs_emg_raw/2), 'bandpass');

% --- Filtres coupe-bande (notch) 60/120/180 Hz ---
notch_freqs = [60, 120, 180];
notch_bw    = 2;  % largeur de bande +-1Hz autour de chaque fréquence

% --- Filtre passe-bas pour l'enveloppe (redressée) ---
env_lowpass = 9;  % Hz, cohérent avec normalisation EMG documentée (9 Hz)
[b_env, a_env] = butter(2, env_lowpass/(Fs_emg_raw/2), 'low');

% --- Paramètres SampleEntropy ---
% resample d'abord (2000Hz -> 100Hz, facteur 20), puis zscore, puis SampEn
Fs_sampen_target = 100;
resample_factor  = round(Fs_emg_raw / Fs_sampen_target);
SampEn_m = 2;
SampEn_r_factor = 0.2;  % r_tol = 0.2 * std(signal z-scoré) = 0.2 par construction

% --- Paramètres CWT (NON CONFIRMÉS — point d'ajustement) ---
Fs_cwt    = 1000;          % sous-échantillonnage cible avant CWT
ond       = 'cmor8-1';
fc_morlet = 1.0;
f_cwt     = 2:2:450;       % grille de fréquences EMG, 2-450 Hz par pas de 2Hz
scale_cwt = Fs_cwt * fc_morlet ./ f_cwt;

fprintf('=== Paramètres ===\n');
fprintf('  Fs_emg_raw=%d Hz | N_bins=%d\n', Fs_emg_raw, N_bins);
fprintf('  Passe-bande : %d-%d Hz | Notch : %s Hz\n', bp_low, bp_high, mat2str(notch_freqs));
fprintf('  SampEn : resample %dHz -> %dHz (facteur %d) | m=%d, r=%.1f*std\n', ...
    Fs_emg_raw, Fs_sampen_target, resample_factor, SampEn_m, SampEn_r_factor);
fprintf('  CWT : Fs_cwt=%d Hz | bande %.0f-%.0f Hz | %s\n\n', ...
    Fs_cwt, f_cwt(1), f_cwt(end), ond);

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
fprintf('Chargement EMG_cut...\n');
load(fullfile(path_data, 'EMG_cut.mat'));   % -> EMG_cut

subjects = fieldnames(EMG_cut);
fprintf('  %d sujets trouvés : %s\n\n', length(subjects), strjoin(subjects, ', '));

%% -----------------------------------------------------------------------
%  STRUCTURE DE RÉSULTATS
%  -----------------------------------------------------------------------
Results_EMG = struct();

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
fprintf('=== Calcul des features EMG ===\n');

for iS = 1:length(subjects)
    subj = subjects{iS};

    % Exclusion P07 (cohérent avec le reste du pipeline)
    if strcmp(subj, 'P07')
        fprintf('  %s : exclu entièrement\n', subj);
        continue;
    end

    conditions = fieldnames(EMG_cut.(subj));
    fprintf('  %s — %s\n', subj, strjoin(conditions, ' + '));

    for iK = 1:length(conditions)
        kbd = conditions{iK};
        fprintf('    [%s] ', kbd);

        sig_raw = EMG_cut.(subj).(kbd);   % [N_frames x 7]
        N_frames = size(sig_raw, 1);

        N_frames_min = Fs_emg_raw * 30;  % au moins 30s de signal
        if N_frames < N_frames_min
            fprintf('Session trop courte (%d frames) — exclue\n', N_frames);
            continue;
        end

        bin_size = floor(N_frames / N_bins);

        % --- Filtrage : passe-bande + notch x3, par canal utilisé ---
        sig_filt = NaN(N_frames, N_ch);
        for iC = 1:N_ch
            x = sig_raw(:, Channels_cols(iC));
            x = filtfilt(b_bp, a_bp, x);
            for iN = 1:length(notch_freqs)
                fn = notch_freqs(iN);
                [b_n, a_n] = butter(2, [fn-notch_bw, fn+notch_bw]/(Fs_emg_raw/2), 'stop');
                x = filtfilt(b_n, a_n, x);
            end
            sig_filt(:, iC) = x;
        end

        % --- Enveloppe : redressement (abs) + passe-bas ---
        sig_env = NaN(N_frames, N_ch);
        for iC = 1:N_ch
            sig_env(:, iC) = filtfilt(b_env, a_env, abs(sig_filt(:, iC)));
        end

        % --- Initialisation des sorties par bin ---
        Activity_bins      = NaN(N_bins, N_ch);
        Mobility_bins       = NaN(N_bins, N_ch);
        SampleEntropy_bins  = NaN(N_bins, N_ch);
        TFR_MedianFreq_bins = NaN(N_bins, N_ch);
        TFR_SpecEnt_bins    = NaN(N_bins, N_ch);

        for iBin = 1:N_bins
            i0 = (iBin-1)*bin_size + 1;
            i1 = min(iBin*bin_size, N_frames);
            if (i1 - i0 + 1) < Fs_emg_raw * 1, continue; end

            for iC = 1:N_ch
                x_env = sig_env(i0:i1, iC);

                % --- Activity : variance de l'enveloppe ---
                Activity_bins(iBin, iC) = var(x_env, 'omitnan');

                % --- Mobility : var(dérivée) / var(enveloppe) ---
                dx = diff(x_env);
                v_env = var(x_env, 'omitnan');
                if v_env > 0
                    Mobility_bins(iBin, iC) = var(dx, 'omitnan') / v_env;
                end

                % --- SampleEntropy : resample -> zscore -> SampEn ---
                try
                    x_sub = resample(x_env, 1, resample_factor);
                    x_sub_z = (x_sub - mean(x_sub)) / std(x_sub);
                    r_tol = SampEn_r_factor * std(x_sub_z);  % = 0.2 par construction
                    SampleEntropy_bins(iBin, iC) = computeSampEn(x_sub_z, SampEn_m, r_tol);
                catch
                    % NaN conservé
                end

                % --- TFR : CWT Morlet sur le signal filtré (pas l'enveloppe) ---
                try
                    x_raw_bin = sig_filt(i0:i1, iC);
                    % Sous-échantillonnage vers Fs_cwt avant CWT (calcul plus rapide)
                    ds_factor = round(Fs_emg_raw / Fs_cwt);
                    x_cwt_in  = resample(x_raw_bin, 1, ds_factor);

                    coef = cwt(x_cwt_in(:), scale_cwt, ond, 'ExtendSignal', 1);
                    TFR  = abs(coef);

                    if any(isnan(TFR(:))), continue; end

                    mf = Compute_Median_Frequency(TFR, f_cwt);
                    TFR_MedianFreq_bins(iBin, iC) = mean(mf(~isnan(mf)));

                    se = Compute_Spectral_Entropy(TFR, f_cwt);
                    TFR_SpecEnt_bins(iBin, iC) = mean(se(~isnan(se)));
                catch
                    % NaN conservé
                end
            end
        end % iBin

        % --- Stockage ---
        for iC = 1:N_ch
            ch = Channels{iC};
            Results_EMG.(subj).(kbd).(['Activity_' ch])            = Activity_bins(:, iC)';
            Results_EMG.(subj).(kbd).(['Mobility_' ch])             = Mobility_bins(:, iC)';
            Results_EMG.(subj).(kbd).(['SampleEntropy_' ch])        = SampleEntropy_bins(:, iC)';
            Results_EMG.(subj).(kbd).(['TFR_MedianFreq_' ch])       = TFR_MedianFreq_bins(:, iC)';
            Results_EMG.(subj).(kbd).(['TFR_SpectralEntropy_' ch])  = TFR_SpecEnt_bins(:, iC)';
        end

        fprintf('OK\n');
    end % conditions
end % sujets

%% -----------------------------------------------------------------------
%  RÉSUMÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Résumé ===\n');
subjects_done = fieldnames(Results_EMG);
N_subj_done   = length(subjects_done);
n_crossover   = 0;
fprintf('  %-6s  %s\n', 'Sujet', 'Conditions');
fprintf('  %s\n', repmat('-', 1, 30));
for iS = 1:N_subj_done
    subj = subjects_done{iS};
    kbds = fieldnames(Results_EMG.(subj));
    fprintf('  %-6s  %s\n', subj, strjoin(kbds, ' + '));
    if length(kbds) == 2, n_crossover = n_crossover + 1; end
end
fprintf('\n  Total : %d sujets | %d crossover (Norm+CS60)\n', N_subj_done, n_crossover);

fprintf('\n  Vérification valeurs valides :\n');
feat_names = {'Activity','Mobility','SampleEntropy','TFR_MedianFreq','TFR_SpectralEntropy'};
fprintf('  %-25s  %s\n', 'Variable', 'Valides/Total');
fprintf('  %s\n', repmat('-', 1, 45));
for iF = 1:length(feat_names)
    for iC = 1:N_ch
        vn = [feat_names{iF} '_' Channels{iC}];
        n_ok = 0; n_tot = 0;
        for iS = 1:N_subj_done
            subj = subjects_done{iS};
            kbds = fieldnames(Results_EMG.(subj));
            for iK = 1:length(kbds)
                v = Results_EMG.(subj).(kbds{iK}).(vn);
                n_ok  = n_ok  + sum(~isnan(v));
                n_tot = n_tot + length(v);
            end
        end
        fprintf('  %-25s  %d/%d\n', vn, n_ok, n_tot);
    end
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'EMG_Features_NewDataset_25vars.mat'), ...
    'Results_EMG', 'Channels', 'Channels_cols', 'N_bins', 'Fs_emg_raw', ...
    'Fs_cwt', 'f_cwt', '-v7.3');
fprintf('\nSauvegardé : EMG_Features_NewDataset_25vars.mat\n');
fprintf('Structure  : Results_EMG.(Sujet).(Condition).(Variable) = [1 x %d]\n', N_bins);

%% =========================================================================
%% FONCTIONS LOCALES
%% =========================================================================

function se = computeSampEn(x, m, r_tol)
% Sample Entropy — algorithme direct, O(N^2)
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