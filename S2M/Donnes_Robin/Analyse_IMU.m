%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES 6 VARIABLES IMU SUR NOUVEAU JEU DE DONNÉES            %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                         %%%%
%%%%                                                                    %%%%
%%%%  BLOC SEGMOD (2 variables) — statistique brute, pas de CWT :      %%%%
%%%%    1. Accel_Mod_Head   — moyenne du module d'accélération tête    %%%%
%%%%    2. Accel_Mod_Hand   — moyenne du module d'accélération main    %%%%
%%%%                                                                    %%%%
%%%%  BLOC GOUBAULT (4 variables) — CWT Morlet cmor8-1 à 60 Hz :      %%%%
%%%%    3. MedianFreq_Accel_Y_Hand       — axe Y uniquement            %%%%
%%%%    4. PeakPower_Accel_Mod_Hand      — somme CWT 3 axes            %%%%
%%%%    5. PeakPower_AngVel_X_Head       — axe X uniquement            %%%%
%%%%    6. SpectralEntropy_AngVel_Mod_Forearm — somme CWT 3 axes       %%%%
%%%%                                                                    %%%%
%%%%  Pipeline CWT conforme à Goubault et al. (2021, 2023)             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_xsens = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\April28\';
path_save  = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

addpath(path_data);  % Compute_Median_Frequency.m et Compute_Spectral_Entropy.m

%% -----------------------------------------------------------------------
%  PARAMÈTRES CWT — conforme pipeline Goubault
%  -----------------------------------------------------------------------
Fs        = 60;
ond       = 'cmor8-1';
fc_morlet = 1.0;
f         = 0.05:0.05:15;
scale_cwt = Fs * fc_morlet ./ f;
N_bins    = 10;

% Indices segments (ordre labels XSens)
% 1=Pelvis, 2=L5, 3=L3, 4=T12, 5=T8, 6=Neck, 7=Head,
% 8=RightShoulder, 9=RightUpperArm, 10=RightForeArm, 11=RightHand
% Indices sensorData (free acceleration — sans gravité)
iSens_Head    = 3;
iSens_Hand    = 7;

% Indices segmentData (angular velocity)
iSeg_Head    = 7;
iSeg_Forearm = 10;

fprintf('=== Paramètres ===\n');
fprintf('  CWT : %s | Fs=%d Hz | f=%.2f-%.2f Hz\n', ond, Fs, f(1), f(end));
fprintf('  Intervalles : %d x %.0f%%\n', N_bins, 100/N_bins);
fprintf('  sensorData : Head=%d | Hand=%d | segmentData : Forearm=%d\n', ...
    3, 7, iSeg_Forearm);

%% -----------------------------------------------------------------------
%  FICHIERS XSENS
%  -----------------------------------------------------------------------
xsens_files = {
    fullfile(path_xsens, 'XSens_Fatigue_1.mat');
    fullfile(path_xsens, 'XSens_Fatigue_2.mat');
    fullfile(path_xsens, 'XSens_Fatigue_3.mat');
    fullfile(path_xsens, 'XSens_Fatigue_4.mat');
    fullfile(path_xsens, 'XSens_Fatigue_5.mat');
};

%% -----------------------------------------------------------------------
%  STRUCTURE DE RÉSULTATS
%  -----------------------------------------------------------------------
Results_IMU = struct();

% Noms des variables par bloc
var_SegMod   = {'Accel_Mod_Head', 'Accel_Mod_Hand'};
var_Goubault = {'MedianFreq_Accel_Y_Hand', 'PeakPower_Accel_Mod_Hand', ...
                'PeakPower_AngVel_X_Head', 'SpectralEntropy_AngVel_Mod_Forearm'};
var_all      = [var_SegMod, var_Goubault];

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
fprintf('=== Calcul des features IMU ===\n');

for iFile = 1:length(xsens_files)
    fpath = xsens_files{iFile};
    if ~exist(fpath, 'file')
        fprintf('  Fichier introuvable : %s — skip\n', fpath);
        continue;
    end

    fprintf('\nChargement : %s\n', fpath);
    tmp = load(fpath);
    XF  = tmp.XSens_Fatigue;
    subjects = fieldnames(XF);

    for iS = 1:length(subjects)
        subj       = subjects{iS};
        conditions = fieldnames(XF.(subj));
        fprintf('  %s — %s\n', subj, strjoin(conditions, ' + '));

        for iK = 1:length(conditions)
            kbd = conditions{iK};
            fprintf('    [%s] ', kbd);

            try
                acc_head    = XF.(subj).(kbd).sensorData(3).sensorFreeAcceleration;
                acc_hand    = XF.(subj).(kbd).sensorData(7).sensorFreeAcceleration;
                angvel_head = XF.(subj).(kbd).segmentData(iSeg_Head).angularVelocity;
                angvel_fore = XF.(subj).(kbd).segmentData(iSeg_Forearm).angularVelocity;
            catch ME
                fprintf('Erreur extraction : %s — skip\n', ME.message);
                continue;
            end

            N_frames = size(acc_head, 1);

            % Exclusion sessions trop courtes (< 50% de la durée médiane ~18400 frames)
            % P07 Norm identifié comme session tronquée (6784 frames = ~2 min au lieu de ~5 min)
            % P07 exclu entièrement (Norm + CS60) car données incomplètes sur une condition
            if strcmp(subj, 'P07')
                fprintf('      P07 exclu entièrement (session Norm tronquée)\n');
                break;  % sort de la boucle conditions pour ce sujet
            end

            N_frames_min = 9000;
            if N_frames < N_frames_min
                fprintf('      Session trop courte (%d frames < %d) — exclue\n', ...
                    N_frames, N_frames_min);
                continue;
            end

            bin_size = floor(N_frames / N_bins);

            % Modules pré-calculés (SegMod)
            mod_acc_head = sqrt(sum(acc_head.^2, 2));
            mod_acc_hand = sqrt(sum(acc_hand.^2, 2));

            % Initialisation
            Accel_Mod_Head                     = NaN(1, N_bins);
            Accel_Mod_Hand                     = NaN(1, N_bins);
            MedianFreq_Accel_Y_Hand            = NaN(1, N_bins);
            PeakPower_Accel_Mod_Hand           = NaN(1, N_bins);
            PeakPower_AngVel_X_Head            = NaN(1, N_bins);
            SpectralEntropy_AngVel_Mod_Forearm = NaN(1, N_bins);

            for iBin = 1:N_bins
                i0 = (iBin-1)*bin_size + 1;
                i1 = min(iBin*bin_size, N_frames);

                if (i1 - i0 + 1) < Fs * 0.5, continue; end

                % ---------------------------------------------------
                % BLOC SEGMOD — moyenne du module brut
                % ---------------------------------------------------
                Accel_Mod_Head(iBin) = mean(mod_acc_head(i0:i1), 'omitnan');
                Accel_Mod_Hand(iBin) = mean(mod_acc_hand(i0:i1), 'omitnan');

                % ---------------------------------------------------
                % BLOC GOUBAULT — CWT Morlet
                % ---------------------------------------------------

                % 1. MedianFreq Accel_Y Hand — axe Y uniquement
                try
                    x = acc_hand(i0:i1, 2); x = zscore(x(:));
                    TFR = abs(cwt(x, scale_cwt, ond, 'ExtendSignal', 1));
                    mf  = Compute_Median_Frequency(TFR, f);
                    MedianFreq_Accel_Y_Hand(iBin) = mean(mf(~isnan(mf)));
                catch; end

                % 2. PeakPower Accel Mod Hand — somme CWT 3 axes
                try
                    ax = acc_hand(i0:i1, 1); ax = zscore(ax(:));
                    ay = acc_hand(i0:i1, 2); ay = zscore(ay(:));
                    az = acc_hand(i0:i1, 3); az = zscore(az(:));
                    TFR = abs(cwt(ax, scale_cwt, ond, 'ExtendSignal', 1)) + ...
                          abs(cwt(ay, scale_cwt, ond, 'ExtendSignal', 1)) + ...
                          abs(cwt(az, scale_cwt, ond, 'ExtendSignal', 1));
                    PeakPower_Accel_Mod_Hand(iBin) = max(mean(TFR, 2));
                catch; end

                % 3. PeakPower AngVel_X Head — axe X uniquement
                try
                    x = angvel_head(i0:i1, 1); x = zscore(x(:));
                    TFR = abs(cwt(x, scale_cwt, ond, 'ExtendSignal', 1));
                    PeakPower_AngVel_X_Head(iBin) = max(mean(TFR, 2));
                catch; end

                % 4. SpectralEntropy AngVel Mod Forearm — somme CWT 3 axes
                try
                    wx = angvel_fore(i0:i1, 1); wx = zscore(wx(:));
                    wy = angvel_fore(i0:i1, 2); wy = zscore(wy(:));
                    wz = angvel_fore(i0:i1, 3); wz = zscore(wz(:));
                    TFR = abs(cwt(wx, scale_cwt, ond, 'ExtendSignal', 1)) + ...
                          abs(cwt(wy, scale_cwt, ond, 'ExtendSignal', 1)) + ...
                          abs(cwt(wz, scale_cwt, ond, 'ExtendSignal', 1));
                    se  = Compute_Spectral_Entropy(TFR, f);
                    SpectralEntropy_AngVel_Mod_Forearm(iBin) = mean(se(~isnan(se)));
                catch; end

            end % iBin

            % Stockage
            Results_IMU.(subj).(kbd).Accel_Mod_Head                     = Accel_Mod_Head;
            Results_IMU.(subj).(kbd).Accel_Mod_Hand                     = Accel_Mod_Hand;
            Results_IMU.(subj).(kbd).MedianFreq_Accel_Y_Hand            = MedianFreq_Accel_Y_Hand;
            Results_IMU.(subj).(kbd).PeakPower_Accel_Mod_Hand           = PeakPower_Accel_Mod_Hand;
            Results_IMU.(subj).(kbd).PeakPower_AngVel_X_Head            = PeakPower_AngVel_X_Head;
            Results_IMU.(subj).(kbd).SpectralEntropy_AngVel_Mod_Forearm = SpectralEntropy_AngVel_Mod_Forearm;

            fprintf('OK\n');

        end % conditions
    end % sujets

    clear XF tmp;
end % fichiers

%% -----------------------------------------------------------------------
%  RÉSUMÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Résumé ===\n');
subjects_done = fieldnames(Results_IMU);
N_subj_done   = length(subjects_done);
n_crossover   = 0;
fprintf('  %-6s  %s\n', 'Sujet', 'Conditions');
fprintf('  %s\n', repmat('-', 1, 30));
for iS = 1:N_subj_done
    subj = subjects_done{iS};
    kbds = fieldnames(Results_IMU.(subj));
    fprintf('  %-6s  %s\n', subj, strjoin(kbds, ' + '));
    if length(kbds) == 2, n_crossover = n_crossover + 1; end
end
fprintf('\n  Total : %d sujets | %d crossover (Norm+CS60)\n', N_subj_done, n_crossover);

% Vérification NaN
fprintf('\n  Vérification valeurs valides :\n');
fprintf('  %-45s  %s\n', 'Variable', 'Valides/Total');
fprintf('  %s\n', repmat('-', 1, 65));
for iV = 1:length(var_all)
    vn = var_all{iV};
    bloc = 'SegMod  ';
    if ismember(vn, var_Goubault), bloc = 'Goubault'; end
    n_ok = 0; n_tot = 0;
    for iS = 1:N_subj_done
        subj = subjects_done{iS};
        kbds = fieldnames(Results_IMU.(subj));
        for iK = 1:length(kbds)
            v = Results_IMU.(subj).(kbds{iK}).(vn);
            n_ok  = n_ok  + sum(~isnan(v));
            n_tot = n_tot + length(v);
        end
    end
    fprintf('  [%s] %-40s  %d/%d\n', bloc, vn, n_ok, n_tot);
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'IMU_Features_NewDataset.mat'), ...
    'Results_IMU', 'var_SegMod', 'var_Goubault', 'var_all', 'N_bins', 'Fs', '-v7.3');
fprintf('\nSauvegardé : IMU_Features_NewDataset.mat\n');
fprintf('Structure  : Results_IMU.(Sujet).(Condition).(Variable) = [1 x %d]\n', N_bins);