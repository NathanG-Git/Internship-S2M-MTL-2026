%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES 10 VARIABLES GOUBAULT SUR NOUVEAU JEU DE DONNÉES      %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                         %%%%
%%%%                                                                    %%%%
%%%%  10 variables (Goubault 2023, Chord task — Table III) :           %%%%
%%%%   1  Hand      MedianFreq       Acceleration Y                    %%%%
%%%%   2  Hand      Percentile90     Acceleration X                    %%%%
%%%%   3  Hand      PeakPower        Module_Acceleration                %%%%
%%%%   4  Hand      SpectralEntropy  Acceleration Y                    %%%%
%%%%   5  Head      PeakPower        Angular_Velocity X                %%%%
%%%%   6  Forearm   Mean             Acceleration Y                    %%%%
%%%%   7  Forearm   SpectralEntropy  Module_Angular_Velocity            %%%%
%%%%   8  Forearm   PeakPower_Freq   Module_Angular_Velocity            %%%%
%%%%   9  Shoulder  Mean             Angular_Velocity X                %%%%
%%%%  10  Shoulder  Mean             Angular_Velocity Y                %%%%
%%%%                                                                    %%%%
%%%%  Segments (mapping Xsens confirmé) :                              %%%%
%%%%    Head=7, Shoulder=8(RightShoulder), Forearm=10(RightForeArm),  %%%%
%%%%    Hand=11(RightHand)                                              %%%%
%%%%                                                                    %%%%
%%%%  angularVelocity confirmée en rad/s (mean(abs)≈1.16, std≈1.66)    %%%%
%%%%                                                                    %%%%
%%%%  CWT : mêmes paramètres que Calcul_variables.m (jeu 1, validé) :  %%%%
%%%%    Fs_cwt=60Hz (=Fs_acq, pas de rééchantillonnage), cmor8-1,      %%%%
%%%%    fc_morlet=1.0, f=0.05:0.05:15 Hz                               %%%%
%%%%                                                                    %%%%
%%%%  Mean : valeur absolue (cohérent avec Visualize_Workload_..._     %%%%
%%%%  Goubault.m — évite l'annulation de signe entre participants)    %%%%
%%%%                                                                    %%%%
%%%%  Découpage : 10 bins directs sur la durée totale (comme SegMod/   %%%%
%%%%  SegAxis — pas de cycles musicaux pour ce dataset)                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_xsens = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\April28\';
path_save  = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
Fs       = 60;
N_bins   = 10;

% Mapping segment -> index Xsens (confirmé par .label)
iSeg_Head     = 7;
iSeg_Shoulder = 8;   % RightShoulder
iSeg_Forearm  = 10;  % RightForeArm
iSeg_Hand     = 11;  % RightHand

% --- Paramètres CWT (identiques à Calcul_variables.m, jeu 1 validé) ---
Fs_cwt    = 60;     % = Fs_acq, pas de rééchantillonnage
ond       = 'cmor8-1';
fc_morlet = 1.0;
f_cwt     = 0.05:0.05:15;
scale_cwt = Fs_cwt * fc_morlet ./ f_cwt;

fprintf('=== Paramètres ===\n');
fprintf('  Fs=%d Hz | N_bins=%d\n', Fs, N_bins);
fprintf('  CWT : Fs_cwt=%d Hz | bande %.2f-%.0f Hz | %s\n\n', ...
    Fs_cwt, f_cwt(1), f_cwt(end), ond);

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
Results_Goubault = struct();

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
fprintf('=== Calcul des 10 variables Goubault ===\n');

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

            % Exclusion P07 (cohérent avec le reste du pipeline)
            if strcmp(subj, 'P07')
                fprintf('exclu entièrement\n');
                break;
            end

            % --- Extraction signaux bruts des 4 segments utiles ---
            try
                acc_head     = XF.(subj).(kbd).segmentData(iSeg_Head).acceleration;
                acc_forearm  = XF.(subj).(kbd).segmentData(iSeg_Forearm).acceleration;
                acc_hand     = XF.(subj).(kbd).segmentData(iSeg_Hand).acceleration;
                angvel_head     = XF.(subj).(kbd).segmentData(iSeg_Head).angularVelocity;
                angvel_forearm  = XF.(subj).(kbd).segmentData(iSeg_Forearm).angularVelocity;
                angvel_shoulder = XF.(subj).(kbd).segmentData(iSeg_Shoulder).angularVelocity;
            catch ME
                fprintf('Erreur extraction : %s — skip\n', ME.message);
                continue;
            end

            N_frames = size(acc_hand, 1);
            N_frames_min = 9000;  % cohérent avec Compute_IMU_NewDataset.m
            if N_frames < N_frames_min
                fprintf('Session trop courte (%d frames) — exclue\n', N_frames);
                continue;
            end

            bin_size = floor(N_frames / N_bins);

            % --- Modules (norme euclidienne 3 axes) ---
            mod_acc_hand    = sqrt(sum(acc_hand.^2, 2));
            mod_angvel_fore = sqrt(sum(angvel_forearm.^2, 2));

            % --- Initialisation des 10 sorties par bin ---
            V1_MedianFreq_Hand_AccY   = NaN(N_bins, 1);
            V2_Prc90_Hand_AccX        = NaN(N_bins, 1);
            V3_PeakPower_Hand_AccMod  = NaN(N_bins, 1);
            V4_SpecEnt_Hand_AccY      = NaN(N_bins, 1);
            V5_PeakPower_Head_AngVelX = NaN(N_bins, 1);
            V6_Mean_Forearm_AccY      = NaN(N_bins, 1);
            V7_SpecEnt_Forearm_AngVelMod    = NaN(N_bins, 1);
            V8_PeakPowerFreq_Forearm_AngVelMod = NaN(N_bins, 1);
            V9_Mean_Shoulder_AngVelX  = NaN(N_bins, 1);
            V10_Mean_Shoulder_AngVelY = NaN(N_bins, 1);

            for iBin = 1:N_bins
                i0 = (iBin-1)*bin_size + 1;
                i1 = min(iBin*bin_size, N_frames);
                if (i1 - i0 + 1) < Fs * 0.5, continue; end

                % --- Variables 2, 6, 9, 10 : stats simples (pas de CWT) ---
                V2_Prc90_Hand_AccX(iBin)       = prctile(acc_hand(i0:i1, 1), 90);
                V6_Mean_Forearm_AccY(iBin)     = mean(abs(acc_forearm(i0:i1, 2)), 'omitnan');
                V9_Mean_Shoulder_AngVelX(iBin) = mean(abs(angvel_shoulder(i0:i1, 1)), 'omitnan');
                V10_Mean_Shoulder_AngVelY(iBin)= mean(abs(angvel_shoulder(i0:i1, 2)), 'omitnan');

                % --- Variables 1,3,4 : CWT sur Acceleration Hand (Y et Module) ---
                try
                    x_accY = acc_hand(i0:i1, 2);
                    coef_y = cwt(x_accY(:), scale_cwt, ond, 'ExtendSignal', 1);
                    TFR_accY = abs(coef_y);

                    if size(TFR_accY,1) >= length(f_cwt) && ~any(isnan(TFR_accY(:)))
                        mf = Compute_Median_Frequency(TFR_accY, f_cwt);
                        V1_MedianFreq_Hand_AccY(iBin) = mean(mf(~isnan(mf)));

                        se = Compute_Spectral_Entropy(TFR_accY, f_cwt);
                        V4_SpecEnt_Hand_AccY(iBin) = mean(se(~isnan(se)));
                    end

                    % Module Acceleration Hand (somme CWT des 3 axes)
                    x_accMod_x = acc_hand(i0:i1, 1);
                    x_accMod_z = acc_hand(i0:i1, 3);
                    coef_x = cwt(x_accMod_x(:), scale_cwt, ond, 'ExtendSignal', 1);
                    coef_z = cwt(x_accMod_z(:), scale_cwt, ond, 'ExtendSignal', 1);
                    TFR_accMod = abs(coef_x) + abs(coef_y) + abs(coef_z);

                    if size(TFR_accMod,1) >= length(f_cwt) && ~any(isnan(TFR_accMod(:)))
                        mean_t = mean(TFR_accMod, 2);
                        [pk_val, ~] = max(mean_t);
                        V3_PeakPower_Hand_AccMod(iBin) = pk_val;
                    end
                catch
                    % NaN conservé
                end

                % --- Variable 5 : CWT sur Angular_Velocity Head, axe X ---
                try
                    x_angvelX_head = angvel_head(i0:i1, 1);
                    coef = cwt(x_angvelX_head(:), scale_cwt, ond, 'ExtendSignal', 1);
                    TFR = abs(coef);
                    if size(TFR,1) >= length(f_cwt) && ~any(isnan(TFR(:)))
                        mean_t = mean(TFR, 2);
                        [pk_val, ~] = max(mean_t);
                        V5_PeakPower_Head_AngVelX(iBin) = pk_val;
                    end
                catch
                end

                % --- Variables 7, 8 : CWT sur Module Angular_Velocity Forearm ---
                try
                    x1 = angvel_forearm(i0:i1, 1);
                    x2 = angvel_forearm(i0:i1, 2);
                    x3 = angvel_forearm(i0:i1, 3);
                    c1 = cwt(x1(:), scale_cwt, ond, 'ExtendSignal', 1);
                    c2 = cwt(x2(:), scale_cwt, ond, 'ExtendSignal', 1);
                    c3 = cwt(x3(:), scale_cwt, ond, 'ExtendSignal', 1);
                    TFR_mod = abs(c1) + abs(c2) + abs(c3);

                    if size(TFR_mod,1) >= length(f_cwt) && ~any(isnan(TFR_mod(:)))
                        se = Compute_Spectral_Entropy(TFR_mod, f_cwt);
                        V7_SpecEnt_Forearm_AngVelMod(iBin) = mean(se(~isnan(se)));

                        mean_t = mean(TFR_mod, 2);
                        [~, pk_idx] = max(mean_t);
                        V8_PeakPowerFreq_Forearm_AngVelMod(iBin) = f_cwt(pk_idx);
                    end
                catch
                end
            end % iBin

            % --- Stockage ---
            Results_Goubault.(subj).(kbd).Hand_MedianFreq_AccelY            = V1_MedianFreq_Hand_AccY';
            Results_Goubault.(subj).(kbd).Hand_Prc90_AccelX                 = V2_Prc90_Hand_AccX';
            Results_Goubault.(subj).(kbd).Hand_PeakPower_AccelModule        = V3_PeakPower_Hand_AccMod';
            Results_Goubault.(subj).(kbd).Hand_SpectralEntropy_AccelY       = V4_SpecEnt_Hand_AccY';
            Results_Goubault.(subj).(kbd).Head_PeakPower_AngVelX            = V5_PeakPower_Head_AngVelX';
            Results_Goubault.(subj).(kbd).Forearm_Mean_AccelY               = V6_Mean_Forearm_AccY';
            Results_Goubault.(subj).(kbd).Forearm_SpectralEntropy_AngVelMod = V7_SpecEnt_Forearm_AngVelMod';
            Results_Goubault.(subj).(kbd).Forearm_PeakPowerFreq_AngVelMod   = V8_PeakPowerFreq_Forearm_AngVelMod';
            Results_Goubault.(subj).(kbd).Shoulder_Mean_AngVelX             = V9_Mean_Shoulder_AngVelX';
            Results_Goubault.(subj).(kbd).Shoulder_Mean_AngVelY             = V10_Mean_Shoulder_AngVelY';

            fprintf('OK\n');
        end % conditions
    end % sujets

    clear XF tmp;
end % fichiers

%% -----------------------------------------------------------------------
%  RÉSUMÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Résumé ===\n');
subjects_done = fieldnames(Results_Goubault);
N_subj_done   = length(subjects_done);
n_crossover   = 0;
fprintf('  %-6s  %s\n', 'Sujet', 'Conditions');
fprintf('  %s\n', repmat('-', 1, 30));
for iS = 1:N_subj_done
    subj = subjects_done{iS};
    kbds = fieldnames(Results_Goubault.(subj));
    fprintf('  %-6s  %s\n', subj, strjoin(kbds, ' + '));
    if length(kbds) == 2, n_crossover = n_crossover + 1; end
end
fprintf('\n  Total : %d sujets | %d crossover (Norm+CS60)\n', N_subj_done, n_crossover);

var_names_goubault = {
    'Hand_MedianFreq_AccelY', 'Hand_Prc90_AccelX', 'Hand_PeakPower_AccelModule', ...
    'Hand_SpectralEntropy_AccelY', 'Head_PeakPower_AngVelX', 'Forearm_Mean_AccelY', ...
    'Forearm_SpectralEntropy_AngVelMod', 'Forearm_PeakPowerFreq_AngVelMod', ...
    'Shoulder_Mean_AngVelX', 'Shoulder_Mean_AngVelY'
};

fprintf('\n  Vérification valeurs valides :\n');
fprintf('  %-38s  %s\n', 'Variable', 'Valides/Total');
fprintf('  %s\n', repmat('-', 1, 55));
for iV = 1:length(var_names_goubault)
    vn = var_names_goubault{iV};
    n_ok = 0; n_tot = 0;
    for iS = 1:N_subj_done
        subj = subjects_done{iS};
        kbds = fieldnames(Results_Goubault.(subj));
        for iK = 1:length(kbds)
            v = Results_Goubault.(subj).(kbds{iK}).(vn);
            n_ok  = n_ok  + sum(~isnan(v));
            n_tot = n_tot + length(v);
        end
    end
    fprintf('  %-38s  %d/%d\n', vn, n_ok, n_tot);
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'Goubault_Features_NewDataset.mat'), ...
    'Results_Goubault', 'N_bins', 'Fs', 'Fs_cwt', 'f_cwt', '-v7.3');
fprintf('\nSauvegardé : Goubault_Features_NewDataset.mat\n');
fprintf('Structure  : Results_Goubault.(Sujet).(Condition).(Variable) = [1 x %d]\n', N_bins);