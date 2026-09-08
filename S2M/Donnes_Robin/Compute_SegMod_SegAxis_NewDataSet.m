%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES VARIABLES SegMod + SegAxis SUR NOUVEAU JEU DE DONNÉES %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                         %%%%
%%%%                                                                    %%%%
%%%%  BLOC SegMod (16 variables) — moyenne du module par bin :        %%%%
%%%%    Accel_L5, Accel_T8, Accel_Head, Accel_Shoulder, Accel_Arm,    %%%%
%%%%    Accel_Forearm, Accel_Hand, Total_Accel (somme des 7 segments) %%%%
%%%%    + idem pour Jerk (8 variables)                                %%%%
%%%%                                                                    %%%%
%%%%  BLOC SegAxis (42 variables) — moyenne du signal par axe par bin :%%%%
%%%%    7 segments x 3 axes (X,Y,Z) x {Accel, Jerk} = 42               %%%%
%%%%                                                                    %%%%
%%%%  Segments (mapping Xsens 23-seg confirmé par label) :             %%%%
%%%%    L5=2, T8=5, Head=7, Shoulder=8(RightShoulder),                %%%%
%%%%    Arm=9(RightUpperArm), Forearm=10(RightForeArm), Hand=11(RightHand)%%%%
%%%%                                                                    %%%%
%%%%  Jerk = diff(Accel,1,1) * Fs  (identique à Calcul_variables.m)    %%%%
%%%%  Module = norme euclidienne des 3 axes (sqrt(x^2+y^2+z^2))        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_xsens = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\April28\';
path_save  = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
Fs       = 60;
N_bins   = 10;

% Mapping segment -> index Xsens (confirmé par .label, ordre 23-segments)
Seg_names = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Seg_idx   = [2, 5, 7, 8, 9, 10, 11];
N_seg     = length(Seg_names);
Axes_lbl  = {'X','Y','Z'};

fprintf('=== Paramètres ===\n');
fprintf('  Fs=%d Hz | N_bins=%d\n', Fs, N_bins);
fprintf('  Segments : %s\n', strjoin(Seg_names, ', '));
fprintf('  Indices  : %s\n\n', mat2str(Seg_idx));

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
Results_SegMod  = struct();
Results_SegAxis = struct();

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
fprintf('=== Calcul des features SegMod + SegAxis ===\n');

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

            % Exclusion P07 (session Norm tronquée — identique à Compute_IMU_NewDataset.m)
            if strcmp(subj, 'P07')
                fprintf('exclu entièrement (session Norm tronquée)\n');
                break;
            end

            % --- Extraction acceleration brute des 7 segments ---
            try
                acc_seg = cell(1, N_seg);
                for iSeg = 1:N_seg
                    acc_seg{iSeg} = XF.(subj).(kbd).segmentData(Seg_idx(iSeg)).acceleration;
                end
            catch ME
                fprintf('Erreur extraction : %s — skip\n', ME.message);
                continue;
            end

            N_frames = size(acc_seg{1}, 1);

            N_frames_min = 9000;
            if N_frames < N_frames_min
                fprintf('Session trop courte (%d frames < %d) — exclue\n', ...
                    N_frames, N_frames_min);
                continue;
            end

            bin_size = floor(N_frames / N_bins);

            % --- Jerk = diff(Accel,1,1) * Fs, complété par une ligne de zéros ---
            jerk_seg = cell(1, N_seg);
            for iSeg = 1:N_seg
                jerk_seg{iSeg} = [diff(acc_seg{iSeg}, 1, 1) * Fs; ...
                                   zeros(1, size(acc_seg{iSeg}, 2))];
            end

            % --- Modules (norme euclidienne 3 axes) ---
            mod_acc  = NaN(N_frames, N_seg);
            mod_jerk = NaN(N_frames, N_seg);
            for iSeg = 1:N_seg
                mod_acc(:, iSeg)  = sqrt(sum(acc_seg{iSeg}.^2, 2));
                mod_jerk(:, iSeg) = sqrt(sum(jerk_seg{iSeg}.^2, 2));
            end

            % --- Initialisation des sorties par bin ---
            SegMod_Accel  = NaN(N_bins, N_seg);
            SegMod_Jerk   = NaN(N_bins, N_seg);
            Total_Accel   = NaN(N_bins, 1);
            Total_Jerk    = NaN(N_bins, 1);
            SegAxis_Accel = NaN(N_bins, N_seg, 3);
            SegAxis_Jerk  = NaN(N_bins, N_seg, 3);

            for iBin = 1:N_bins
                i0 = (iBin-1)*bin_size + 1;
                i1 = min(iBin*bin_size, N_frames);
                if (i1 - i0 + 1) < Fs * 0.5, continue; end

                % SegMod : moyenne du module par segment
                SegMod_Accel(iBin, :) = mean(mod_acc(i0:i1, :), 1, 'omitnan');
                SegMod_Jerk(iBin, :)  = mean(mod_jerk(i0:i1, :), 1, 'omitnan');

                % Total : somme du module sur les 7 segments
                % (pas une moyenne — confirmé par Nathan)
                Total_Accel(iBin) = sum(SegMod_Accel(iBin, :), 'omitnan');
                Total_Jerk(iBin)  = sum(SegMod_Jerk(iBin, :), 'omitnan');

                % SegAxis : moyenne du signal brut par segment et par axe
                for iSeg = 1:N_seg
                    for iAx = 1:3
                        SegAxis_Accel(iBin, iSeg, iAx) = ...
                            mean(acc_seg{iSeg}(i0:i1, iAx), 'omitnan');
                        SegAxis_Jerk(iBin, iSeg, iAx) = ...
                            mean(jerk_seg{iSeg}(i0:i1, iAx), 'omitnan');
                    end
                end
            end % iBin

            % --- Stockage SegMod ---
            for iSeg = 1:N_seg
                Results_SegMod.(subj).(kbd).(['Accel_' Seg_names{iSeg}]) = SegMod_Accel(:, iSeg)';
                Results_SegMod.(subj).(kbd).(['Jerk_'  Seg_names{iSeg}]) = SegMod_Jerk(:, iSeg)';
            end
            Results_SegMod.(subj).(kbd).Total_Accel = Total_Accel';
            Results_SegMod.(subj).(kbd).Total_Jerk  = Total_Jerk';

            % --- Stockage SegAxis ---
            for iSeg = 1:N_seg
                for iAx = 1:3
                    vn = ['Accel_Axe' Axes_lbl{iAx} '_' Seg_names{iSeg}];
                    Results_SegAxis.(subj).(kbd).(vn) = SegAxis_Accel(:, iSeg, iAx)';
                    vn = ['Jerk_Axe' Axes_lbl{iAx} '_' Seg_names{iSeg}];
                    Results_SegAxis.(subj).(kbd).(vn) = SegAxis_Jerk(:, iSeg, iAx)';
                end
            end

            fprintf('OK\n');

        end % conditions
    end % sujets

    clear XF tmp;
end % fichiers

%% -----------------------------------------------------------------------
%  RÉSUMÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Résumé ===\n');
subjects_done = fieldnames(Results_SegMod);
N_subj_done   = length(subjects_done);
n_crossover   = 0;
fprintf('  %-6s  %s\n', 'Sujet', 'Conditions');
fprintf('  %s\n', repmat('-', 1, 30));
for iS = 1:N_subj_done
    subj = subjects_done{iS};
    kbds = fieldnames(Results_SegMod.(subj));
    fprintf('  %-6s  %s\n', subj, strjoin(kbds, ' + '));
    if length(kbds) == 2, n_crossover = n_crossover + 1; end
end
fprintf('\n  Total : %d sujets | %d crossover (Norm+CS60)\n', N_subj_done, n_crossover);

% Vérification NaN — SegMod (16 variables)
fprintf('\n  Vérification valeurs valides — SegMod :\n');
var_names_SegMod = {};
for iSeg = 1:N_seg
    var_names_SegMod{end+1} = ['Accel_' Seg_names{iSeg}];
end
for iSeg = 1:N_seg
    var_names_SegMod{end+1} = ['Jerk_' Seg_names{iSeg}];
end
var_names_SegMod{end+1} = 'Total_Accel';
var_names_SegMod{end+1} = 'Total_Jerk';

fprintf('  %-20s  %s\n', 'Variable', 'Valides/Total');
fprintf('  %s\n', repmat('-', 1, 40));
for iV = 1:length(var_names_SegMod)
    vn = var_names_SegMod{iV};
    n_ok = 0; n_tot = 0;
    for iS = 1:N_subj_done
        subj = subjects_done{iS};
        kbds = fieldnames(Results_SegMod.(subj));
        for iK = 1:length(kbds)
            v = Results_SegMod.(subj).(kbds{iK}).(vn);
            n_ok  = n_ok  + sum(~isnan(v));
            n_tot = n_tot + length(v);
        end
    end
    fprintf('  %-20s  %d/%d\n', vn, n_ok, n_tot);
end

% Vérification NaN — SegAxis (42 variables, résumé compact)
fprintf('\n  Vérification valeurs valides — SegAxis (résumé) :\n');
n_ok_tot = 0; n_tot_tot = 0;
for iSeg = 1:N_seg
    for iAx = 1:3
        for sigtype = {'Accel','Jerk'}
            vn = [sigtype{1} '_Axe' Axes_lbl{iAx} '_' Seg_names{iSeg}];
            for iS = 1:N_subj_done
                subj = subjects_done{iS};
                kbds = fieldnames(Results_SegAxis.(subj));
                for iK = 1:length(kbds)
                    v = Results_SegAxis.(subj).(kbds{iK}).(vn);
                    n_ok_tot  = n_ok_tot  + sum(~isnan(v));
                    n_tot_tot = n_tot_tot + length(v);
                end
            end
        end
    end
end
fprintf('  Total SegAxis (42 vars) : %d/%d valides\n', n_ok_tot, n_tot_tot);

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'SegMod_SegAxis_Features_NewDataset.mat'), ...
    'Results_SegMod', 'Results_SegAxis', 'Seg_names', 'Seg_idx', 'Axes_lbl', ...
    'N_bins', 'Fs', '-v7.3');
fprintf('\nSauvegardé : SegMod_SegAxis_Features_NewDataset.mat\n');
fprintf('Structure  : Results_SegMod.(Sujet).(Condition).(Variable)  = [1 x %d]\n', N_bins);
fprintf('             Results_SegAxis.(Sujet).(Condition).(Variable) = [1 x %d]\n', N_bins);