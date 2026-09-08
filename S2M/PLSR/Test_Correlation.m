%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  TEST DE CORRÉLATION — 17 sous-groupes                           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

path_plsr = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\PLSR\';
path_data = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
F_sp  = load(fullfile(path_plsr, 'SpectralFeatures_Ergo.mat'));
F_emg = load(fullfile(path_data, 'EMG_Features_NewDataset_25vars.mat'));
fprintf('Données chargées\n\n');

    function col = collect_var(F_sp, F_emg, vname)
        col = [];
        subjects = fieldnames(F_sp.Results_Spectral);
        for iS = 1:numel(subjects)
            subj = subjects{iS};
            conds = fieldnames(F_sp.Results_Spectral.(subj));
            for iK = 1:numel(conds)
                kbd = conds{iK};
                if isfield(F_sp.Results_Spectral.(subj).(kbd), vname)
                    col = [col; F_sp.Results_Spectral.(subj).(kbd).(vname)(:)];
                elseif isfield(F_emg.Results_EMG, subj) && ...
                       isfield(F_emg.Results_EMG.(subj), kbd) && ...
                       isfield(F_emg.Results_EMG.(subj).(kbd), vname)
                    col = [col; F_emg.Results_EMG.(subj).(kbd).(vname)(:)];
                else
                    col = [col; NaN(10,1)];
                end
            end
        end
    end

    function test_group(F_sp, F_emg, num, gname, vars, thresh)
        fprintf('\n══ Sous-groupe %d — %s\n', num, gname);
        X = []; valid = {};
        for iV = 1:numel(vars)
            col = collect_var(F_sp, F_emg, vars{iV});
            if sum(~isnan(col)) > 20
                X = [X, col]; valid{end+1} = vars{iV};
            else
                fprintf('   [absent] %s\n', vars{iV});
            end
        end
        if numel(valid) < 2
            fprintf('   < 2 variables valides\n'); return;
        end
        R = corr(X,'rows','pairwise');
        n = numel(valid);
        % Matrice
        maxlen = max(cellfun(@length,valid));
        fmt = sprintf('  %%-%ds', maxlen);
        fprintf(fmt,'');
        for j=1:n, fprintf('  %6s',valid{j}(max(1,end-5):end)); end
        fprintf('\n');
        for i=1:n
            fprintf(fmt, valid{i});
            for j=1:n
                if i==j, fprintf('    1.00');
                else,     fprintf('  %+6.2f', R(i,j)); end
            end
            fprintf('\n');
        end
        % Paires > thresh
        found = false;
        for i=1:n
            for j=i+1:n
                if abs(R(i,j)) >= thresh
                    fprintf('  *** r=%+.3f  %s  ↔  %s\n', R(i,j), valid{i}, valid{j});
                    found = true;
                end
            end
        end
        if ~found, fprintf('  (aucune paire r > %.2f)\n', thresh); end
    end

thresh = 0.70;

%% 1. TÊTE
test_group(F_sp,F_emg, 1,'TÊTE 🟡 — features spectrales AngVel + Accel',{
    'AngVel_Mod_MedianFreq_Head',
    'AngVel_Mod_PeakPower_Head',
    'Accel_Mod_Mean_Head',
    'AngVel_X_PeakPower_Head',
    'Accel_X_SpectralEntropy_Head'},thresh);

test_group(F_sp,F_emg, 2,'TÊTE 🔵 — tendance centrale AngVel',{
    'AngVel_Z_Mean_Head',
    'AngVel_Y_Median_Head'},thresh);

%% 2. ÉPAULE
test_group(F_sp,F_emg, 3,'ÉPAULE ⚪ sg1 — percentiles AngVelZ M1',{
    'AngVel_Z_p75_Shoulder',
    'AngVel_Z_p10_Shoulder'},thresh);

test_group(F_sp,F_emg, 4,'ÉPAULE ⚪ sg2 — Mean AngVelX/Y Goubault Chord',{
    'AngVel_X_Mean_Shoulder',
    'AngVel_Y_Mean_Shoulder'},thresh);

test_group(F_sp,F_emg, 5,'ÉPAULE ⚪ sg3 — PeakPowerFreq AccelMag + SpEnt JerkZ Goubault Digital',{
    'Accel_Mod_PeakPowerFreq_Shoulder',
    'Jerk_Z_SpectralEntropy_Shoulder'},thresh);

test_group(F_sp,F_emg, 6,'ÉPAULE 🔵 — EMG DeltAnt Activity vs TFR SpEnt',{
    'Activity_DeltAnt',
    'TFR_SpectralEntropy_DeltAnt'},thresh);

%% 3. TRONC
test_group(F_sp,F_emg, 7,'TRONC 🟢 — Power above4Hz AngVel T8',{
    'AngVel_X_Power_above4Hz_T8',
    'AngVel_Mod_Power_above4Hz_T8'},thresh);
% Trunk_PeakPower_AccelZ absent (pas de segment Trunk dans nos données)

%% 4. BRAS
test_group(F_sp,F_emg, 8,'BRAS 🟢 — PeakPowerFreq + PeakPower AngVel Arm',{
    'AngVel_Mod_PeakPowerFreq_Arm',
    'AngVel_X_PeakPower_Arm'},thresh);

test_group(F_sp,F_emg, 9,'BRAS 🔵 — SampleEntropy Triceps + Max JerkY Arm + TFR MedianFreq Triceps + SpEnt AccelMag Arm',{
    'SampleEntropy_Triceps',
    'Jerk_Y_Max_Arm',
    'TFR_MedianFreq_Triceps',
    'Accel_Mod_SpectralEntropy_Arm'},thresh);

%% 5. AVANT-BRAS
test_group(F_sp,F_emg,10,'FOREARM 🟢 — SpEnt + PeakPowerFreq AngVelMagnitude',{
    'AngVel_Mod_SpectralEntropy_Forearm',
    'AngVel_Mod_PeakPowerFreq_Forearm'},thresh);

test_group(F_sp,F_emg,11,'FOREARM 🔵 — PeakPower vs Power_above4Hz AccelX',{
    'Accel_X_PeakPower_Forearm',
    'Accel_X_Power_above4Hz_Forearm'},thresh);

test_group(F_sp,F_emg,12,'FOREARM 🟠 — Power_above4Hz AccelY + JerkY',{
    'Accel_Y_Power_above4Hz_Forearm',
    'Jerk_Y_Power_above4Hz_Forearm'},thresh);

test_group(F_sp,F_emg,13,'FOREARM 🟣 — Median AccelZ (même variable M1+M2)',{
    'Accel_Z_Median_Forearm'},thresh);

%% 6. MAIN
test_group(F_sp,F_emg,14,'HAND 🟡 — percentiles/Std JerkY',{
    'Jerk_Y_p90_Hand',
    'Jerk_Y_Std_Hand',
    'Jerk_Y_p75_Hand'},thresh);

test_group(F_sp,F_emg,15,'HAND ⚪ — features Jerk Z + AccelZ M3',{
    'Jerk_Z_Power_above4Hz_Hand',
    'Jerk_Z_TotalPower_Hand',
    'Jerk_X_PeakPower_Hand',
    'Jerk_Z_Std_Hand',
    'Accel_Z_Power_above4Hz_Hand',
    'Jerk_Z_p10_Hand'},thresh);

test_group(F_sp,F_emg,16,'HAND 🔵 — Mean vs PeakPower AccelMagnitude',{
    'Accel_Mod_Mean_Hand',
    'Accel_Mod_PeakPower_Hand'},thresh);

test_group(F_sp,F_emg,17,'HAND 🟣 — MedianFreq vs SpEnt AccelY',{
    'Accel_Y_MedianFreq_Hand',
    'Accel_Y_SpectralEntropy_Hand'},thresh);

fprintf('\n=== TERMINÉ ===\n');