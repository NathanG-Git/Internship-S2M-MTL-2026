%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  ANALYSE DE VALIDATION — 9 VARIABLES IMU + EMG                    %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                         %%%%
%%%%                                                                    %%%%
%%%%  Pour chaque variable :                                           %%%%
%%%%    LMM : Variable ~ Temps + Condition + Temps:Condition           %%%%
%%%%          + (1|Participant)                                         %%%%
%%%%    Figures : trajectoires individuelles Normal + CS60             %%%%
%%%%    Tableau : betas, p-valeurs, R²m, R²c                          %%%%
%%%%                                                                    %%%%
%%%%  Blocs :                                                           %%%%
%%%%    SegMod   (2 var) : Accel_Mod_Head, Accel_Mod_Hand             %%%%
%%%%    Goubault (4 var) : MedianFreq, PeakPower x2, SpectralEntropy  %%%%
%%%%    EMG      (3 var) : MedianFreq_Triceps, SpectEnt_Deltoid,      %%%%
%%%%                       SampleEntropy_Biceps                        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_save = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
fprintf('Chargement des données...\n');
load(fullfile(path_data, 'IMU_Features_NewDataset.mat'));   % Results_IMU
load(fullfile(path_data, 'EMG_Features_NewDataset.mat'));   % Results (EMG)
fprintf('  OK\n\n');

%% -----------------------------------------------------------------------
%  EXCLUSION P07 — session Normal tronquée (6784 frames au lieu de ~18400)
%  Exclu entièrement : IMU + EMG + RPE, Normal + CS60
%  -----------------------------------------------------------------------
if isfield(Results_IMU, 'P07')
    Results_IMU = rmfield(Results_IMU, 'P07');
    fprintf('P07 retiré de Results_IMU\n');
end
if isfield(Results, 'P07')
    Results = rmfield(Results, 'P07');
    fprintf('P07 retiré de Results (EMG)\n');
end

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
N_bins  = 10;
t_pct   = linspace(5, 95, N_bins);  % abscisse : % de session

% Sujets crossover uniquement pour les figures et LMM Condition — P07 exclu
Subject_Both = {'P02','P03','P04','P05','P08','P10','P12',...
                'P15','P17','P18','P20','P21','P24'};
% Tous les sujets (pour effet Temps) — P07 déjà retiré de Results_IMU
Subject_All_IMU = fieldnames(Results_IMU);
% Sujets grandes mains (LD) — Normal uniquement, pas de CS60 (cf.
% Analyse_RPE_SDvsLD.m) — sous-ensemble de Subject_All_IMU, complémentaire
% de Subject_Both (qui correspond exactement au groupe SD/petites mains)
Subject_LD = {'P01','P06','P09','P11','P13','P14','P16','P19','P22','P23'};

% Couleurs
col_norm  = [0.20 0.40 0.70];   % bleu  — Normal, N=13 crossover (= SD, petites mains)
col_ergo  = [0.15 0.60 0.35];   % vert  — CS60, N=13 crossover
col_all23 = [0.90 0.70 0.00];   % jaune (gold) — Normal, N=23 (tous sujets)
col_ld    = [0.90 0.40 0.60];   % rose  — Normal, N=10 (LD, grandes mains)

% Définition des blocs
blocs = struct();

blocs(1).name     = 'SegMod';
blocs(1).vars     = {'Accel_Mod_Head', 'Accel_Mod_Hand'};
blocs(1).labels   = {'Accel Mod Tête (m/s²)', 'Accel Mod Main (m/s²)'};
blocs(1).source   = 'IMU';

blocs(2).name     = 'Goubault';
blocs(2).vars     = {'MedianFreq_Accel_Y_Hand', 'PeakPower_Accel_Mod_Hand', ...
                     'PeakPower_AngVel_X_Head', 'SpectralEntropy_AngVel_Mod_Forearm'};
blocs(2).labels   = {'MedianFreq Accel Y Main (Hz)', 'PeakPower Accel Mod Main', ...
                     'PeakPower AngVel X Tête', 'SpectralEntropy AngVel Mod Avant-bras'};
blocs(2).source   = 'IMU';

blocs(3).name     = 'EMG';
blocs(3).vars     = {'TFR_MedianFreq_Triceps', 'TFR_SpectralEntropy_Deltoid', ...
                     'SampleEntropy_Biceps'};
blocs(3).labels   = {'TFR MedianFreq Triceps (Hz)', 'TFR SpectralEntropy Deltoïde', ...
                     'SampleEntropy Biceps'};
blocs(3).source   = 'EMG';

%% -----------------------------------------------------------------------
%  LMM PAR VARIABLE
%  -----------------------------------------------------------------------
fprintf('=== LMM : Variable ~ Temps + Condition + Temps:Condition + (1|Participant) ===\n');
fprintf('    Variables z-scorées | Temps en %% de session\n');
fprintf('    Effet Temps : tous sujets Norm (N=24)\n');
fprintf('    Effet Condition + Interaction : sujets crossover (N=%d)\n\n', length(Subject_Both));

LMM_Results = struct();

for ib = 1:length(blocs)
    bloc   = blocs(ib);
    N_vars = length(bloc.vars);
    fprintf('--- Bloc %s ---\n', bloc.name);

    for iV = 1:N_vars
        vn = bloc.vars{iV};

        % ---------------------------------------------------------------
        % Collecter TOUTES les données pour z-score global
        % ---------------------------------------------------------------
        all_vals = [];
        for iS = 1:length(Subject_All_IMU)
            subj = Subject_All_IMU{iS};
            try
                if strcmp(bloc.source, 'IMU')
                    kbds = fieldnames(Results_IMU.(subj));
                else
                    kbds = fieldnames(Results.(subj));
                end
                for iK = 1:length(kbds)
                    if strcmp(bloc.source, 'IMU')
                        v = Results_IMU.(subj).(kbds{iK}).(vn);
                    else
                        v = Results.(subj).(kbds{iK}).(vn);
                    end
                    all_vals = [all_vals; v(:)];
                end
            catch; end
        end
        mu_z = mean(all_vals, 'omitnan');
        sd_z = std(all_vals,  0, 'omitnan');
        if sd_z == 0, sd_z = 1; end

        % ---------------------------------------------------------------
        % Construire tableau long avec z-score
        % Condition Norm : tous les sujets disponibles
        % Condition CS60 : sujets crossover uniquement
        % ---------------------------------------------------------------
        Y_lmm = []; T_lmm = []; C_lmm = []; S_lmm = [];
        iSubj = 0;

        % --- Tous sujets en Norm ---
        for iS = 1:length(Subject_All_IMU)
            subj = Subject_All_IMU{iS};
            iSubj = iSubj + 1;
            try
                if strcmp(bloc.source, 'IMU')
                    y_n = Results_IMU.(subj).Norm.(vn)(:);
                else
                    y_n = Results.(subj).Norm.(vn)(:);
                end
                y_n_z = (y_n - mu_z) / sd_z;
                ok_n  = ~isnan(y_n_z);
                t_vec = t_pct(:);
                Y_lmm = [Y_lmm; y_n_z(ok_n)];
                T_lmm = [T_lmm; t_vec(ok_n)];
                C_lmm = [C_lmm; zeros(sum(ok_n),1)];   % 0 = Norm
                S_lmm = [S_lmm; iSubj*ones(sum(ok_n),1)];
            catch; end
        end

        % --- Sujets crossover en CS60 ---
        iSubj_cross = 0;
        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            iSubj_cross = iSubj_cross + 1;
            % Retrouver l'indice global du sujet
            idx_global = find(strcmp(Subject_All_IMU, subj));
            if isempty(idx_global), continue; end
            try
                if strcmp(bloc.source, 'IMU')
                    y_e = Results_IMU.(subj).CS60.(vn)(:);
                else
                    y_e = Results.(subj).CS60.(vn)(:);
                end
                y_e_z = (y_e - mu_z) / sd_z;
                ok_e  = ~isnan(y_e_z);
                t_vec = t_pct(:);
                Y_lmm = [Y_lmm; y_e_z(ok_e)];
                T_lmm = [T_lmm; t_vec(ok_e)];
                C_lmm = [C_lmm; ones(sum(ok_e),1)];    % 1 = CS60
                S_lmm = [S_lmm; idx_global*ones(sum(ok_e),1)];
            catch; end
        end

        if length(Y_lmm) < 20
            fprintf('  %-45s : pas assez de données\n', vn);
            continue;
        end

        % Ajustement LMM
        T_tbl = table(Y_lmm, T_lmm, categorical(C_lmm), categorical(S_lmm), ...
            'VariableNames', {'Y','Temps','Condition','Participant'});

        try
            lme = fitlme(T_tbl, 'Y ~ Temps + Condition + Temps:Condition + (1|Participant)');
            [~,~,FE] = fixedEffects(lme, 'DFMethod','satterthwaite');
            [R2m, R2c] = computeNakagawaR2(lme);

            LMM_Results.(vn).lme   = lme;
            LMM_Results.(vn).FE    = FE;
            LMM_Results.(vn).R2m   = R2m;
            LMM_Results.(vn).R2c   = R2c;
            LMM_Results.(vn).bloc  = bloc.name;
            LMM_Results.(vn).mu_z  = mu_z;
            LMM_Results.(vn).sd_z  = sd_z;

            fprintf('  %-45s  R²m=%.3f R²c=%.3f\n', vn, R2m, R2c);
            for ir = 1:height(FE)
                fprintf('    %-35s  β=%+7.4f  p=%s\n', ...
                    FE.Name{ir}, FE.Estimate(ir), p2star(FE.pValue(ir)));
            end
            fprintf('\n');

        catch ME
            fprintf('  %-45s : LMM échoué — %s\n', vn, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
%  LMM AVEC EFFET D'ORDRE — sujets crossover uniquement (N=13)
%  -----------------------------------------------------------------------
%  Y ~ Temps + Condition + Order + Temps:Condition + Temps:Order +
%      Condition:Order + (1|Participant)
%
%  Order = condition jouée en 1ère position (cf. tableau ordre de passage,
%  fourni par l'utilisateur), facteur catégoriel fixe par sujet. N'existe
%  que pour les 13 crossover (Subject_Both) -> ce modèle est ajusté sur ce
%  seul sous-ensemble (pas sur N=23/24 comme le LMM principal ci-dessus),
%  ce qui est nécessaire pour un LRT valide (modèles comparés doivent être
%  ajustés sur EXACTEMENT les mêmes données).
%
%  (1|Participant) conservé : Order ne varie qu'entre sujets (18 obs par
%  sujet partagent la même valeur) -> le retirer créerait une
%  pseudo-réplication massive et des p anti-conservatrices (même
%  mécanisme que l'autocorrélation intra-sujet débusquée dans le GAMM).
%
%  Test retenu : LRT (rapport de vraisemblance) entre modèle réduit
%  (sans Order) et modèle complet -> UN seul test à 3 degrés de liberté
%  pour "l'ordre joue-t-il un rôle quelconque", plutôt que 3 coefficients
%  individuels mal estimés sur seulement 5 sujets Normal-1er vs 8 CS60-1er.
%  -----------------------------------------------------------------------
fprintf('\n=== LMM avec effet d''ordre (LRT, N=13 crossover) — Y ~ Temps + Condition + Order + Temps:Condition + Temps:Order + Condition:Order + (1|Participant) ===\n');
fprintf('    Order : Normal en 1er (N=5) vs CS60 en 1er (N=8)\n\n');

Subject_Order = struct(...
    'P02','Norm', 'P05','Norm', 'P15','Norm', 'P17','Norm', 'P20','Norm', ...
    'P03','CS60', 'P04','CS60', 'P08','CS60', 'P10','CS60', 'P12','CS60', ...
    'P18','CS60', 'P21','CS60', 'P24','CS60');

LMM_Order_Results = struct();

for ib = 1:length(blocs)
    bloc = blocs(ib);
    fprintf('--- Bloc %s ---\n', bloc.name);

    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end
        mu_z = LMM_Results.(vn).mu_z;
        sd_z = LMM_Results.(vn).sd_z;

        Y_lmm = []; T_lmm = []; C_lmm = []; O_lmm = {}; S_lmm = [];

        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            ord  = Subject_Order.(subj);

            try
                if strcmp(bloc.source, 'IMU')
                    y_n = Results_IMU.(subj).Norm.(vn)(:);
                    y_e = Results_IMU.(subj).CS60.(vn)(:);
                else
                    y_n = Results.(subj).Norm.(vn)(:);
                    y_e = Results.(subj).CS60.(vn)(:);
                end
            catch
                continue;
            end

            y_n_z = (y_n - mu_z) / sd_z;
            y_e_z = (y_e - mu_z) / sd_z;
            ok_n  = ~isnan(y_n_z);
            ok_e  = ~isnan(y_e_z);
            t_vec = t_pct(:);
            n_n   = sum(ok_n);
            n_e   = sum(ok_e);

            Y_lmm = [Y_lmm; y_n_z(ok_n); y_e_z(ok_e)];
            T_lmm = [T_lmm; t_vec(ok_n); t_vec(ok_e)];
            C_lmm = [C_lmm; zeros(n_n,1); ones(n_e,1)];      % 0=Norm, 1=CS60
            O_lmm = [O_lmm; repmat({ord}, n_n+n_e, 1)];
            S_lmm = [S_lmm; iS*ones(n_n+n_e,1)];
        end

        if length(Y_lmm) < 20
            fprintf('  %-45s : pas assez de données\n', vn);
            continue;
        end

        T_tbl = table(Y_lmm, T_lmm, categorical(C_lmm), categorical(O_lmm), categorical(S_lmm), ...
            'VariableNames', {'Y','Temps','Condition','Order','Participant'});

        try
            % FitMethod ML obligatoire : un LRT entre modèles différant par
            % leurs effets fixes n'est valide qu'en ML, pas en REML (défaut).
            lme_reduced = fitlme(T_tbl, ...
                'Y ~ Temps + Condition + Temps:Condition + (1|Participant)', ...
                'FitMethod', 'ML');
            lme_full = fitlme(T_tbl, ...
                'Y ~ Temps + Condition + Order + Temps:Condition + Temps:Order + Condition:Order + (1|Participant)', ...
                'FitMethod', 'ML');

            lrt_tbl = compare(lme_reduced, lme_full);
            p_LRT   = lrt_tbl.pValue(2);

            [~,~,FE] = fixedEffects(lme_full, 'DFMethod', 'satterthwaite');

            LMM_Order_Results.(vn).lme_reduced = lme_reduced;
            LMM_Order_Results.(vn).lme_full    = lme_full;
            LMM_Order_Results.(vn).FE          = FE;
            LMM_Order_Results.(vn).p_LRT       = p_LRT;
            LMM_Order_Results.(vn).bloc        = bloc.name;

            fprintf('  %-45s  LRT(Order, 3 df) p=%.4f%s\n', vn, p_LRT, p2star(p_LRT));
            for ir = 1:height(FE)
                if contains(FE.Name{ir}, 'Order')
                    fprintf('    %-35s  β=%+7.4f  p=%s\n', ...
                        FE.Name{ir}, FE.Estimate(ir), p2star(FE.pValue(ir)));
                end
            end
            fprintf('\n');

        catch ME
            fprintf('  %-45s : LMM ordre échoué — %s\n', vn, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
%  TABLEAU RÉCAPITULATIF — EFFET D'ORDRE (LRT)
%  -----------------------------------------------------------------------
fprintf('\n=== TABLEAU RÉCAPITULATIF — Effet d''ordre (LRT, N=13) ===\n');
fprintf('%-45s  %-10s  %-12s  %-10s  %-10s  %-10s\n', ...
    'Variable', 'Bloc', 'p_LRT(3df)', 'β_Order', 'β_T:Order', 'β_C:Order');
fprintf('%s\n', repmat('-', 1, 105));

for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Order_Results, vn), continue; end

        FE    = LMM_Order_Results.(vn).FE;
        p_LRT = LMM_Order_Results.(vn).p_LRT;

        idx_o  = find(contains(FE.Name, 'Order') & ~contains(FE.Name, ':'));
        idx_to = find(contains(FE.Name, 'Temps') & contains(FE.Name, 'Order') & contains(FE.Name, ':'));
        idx_co = find(contains(FE.Name, 'Condition') & contains(FE.Name, 'Order') & contains(FE.Name, ':'));

        beta_o = NaN; beta_to = NaN; beta_co = NaN;
        if ~isempty(idx_o),  beta_o  = FE.Estimate(idx_o(1));  end
        if ~isempty(idx_to), beta_to = FE.Estimate(idx_to(1)); end
        if ~isempty(idx_co), beta_co = FE.Estimate(idx_co(1)); end

        fprintf('%-45s  %-10s  %.4f%s       %+7.4f   %+7.4f   %+7.4f\n', ...
            vn, bloc.name, p_LRT, p2star(p_LRT), beta_o, beta_to, beta_co);
    end
end

%% -----------------------------------------------------------------------
%  COMPARAISON SD (petites mains) vs LD (grandes mains), Normal uniquement
%  -----------------------------------------------------------------------
%  Variable ~ Temps + Group + Temps:Group + (1|Participant)
%
%  Deux variantes, imprimées à la suite pour comparaison directe :
%   - Complète   : SD = Subject_Both (N=13), ordre non contrôlé.
%   - Restreinte : SD = sous-ensemble de Subject_Both ayant joué Normal en
%                  1er (N=5, cf. Subject_Order) -> élimine l'asymétrie de
%                  position de session avec LD (qui n'a qu'une seule
%                  session, donc "position 1" par construction). Motivée
%                  par l'effet d'ordre démontré juste au-dessus (6/9
%                  variables confondues par l'ordre).
%
%  LRT à 2 df (Group + Temps:Group) entre le modèle réduit (Temps seul) et
%  le modèle complet, par variable et par variante.
%  -----------------------------------------------------------------------
fprintf('\n=== Comparaison SD vs LD (Normal) — Variable ~ Temps + Group + Temps:Group + (1|Participant) ===\n');

is_norm1er = cellfun(@(s) strcmp(Subject_Order.(s), 'Norm'), Subject_Both);
Subject_SD_Norm1er = Subject_Both(is_norm1er);
fprintf('  SD complet (Subject_Both) : N=%d | SD Normal-en-1er (restreint) : N=%d | LD : N=%d\n\n', ...
    length(Subject_Both), length(Subject_SD_Norm1er), length(Subject_LD));

group_variants = {Subject_Both, 'Complet (SD N=13 vs LD N=10)'; ...
                   Subject_SD_Norm1er, 'Restreint (SD Normal-1er N=5 vs LD N=10)'};

LMM_SDvsLD_Results = struct();

for ib = 1:length(blocs)
    bloc = blocs(ib);
    fprintf('--- Bloc %s ---\n', bloc.name);

    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end
        mu_z = LMM_Results.(vn).mu_z;
        sd_z = LMM_Results.(vn).sd_z;

        for iVar = 1:size(group_variants, 1)
            Subj_SD_this = group_variants{iVar, 1};
            variant_lbl  = group_variants{iVar, 2};

            Y_g = []; T_g = []; G_g = []; S_g = [];
            iSubj = 0;

            for iS = 1:length(Subj_SD_this)
                subj  = Subj_SD_this{iS};
                iSubj = iSubj + 1;
                try
                    if strcmp(bloc.source, 'IMU')
                        y = Results_IMU.(subj).Norm.(vn)(:);
                    else
                        y = Results.(subj).Norm.(vn)(:);
                    end
                catch
                    continue;
                end
                y_z   = (y - mu_z) / sd_z;
                ok    = ~isnan(y_z);
                t_vec = t_pct(:);
                Y_g = [Y_g; y_z(ok)];
                T_g = [T_g; t_vec(ok)];
                G_g = [G_g; zeros(sum(ok),1)];   % 0 = SD
                S_g = [S_g; iSubj*ones(sum(ok),1)];
            end

            for iS = 1:length(Subject_LD)
                subj  = Subject_LD{iS};
                iSubj = iSubj + 1;
                try
                    if strcmp(bloc.source, 'IMU')
                        y = Results_IMU.(subj).Norm.(vn)(:);
                    else
                        y = Results.(subj).Norm.(vn)(:);
                    end
                catch
                    continue;
                end
                y_z   = (y - mu_z) / sd_z;
                ok    = ~isnan(y_z);
                t_vec = t_pct(:);
                Y_g = [Y_g; y_z(ok)];
                T_g = [T_g; t_vec(ok)];
                G_g = [G_g; ones(sum(ok),1)];    % 1 = LD
                S_g = [S_g; iSubj*ones(sum(ok),1)];
            end

            if length(Y_g) < 20
                fprintf('  [%s] %-40s : pas assez de données\n', variant_lbl, vn);
                continue;
            end

            T_tbl_g = table(Y_g, T_g, categorical(G_g), categorical(S_g), ...
                'VariableNames', {'Y','Temps','Group','Participant'});

            try
                lme_g_reduced = fitlme(T_tbl_g, 'Y ~ Temps + (1|Participant)', 'FitMethod','ML');
                lme_g_full    = fitlme(T_tbl_g, 'Y ~ Temps + Group + Temps:Group + (1|Participant)', 'FitMethod','ML');

                lrt_g   = compare(lme_g_reduced, lme_g_full);
                p_LRT_g = lrt_g.pValue(2);

                [~,~,FE_g] = fixedEffects(lme_g_full, 'DFMethod','satterthwaite');

                idx_grp  = find(strcmp(FE_g.Name, 'Group_1'));
                idx_tgrp = find(contains(FE_g.Name, 'Temps') & contains(FE_g.Name, 'Group') & contains(FE_g.Name, ':'));
                beta_grp = NaN; p_grp = NaN; beta_tgrp = NaN; p_tgrp = NaN;
                if ~isempty(idx_grp),  beta_grp  = FE_g.Estimate(idx_grp);  p_grp  = FE_g.pValue(idx_grp);  end
                if ~isempty(idx_tgrp), beta_tgrp = FE_g.Estimate(idx_tgrp); p_tgrp = FE_g.pValue(idx_tgrp); end

                field_v = sprintf('variant%d', iVar);
                LMM_SDvsLD_Results.(vn).(field_v).label = variant_lbl;
                LMM_SDvsLD_Results.(vn).(field_v).p_LRT = p_LRT_g;
                LMM_SDvsLD_Results.(vn).(field_v).FE    = FE_g;
                LMM_SDvsLD_Results.(vn).(field_v).N_obs = height(T_tbl_g);

                fprintf('  [%s] %-40s  LRT(Group,2df) p=%.4f%s | Group β=%+7.4f p=%.4f%s | Temps:Group β=%+7.4f p=%.4f%s | N=%d\n', ...
                    variant_lbl, vn, p_LRT_g, p2star(p_LRT_g), beta_grp, p_grp, p2star(p_grp), ...
                    beta_tgrp, p_tgrp, p2star(p_tgrp), height(T_tbl_g));

            catch ME
                fprintf('  [%s] %-40s : modèle échoué — %s\n', variant_lbl, vn, ME.message);
            end
        end
    end
    fprintf('\n');
end

%% -----------------------------------------------------------------------
%  TABLEAU RÉCAPITULATIF — SD vs LD, complet vs restreint par ordre
%  -----------------------------------------------------------------------
fprintf('\n=== TABLEAU RÉCAPITULATIF — SD vs LD (Normal), complet vs restreint ===\n');
fprintf('%-45s  %-10s  %-16s  %-16s\n', 'Variable', 'Bloc', 'LRT complet', 'LRT restreint');
fprintf('%s\n', repmat('-', 1, 95));
for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_SDvsLD_Results, vn), continue; end
        p1 = NaN; p2 = NaN;
        if isfield(LMM_SDvsLD_Results.(vn), 'variant1'), p1 = LMM_SDvsLD_Results.(vn).variant1.p_LRT; end
        if isfield(LMM_SDvsLD_Results.(vn), 'variant2'), p2 = LMM_SDvsLD_Results.(vn).variant2.p_LRT; end
        fprintf('%-45s  %-10s  %.4f%-8s  %.4f%-8s\n', ...
            vn, bloc.name, p1, p2star(p1), p2, p2star(p2));
    end
end


%  -----------------------------------------------------------------------
fprintf('=== Génération des figures ===\n');

% Subject_Both et Subject_All_IMU déjà nettoyés de P07 en début de script
N_cross = length(Subject_Both);

for ib = 1:length(blocs)
    bloc   = blocs(ib);
    N_vars = length(bloc.vars);

    fig = figure('Name', sprintf('Trajectoires_%s', bloc.name), ...
        'NumberTitle','off', ...
        'Position', [50 50 350*N_vars 380], ...
        'Color','white');

    for iV = 1:N_vars
        vn  = bloc.vars{iV};
        lbl = bloc.labels{iV};

        % Collecter matrices [N_cross x N_bins]
        mat_n = NaN(N_cross, N_bins);
        mat_e = NaN(N_cross, N_bins);
        for iS = 1:N_cross
            subj = Subject_Both{iS};
            try
                if strcmp(bloc.source, 'IMU')
                    mat_n(iS,:) = Results_IMU.(subj).Norm.(vn);
                    mat_e(iS,:) = Results_IMU.(subj).CS60.(vn);
                else
                    mat_n(iS,:) = Results.(subj).Norm.(vn);
                    mat_e(iS,:) = Results.(subj).CS60.(vn);
                end
            catch; end
        end

        % Troisième série : TOUS les sujets Normal (N=23), même logique que
        % pour Fig_Index_* et Fig_RPE_Scatter_*
        mat_all23 = NaN(length(Subject_All_IMU), N_bins);
        for iS = 1:length(Subject_All_IMU)
            subj = Subject_All_IMU{iS};
            try
                if strcmp(bloc.source, 'IMU')
                    mat_all23(iS,:) = Results_IMU.(subj).Norm.(vn);
                else
                    mat_all23(iS,:) = Results.(subj).Norm.(vn);
                end
            catch; end
        end

        % Quatrième série : sujets LD (grandes mains, N=10), Normal uniquement
        mat_ld = NaN(length(Subject_LD), N_bins);
        for iS = 1:length(Subject_LD)
            subj = Subject_LD{iS};
            try
                if strcmp(bloc.source, 'IMU')
                    mat_ld(iS,:) = Results_IMU.(subj).Norm.(vn);
                else
                    mat_ld(iS,:) = Results.(subj).Norm.(vn);
                end
            catch; end
        end

        mn_n     = mean(mat_n, 1, 'omitnan');
        mn_e     = mean(mat_e, 1, 'omitnan');
        mn_all23 = mean(mat_all23, 1, 'omitnan');
        mn_ld    = mean(mat_ld, 1, 'omitnan');

        % Normalisation par rapport à la première intervalle
        % Soustraction de la moyenne du premier intervalle → les deux
        % conditions partent de 0, ce qui permet de comparer l'évolution
        % indépendamment du niveau de base
        baseline_n     = mn_n(1);
        baseline_e     = mn_e(1);
        baseline_all23 = mn_all23(1);
        baseline_ld    = mn_ld(1);
        mn_n_norm      = mn_n - baseline_n;
        mn_e_norm      = mn_e - baseline_e;
        mn_all23_norm  = mn_all23 - baseline_all23;
        mn_ld_norm     = mn_ld - baseline_ld;

        subplot(1, N_vars, iV); hold on;

        % Courbe Norm — bleu
        plot(t_pct, mn_n_norm, 'o-', ...
            'Color', col_norm, 'LineWidth', 1.3, ...
            'MarkerSize', 4, 'MarkerFaceColor', col_norm);

        % Courbe CS60 — vert
        plot(t_pct, mn_e_norm, 's-', ...
            'Color', col_ergo, 'LineWidth', 1.3, ...
            'MarkerSize', 4, 'MarkerFaceColor', col_ergo);

        % Courbe Normal N=23 — jaune
        plot(t_pct, mn_all23_norm, '^-', ...
            'Color', col_all23, 'LineWidth', 1.3, ...
            'MarkerSize', 4, 'MarkerFaceColor', col_all23);

        % Courbe Normal LD (grandes mains, N=10) — rose
        plot(t_pct, mn_ld_norm, 'o-', ...
            'Color', col_ld, 'LineWidth', 1.3, ...
            'MarkerSize', 4, 'MarkerFaceColor', col_ld);

        % Ligne de référence à zéro
        yline(0, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);

        % Annotations LMM
        if isfield(LMM_Results, vn)
            FE  = LMM_Results.(vn).FE;
            R2m = LMM_Results.(vn).R2m;
            R2c = LMM_Results.(vn).R2c;

            idx_t = find(strcmp(FE.Name, 'Temps'));
            idx_i = find(contains(FE.Name, 'Temps:Condition'));

            p_t = NaN; p_i = NaN;
            if ~isempty(idx_t), p_t = FE.pValue(idx_t); end
            if ~isempty(idx_i), p_i = FE.pValue(idx_i); end

            y_all = [mn_n_norm mn_e_norm mn_all23_norm mn_ld_norm];
            y_min = min(y_all, [], 'omitnan');
            y_max = max(y_all, [], 'omitnan');
            y_rng = y_max - y_min;
            if y_rng == 0, y_rng = 0.01; end

            if isnan(R2c)
                r2_str = sprintf('R²m=%.3f', R2m);
            else
                r2_str = sprintf('R²m=%.3f R²c=%.3f', R2m, R2c);
            end
            txt = sprintf('Temps %s | Inter %s | %s', ...
                p2star_full(p_t), p2star_full(p_i), r2_str);
            text(5, y_min - y_rng*0.12, txt, ...
                'FontSize', 7, 'Color', [0.3 0.3 0.3], ...
                'VerticalAlignment', 'top');
        end

        ylabel(sprintf('Δ %s', lbl), 'FontSize', 9);
        xlabel('% session', 'FontSize', 9);
        title(strrep(vn, '_', ' '), 'FontSize', 9, 'FontWeight', 'bold');
        xlim([0 100]); grid on; box on;

        if iV == 1
            legend({'Normal (N=13)', 'CS60 (N=13)', 'Normal (N=23)', 'Normal LD (N=10)'}, ...
                'Location', 'best', 'FontSize', 7);
        end
    end

    sgtitle(sprintf('Bloc %s — Normal (N=%d) vs CS60 (N=%d) vs Normal total (N=%d) vs Normal LD (N=%d)', ...
        bloc.name, N_cross, N_cross, length(Subject_All_IMU), length(Subject_LD)), ...
        'FontWeight', 'bold', 'FontSize', 11);

    fname = fullfile(path_save, sprintf('Fig_Trajectoires_%s.png', bloc.name));
    saveas(fig, fname);
    fprintf('  Fig_Trajectoires_%s.png sauvegardée\n', bloc.name);
    close(fig);
end

%% -----------------------------------------------------------------------
%  TABLEAU RÉCAPITULATIF LMM
%  -----------------------------------------------------------------------
fprintf('\n=== TABLEAU RÉCAPITULATIF LMM ===\n');
fprintf('%-45s  %-10s  %-10s  %-8s  %-10s  %-8s  %-10s  %-8s  %-12s\n', ...
    'Variable', 'Bloc', 'β_Temps', 'p_Temps', 'β_Cond', 'p_Cond', 'β_Inter', 'p_Inter', 'R²m / R²c');
fprintf('%s\n', repmat('-', 1, 140));

for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end

        FE  = LMM_Results.(vn).FE;
        R2m = LMM_Results.(vn).R2m;
        R2c = LMM_Results.(vn).R2c;

        % Extraire beta et p
        idx_t = find(strcmp(FE.Name, 'Temps'));
        idx_c = find(contains(FE.Name, 'Condition_') & ~contains(FE.Name, 'Temps'));
        idx_i = find(contains(FE.Name, 'Temps:Condition'));

        beta_t = NaN; p_t = NaN;
        beta_c = NaN; p_c = NaN;
        beta_i = NaN; p_i = NaN;

        if ~isempty(idx_t), beta_t = FE.Estimate(idx_t); p_t = FE.pValue(idx_t); end
        if ~isempty(idx_c), beta_c = FE.Estimate(idx_c); p_c = FE.pValue(idx_c); end
        if ~isempty(idx_i), beta_i = FE.Estimate(idx_i); p_i = FE.pValue(idx_i); end

        r2c_str = 'NaN';
        if ~isnan(R2c), r2c_str = sprintf('%.3f', R2c); end

        fprintf('%-45s  %-10s  %+7.4f%s  %-8s  %+7.4f%s  %-8s  %+7.4f%s  %-8s  %.3f / %s\n', ...
            vn, bloc.name, ...
            beta_t, p2star(p_t), sprintf('%.3f', p_t), ...
            beta_c, p2star(p_c), sprintf('%.3f', p_c), ...
            beta_i, p2star(p_i), sprintf('%.3f', p_i), ...
            R2m, r2c_str);
    end
end

%% -----------------------------------------------------------------------
%  INDEX DE BLOC — NON PONDÉRÉ ET PONDÉRÉ, SIGNES GOUBAULT ET OBSERVÉS
%  -----------------------------------------------------------------------
fprintf('\n=== Calcul des index de bloc ===\n');

% -------------------------------------------------------------------
% SIGNES GOUBAULT — direction attendue d'après Goubault G1 (forest plots)
% +1 : variable augmente avec la fatigue → z-score tel quel
% -1 : variable diminue avec la fatigue → inverser le z-score
% -------------------------------------------------------------------
sign_goubault = struct(...
    'Accel_Mod_Head',                    +1, ...  % β=+1.09
    'Accel_Mod_Hand',                    -1, ...  % β=-1.12
    'MedianFreq_Accel_Y_Hand',           +1, ...  % β=+0.96
    'PeakPower_Accel_Mod_Hand',          -1, ...  % β=-0.98
    'PeakPower_AngVel_X_Head',           +1, ...  % β=+0.56
    'SpectralEntropy_AngVel_Mod_Forearm',+1, ...  % β=+0.90
    'TFR_MedianFreq_Triceps',            -1, ...  % β=-0.58
    'TFR_SpectralEntropy_Deltoid',       -1, ...  % β=-0.59
    'SampleEntropy_Biceps',              -1);      % β=-0.88

% -------------------------------------------------------------------
% SIGNES OBSERVÉS — direction réelle du β_Temps dans ce dataset
% Déterminés automatiquement à partir de LMM_Results (déjà calculés)
% +1 si β_Temps > 0, -1 si β_Temps < 0
% -------------------------------------------------------------------
all_vars_idx = {'Accel_Mod_Head','Accel_Mod_Hand','MedianFreq_Accel_Y_Hand', ...
    'PeakPower_Accel_Mod_Hand','PeakPower_AngVel_X_Head', ...
    'SpectralEntropy_AngVel_Mod_Forearm','TFR_MedianFreq_Triceps', ...
    'TFR_SpectralEntropy_Deltoid','SampleEntropy_Biceps'};

sign_observed = struct();
fprintf('Signes observés (β_Temps de ce dataset) :\n');
for iV = 1:length(all_vars_idx)
    vn = all_vars_idx{iV};
    if isfield(LMM_Results, vn)
        FE = LMM_Results.(vn).FE;
        idx_t = find(strcmp(FE.Name, 'Temps'));
        if ~isempty(idx_t)
            b = FE.Estimate(idx_t);
            s = sign(b); if s == 0, s = 1; end
            sign_observed.(vn) = s;
            match = '';
            if s ~= sign_goubault.(vn), match = '  <-- diffère de Goubault'; end
            fprintf('  %-40s  %+d%s\n', vn, s, match);
        else
            sign_observed.(vn) = sign_goubault.(vn);
        end
    end
end

%% -----------------------------------------------------------------------
%  SENS DE LA PENTE (β_Temps) — SD vs LD vs SIGNE G1 HARDCODÉ
%  -----------------------------------------------------------------------
%  Extraction gratuite depuis LMM_SDvsLD_Results (déjà ajusté plus haut,
%  aucun nouveau modèle ici) : pour la variante Restreinte (SD Normal-en-
%  1er N=5 vs LD N=10, contrôlée pour l'ordre — préférée ici à la variante
%  Complète pour la même raison méthodologique que partout ailleurs dans
%  cette phase), on lit directement dans la table FE :
%    beta_Temps_SD = coefficient 'Temps'            (Group=0 -> SD = référence)
%    beta_Temps_LD = beta_Temps_SD + beta 'Temps:Group'
%  Comparé au signe hardcodé sign_goubault (thèse M1, G1 uniquement — G2
%  non disponible pour ces 9 variables, cf. discussion précédente).
%  -----------------------------------------------------------------------
fprintf('\n=== Sens de la pente (β_Temps) — SD vs LD (restreint par ordre) vs signe G1 ===\n');
fprintf('%-40s  %-8s  %-10s  %-10s  %-12s  %-12s  %-10s\n', ...
    'Variable', 'Signe G1', 'Pente SD', 'Pente LD', 'G1 vs SD', 'G1 vs LD', 'SD vs LD');
fprintf('%s\n', repmat('-', 1, 105));

Slope_SDvsLD_vs_G1 = struct();

for iV = 1:length(all_vars_idx)
    vn = all_vars_idx{iV};
    if ~isfield(LMM_SDvsLD_Results, vn) || ~isfield(LMM_SDvsLD_Results.(vn), 'variant2')
        continue;
    end

    FE_r = LMM_SDvsLD_Results.(vn).variant2.FE;   % variante Restreinte
    idx_t  = find(strcmp(FE_r.Name, 'Temps'));
    idx_tg = find(contains(FE_r.Name, 'Temps') & contains(FE_r.Name, 'Group') & contains(FE_r.Name, ':'));
    if isempty(idx_t) || isempty(idx_tg), continue; end

    beta_SD = FE_r.Estimate(idx_t);
    beta_LD = beta_SD + FE_r.Estimate(idx_tg);

    sign_SD = sign(beta_SD); if sign_SD == 0, sign_SD = 1; end
    sign_LD = sign(beta_LD); if sign_LD == 0, sign_LD = 1; end
    sign_G1 = sign_goubault.(vn);

    match_G1_SD = 'identique';  if sign_G1 ~= sign_SD, match_G1_SD = 'DIFFÈRE'; end
    match_G1_LD = 'identique';  if sign_G1 ~= sign_LD, match_G1_LD = 'DIFFÈRE'; end
    match_SD_LD = 'identique';  if sign_SD ~= sign_LD, match_SD_LD = 'DIFFÈRE'; end

    Slope_SDvsLD_vs_G1.(vn).beta_SD = beta_SD;
    Slope_SDvsLD_vs_G1.(vn).beta_LD = beta_LD;
    Slope_SDvsLD_vs_G1.(vn).sign_SD = sign_SD;
    Slope_SDvsLD_vs_G1.(vn).sign_LD = sign_LD;
    Slope_SDvsLD_vs_G1.(vn).sign_G1 = sign_G1;

    fprintf('%-40s  %+8d  %+7.4f   %+7.4f   %-12s  %-12s  %-10s\n', ...
        vn, sign_G1, beta_SD, beta_LD, match_G1_SD, match_G1_LD, match_SD_LD);
end

% -------------------------------------------------------------------
% POIDS — 3 schémas de pondération basés sur le pipeline Goubault G1
%   beta : |β_std| issu du LMM Goubault (force d'effet)
%   taux : fréquence de sélection sur les 100 runs (robustesse statistique)
%   bxt  : β_std x fréquence (force ET robustesse combinées)
% -------------------------------------------------------------------
weights_beta = struct(...
    'Accel_Mod_Head',                    1.093, ...
    'Accel_Mod_Hand',                    1.120, ...
    'MedianFreq_Accel_Y_Hand',           0.956, ...
    'PeakPower_Accel_Mod_Hand',          0.984, ...
    'PeakPower_AngVel_X_Head',           0.555, ...
    'SpectralEntropy_AngVel_Mod_Forearm',0.897, ...
    'TFR_MedianFreq_Triceps',            0.579, ...
    'TFR_SpectralEntropy_Deltoid',       0.591, ...
    'SampleEntropy_Biceps',              0.879);

% Fréquence de sélection sur les 100 runs (en proportion 0-1)
weights_taux = struct(...
    'Accel_Mod_Head',                    0.80, ...
    'Accel_Mod_Hand',                    0.82, ...
    'MedianFreq_Accel_Y_Hand',           1.00, ...
    'PeakPower_Accel_Mod_Hand',          1.00, ...
    'PeakPower_AngVel_X_Head',           1.00, ...
    'SpectralEntropy_AngVel_Mod_Forearm',1.00, ...
    'TFR_MedianFreq_Triceps',            0.97, ...
    'TFR_SpectralEntropy_Deltoid',       0.92, ...
    'SampleEntropy_Biceps',              0.68);

% β x taux — combinaison force d'effet et robustesse
weights_bxt = struct();
wb_fields = fieldnames(weights_beta);
for iF = 1:length(wb_fields)
    fn = wb_fields{iF};
    weights_bxt.(fn) = weights_beta.(fn) * weights_taux.(fn);
end

% -------------------------------------------------------------------
% DÉFINITION DES BLOCS POUR L'INDEX
% -------------------------------------------------------------------
idx_blocs = struct();

idx_blocs(1).name   = 'SegMod';
idx_blocs(1).vars   = {'Accel_Mod_Head', 'Accel_Mod_Hand'};
idx_blocs(1).source = 'IMU';

idx_blocs(2).name   = 'Goubault';
idx_blocs(2).vars   = {'MedianFreq_Accel_Y_Hand', 'PeakPower_Accel_Mod_Hand', ...
                       'PeakPower_AngVel_X_Head', 'SpectralEntropy_AngVel_Mod_Forearm'};
idx_blocs(2).source = 'IMU';

idx_blocs(3).name   = 'EMG';
idx_blocs(3).vars   = {'TFR_MedianFreq_Triceps', 'TFR_SpectralEntropy_Deltoid', ...
                       'SampleEntropy_Biceps'};
idx_blocs(3).source = 'EMG';

% -------------------------------------------------------------------
% CALCUL DES VARIANTES D'INDEX PAR SUJET, CONDITION, INTERVALLE
%   Pour chaque jeu de signes (goubault | observed) :
%     np  : non pondéré
%     pB  : pondéré par β_std Goubault
%     pT  : pondéré par taux de sélection (100 runs)
%     pBT : pondéré par β_std x taux
% -------------------------------------------------------------------
Index_Results = struct();

sign_sets  = {sign_goubault, sign_observed};
sign_names = {'goubault', 'observed'};

weight_sets  = {[], weights_beta, weights_taux, weights_bxt};  % [] = non pondéré
weight_names = {'np', 'pB', 'pT', 'pBT'};

for ib = 1:length(idx_blocs)
    bloc  = idx_blocs(ib);
    bname = bloc.name;
    vars  = bloc.vars;
    N_v   = length(vars);

    fprintf('\n--- Index %s (%d variables) ---\n', bname, N_v);

    % Liste (sujet, condition) à traiter : TOUS les sujets en Normal (N=23,
    % Subject_All_IMU) + les sujets crossover en CS60 (N=13, Subject_Both).
    % Permet d'avoir, en plus des deux courbes crossover existantes, une
    % troisième série "Normal, N=23" (comportement attendu : similaire mais
    % amplifié par rapport à "Normal, N=13", avec plus de puissance).
    subj_kbd_list = {};
    for iS = 1:length(Subject_All_IMU)
        subj_kbd_list(end+1,:) = {Subject_All_IMU{iS}, 'Norm'};
    end
    for iS = 1:length(Subject_Both)
        subj_kbd_list(end+1,:) = {Subject_Both{iS}, 'CS60'};
    end

    for iSK = 1:size(subj_kbd_list, 1)
        subj = subj_kbd_list{iSK, 1};
        kbd  = subj_kbd_list{iSK, 2};

        for iSign = 1:2
            sign_set  = sign_sets{iSign};
            sign_name = sign_names{iSign};

            % Construire la matrice de z-scores signés une seule fois
            Z_mat = NaN(N_bins, N_v);
            for iV = 1:N_v
                vn = vars{iV};
                if ~isfield(LMM_Results, vn), continue; end
                try
                    if strcmp(bloc.source, 'IMU')
                        y = Results_IMU.(subj).(kbd).(vn)(:);
                    else
                        y = Results.(subj).(kbd).(vn)(:);
                    end
                catch; continue; end

                mu_z = LMM_Results.(vn).mu_z;
                sd_z = LMM_Results.(vn).sd_z;
                z    = (y - mu_z) / sd_z;
                Z_mat(:, iV) = z * sign_set.(vn);
            end

            % Pour chaque schéma de pondération
            for iW = 1:length(weight_sets)
                w_set  = weight_sets{iW};
                w_name = weight_names{iW};

                if isempty(w_set)
                    idx_val = mean(Z_mat, 2, 'omitnan');
                else
                    Z_pond = NaN(N_bins, N_v);
                    w_sum  = 0;
                    for iV = 1:N_v
                        vn = vars{iV};
                        if ~isfield(w_set, vn), continue; end
                        w = w_set.(vn);
                        Z_pond(:, iV) = Z_mat(:, iV) * w;
                        w_sum = w_sum + w;
                    end
                    if w_sum > 0
                        idx_val = sum(Z_pond, 2, 'omitnan') / w_sum;
                    else
                        idx_val = mean(Z_mat, 2, 'omitnan');
                    end
                end

                field_name = sprintf('%s_%s', sign_name, w_name);
                Index_Results.(subj).(kbd).(bname).(field_name) = idx_val';
            end
        end
    end
end

% -------------------------------------------------------------------
% LMM SUR TOUTES LES VARIANTES D'INDEX (8 par bloc : 2 signes x 4 poids)
% -------------------------------------------------------------------
fprintf('\n=== LMM sur les index de bloc — Index ~ Temps + Condition + Temps:Condition + (1|Participant) ===\n');
LMM_Index = struct();

variant_types = {};
variant_lbls  = {};
sign_lbls_map = struct('goubault','Goubault', 'observed','Observé');
weight_lbls_map = struct('np','Non pondéré', 'pB','Pondéré β', ...
                          'pT','Pondéré taux', 'pBT','Pondéré β×taux');
for iSign = 1:2
    for iW = 1:length(weight_names)
        variant_types{end+1} = sprintf('%s_%s', sign_names{iSign}, weight_names{iW});
        variant_lbls{end+1}  = sprintf('%s — %s', ...
            sign_lbls_map.(sign_names{iSign}), weight_lbls_map.(weight_names{iW}));
    end
end

for ib = 1:length(idx_blocs)
    bname = idx_blocs(ib).name;

    for iType = 1:length(variant_types)
        type_name = variant_types{iType};
        type_lbl  = variant_lbls{iType};

        fprintf('\n--- %s | %s ---\n', bname, type_lbl);

        Y_lmm = []; T_lmm = []; C_lmm = []; S_lmm = [];
        iSubj = 0;

        for iS = 1:length(Subject_Both)
            subj  = Subject_Both{iS};
            iSubj = iSubj + 1;

            try
                y_n = Index_Results.(subj).Norm.(bname).(type_name)(:);
                y_e = Index_Results.(subj).CS60.(bname).(type_name)(:);
            catch; continue; end

            t_vec = t_pct(:);
            ok_n  = ~isnan(y_n); ok_e = ~isnan(y_e);

            Y_lmm = [Y_lmm; y_n(ok_n); y_e(ok_e)];
            T_lmm = [T_lmm; t_vec(ok_n); t_vec(ok_e)];
            C_lmm = [C_lmm; zeros(sum(ok_n),1); ones(sum(ok_e),1)];
            S_lmm = [S_lmm; iSubj*ones(sum(ok_n),1); iSubj*ones(sum(ok_e),1)];
        end

        if length(Y_lmm) < 20, continue; end

        T_tbl = table(Y_lmm, T_lmm, categorical(C_lmm), categorical(S_lmm), ...
            'VariableNames', {'Y','Temps','Condition','Participant'});

        try
            lme = fitlme(T_tbl, 'Y ~ Temps + Condition + Temps:Condition + (1|Participant)');
            [~,~,FE] = fixedEffects(lme, 'DFMethod','satterthwaite');
            [R2m, R2c] = computeNakagawaR2(lme);

            LMM_Index.(bname).(type_name).lme  = lme;
            LMM_Index.(bname).(type_name).FE   = FE;
            LMM_Index.(bname).(type_name).R2m  = R2m;
            LMM_Index.(bname).(type_name).R2c  = R2c;

            fprintf('  R²m=%.3f R²c=%.3f\n', R2m, R2c);
            for ir = 1:height(FE)
                fprintf('  %-35s  β=%+7.4f  p=%s\n', ...
                    FE.Name{ir}, FE.Estimate(ir), p2star(FE.pValue(ir)));
            end
        catch ME
            fprintf('  LMM échoué : %s\n', ME.message);
        end
    end
end

% -------------------------------------------------------------------
% FIGURES INDEX — 4 SOUS-GRAPHIQUES PAR BLOC
%   Haut   : signes Goubault   (non pondéré | pondéré)
%   Bas    : signes observés   (non pondéré | pondéré)
% -------------------------------------------------------------------
fprintf('\n=== Figures index de bloc ===\n');

for ib = 1:length(idx_blocs)
    bname = idx_blocs(ib).name;

    fig = figure('Name', sprintf('Index_%s', bname), ...
        'NumberTitle','off', 'Position',[30 30 1700 750], 'Color','white');

    for iType = 1:length(variant_types)
        type_name = variant_types{iType};
        type_lbl  = variant_lbls{iType};

        mat_n = NaN(length(Subject_Both), N_bins);
        mat_e = NaN(length(Subject_Both), N_bins);
        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            try
                mat_n(iS,:) = Index_Results.(subj).Norm.(bname).(type_name);
                mat_e(iS,:) = Index_Results.(subj).CS60.(bname).(type_name);
            catch; end
        end

        % Troisième série : TOUS les sujets Normal (N=23), pas seulement
        % les crossover -> comportement attendu : similaire mais amplifié
        % par rapport à la courbe Normal N=13 (plus de puissance, même
        % condition, cf. demande explicite de comparaison N=13 vs N=23)
        mat_all23 = NaN(length(Subject_All_IMU), N_bins);
        for iS = 1:length(Subject_All_IMU)
            subj = Subject_All_IMU{iS};
            try
                mat_all23(iS,:) = Index_Results.(subj).Norm.(bname).(type_name);
            catch; end
        end

        % Quatrième série : sujets LD (grandes mains, N=10), Normal uniquement
        mat_ld = NaN(length(Subject_LD), N_bins);
        for iS = 1:length(Subject_LD)
            subj = Subject_LD{iS};
            try
                mat_ld(iS,:) = Index_Results.(subj).Norm.(bname).(type_name);
            catch; end
        end

        mn_n     = mean(mat_n, 1, 'omitnan');
        mn_e     = mean(mat_e, 1, 'omitnan');
        mn_all23 = mean(mat_all23, 1, 'omitnan');
        mn_ld    = mean(mat_ld, 1, 'omitnan');

        % Normalisation baseline → partir de 0
        mn_n     = mn_n - mn_n(1);
        mn_e     = mn_e - mn_e(1);
        mn_all23 = mn_all23 - mn_all23(1);
        mn_ld    = mn_ld - mn_ld(1);

        subplot(2, 4, iType); hold on;

        plot(t_pct, mn_n, 'o-', 'Color',col_norm, 'LineWidth',1.1, ...
            'MarkerSize',3.5, 'MarkerFaceColor',col_norm);
        plot(t_pct, mn_e, 's-', 'Color',col_ergo, 'LineWidth',1.1, ...
            'MarkerSize',3.5, 'MarkerFaceColor',col_ergo);
        plot(t_pct, mn_all23, '^-', 'Color',col_all23, 'LineWidth',1.1, ...
            'MarkerSize',3.5, 'MarkerFaceColor',col_all23);
        plot(t_pct, mn_ld, 'o-', 'Color',col_ld, 'LineWidth',1.1, ...
            'MarkerSize',3.5, 'MarkerFaceColor',col_ld);
        yline(0, '--k', 'LineWidth',0.8, 'Alpha',0.4);

        % Annotation LMM
        if isfield(LMM_Index, bname) && isfield(LMM_Index.(bname), type_name)
            FE  = LMM_Index.(bname).(type_name).FE;
            R2m = LMM_Index.(bname).(type_name).R2m;
            R2c = LMM_Index.(bname).(type_name).R2c;

            idx_t = find(strcmp(FE.Name,'Temps'));
            idx_i = find(contains(FE.Name,'Temps:Condition'));
            p_t = NaN; p_i = NaN;
            if ~isempty(idx_t), p_t = FE.pValue(idx_t); end
            if ~isempty(idx_i), p_i = FE.pValue(idx_i); end

            y_all = [mn_n mn_e mn_all23 mn_ld];
            y_min = min(y_all,[],'omitnan');
            y_max = max(y_all,[],'omitnan');
            y_rng = max(y_max - y_min, 0.01);

            txt = sprintf('T %s | I %s\nR²m=%.3f R²c=%.3f', ...
                p2star_full(p_t), p2star_full(p_i), R2m, R2c);
            text(5, y_min - y_rng*0.15, txt, 'FontSize',6, ...
                'Color',[0.3 0.3 0.3], 'VerticalAlignment','top');
        end

        if iType == 1
            legend({'Normal (N=13)','CS60 (N=13)','Normal (N=23)','Normal LD (N=10)'}, ...
                'Location','best', 'FontSize',6.5);
        end
        ylabel('Δ Index (z)', 'FontSize',8);
        xlabel('% session', 'FontSize',8);
        title(type_lbl, 'FontSize',8, 'FontWeight','bold');
        xlim([0 100]); grid on; box on;
    end

    sgtitle(sprintf('Index de bloc %s — Normal (N=%d) vs CS60 (N=%d) vs Normal total (N=%d) vs Normal LD (N=%d) — Goubault (haut) vs Observé (bas)', ...
        bname, length(Subject_Both), length(Subject_Both), length(Subject_All_IMU), length(Subject_LD)), 'FontWeight','bold', 'FontSize',12);

    fname = fullfile(path_save, sprintf('Fig_Index_%s.png', bname));
    saveas(fig, fname);
    fprintf('  Fig_Index_%s.png sauvegardée\n', bname);
    close(fig);
end

%% -----------------------------------------------------------------------
%  DONNÉES RPE — saisies manuellement, P07 exclu
%  (reconstruit ici : cette section avait été supprimée du pipeline ;
%  mêmes valeurs que Analyse_RPE.m, dupliquées pour générer directement
%  les scatter dérive<->ΔRPE avec la 3e série Normal N=23)
%  -----------------------------------------------------------------------
fprintf('\n=== Chargement RPE pour scatter dérive <-> RPE ===\n');
RPE_data = {
    'P01',  0, 30, NaN, NaN;
    'P02', 30, 75,  23,  35;
    'P03', 25, 85,  30,  70;
    'P04', 15, 23,  13,  25;
    'P05', 35, 75,  45,  72;
    'P06', 10, 55, NaN, NaN;
    'P08',  7, 18,   7,  23;
    'P09',  5, 60, NaN, NaN;
    'P10',  7, 30,   4,  25;
    'P11',  7, 60, NaN, NaN;
    'P12', 13, 95,  13, 105;
    'P13', 17, 55, NaN, NaN;
    'P14', 22, 28, NaN, NaN;
    'P15',  5, 70,   2,  55;
    'P16', 20, 60, NaN, NaN;
    'P17', 13, 60,  15,  50;
    'P18', 13, 60,  15,  50;
    'P19',  0, 25, NaN, NaN;
    'P20', 15, 55,  15,  40;
    'P21', 10, 50,  10,  24;
    'P22',  1, 20, NaN, NaN;
    'P23',  3, 35, NaN, NaN;
    'P24', 10, 25,   7,  16;
};

RPE = struct();
for i = 1:size(RPE_data,1)
    subj = RPE_data{i,1};
    RPE.(subj).Norm.debut = RPE_data{i,2};
    RPE.(subj).Norm.fin   = RPE_data{i,3};
    RPE.(subj).Norm.delta = RPE_data{i,3} - RPE_data{i,2};
    if ~isnan(RPE_data{i,4})
        RPE.(subj).CS60.debut = RPE_data{i,4};
        RPE.(subj).CS60.fin   = RPE_data{i,5};
        RPE.(subj).CS60.delta = RPE_data{i,5} - RPE_data{i,4};
    end
end
fprintf('  RPE chargé pour %d sujets\n', size(RPE_data,1));

%% -----------------------------------------------------------------------
%  LMM SUR LE RPE — EFFET D'ORDRE (LRT, N=13 crossover, N_obs=26)
%  -----------------------------------------------------------------------
%  Delta_RPE ~ Condition + Order + Condition:Order + (1|Participant)
%
%  Même logique que le LMM d'ordre sur les biomarqueurs ci-dessus, mais sur
%  Delta_RPE (RPE_fin - RPE_début) plutôt que sur la trajectoire à 10 bins
%  -> pas de terme Temps ici (un seul point par sujet/condition), donc pas
%  d'interaction à 3 voies à estimer (qui serait inestimable sur N=13).
%  Condition:Order reste cependant directement comparable en esprit au
%  Temps:Condition:Order qu'on aurait eu sur les biomarqueurs : "l'écart
%  CS60-Normal sur le RPE dépend-il de l'ordre de passage ?"
%
%  Order n'est défini que pour les 13 crossover (Subject_Order, déjà
%  construit plus haut dans ce script). (1|Participant) conservé pour la
%  même raison qu'avant (Condition ne varie qu'intra-sujet, mais Order
%  reste fixe par sujet -> pseudo-réplication si on l'enlève).
%  -----------------------------------------------------------------------
fprintf('\n=== LMM RPE avec effet d''ordre (LRT, N=13 crossover, N_obs=26) — DeltaRPE ~ Condition + Order + Condition:Order + (1|Participant) ===\n');

Y_rpe = []; C_rpe = []; O_rpe = []; S_rpe = [];
for iS = 1:length(Subject_Both)
    subj = Subject_Both{iS};
    if ~isfield(RPE, subj), continue; end
    ord = Subject_Order.(subj);

    if isfield(RPE.(subj), 'Norm')
        Y_rpe = [Y_rpe; RPE.(subj).Norm.delta];
        C_rpe = [C_rpe; 0];                 % 0 = Norm
        O_rpe = [O_rpe; {ord}];
        S_rpe = [S_rpe; iS];
    end
    if isfield(RPE.(subj), 'CS60')
        Y_rpe = [Y_rpe; RPE.(subj).CS60.delta];
        C_rpe = [C_rpe; 1];                 % 1 = CS60
        O_rpe = [O_rpe; {ord}];
        S_rpe = [S_rpe; iS];
    end
end

T_rpe = table(Y_rpe, categorical(C_rpe), categorical(O_rpe), categorical(S_rpe), ...
    'VariableNames', {'DeltaRPE','Condition','Order','Participant'});

fprintf('  N_obs=%d (attendu 26)\n', height(T_rpe));

try
    lme_rpe_reduced = fitlme(T_rpe, 'DeltaRPE ~ Condition + (1|Participant)', 'FitMethod','ML');
    lme_rpe_full    = fitlme(T_rpe, 'DeltaRPE ~ Condition + Order + Condition:Order + (1|Participant)', 'FitMethod','ML');

    lrt_rpe   = compare(lme_rpe_reduced, lme_rpe_full);
    p_LRT_rpe = lrt_rpe.pValue(2);

    [~,~,FE_rpe] = fixedEffects(lme_rpe_full, 'DFMethod','satterthwaite');

    fprintf('  LRT(Order, 2 df) p=%.4f%s\n', p_LRT_rpe, p2star(p_LRT_rpe));
    for ir = 1:height(FE_rpe)
        fprintf('    %-35s  β=%+7.4f  p=%.4f%s\n', ...
            FE_rpe.Name{ir}, FE_rpe.Estimate(ir), FE_rpe.pValue(ir), p2star(FE_rpe.pValue(ir)));
    end

    LMM_RPE_Order_Results.lme_reduced = lme_rpe_reduced;
    LMM_RPE_Order_Results.lme_full    = lme_rpe_full;
    LMM_RPE_Order_Results.FE          = FE_rpe;
    LMM_RPE_Order_Results.p_LRT       = p_LRT_rpe;
    LMM_RPE_Order_Results.N_obs       = height(T_rpe);

catch ME
    fprintf('  LMM RPE ordre échoué : %s\n', ME.message);
    LMM_RPE_Order_Results = struct();
end

%% -----------------------------------------------------------------------
%  MODÈLE A — LE RPE0 DÉPEND-IL DE L'ORDRE DE PASSAGE ? (N=13, N_obs=26)
%  -----------------------------------------------------------------------
%  RPE_début ~ Condition + Order + Condition:Order + (1|Participant)
%  Même structure que le modèle ΔRPE ci-dessus, mais sur RPE_début plutôt
%  que sur ΔRPE -> teste si le niveau de DÉPART (pas la progression)
%  diffère selon qui a joué quoi en premier.
%  -----------------------------------------------------------------------
fprintf('\n=== LMM RPE_début avec effet d''ordre (LRT, N=13 crossover, N_obs=26) — RPE_debut ~ Condition + Order + Condition:Order + (1|Participant) ===\n');

Y0_rpe = []; C0_rpe = []; O0_rpe = {}; S0_rpe = [];
for iS = 1:length(Subject_Both)
    subj = Subject_Both{iS};
    if ~isfield(RPE, subj), continue; end
    ord = Subject_Order.(subj);

    if isfield(RPE.(subj), 'Norm')
        Y0_rpe = [Y0_rpe; RPE.(subj).Norm.debut];
        C0_rpe = [C0_rpe; 0];
        O0_rpe = [O0_rpe; {ord}];
        S0_rpe = [S0_rpe; iS];
    end
    if isfield(RPE.(subj), 'CS60')
        Y0_rpe = [Y0_rpe; RPE.(subj).CS60.debut];
        C0_rpe = [C0_rpe; 1];
        O0_rpe = [O0_rpe; {ord}];
        S0_rpe = [S0_rpe; iS];
    end
end

T0_rpe = table(Y0_rpe, categorical(C0_rpe), categorical(O0_rpe), categorical(S0_rpe), ...
    'VariableNames', {'RPE_debut','Condition','Order','Participant'});
fprintf('  N_obs=%d (attendu 26)\n', height(T0_rpe));

try
    lme0_reduced = fitlme(T0_rpe, 'RPE_debut ~ Condition + (1|Participant)', 'FitMethod','ML');
    lme0_full    = fitlme(T0_rpe, 'RPE_debut ~ Condition + Order + Condition:Order + (1|Participant)', 'FitMethod','ML');

    lrt0   = compare(lme0_reduced, lme0_full);
    p_LRT0 = lrt0.pValue(2);

    [~,~,FE0] = fixedEffects(lme0_full, 'DFMethod','satterthwaite');

    fprintf('  LRT(Order, 2 df) p=%.4f%s\n', p_LRT0, p2star(p_LRT0));
    for ir = 1:height(FE0)
        fprintf('    %-35s  β=%+7.4f  p=%.4f%s\n', ...
            FE0.Name{ir}, FE0.Estimate(ir), FE0.pValue(ir), p2star(FE0.pValue(ir)));
    end

    LMM_RPE0_Order_Results.lme_reduced = lme0_reduced;
    LMM_RPE0_Order_Results.lme_full    = lme0_full;
    LMM_RPE0_Order_Results.FE          = FE0;
    LMM_RPE0_Order_Results.p_LRT       = p_LRT0;
    LMM_RPE0_Order_Results.N_obs       = height(T0_rpe);

catch ME
    fprintf('  LMM RPE_début ordre échoué : %s\n', ME.message);
    LMM_RPE0_Order_Results = struct();
end

%% -----------------------------------------------------------------------
%  MODÈLE B — LE RPE0 INFLUENCE-T-IL LE RPE FINAL ? (tous sujets, N_obs<=36)
%  -----------------------------------------------------------------------
%  RPE_fin ~ RPE_début + Condition + (1|Participant)
%
%  ATTENTION : ce n'est PAS le même test que corr(RPE_début, ΔRPE) fait
%  précédemment (Analyse_RPE_SDvsLD.m, ρ≈0,09-0,18, ns). ΔRPE = RPE_fin −
%  RPE_début est mécaniquement défavorable à cette question : si
%  RPE_fin ≈ RPE_début + constante + bruit (relation forte, pente ≈1),
%  alors corr(RPE_début, ΔRPE) ≈ 0 MÊME QUAND RPE_début explique presque
%  toute la variance de RPE_fin. La régression directe RPE_fin ~ RPE_début
%  est le bon test pour répondre à "RPE0 influence-t-il RPE_final ?".
%
%  Échantillon plus large que les modèles d'ordre ci-dessus (23 sujets en
%  Normal + 13 en CS60, jusqu'à 36 obs) : cette question ne nécessite pas
%  de connaître l'ordre de passage, donc pas de restriction aux crossover.
%  -----------------------------------------------------------------------
fprintf('\n=== LMM RPE_fin ~ RPE_début + Condition (influence de RPE0, N_obs<=36) ===\n');

Yb_rpe = []; Xb_rpe = []; Cb_rpe = []; Sb_rpe = [];
for iS = 1:length(Subject_All_IMU)
    subj = Subject_All_IMU{iS};
    if ~isfield(RPE, subj), continue; end

    if isfield(RPE.(subj), 'Norm')
        Yb_rpe = [Yb_rpe; RPE.(subj).Norm.fin];
        Xb_rpe = [Xb_rpe; RPE.(subj).Norm.debut];
        Cb_rpe = [Cb_rpe; 0];
        Sb_rpe = [Sb_rpe; iS];
    end
    if isfield(RPE.(subj), 'CS60')
        Yb_rpe = [Yb_rpe; RPE.(subj).CS60.fin];
        Xb_rpe = [Xb_rpe; RPE.(subj).CS60.debut];
        Cb_rpe = [Cb_rpe; 1];
        Sb_rpe = [Sb_rpe; iS];
    end
end

Tb_rpe = table(Yb_rpe, Xb_rpe, categorical(Cb_rpe), categorical(Sb_rpe), ...
    'VariableNames', {'RPE_fin','RPE_debut','Condition','Participant'});
fprintf('  N_obs=%d\n', height(Tb_rpe));

try
    lmeB_reduced = fitlme(Tb_rpe, 'RPE_fin ~ Condition + (1|Participant)', 'FitMethod','ML');
    lmeB_full    = fitlme(Tb_rpe, 'RPE_fin ~ RPE_debut + Condition + (1|Participant)', 'FitMethod','ML');

    lrtB   = compare(lmeB_reduced, lmeB_full);
    p_LRTB = lrtB.pValue(2);

    [~,~,FEB]      = fixedEffects(lmeB_full, 'DFMethod','satterthwaite');
    [R2m_B, R2c_B] = computeNakagawaR2(lmeB_full);

    fprintf('  LRT(RPE_debut, 1 df) p=%.4f%s | R²m=%.3f R²c=%.3f\n', p_LRTB, p2star(p_LRTB), R2m_B, R2c_B);
    for ir = 1:height(FEB)
        fprintf('    %-35s  β=%+7.4f  p=%.4f%s\n', ...
            FEB.Name{ir}, FEB.Estimate(ir), FEB.pValue(ir), p2star(FEB.pValue(ir)));
    end

    LMM_RPE0_Influence_Results.lme_reduced = lmeB_reduced;
    LMM_RPE0_Influence_Results.lme_full    = lmeB_full;
    LMM_RPE0_Influence_Results.FE          = FEB;
    LMM_RPE0_Influence_Results.p_LRT       = p_LRTB;
    LMM_RPE0_Influence_Results.R2m         = R2m_B;
    LMM_RPE0_Influence_Results.R2c         = R2c_B;
    LMM_RPE0_Influence_Results.N_obs       = height(Tb_rpe);

catch ME
    fprintf('  LMM RPE_fin~RPE_debut échoué : %s\n', ME.message);
    LMM_RPE0_Influence_Results = struct();
end

%% -----------------------------------------------------------------------
%  MODÈLE C — RPE0 PRÉDIT-IL LA PENTE (ΔRPE), PAS SEULEMENT LE NIVEAU FINAL ?
%  -----------------------------------------------------------------------
%  ΔRPE ~ RPE_début + Condition + (1|Participant)
%
%  Distinction avec le Modèle B : RPE_fin ~ RPE_début (β1) ne dit PAS si
%  un sujet qui part avec un RPE0 élevé augmente AUSSI plus que les autres
%  (pente) ou s'il reste juste décalé vers le haut d'un montant constant
%  (pur effet de niveau). En écrivant ΔRPE = RPE_fin - RPE_début, le même
%  modèle se reformule en :
%     ΔRPE ~ (β1-1)*RPE_début + Condition + ...
%  Le test de significativité sur le coefficient de RPE_début ICI répond
%  directement à la question pente vs niveau (β1=1 dans le Modèle B
%  correspondrait à un coefficient nul ici) :
%    - coefficient non significatif -> pur effet de niveau (décalage
%      constant, RPE0 ne change pas la dynamique de fatigue elle-même).
%    - coefficient significatif -> RPE0 module aussi l'AMPLEUR de la
%      fatigue accumulée, pas seulement le point de départ.
%  Même échantillon que le Modèle B (jusqu'à N_obs=36).
%  -----------------------------------------------------------------------
fprintf('\n=== LMM Delta_RPE ~ RPE_début + Condition (pente vs niveau, N_obs<=36) — DeltaRPE ~ RPE_debut + Condition + (1|Participant) ===\n');

Yc_rpe = Yb_rpe - Xb_rpe;   % DeltaRPE = RPE_fin - RPE_debut, mêmes lignes/sujets que le Modèle B
Tc_rpe = table(Yc_rpe, Xb_rpe, categorical(Cb_rpe), categorical(Sb_rpe), ...
    'VariableNames', {'DeltaRPE','RPE_debut','Condition','Participant'});
fprintf('  N_obs=%d\n', height(Tc_rpe));

try
    lmeC_reduced = fitlme(Tc_rpe, 'DeltaRPE ~ Condition + (1|Participant)', 'FitMethod','ML');
    lmeC_full    = fitlme(Tc_rpe, 'DeltaRPE ~ RPE_debut + Condition + (1|Participant)', 'FitMethod','ML');

    lrtC   = compare(lmeC_reduced, lmeC_full);
    p_LRTC = lrtC.pValue(2);

    [~,~,FEC]      = fixedEffects(lmeC_full, 'DFMethod','satterthwaite');
    [R2m_C, R2c_C] = computeNakagawaR2(lmeC_full);

    fprintf('  LRT(RPE_debut, 1 df) p=%.4f%s | R²m=%.3f R²c=%.3f\n', p_LRTC, p2star(p_LRTC), R2m_C, R2c_C);
    for ir = 1:height(FEC)
        fprintf('    %-35s  β=%+7.4f  p=%.4f%s\n', ...
            FEC.Name{ir}, FEC.Estimate(ir), FEC.pValue(ir), p2star(FEC.pValue(ir)));
    end

    LMM_RPE0_Slope_Results.lme_reduced = lmeC_reduced;
    LMM_RPE0_Slope_Results.lme_full    = lmeC_full;
    LMM_RPE0_Slope_Results.FE          = FEC;
    LMM_RPE0_Slope_Results.p_LRT       = p_LRTC;
    LMM_RPE0_Slope_Results.R2m         = R2m_C;
    LMM_RPE0_Slope_Results.R2c         = R2c_C;
    LMM_RPE0_Slope_Results.N_obs       = height(Tc_rpe);

catch ME
    fprintf('  LMM Delta_RPE~RPE_debut échoué : %s\n', ME.message);
    LMM_RPE0_Slope_Results = struct();
end

%% -----------------------------------------------------------------------
%  DÉRIVE (intervalle 10 - intervalle 1, z-scorée) — variables individuelles
%  + index de bloc. Calculée pour TOUS les sujets en Norm (couvre à la fois
%  la série Normal N=13 crossover et la série Normal N=23 totale) et pour
%  les crossover en CS60 (N=13).
%  -----------------------------------------------------------------------
Derive_Results = struct();  % Derive_Results.(subj).(kbd).(vn) = scalaire

% --- Variables individuelles ---
for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end
        mu_z = LMM_Results.(vn).mu_z;
        sd_z = LMM_Results.(vn).sd_z;

        for iS = 1:length(Subject_All_IMU)
            subj = Subject_All_IMU{iS};
            try
                if strcmp(bloc.source,'IMU')
                    y = Results_IMU.(subj).Norm.(vn)(:);
                else
                    y = Results.(subj).Norm.(vn)(:);
                end
                z = (y - mu_z) / sd_z;
                Derive_Results.(subj).Norm.(vn) = z(end) - z(1);
            catch
                Derive_Results.(subj).Norm.(vn) = NaN;
            end
        end

        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            try
                if strcmp(bloc.source,'IMU')
                    y = Results_IMU.(subj).CS60.(vn)(:);
                else
                    y = Results.(subj).CS60.(vn)(:);
                end
                z = (y - mu_z) / sd_z;
                Derive_Results.(subj).CS60.(vn) = z(end) - z(1);
            catch
                Derive_Results.(subj).CS60.(vn) = NaN;
            end
        end
    end
end

% --- Index de bloc (réutilise Index_Results déjà calculé plus haut,
%     variante observed_np, qui couvre déjà Normal N=23 + CS60 N=13) ---
index_variant_used = 'observed_np';
fprintf('  Variante d''index utilisée pour la dérive : %s\n', index_variant_used);
for ib = 1:length(blocs)
    bname = blocs(ib).name;
    for iS = 1:length(Subject_All_IMU)
        subj = Subject_All_IMU{iS};
        try
            v = Index_Results.(subj).Norm.(bname).(index_variant_used);
            Derive_Results.(subj).Norm.(['Index_' bname]) = v(end) - v(1);
        catch
            Derive_Results.(subj).Norm.(['Index_' bname]) = NaN;
        end
    end
    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        try
            v = Index_Results.(subj).CS60.(bname).(index_variant_used);
            Derive_Results.(subj).CS60.(['Index_' bname]) = v(end) - v(1);
        catch
            Derive_Results.(subj).CS60.(['Index_' bname]) = NaN;
        end
    end
end

%% -----------------------------------------------------------------------
%  CORRÉLATIONS PEARSON : Dérive <-> ΔRPE, pour les 3 séries
%  (Normal N=13 crossover | CS60 N=13 | Normal N=23 total)
%  -----------------------------------------------------------------------
fprintf('\n=== Corrélations Pearson : Dérive <-> ΔRPE (4 séries) ===\n');

RPE_Scatter_Corr = struct();

for ib = 1:length(blocs)
    bname    = blocs(ib).name;
    vars_fig = [blocs(ib).vars, {['Index_' bname]}];

    for iV = 1:length(vars_fig)
        vn = vars_fig{iV};

        % --- Normal N=13 (crossover) ---
        x13 = []; y13 = [];
        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),'Norm') ...
                    && isfield(Derive_Results.(subj).Norm,vn) && isfield(RPE,subj)
                d = Derive_Results.(subj).Norm.(vn);
                if ~isnan(d), x13(end+1) = d; y13(end+1) = RPE.(subj).Norm.delta; end %#ok<AGROW>
            end
        end

        % --- CS60 N=13 ---
        xcs = []; ycs = [];
        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),'CS60') ...
                    && isfield(Derive_Results.(subj).CS60,vn) && isfield(RPE,subj) && isfield(RPE.(subj),'CS60')
                d = Derive_Results.(subj).CS60.(vn);
                if ~isnan(d), xcs(end+1) = d; ycs(end+1) = RPE.(subj).CS60.delta; end %#ok<AGROW>
            end
        end

        % --- Normal N=23 (tous sujets) ---
        x23 = []; y23 = [];
        for iS = 1:length(Subject_All_IMU)
            subj = Subject_All_IMU{iS};
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),'Norm') ...
                    && isfield(Derive_Results.(subj).Norm,vn) && isfield(RPE,subj)
                d = Derive_Results.(subj).Norm.(vn);
                if ~isnan(d), x23(end+1) = d; y23(end+1) = RPE.(subj).Norm.delta; end %#ok<AGROW>
            end
        end

        % --- Normal LD (grandes mains, N=10) ---
        xld = []; yld = [];
        for iS = 1:length(Subject_LD)
            subj = Subject_LD{iS};
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),'Norm') ...
                    && isfield(Derive_Results.(subj).Norm,vn) && isfield(RPE,subj)
                d = Derive_Results.(subj).Norm.(vn);
                if ~isnan(d), xld(end+1) = d; yld(end+1) = RPE.(subj).Norm.delta; end %#ok<AGROW>
            end
        end

        r13=NaN; p13=NaN; rcs=NaN; pcs=NaN; r23=NaN; p23=NaN; rld=NaN; pld=NaN;
        if length(x13) >= 4, [r13,p13] = corr(x13(:), y13(:), 'Type','Pearson'); end
        if length(xcs) >= 4, [rcs,pcs] = corr(xcs(:), ycs(:), 'Type','Pearson'); end
        if length(x23) >= 4, [r23,p23] = corr(x23(:), y23(:), 'Type','Pearson'); end
        if length(xld) >= 4, [rld,pld] = corr(xld(:), yld(:), 'Type','Pearson'); end

        RPE_Scatter_Corr.(vn).x13=x13; RPE_Scatter_Corr.(vn).y13=y13; RPE_Scatter_Corr.(vn).r13=r13; RPE_Scatter_Corr.(vn).p13=p13;
        RPE_Scatter_Corr.(vn).xcs=xcs; RPE_Scatter_Corr.(vn).ycs=ycs; RPE_Scatter_Corr.(vn).rcs=rcs; RPE_Scatter_Corr.(vn).pcs=pcs;
        RPE_Scatter_Corr.(vn).x23=x23; RPE_Scatter_Corr.(vn).y23=y23; RPE_Scatter_Corr.(vn).r23=r23; RPE_Scatter_Corr.(vn).p23=p23;
        RPE_Scatter_Corr.(vn).xld=xld; RPE_Scatter_Corr.(vn).yld=yld; RPE_Scatter_Corr.(vn).rld=rld; RPE_Scatter_Corr.(vn).pld=pld;

        fprintf('  %-30s  Norm13 r=%+.2f(%s) | CS60 r=%+.2f(%s) | Norm23 r=%+.2f(%s) | LD r=%+.2f(%s)\n', ...
            vn, r13, p2star_full(p13), rcs, p2star_full(pcs), r23, p2star_full(p23), rld, p2star_full(pld));
    end
end

%% -----------------------------------------------------------------------
%  MODÈLE D — LE LIEN DÉRIVE <-> ΔRPE DÉPEND-IL DE L'ORDRE DE PASSAGE ?
%  -----------------------------------------------------------------------
%  DeltaRPE ~ Derive + Condition + Order + Derive:Order + (1|Participant)
%
%  Test le plus direct de la question posée : pas une comparaison indirecte
%  entre deux LRT séparés (effet d'ordre sur les biomarqueurs VS effet
%  d'ordre sur le RPE, faite au tour précédent), mais l'interaction
%  Derive:Order PRISE DANS LE MÊME modèle que le lien biomarqueur->RPE
%  lui-même : "la pente Dérive->ΔRPE diffère-t-elle selon l'ordre ?"
%
%  N=13 crossover (Order non défini pour les 10 LD), jusqu'à N_obs=26 par
%  variable (mêmes sujets/lignes que les Modèles A/B/C ci-dessus, mais
%  Dérive remplace RPE_début comme prédicteur). LRT à 2 df (Order +
%  Derive:Order) contre le modèle réduit Derive + Condition seuls.
%  -----------------------------------------------------------------------
fprintf('\n=== Modèle D : lien Dérive<->ΔRPE selon l''ordre (LRT, N=13 crossover) — DeltaRPE ~ Derive + Condition + Order + Derive:Order + (1|Participant) ===\n');

LMM_DeriveOrder_Results = struct();

for ib = 1:length(blocs)
    bname    = blocs(ib).name;
    vars_fig = [blocs(ib).vars, {['Index_' bname]}];

    for iV = 1:length(vars_fig)
        vn = vars_fig{iV};

        Yd = []; Xd = []; Cd = []; Od = {}; Sd = [];
        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            ord  = Subject_Order.(subj);

            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),'Norm') ...
                    && isfield(Derive_Results.(subj).Norm,vn) && isfield(RPE,subj)
                d = Derive_Results.(subj).Norm.(vn);
                if ~isnan(d)
                    Yd = [Yd; RPE.(subj).Norm.delta]; Xd = [Xd; d];
                    Cd = [Cd; 0]; Od = [Od; {ord}]; Sd = [Sd; iS];
                end
            end
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),'CS60') ...
                    && isfield(Derive_Results.(subj).CS60,vn) && isfield(RPE,subj) && isfield(RPE.(subj),'CS60')
                d = Derive_Results.(subj).CS60.(vn);
                if ~isnan(d)
                    Yd = [Yd; RPE.(subj).CS60.delta]; Xd = [Xd; d];
                    Cd = [Cd; 1]; Od = [Od; {ord}]; Sd = [Sd; iS];
                end
            end
        end

        if length(Yd) < 15
            fprintf('  %-30s : pas assez de données (N=%d)\n', vn, length(Yd));
            continue;
        end

        Td = table(Yd, Xd, categorical(Cd), categorical(Od), categorical(Sd), ...
            'VariableNames', {'DeltaRPE','Derive','Condition','Order','Participant'});

        try
            lmeD_reduced = fitlme(Td, 'DeltaRPE ~ Derive + Condition + (1|Participant)', 'FitMethod','ML');
            lmeD_full    = fitlme(Td, 'DeltaRPE ~ Derive + Condition + Order + Derive:Order + (1|Participant)', 'FitMethod','ML');

            lrtD   = compare(lmeD_reduced, lmeD_full);
            p_LRTD = lrtD.pValue(2);

            [~,~,FED] = fixedEffects(lmeD_full, 'DFMethod','satterthwaite');

            idx_deriveOrder = find(contains(FED.Name, 'Derive') & contains(FED.Name, 'Order') & contains(FED.Name, ':'));
            beta_do = NaN; p_do = NaN;
            if ~isempty(idx_deriveOrder)
                beta_do = FED.Estimate(idx_deriveOrder(1));
                p_do    = FED.pValue(idx_deriveOrder(1));
            end

            LMM_DeriveOrder_Results.(vn).FE      = FED;
            LMM_DeriveOrder_Results.(vn).p_LRT   = p_LRTD;
            LMM_DeriveOrder_Results.(vn).beta_DO = beta_do;
            LMM_DeriveOrder_Results.(vn).p_DO    = p_do;
            LMM_DeriveOrder_Results.(vn).N_obs   = height(Td);

            fprintf('  %-30s  LRT(Order, 2 df) p=%.4f%s | Derive:Order β=%+7.4f p=%.4f%s | N=%d\n', ...
                vn, p_LRTD, p2star(p_LRTD), beta_do, p_do, p2star(p_do), height(Td));

        catch ME
            fprintf('  %-30s : modèle échoué — %s\n', vn, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
%  FIGURES SCATTER — Dérive <-> ΔRPE, par bloc (4 séries)
%  -----------------------------------------------------------------------
fprintf('\n=== Génération des figures scatter dérive <-> RPE ===\n');

for ib = 1:length(blocs)
    bname    = blocs(ib).name;
    vars_fig = [blocs(ib).vars, {['Index_' bname]}];
    N_v      = length(vars_fig);

    fig = figure('Name', sprintf('RPE_Scatter_%s', bname), 'NumberTitle','off', ...
        'Position',[30 30 380*N_v 430], 'Color','white');

    for iV = 1:N_v
        vn = vars_fig{iV};
        C  = RPE_Scatter_Corr.(vn);

        subplot(1, N_v, iV); hold on;

        scatter(C.x13, C.y13, 28, col_norm,  'o', 'filled');
        scatter(C.xcs, C.ycs, 28, col_ergo,  's', 'filled');
        scatter(C.x23, C.y23, 28, col_all23, '^', 'filled');
        scatter(C.xld, C.yld, 28, col_ld,    'o', 'filled');

        if length(C.x13) >= 2
            pf = polyfit(C.x13, C.y13, 1);
            xr = [min(C.x13) max(C.x13)];
            plot(xr, polyval(pf,xr), '-', 'Color', col_norm, 'LineWidth', 1.0);
        end
        if length(C.xcs) >= 2
            pf = polyfit(C.xcs, C.ycs, 1);
            xr = [min(C.xcs) max(C.xcs)];
            plot(xr, polyval(pf,xr), '-', 'Color', col_ergo, 'LineWidth', 1.0);
        end
        if length(C.x23) >= 2
            pf = polyfit(C.x23, C.y23, 1);
            xr = [min(C.x23) max(C.x23)];
            plot(xr, polyval(pf,xr), '-', 'Color', col_all23, 'LineWidth', 1.0);
        end
        if length(C.xld) >= 2
            pf = polyfit(C.xld, C.yld, 1);
            xr = [min(C.xld) max(C.xld)];
            plot(xr, polyval(pf,xr), '-', 'Color', col_ld, 'LineWidth', 1.0);
        end

        title(strrep(vn,'_',' '), 'FontSize',9, 'FontWeight','bold');
        xlabel('Dérive (\Delta z-score)', 'FontSize',8);
        ylabel('\Delta RPE', 'FontSize',8);

        txt = sprintf('Norm13 r=%+.2f(%s) | CS60 r=%+.2f(%s)\nNorm23 r=%+.2f(%s) | LD r=%+.2f(%s)', ...
            C.r13, p2star_full(C.p13), C.rcs, p2star_full(C.pcs), C.r23, p2star_full(C.p23), C.rld, p2star_full(C.pld));
        text(0.02, 0.02, txt, 'Units','normalized', 'FontSize',6, ...
            'Color',[0.3 0.3 0.3], 'VerticalAlignment','bottom');

        if iV == 1
            legend({'Normal (N=13)','CS60 (N=13)','Normal (N=23)','Normal LD (N=10)'}, ...
                'Location','best','FontSize',6.5);
        end
        grid on; box on;
    end

    sgtitle(sprintf('Corrélation dérive %s \\leftrightarrow \\Delta RPE — Normal (N=13) vs CS60 (N=13) vs Normal (N=23) vs Normal LD (N=10)', bname), ...
        'FontWeight','bold', 'FontSize',12);

    fname = fullfile(path_save, sprintf('Fig_RPE_Scatter_%s.png', bname));
    saveas(fig, fname);
    fprintf('  Fig_RPE_Scatter_%s.png sauvegardée\n', bname);
    close(fig);
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'LMM_Validation_Results.mat'), ...
    'LMM_Results', 'LMM_Index', 'Index_Results', 'RPE', 'Derive_Results', 'RPE_Scatter_Corr', ...
    'LMM_Order_Results', 'LMM_RPE_Order_Results', 'LMM_RPE0_Order_Results', 'LMM_RPE0_Influence_Results', ...
    'LMM_RPE0_Slope_Results', 'LMM_DeriveOrder_Results', 'LMM_SDvsLD_Results', 'Slope_SDvsLD_vs_G1');
fprintf('\nRésultats sauvegardés : LMM_Validation_Results.mat\n');

%% -----------------------------------------------------------------------
%  FIGURE TABLEAU RÉCAPITULATIF — INDEX DE BLOC (4 variantes)
%  -----------------------------------------------------------------------
fprintf('Génération figure tableau récapitulatif index...\n');

% Collecter les lignes : une par bloc x variante
idx_rows = {};
for ib = 1:length(idx_blocs)
    bname = idx_blocs(ib).name;
    for iType = 1:length(variant_types)
        type_name = variant_types{iType};
        type_lbl  = variant_lbls{iType};
        if ~isfield(LMM_Index, bname) || ~isfield(LMM_Index.(bname), type_name)
            continue;
        end
        FE  = LMM_Index.(bname).(type_name).FE;
        R2m = LMM_Index.(bname).(type_name).R2m;
        R2c = LMM_Index.(bname).(type_name).R2c;

        idx_t = find(strcmp(FE.Name, 'Temps'));
        idx_c = find(contains(FE.Name,'Condition_') & ~contains(FE.Name,'Temps'));
        idx_i = find(contains(FE.Name,'Temps:Condition'));

        beta_t=NaN; p_t=NaN; beta_c=NaN; p_c=NaN; beta_i=NaN; p_i=NaN;
        if ~isempty(idx_t), beta_t=FE.Estimate(idx_t); p_t=FE.pValue(idx_t); end
        if ~isempty(idx_c), beta_c=FE.Estimate(idx_c); p_c=FE.pValue(idx_c); end
        if ~isempty(idx_i), beta_i=FE.Estimate(idx_i); p_i=FE.pValue(idx_i); end

        idx_rows{end+1} = {bname, type_lbl, beta_t, p_t, beta_c, p_c, beta_i, p_i, R2m, R2c};
    end
end

N_rows = length(idx_rows);
headers = {'Bloc','Variante','β Temps','p Temps','β Cond','p Cond','β Inter','p Inter','R²m','R²c'};
N_cols  = length(headers);

fig_w = 1400; fig_h = 80 + N_rows*32 + 70;
fig_t2 = figure('Name','Tableau_LMM_Index','NumberTitle','off',...
    'Position',[50 50 fig_w fig_h],'Color','white');
ax = axes('Position',[0 0 1 1],'Visible','off','XLim',[0 1],'YLim',[0 1]);
hold on;

try

col_hdr  = [0.15 0.25 0.45];
col_seg  = [0.85 0.92 0.98];
col_goub = [0.88 0.97 0.88];
col_emg  = [0.97 0.93 0.82];
col_blocs_idx = struct('SegMod',col_seg,'Goubault',col_goub,'EMG',col_emg);

top   = 0.94;
row_h = 0.80 / (N_rows + 1);
col_x = [0.01 0.14 0.28 0.37 0.46 0.55 0.64 0.73 0.84 0.92];
col_w = [0.12 0.14 0.08 0.08 0.08 0.08 0.08 0.08 0.07 0.07];

% Titre
rectangle('Position',[0 top+0.01 1 0.06],'FaceColor',col_hdr,'EdgeColor','none');
text(0.5, top+0.04, 'Tableau récapitulatif — Index de bloc (signes Goubault/Observés × Non pondéré/β/taux/β×taux, N=13)', ...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',10,'FontWeight','bold','Color','white');

% En-têtes
hdr_y = top - row_h*0.5;
rectangle('Position',[0 top-row_h 1 row_h],'FaceColor',[0.75 0.80 0.88],'EdgeColor','none');
for ic = 1:N_cols
    text(col_x(ic)+col_w(ic)/2, hdr_y, headers{ic}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',8,'FontWeight','bold','Color',[0.1 0.1 0.3]);
end
line([0 1],[top-row_h top-row_h],'Color',[0.5 0.5 0.6],'LineWidth',1);

% Lignes
prev_bloc = '';
for ir = 1:N_rows
    r     = idx_rows{ir};
    y_row = top - row_h*(ir+1);
    y_txt = y_row + row_h*0.5;

    bn = r{1};
    bg = [0.97 0.97 0.97];
    if isfield(col_blocs_idx, bn)
        bg = col_blocs_idx.(bn);
    end
    if mod(ir,2)==0, bg = bg * 0.96; end
    rectangle('Position',[0 y_row 1 row_h],'FaceColor',bg,'EdgeColor','none');

    % Colonne Bloc — affichée seulement au changement de bloc
    if ~strcmp(bn, prev_bloc)
        text(col_x(1)+0.005, y_txt, bn, ...
            'VerticalAlignment','middle','FontSize',7.5,'Color',[0.1 0.1 0.1],'FontWeight','bold');
        prev_bloc = bn;
        line([0 1],[y_row+row_h y_row+row_h],'Color',[0.4 0.4 0.5],'LineWidth',0.8);
    end

    % Colonne Variante
    text(col_x(2)+0.005, y_txt, r{2}, ...
        'VerticalAlignment','middle','FontSize',6.5,'Color',[0.2 0.2 0.2]);

    % β et p : Temps, Cond, Inter
    pairs   = {r{3},r{4}; r{5},r{6}; r{7},r{8}};
    col_idx_loc = [3 5 7];
    for ip = 1:3
        beta_v = pairs{ip,1};
        p_v    = pairs{ip,2};
        ic     = col_idx_loc(ip);

        if ~isnan(beta_v)
            col_b = [0.1 0.1 0.1];
            if beta_v > 0, col_b = [0.10 0.40 0.10]; end
            if beta_v < 0, col_b = [0.55 0.10 0.10]; end
            fw = 'normal';
            if ~isnan(p_v) && p_v < 0.05, fw = 'bold'; end
            text(col_x(ic)+col_w(ic)/2, y_txt, sprintf('%+.4f', beta_v), ...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontSize',7.5,'Color',col_b,'FontWeight',fw);
        end

        if ~isnan(p_v)
            stars = p2star(p_v);
            if ~isempty(stars)
                p_str = sprintf('%.3f%s', p_v, stars);
                p_col = [0.7 0.1 0.1];
            else
                p_str = sprintf('%.3f', p_v);
                p_col = [0.5 0.5 0.5];
            end
            text(col_x(ic+1)+col_w(ic+1)/2, y_txt, p_str, ...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontSize',7,'Color',p_col);
        end
    end

    % R²m / R²c
    R2m_v = r{9}; R2c_v = r{10};
    if ~isnan(R2m_v)
        text(col_x(9)+col_w(9)/2, y_txt, sprintf('%.3f', R2m_v), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7.5);
    end
    if ~isnan(R2c_v)
        text(col_x(10)+col_w(10)/2, y_txt, sprintf('%.3f', R2c_v), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7.5);
    end

    line([0 1],[y_row y_row],'Color',[0.85 0.85 0.88],'LineWidth',0.4);
end

line([0 1],[top-row_h*(N_rows+1) top-row_h*(N_rows+1)],'Color',[0.4 0.4 0.5],'LineWidth',1);

% Légende
leg_y = top - row_h*(N_rows+1) - 0.05;
text(0.01, leg_y, 'β : ', 'FontSize',7,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle');
text(0.05, leg_y, '■ positif', 'FontSize',7,'Color',[0.10 0.40 0.10],'FontWeight','bold','VerticalAlignment','middle');
text(0.14, leg_y, '■ négatif', 'FontSize',7,'Color',[0.55 0.10 0.10],'FontWeight','bold','VerticalAlignment','middle');
text(0.24, leg_y, '| p en gras si significatif (p<0.05)', 'FontSize',7,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle');
text(0.62, leg_y, '| Fond : ', 'FontSize',7,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle');
text(0.68, leg_y, '■ SegMod', 'FontSize',7,'Color',col_seg*0.6,'FontWeight','bold','VerticalAlignment','middle');
text(0.77, leg_y, '■ Goubault', 'FontSize',7,'Color',col_goub*0.5,'FontWeight','bold','VerticalAlignment','middle');
text(0.87, leg_y, '■ EMG', 'FontSize',7,'Color',col_emg*0.6,'FontWeight','bold','VerticalAlignment','middle');

hold off;
fname = fullfile(path_save, 'Fig_Tableau_LMM_Index.png');
exportgraphics(fig_t2, fname, 'BackgroundColor','white', 'Resolution',200);
fprintf('  Fig_Tableau_LMM_Index.png sauvegardée avec succès\n');

catch ME_tab
    fprintf('  ERREUR lors de la génération du tableau index : %s\n', ME_tab.message);
    fprintf('  Ligne : %d\n', ME_tab.stack(1).line);
end

close(fig_t2);

%% -----------------------------------------------------------------------
%  FIGURE TABLEAU RÉCAPITULATIF
%  -----------------------------------------------------------------------
fprintf('Génération figure tableau récapitulatif...\n');

% Collecter toutes les lignes
rows = {};
for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end
        FE  = LMM_Results.(vn).FE;
        R2m = LMM_Results.(vn).R2m;
        R2c = LMM_Results.(vn).R2c;

        idx_t = find(strcmp(FE.Name, 'Temps'));
        idx_c = find(contains(FE.Name,'Condition_') & ~contains(FE.Name,'Temps'));
        idx_i = find(contains(FE.Name,'Temps:Condition'));

        beta_t=NaN; p_t=NaN; beta_c=NaN; p_c=NaN; beta_i=NaN; p_i=NaN;
        if ~isempty(idx_t), beta_t=FE.Estimate(idx_t); p_t=FE.pValue(idx_t); end
        if ~isempty(idx_c), beta_c=FE.Estimate(idx_c); p_c=FE.pValue(idx_c); end
        if ~isempty(idx_i), beta_i=FE.Estimate(idx_i); p_i=FE.pValue(idx_i); end

        % Nom court
        vn_short = strrep(vn,'_',' ');
        vn_short = strrep(vn_short,'TFR ','');
        vn_short = strrep(vn_short,'SpectralEntropy','SpEnt');
        vn_short = strrep(vn_short,'SampleEntropy','SampEnt');
        vn_short = strrep(vn_short,'MedianFreq','MedFreq');
        vn_short = strrep(vn_short,'PeakPower','PkPow');
        vn_short = strrep(vn_short,'AngVel','AngV');
        vn_short = strrep(vn_short,'Accel','Acc');
        vn_short = strrep(vn_short,'Forearm','Fore');

        rows{end+1} = {vn_short, bloc.name, ...
            beta_t, p_t, beta_c, p_c, beta_i, p_i, R2m, R2c};
    end
end

N_rows = length(rows);
headers = {'Variable','Bloc','β Temps','p Temps','β Cond','p Cond','β Inter','p Inter','R²m','R²c'};
N_cols  = length(headers);

% Dimensions figure
fig_w = 1400; fig_h = 80 + N_rows*32 + 60;
fig_t = figure('Name','Tableau_LMM','NumberTitle','off',...
    'Position',[50 50 fig_w fig_h],'Color','white');
ax = axes('Position',[0 0 1 1],'Visible','off','XLim',[0 1],'YLim',[0 1]);
hold on;

% Couleurs
col_hdr  = [0.15 0.25 0.45];
col_seg  = [0.85 0.92 0.98];
col_goub = [0.88 0.97 0.88];
col_emg  = [0.97 0.93 0.82];
col_blocs = struct('SegMod',col_seg,'Goubault',col_goub,'EMG',col_emg);

% Géométrie
top   = 0.94;
row_h = 0.80 / (N_rows + 1);
col_x = [0.01 0.22 0.31 0.40 0.49 0.58 0.67 0.76 0.85 0.92];
col_w = [0.20 0.08 0.08 0.08 0.08 0.08 0.08 0.08 0.07 0.07];

% Titre
rectangle('Position',[0 top+0.01 1 0.06],'FaceColor',col_hdr,'EdgeColor','none');
text(0.5, top+0.04, 'Tableau récapitulatif LMM — Validation externe (N=13 crossover)', ...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',11,'FontWeight','bold','Color','white');

% En-têtes colonnes
hdr_y = top - row_h*0.5;
rectangle('Position',[0 top-row_h 1 row_h],'FaceColor',[0.75 0.80 0.88],'EdgeColor','none');
for ic = 1:N_cols
    text(col_x(ic)+col_w(ic)/2, hdr_y, headers{ic}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',8,'FontWeight','bold','Color',[0.1 0.1 0.3]);
end
line([0 1],[top-row_h top-row_h],'Color',[0.5 0.5 0.6],'LineWidth',1);

% Lignes de données
for ir = 1:N_rows
    r     = rows{ir};
    y_row = top - row_h*(ir+1);
    y_txt = y_row + row_h*0.5;

    % Couleur de fond par bloc
    bn = r{2};
    if isfield(col_blocs, bn)
        bg = col_blocs.(bn);
    else
        bg = [0.97 0.97 0.97];
    end
    if mod(ir,2)==0, bg = bg * 0.96; end
    rectangle('Position',[0 y_row 1 row_h],'FaceColor',bg,'EdgeColor','none');

    % Colonne 1 : Variable
    text(col_x(1)+0.005, y_txt, r{1}, ...
        'VerticalAlignment','middle','FontSize',7.5,'Color',[0.1 0.1 0.1],'FontWeight','bold');

    % Colonne 2 : Bloc
    text(col_x(2)+col_w(2)/2, y_txt, r{2}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7,'Color',[0.3 0.3 0.3]);

    % Colonnes β et p : Temps, Condition, Interaction
    pairs = {r{3},r{4}; r{5},r{6}; r{7},r{8}};
    col_idx = [3 5 7];
    for ip = 1:3
        beta_v = pairs{ip,1};
        p_v    = pairs{ip,2};
        ic     = col_idx(ip);

        % β
        if ~isnan(beta_v)
            col_b = [0.1 0.1 0.1];
            if beta_v > 0, col_b = [0.10 0.40 0.10]; end
            if beta_v < 0, col_b = [0.55 0.10 0.10]; end
            fw = 'normal';
            if ~isnan(p_v) && p_v < 0.05, fw = 'bold'; end
            text(col_x(ic)+col_w(ic)/2, y_txt, sprintf('%+.4f', beta_v), ...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontSize',7.5,'Color',col_b,'FontWeight',fw);
        end

        % p
        if ~isnan(p_v)
            stars = p2star(p_v);
            if ~isempty(stars)
                p_str = sprintf('%.3f%s', p_v, stars);
                p_col = [0.7 0.1 0.1];
            else
                p_str = sprintf('%.3f', p_v);
                p_col = [0.5 0.5 0.5];
            end
            text(col_x(ic+1)+col_w(ic+1)/2, y_txt, p_str, ...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontSize',7,'Color',p_col);
        end
    end

    % R²m et R²c
    R2m_v = r{9}; R2c_v = r{10};
    if ~isnan(R2m_v)
        text(col_x(9)+col_w(9)/2, y_txt, sprintf('%.3f', R2m_v), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7.5);
    end
    if ~isnan(R2c_v)
        text(col_x(10)+col_w(10)/2, y_txt, sprintf('%.3f', R2c_v), ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7.5);
    end

    line([0 1],[y_row y_row],'Color',[0.85 0.85 0.88],'LineWidth',0.4);
end

% Ligne finale
line([0 1],[top-row_h*(N_rows+1) top-row_h*(N_rows+1)],'Color',[0.4 0.4 0.5],'LineWidth',1);

% Légende couleurs β
leg_y = top - row_h*(N_rows+1) - 0.04;
text(0.01, leg_y, 'β : ', 'FontSize',7,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle');
text(0.05, leg_y, '■ positif', 'FontSize',7,'Color',[0.10 0.40 0.10],'FontWeight','bold','VerticalAlignment','middle');
text(0.14, leg_y, '■ négatif', 'FontSize',7,'Color',[0.55 0.10 0.10],'FontWeight','bold','VerticalAlignment','middle');
text(0.24, leg_y, '| p en gras si significatif (p<0.05)', 'FontSize',7,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle');
text(0.60, leg_y, '| Fond : ', 'FontSize',7,'Color',[0.3 0.3 0.3],'VerticalAlignment','middle');
text(0.65, leg_y, '■ SegMod', 'FontSize',7,'Color',col_seg*0.6,'FontWeight','bold','VerticalAlignment','middle');
text(0.73, leg_y, '■ Goubault', 'FontSize',7,'Color',col_goub*0.5,'FontWeight','bold','VerticalAlignment','middle');
text(0.82, leg_y, '■ EMG', 'FontSize',7,'Color',col_emg*0.6,'FontWeight','bold','VerticalAlignment','middle');

hold off;
fname = fullfile(path_save, 'Fig_Tableau_LMM.png');
exportgraphics(fig_t, fname, 'BackgroundColor','white', 'Resolution',200);
close(fig_t);
fprintf('  Fig_Tableau_LMM.png sauvegardée\n');

fprintf('\n=== TERMINÉ ===\n');

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------
function sig = p2star(p)
    if isnan(p),       sig = '';
    elseif p < 0.001,  sig = '***';
    elseif p < 0.01,   sig = '**';
    elseif p < 0.05,   sig = '*';
    else,              sig = '';
    end
end

function str = p2star_full(p)
    if isnan(p)
        str = 'p=N/A';
    else
        str = sprintf('p=%.3f%s', p, p2star(p));
    end
end

function [R2m, R2c] = computeNakagawaR2(lme)
    try
        % Variance des effets fixes
        fe   = fixedEffects(lme);
        X    = designMatrix(lme, 'Fixed');
        varF = var(X * fe, 1, 'omitnan');

        % Variance des effets aléatoires
        vc   = lme.covarianceParameters;
        varR = 0;
        for ir = 1:length(vc)
            v = vc{ir};
            if isnumeric(v) && numel(v) == 1
                varR = varR + v;
            elseif isnumeric(v)
                varR = varR + sum(diag(v));
            end
        end

        % Variance résiduelle
        sigma2 = lme.MSE;
        denom  = varF + varR + sigma2;

        if denom <= 0 || isnan(denom)
            R2m = NaN; R2c = NaN;
        else
            R2m = varF / denom;
            R2c = (varF + varR) / denom;
        end
    catch
        R2m = NaN; R2c = NaN;
    end
end