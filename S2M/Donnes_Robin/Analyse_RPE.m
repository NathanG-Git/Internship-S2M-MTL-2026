%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  ANALYSE RPE — CORRÉLATION AVEC DÉRIVE DES BIOMARQUEURS            %%%%
%%%%  Piano Normal vs Piano Ergonomique (CS60)                          %%%%
%%%%                                                                     %%%%
%%%%  Pour chaque variable et chaque index de bloc, par condition :    %%%%
%%%%    Corrélation de Pearson  : Dérive ↔ ΔRPE                       %%%%
%%%%    Corrélation de Spearman : Dérive ↔ ΔRPE                       %%%%
%%%%                                                                     %%%%
%%%%  Dérive = valeur intervalle 10 − valeur intervalle 1 (z-scorée)   %%%%
%%%%  ΔRPE   = RPE_fin − RPE_début                                     %%%%
%%%%                                                                     %%%%
%%%%  Pas d'ajustement pour RPE_début à ce stade (à voir ultérieurement)%%%%
%%%%  P07 exclu entièrement (Normal + CS60), comme pour le reste du    %%%%
%%%%  pipeline de validation                                            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_save = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

%% -----------------------------------------------------------------------
%  CHARGEMENT DES RÉSULTATS DE VALIDATION DÉJÀ CALCULÉS
%  -----------------------------------------------------------------------
fprintf('Chargement des résultats de validation...\n');
load(fullfile(path_data, 'LMM_Validation_Results.mat'));  % LMM_Results, LMM_Index, Index_Results
load(fullfile(path_data, 'IMU_Features_NewDataset.mat')); % Results_IMU
load(fullfile(path_data, 'EMG_Features_NewDataset.mat')); % Results (EMG)
fprintf('  OK\n\n');

% Exclusion P07 (cohérence avec le reste du pipeline)
if isfield(Results_IMU, 'P07'), Results_IMU = rmfield(Results_IMU, 'P07'); end
if isfield(Results, 'P07'),     Results     = rmfield(Results, 'P07');     end

%% -----------------------------------------------------------------------
%  DONNÉES RPE — saisies manuellement, P07 exclu
%  -----------------------------------------------------------------------
% Format : {Sujet, RPE_debut_Norm, RPE_fin_Norm, RPE_debut_CS60, RPE_fin_CS60}
% NaN si le sujet n'a pas la condition CS60
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

N_subj_rpe = size(RPE_data, 1);
RPE = struct();
for i = 1:N_subj_rpe
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

fprintf('RPE chargé pour %d sujets (P07 exclu)\n\n', N_subj_rpe);

%% -----------------------------------------------------------------------
%  SUJETS CROSSOVER (Normal + CS60), P07 exclu
%  -----------------------------------------------------------------------
Subject_Both = {'P02','P03','P04','P05','P08','P10','P12',...
                'P15','P17','P18','P20','P21','P24'};
N_bins = 10;

%% -----------------------------------------------------------------------
%  DÉFINITION DES BLOCS (variables + index)
%  -----------------------------------------------------------------------
blocs = struct();
blocs(1).name   = 'SegMod';
blocs(1).vars   = {'Accel_Mod_Head', 'Accel_Mod_Hand'};
blocs(1).source = 'IMU';

blocs(2).name   = 'Goubault';
blocs(2).vars   = {'MedianFreq_Accel_Y_Hand', 'PeakPower_Accel_Mod_Hand', ...
                   'PeakPower_AngVel_X_Head', 'SpectralEntropy_AngVel_Mod_Forearm'};
blocs(2).source = 'IMU';

blocs(3).name   = 'EMG';
blocs(3).vars   = {'TFR_MedianFreq_Triceps', 'TFR_SpectralEntropy_Deltoid', ...
                   'SampleEntropy_Biceps'};
blocs(3).source = 'EMG';

%% -----------------------------------------------------------------------
%  CALCUL DE LA DÉRIVE (intervalle 10 − intervalle 1, z-scorée) PAR VARIABLE
%  -----------------------------------------------------------------------
fprintf('=== Calcul des dérives (Δ z-score, intervalle 10 - intervalle 1) ===\n');

all_vars = {};
for ib = 1:length(blocs)
    all_vars = [all_vars, blocs(ib).vars];
end

Derive_Results = struct();  % Derive_Results.(subj).(kbd).(vn) = scalaire

for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end
        mu_z = LMM_Results.(vn).mu_z;
        sd_z = LMM_Results.(vn).sd_z;

        for iS = 1:length(Subject_Both)
            subj = Subject_Both{iS};
            kbds = {'Norm','CS60'};
            for iK = 1:length(kbds)
                kbd = kbds{iK};
                try
                    if strcmp(bloc.source,'IMU')
                        y = Results_IMU.(subj).(kbd).(vn)(:);
                    else
                        y = Results.(subj).(kbd).(vn)(:);
                    end
                    z = (y - mu_z) / sd_z;
                    Derive_Results.(subj).(kbd).(vn) = z(end) - z(1);
                catch
                    Derive_Results.(subj).(kbd).(vn) = NaN;
                end
            end
        end
    end
end

%% -----------------------------------------------------------------------
%  DÉRIVE DES INDEX DE BLOC (utilise la variante observed_np par défaut)
%  -----------------------------------------------------------------------
index_variant_used = 'observed_np';  % variante de référence pour le RPE
fprintf('Variante d''index utilisée pour la dérive : %s\n\n', index_variant_used);

for ib = 1:length(blocs)
    bname = blocs(ib).name;
    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        kbds = {'Norm','CS60'};
        for iK = 1:length(kbds)
            kbd = kbds{iK};
            try
                v = Index_Results.(subj).(kbd).(bname).(index_variant_used);
                Derive_Results.(subj).(kbd).(['Index_' bname]) = v(end) - v(1);
            catch
                Derive_Results.(subj).(kbd).(['Index_' bname]) = NaN;
            end
        end
    end
end

%% -----------------------------------------------------------------------
%  CORRÉLATIONS PEARSON ET SPEARMAN — DÉRIVE vs ΔRPE, PAR CONDITION
%  -----------------------------------------------------------------------
fprintf('=== Corrélations Pearson et Spearman : Dérive ↔ ΔRPE ===\n\n');

vars_to_test = [all_vars, {'Index_SegMod','Index_Goubault','Index_EMG'}];

RPE_Corr_Results = struct();

for iV = 1:length(vars_to_test)
    vn = vars_to_test{iV};

    Y_drpe = []; X_drv = []; X_cond = [];

    for iS = 1:length(Subject_Both)
        subj = Subject_Both{iS};
        if ~isfield(RPE, subj), continue; end

        % Normal
        if isfield(RPE.(subj),'Norm') && isfield(Derive_Results,subj) ...
                && isfield(Derive_Results.(subj),'Norm') && isfield(Derive_Results.(subj).Norm, vn)
            drv = Derive_Results.(subj).Norm.(vn);
            if ~isnan(drv) && ~isnan(RPE.(subj).Norm.delta)
                Y_drpe = [Y_drpe; RPE.(subj).Norm.delta];
                X_drv  = [X_drv; drv];
                X_cond = [X_cond; 0];
            end
        end

        % CS60
        if isfield(RPE.(subj),'CS60') && isfield(Derive_Results,subj) ...
                && isfield(Derive_Results.(subj),'CS60') && isfield(Derive_Results.(subj).CS60, vn)
            drv = Derive_Results.(subj).CS60.(vn);
            if ~isnan(drv) && ~isnan(RPE.(subj).CS60.delta)
                Y_drpe = [Y_drpe; RPE.(subj).CS60.delta];
                X_drv  = [X_drv; drv];
                X_cond = [X_cond; 1];
            end
        end
    end

    if length(Y_drpe) < 8
        fprintf('  %-30s : pas assez de données (N=%d)\n', vn, length(Y_drpe));
        continue;
    end

    idx_n = X_cond == 0;
    idx_e = X_cond == 1;

    % --- Pearson ---
    r_n=NaN; p_rn=NaN; r_e=NaN; p_re=NaN;
    if sum(idx_n) >= 4, [r_n, p_rn] = corr(X_drv(idx_n), Y_drpe(idx_n), 'Type','Pearson'); end
    if sum(idx_e) >= 4, [r_e, p_re] = corr(X_drv(idx_e), Y_drpe(idx_e), 'Type','Pearson'); end

    % --- Spearman ---
    rho_n=NaN; p_rhon=NaN; rho_e=NaN; p_rhoe=NaN;
    if sum(idx_n) >= 4, [rho_n, p_rhon] = corr(X_drv(idx_n), Y_drpe(idx_n), 'Type','Spearman'); end
    if sum(idx_e) >= 4, [rho_e, p_rhoe] = corr(X_drv(idx_e), Y_drpe(idx_e), 'Type','Spearman'); end

    RPE_Corr_Results.(vn).r_norm   = r_n;   RPE_Corr_Results.(vn).p_norm   = p_rn;
    RPE_Corr_Results.(vn).r_cs60   = r_e;   RPE_Corr_Results.(vn).p_cs60   = p_re;
    RPE_Corr_Results.(vn).rho_norm = rho_n; RPE_Corr_Results.(vn).prho_norm = p_rhon;
    RPE_Corr_Results.(vn).rho_cs60 = rho_e; RPE_Corr_Results.(vn).prho_cs60 = p_rhoe;
    RPE_Corr_Results.(vn).N_norm   = sum(idx_n);
    RPE_Corr_Results.(vn).N_cs60   = sum(idx_e);

    fprintf('--- %s ---\n', vn);
    fprintf('  Normal (N=%d) : Pearson r=%+.3f (p=%.3f%s) | Spearman ρ=%+.3f (p=%.3f%s)\n', ...
        sum(idx_n), r_n, p_rn, p2star(p_rn), rho_n, p_rhon, p2star(p_rhon));
    fprintf('  CS60   (N=%d) : Pearson r=%+.3f (p=%.3f%s) | Spearman ρ=%+.3f (p=%.3f%s)\n\n', ...
        sum(idx_e), r_e, p_re, p2star(p_re), rho_e, p_rhoe, p2star(p_rhoe));
end

%% -----------------------------------------------------------------------
%  FOREST PLOT — CORRÉLATIONS DÉRIVE ↔ ΔRPE, TOUTES VARIABLES
%  -----------------------------------------------------------------------
fprintf('=== Génération du forest plot des corrélations ===\n');

col_norm = [0.20 0.40 0.70];
col_ergo = [0.15 0.60 0.35];

% Ordre d'affichage : variables par bloc puis index, dans l'ordre inverse
% pour que SegMod apparaisse en haut du graphique
plot_order = {'Accel_Mod_Head','Accel_Mod_Hand','Index_SegMod', ...
              'MedianFreq_Accel_Y_Hand','PeakPower_Accel_Mod_Hand', ...
              'PeakPower_AngVel_X_Head','SpectralEntropy_AngVel_Mod_Forearm','Index_Goubault', ...
              'TFR_MedianFreq_Triceps','TFR_SpectralEntropy_Deltoid','SampleEntropy_Biceps','Index_EMG'};
plot_order = flip(plot_order);  % pour affichage de haut en bas dans l'ordre naturel

% Labels courts pour l'affichage
short_lbl = struct(...
    'Accel_Mod_Head','Accel Mod Head', ...
    'Accel_Mod_Hand','Accel Mod Hand', ...
    'Index_SegMod','Index SegMod', ...
    'MedianFreq_Accel_Y_Hand','MedFreq Accel Y Hand', ...
    'PeakPower_Accel_Mod_Hand','PkPow Accel Mod Hand', ...
    'PeakPower_AngVel_X_Head','PkPow AngV X Head', ...
    'SpectralEntropy_AngVel_Mod_Forearm','SpEnt AngV Mod Fore', ...
    'Index_Goubault','Index Goubault', ...
    'TFR_MedianFreq_Triceps','MedFreq Triceps', ...
    'TFR_SpectralEntropy_Deltoid','SpEnt Deltoid', ...
    'SampleEntropy_Biceps','SampEnt Biceps', ...
    'Index_EMG','Index EMG');

% Fisher z-transform pour IC95% des corrélations : z = atanh(r), SE = 1/sqrt(N-3)
N_v = length(plot_order);
y_pos_n = (1:N_v) + 0.15;
y_pos_e = (1:N_v) - 0.15;

fig = figure('Name','Forest_RPE_Correlations','NumberTitle','off', ...
    'Position',[30 30 1500 700], 'Color','white');

method_names = {'Pearson','Spearman'};
field_r   = {'r_norm','rho_norm'};
field_re  = {'r_cs60','rho_cs60'};

for iM = 1:2
    subplot(1,2,iM); hold on;
    mname = method_names{iM};

    for iV = 1:N_v
        vn = plot_order{iV};
        if ~isfield(RPE_Corr_Results, vn), continue; end

        r_n = RPE_Corr_Results.(vn).(field_r{iM});  N_n = RPE_Corr_Results.(vn).N_norm;
        r_e = RPE_Corr_Results.(vn).(field_re{iM}); N_e = RPE_Corr_Results.(vn).N_cs60;

        % IC95% via Fisher z (approximation valable aussi pour Spearman avec N modéré)
        if ~isnan(r_n) && N_n > 3
            z_n = atanh(max(min(r_n,0.999),-0.999)); se_n = 1/sqrt(N_n-3);
            ci_n = tanh([z_n-1.96*se_n, z_n+1.96*se_n]);
        else
            ci_n = [NaN NaN];
        end
        if ~isnan(r_e) && N_e > 3
            z_e = atanh(max(min(r_e,0.999),-0.999)); se_e = 1/sqrt(N_e-3);
            ci_e = tanh([z_e-1.96*se_e, z_e+1.96*se_e]);
        else
            ci_e = [NaN NaN];
        end

        if ~any(isnan(ci_n))
            plot(ci_n, [y_pos_n(iV) y_pos_n(iV)], '-', 'Color', col_norm, 'LineWidth', 2);
        end
        if ~any(isnan(ci_e))
            plot(ci_e, [y_pos_e(iV) y_pos_e(iV)], '-', 'Color', col_ergo, 'LineWidth', 2);
        end

        if ~isnan(r_n)
            plot(r_n, y_pos_n(iV), 'o', 'MarkerSize',8, 'MarkerFaceColor',col_norm, ...
                'MarkerEdgeColor','white', 'LineWidth',1);
        end
        if ~isnan(r_e)
            plot(r_e, y_pos_e(iV), 's', 'MarkerSize',8, 'MarkerFaceColor',col_ergo, ...
                'MarkerEdgeColor','white', 'LineWidth',1);
        end
    end

    xline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);

    yticks(1:N_v);
    ylabels = cell(N_v,1);
    for iV = 1:N_v
        ylabels{iV} = short_lbl.(plot_order{iV});
    end
    yticklabels(ylabels);
    ylim([0.5 N_v+0.5]);
    xlim([-1 1]);
    xlabel(sprintf('Corrélation %s (dérive ↔ ΔRPE)', mname), 'FontSize',9);
    title(mname, 'FontSize',11, 'FontWeight','bold');

    yline(3.5, ':', 'Color',[0.6 0.6 0.6]);
    yline(8.5, ':', 'Color',[0.6 0.6 0.6]);

    if iM == 1
        plot(NaN,NaN,'o','MarkerSize',8,'MarkerFaceColor',col_norm,'MarkerEdgeColor','white','LineWidth',1,'DisplayName','Normal');
        plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor',col_ergo,'MarkerEdgeColor','white','LineWidth',1,'DisplayName','CS60');
        legend('Location','southoutside','Orientation','horizontal','FontSize',9);
    end

    grid on; box on;
    set(gca,'FontSize',8);
end

sgtitle('Corrélations dérive ↔ ΔRPE — Normal vs CS60 (IC95% Fisher z)', ...
    'FontSize',12, 'FontWeight','bold');

fname = fullfile(path_save, 'Fig_Forest_RPE_Correlations.png');
saveas(fig, fname);
fprintf('  Fig_Forest_RPE_Correlations.png sauvegardée\n');
close(fig);

%% -----------------------------------------------------------------------
%  TABLEAU RÉCAPITULATIF TEXTE
%  -----------------------------------------------------------------------
fprintf('\n=== TABLEAU RÉCAPITULATIF RPE — Pearson et Spearman ===\n');
fprintf('%-30s  %-14s  %-14s  %-14s  %-14s\n', ...
    'Variable', 'Pearson Norm', 'Pearson CS60', 'Spearman Norm', 'Spearman CS60');
fprintf('%s\n', repmat('-',1,95));

for iV = 1:length(vars_to_test)
    vn = vars_to_test{iV};
    if ~isfield(RPE_Corr_Results, vn), continue; end

    r_n   = RPE_Corr_Results.(vn).r_norm;   p_n   = RPE_Corr_Results.(vn).p_norm;
    r_e   = RPE_Corr_Results.(vn).r_cs60;   p_e   = RPE_Corr_Results.(vn).p_cs60;
    rho_n = RPE_Corr_Results.(vn).rho_norm; prn   = RPE_Corr_Results.(vn).prho_norm;
    rho_e = RPE_Corr_Results.(vn).rho_cs60; pre   = RPE_Corr_Results.(vn).prho_cs60;

    fprintf('%-30s  %+.3f (%s)    %+.3f (%s)    %+.3f (%s)    %+.3f (%s)\n', ...
        vn, r_n, p2star(p_n), r_e, p2star(p_e), rho_n, p2star(prn), rho_e, p2star(pre));
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'RPE_Analysis_Results.mat'), ...
    'RPE', 'Derive_Results', 'RPE_Corr_Results');
fprintf('\nRésultats RPE sauvegardés : RPE_Analysis_Results.mat\n');
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