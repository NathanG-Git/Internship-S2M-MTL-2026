%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  ANALYSE RPE PAR SOUS-GROUPE — SD (petites mains) vs LD (grandes) %%%%
%%%%  Piano Normal uniquement                                          %%%%
%%%%                                                                     %%%%
%%%%  Hypothèse : petites mains = plus de contrainte mécanique         %%%%
%%%%  relative sur clavier standard = fatigue plus rapide/marquée      %%%%
%%%%  (analogue à G1/G2 dans Goubault et al.)                          %%%%
%%%%                                                                     %%%%
%%%%  Étape 1 : Mann-Whitney sur RPE_début, RPE_fin, ΔRPE entre groupes%%%%
%%%%  Étape 2 : Spearman Dérive ↔ ΔRPE, séparément par groupe         %%%%
%%%%                                                                     %%%%
%%%%  P07 exclu (déjà absent des deux listes ci-dessous)               %%%%
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
load(fullfile(path_data, 'LMM_Validation_Results.mat'));  % LMM_Results, Index_Results
load(fullfile(path_data, 'IMU_Features_NewDataset.mat')); % Results_IMU
load(fullfile(path_data, 'EMG_Features_NewDataset.mat')); % Results (EMG)
fprintf('  OK\n\n');

if isfield(Results_IMU, 'P07'), Results_IMU = rmfield(Results_IMU, 'P07'); end
if isfield(Results, 'P07'),     Results     = rmfield(Results, 'P07');     end

%% -----------------------------------------------------------------------
%  GROUPES SD (petites mains) ET LD (grandes mains) — condition Normal
%  -----------------------------------------------------------------------
SD_data = {
    'P02', 30, 75;
    'P03', 25, 85;
    'P04', 15, 23;
    'P05', 35, 75;
    'P08',  7, 18;
    'P10',  7, 30;
    'P12', 13, 95;
    'P15',  5, 70;
    'P17', 13, 60;
    'P18', 13, 60;
    'P20', 15, 55;
    'P21', 10, 50;
    'P24', 10, 25;
};

LD_data = {
    'P01',  0, 30;
    'P06', 10, 55;
    'P09',  5, 60;
    'P11',  7, 60;
    'P13', 17, 55;
    'P14', 22, 28;
    'P16', 20, 60;
    'P19',  0, 25;
    'P22',  1, 20;
    'P23',  3, 35;
};

Subject_SD = SD_data(:,1)';
Subject_LD = LD_data(:,1)';

RPE = struct();
for i = 1:size(SD_data,1)
    subj = SD_data{i,1};
    RPE.(subj).debut = SD_data{i,2};
    RPE.(subj).fin   = SD_data{i,3};
    RPE.(subj).delta = SD_data{i,3} - SD_data{i,2};
    RPE.(subj).group = 'SD';
end
for i = 1:size(LD_data,1)
    subj = LD_data{i,1};
    RPE.(subj).debut = LD_data{i,2};
    RPE.(subj).fin   = LD_data{i,3};
    RPE.(subj).delta = LD_data{i,3} - LD_data{i,2};
    RPE.(subj).group = 'LD';
end

fprintf('SD (petites mains) : N=%d | LD (grandes mains) : N=%d\n\n', ...
    length(Subject_SD), length(Subject_LD));

%% -----------------------------------------------------------------------
%  EXCLUSION DES SUJETS IDENTIFIÉS COMME PERTURBATEURS DE LA CORRÉLATION
%  (RPE atypique par rapport à leur dérive biomécanique, identifiés
%   visuellement sur les figures scatter triées par dérive)
%  -----------------------------------------------------------------------
exclude_outliers = false;   % mettre à false pour revenir à l'analyse complète
Outliers_SD = {'P04','P08','P12','P15'};
Outliers_LD = {'P14'};

if exclude_outliers
    fprintf('=== EXCLUSION DES OUTLIERS ===\n');
    fprintf('  SD : %s\n', strjoin(Outliers_SD, ', '));
    fprintf('  LD : %s\n\n', strjoin(Outliers_LD, ', '));
    Subject_SD = Subject_SD(~ismember(Subject_SD, Outliers_SD));
    Subject_LD = Subject_LD(~ismember(Subject_LD, Outliers_LD));
    fprintf('  Nouveaux effectifs : SD N=%d | LD N=%d\n\n', length(Subject_SD), length(Subject_LD));
end

%% -----------------------------------------------------------------------
%  ÉTAPE 1 — MANN-WHITNEY SUR RPE_DÉBUT, RPE_FIN, ΔRPE
%  -----------------------------------------------------------------------
fprintf('=== ÉTAPE 1 : Mann-Whitney SD vs LD ===\n\n');

rpe_debut_SD = cellfun(@(s) RPE.(s).debut, Subject_SD)';
rpe_fin_SD   = cellfun(@(s) RPE.(s).fin,   Subject_SD)';
rpe_delta_SD = cellfun(@(s) RPE.(s).delta, Subject_SD)';

rpe_debut_LD = cellfun(@(s) RPE.(s).debut, Subject_LD)';
rpe_fin_LD   = cellfun(@(s) RPE.(s).fin,   Subject_LD)';
rpe_delta_LD = cellfun(@(s) RPE.(s).delta, Subject_LD)';

[p_debut, ~, stats_debut] = ranksum(rpe_debut_SD, rpe_debut_LD);
[p_fin,   ~, stats_fin]   = ranksum(rpe_fin_SD,   rpe_fin_LD);
[p_delta, ~, stats_delta] = ranksum(rpe_delta_SD, rpe_delta_LD);

fprintf('%-15s  %-8s  %-8s  %-10s  %-10s\n', 'Mesure','Med SD','Med LD','U-stat','p (MW)');
fprintf('%s\n', repmat('-',1,60));
fprintf('%-15s  %-8.1f  %-8.1f  %-10.1f  %.4f%s\n', 'RPE_début', ...
    median(rpe_debut_SD), median(rpe_debut_LD), stats_debut.ranksum, p_debut, p2star(p_debut));
fprintf('%-15s  %-8.1f  %-8.1f  %-10.1f  %.4f%s\n', 'RPE_fin', ...
    median(rpe_fin_SD), median(rpe_fin_LD), stats_fin.ranksum, p_fin, p2star(p_fin));
fprintf('%-15s  %-8.1f  %-8.1f  %-10.1f  %.4f%s\n', 'ΔRPE', ...
    median(rpe_delta_SD), median(rpe_delta_LD), stats_delta.ranksum, p_delta, p2star(p_delta));
fprintf('\n');

%% -----------------------------------------------------------------------
%  BOXPLOT RPE PAR GROUPE
%  -----------------------------------------------------------------------
fig0 = figure('Name','RPE_SD_vs_LD','NumberTitle','off', ...
    'Position',[50 50 900 350], 'Color','white');

measures = {rpe_debut_SD, rpe_debut_LD, 'RPE début', p_debut; ...
            rpe_fin_SD,   rpe_fin_LD,   'RPE fin',   p_fin; ...
            rpe_delta_SD, rpe_delta_LD, 'ΔRPE',      p_delta};

for im = 1:3
    subplot(1,3,im); hold on;
    g_sd = measures{im,1}; g_ld = measures{im,2};
    grp_lbl = [repmat({'SD'},length(g_sd),1); repmat({'LD'},length(g_ld),1)];
    boxplot([g_sd; g_ld], grp_lbl, 'Colors', [0.20 0.40 0.70; 0.80 0.40 0.10]);
    title(sprintf('%s (p=%.3f%s)', measures{im,3}, measures{im,4}, p2star(measures{im,4})), ...
        'FontSize',10,'FontWeight','bold');
    ylabel(measures{im,3});
    grid on; box on;
end
title_suffix = '';
if exclude_outliers
    title_suffix = ' — outliers exclus';
end
sgtitle(sprintf('RPE — SD (petites mains) vs LD (grandes mains), Piano Normal%s', title_suffix), ...
    'FontWeight','bold','FontSize',12);

fname0 = fullfile(path_save, 'Fig_RPE_SD_vs_LD.png');
saveas(fig0, fname0);
fprintf('Fig_RPE_SD_vs_LD.png sauvegardée\n\n');
close(fig0);

%% -----------------------------------------------------------------------
%  DÉFINITION DES BLOCS ET DÉRIVE (intervalle 10 - intervalle 1, z-scorée)
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

all_vars = {};
for ib = 1:length(blocs)
    all_vars = [all_vars, blocs(ib).vars];
end

all_subj = [Subject_SD, Subject_LD];
N_bins_idx = 10;

for ib = 1:length(blocs)
    bloc = blocs(ib);
    for iV = 1:length(bloc.vars)
        vn = bloc.vars{iV};
        if ~isfield(LMM_Results, vn), continue; end
        mu_z = LMM_Results.(vn).mu_z;
        sd_z = LMM_Results.(vn).sd_z;

        for iS = 1:length(all_subj)
            subj = all_subj{iS};
            try
                if strcmp(bloc.source,'IMU')
                    y = Results_IMU.(subj).Norm.(vn)(:);
                else
                    y = Results.(subj).Norm.(vn)(:);
                end
                z = (y - mu_z) / sd_z;
                Derive_Results.(subj).(vn) = z(end) - z(1);
            catch
                Derive_Results.(subj).(vn) = NaN;
            end
        end
    end
end

% Dérive des index — recalculée directement depuis les variables brutes
% (Index_Results ne contient que les 13 sujets crossover de la validation
%  principale ; ici on a besoin de tous les sujets SD+LD en Normal)
%
% Index = moyenne des z-scores signés (signes "observed", non pondéré)
% pour rester cohérent avec la variante utilisée par défaut ailleurs

% Signes observés — recalculés à partir du signe de β_Temps dans LMM_Results
sign_observed = struct();
for ib = 1:length(blocs)
    for iV = 1:length(blocs(ib).vars)
        vn = blocs(ib).vars{iV};
        if isfield(LMM_Results, vn)
            FE = LMM_Results.(vn).FE;
            idx_t = find(strcmp(FE.Name, 'Temps'));
            if ~isempty(idx_t)
                b = FE.Estimate(idx_t);
                s = sign(b); if s == 0, s = 1; end
                sign_observed.(vn) = s;
            else
                sign_observed.(vn) = 1;
            end
        end
    end
end

for ib = 1:length(blocs)
    bloc  = blocs(ib);
    bname = bloc.name;
    vars  = bloc.vars;
    N_v   = length(vars);

    for iS = 1:length(all_subj)
        subj = all_subj{iS};

        Z_mat = NaN(N_bins_idx, N_v);
        for iV = 1:N_v
            vn = vars{iV};
            if ~isfield(LMM_Results, vn), continue; end
            try
                if strcmp(bloc.source,'IMU')
                    y = Results_IMU.(subj).Norm.(vn)(:);
                else
                    y = Results.(subj).Norm.(vn)(:);
                end
                mu_z = LMM_Results.(vn).mu_z;
                sd_z = LMM_Results.(vn).sd_z;
                z    = (y - mu_z) / sd_z;
                Z_mat(:,iV) = z * sign_observed.(vn);
            catch; end
        end

        idx_val = mean(Z_mat, 2, 'omitnan');
        if all(isnan(idx_val))
            Derive_Results.(subj).(['Index_' bname]) = NaN;
        else
            Derive_Results.(subj).(['Index_' bname]) = idx_val(end) - idx_val(1);
        end
    end
end

vars_to_test = [all_vars, {'Index_SegMod','Index_Goubault','Index_EMG'}];

%% -----------------------------------------------------------------------
%  ÉTAPE 2 — SPEARMAN PAR GROUPE (SD et LD séparément)
%  -----------------------------------------------------------------------
fprintf('=== ÉTAPE 2 : Spearman Dérive ↔ ΔRPE, par groupe (SD / LD) ===\n\n');

RPE_Corr_Groups = struct();

for iV = 1:length(vars_to_test)
    vn = vars_to_test{iV};

    X_sd=[]; Y_sd=[];
    for iS = 1:length(Subject_SD)
        subj = Subject_SD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn)
            drv = Derive_Results.(subj).(vn);
            if ~isnan(drv)
                X_sd(end+1) = drv;
                Y_sd(end+1) = RPE.(subj).delta;
            end
        end
    end

    X_ld=[]; Y_ld=[];
    for iS = 1:length(Subject_LD)
        subj = Subject_LD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn)
            drv = Derive_Results.(subj).(vn);
            if ~isnan(drv)
                X_ld(end+1) = drv;
                Y_ld(end+1) = RPE.(subj).delta;
            end
        end
    end

    rho_sd=NaN; p_sd=NaN; rho_ld=NaN; p_ld=NaN;
    if length(X_sd) >= 4, [rho_sd, p_sd] = corr(X_sd(:), Y_sd(:), 'Type','Spearman'); end
    if length(X_ld) >= 4, [rho_ld, p_ld] = corr(X_ld(:), Y_ld(:), 'Type','Spearman'); end

    RPE_Corr_Groups.(vn).rho_sd = rho_sd; RPE_Corr_Groups.(vn).p_sd = p_sd; RPE_Corr_Groups.(vn).N_sd = length(X_sd);
    RPE_Corr_Groups.(vn).rho_ld = rho_ld; RPE_Corr_Groups.(vn).p_ld = p_ld; RPE_Corr_Groups.(vn).N_ld = length(X_ld);

    fprintf('%-35s  SD (N=%d) ρ=%+.3f (%s)  |  LD (N=%d) ρ=%+.3f (%s)\n', ...
        vn, length(X_sd), rho_sd, p2star_full(p_sd), length(X_ld), rho_ld, p2star_full(p_ld));

    % --- Détail par sujet (pour identifier les points influents) ---
    fprintf('    SD : ');
    iSubj_sd_used = 0;
    for iS = 1:length(Subject_SD)
        subj = Subject_SD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn) ...
                && ~isnan(Derive_Results.(subj).(vn))
            iSubj_sd_used = iSubj_sd_used + 1;
            fprintf('%s(Δz=%+.2f,ΔRPE=%d) ', subj, Derive_Results.(subj).(vn), RPE.(subj).delta);
        end
    end
    fprintf('\n    LD : ');
    for iS = 1:length(Subject_LD)
        subj = Subject_LD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn) ...
                && ~isnan(Derive_Results.(subj).(vn))
            fprintf('%s(Δz=%+.2f,ΔRPE=%d) ', subj, Derive_Results.(subj).(vn), RPE.(subj).delta);
        end
    end
    fprintf('\n\n');
end

%% -----------------------------------------------------------------------
%  FOREST PLOT — SPEARMAN PAR GROUPE
%  -----------------------------------------------------------------------
fprintf('\n=== Génération du forest plot SD vs LD ===\n');

col_sd = [0.20 0.40 0.70];
col_ld = [0.80 0.40 0.10];

plot_order = flip({'Accel_Mod_Head','Accel_Mod_Hand','Index_SegMod', ...
              'MedianFreq_Accel_Y_Hand','PeakPower_Accel_Mod_Hand', ...
              'PeakPower_AngVel_X_Head','SpectralEntropy_AngVel_Mod_Forearm','Index_Goubault', ...
              'TFR_MedianFreq_Triceps','TFR_SpectralEntropy_Deltoid','SampleEntropy_Biceps','Index_EMG'});

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

N_v = length(plot_order);
y_pos_sd = (1:N_v) + 0.15;
y_pos_ld = (1:N_v) - 0.15;

fig = figure('Name','Forest_RPE_SD_LD','NumberTitle','off', ...
    'Position',[50 30 900 700], 'Color','white');
hold on;

for iV = 1:N_v
    vn = plot_order{iV};
    if ~isfield(RPE_Corr_Groups, vn), continue; end

    rho_sd = RPE_Corr_Groups.(vn).rho_sd; N_sd = RPE_Corr_Groups.(vn).N_sd;
    rho_ld = RPE_Corr_Groups.(vn).rho_ld; N_ld = RPE_Corr_Groups.(vn).N_ld;

    if ~isnan(rho_sd) && N_sd > 3
        z_sd = atanh(max(min(rho_sd,0.999),-0.999)); se_sd = 1/sqrt(N_sd-3);
        ci_sd = tanh([z_sd-1.96*se_sd, z_sd+1.96*se_sd]);
    else
        ci_sd = [NaN NaN];
    end
    if ~isnan(rho_ld) && N_ld > 3
        z_ld = atanh(max(min(rho_ld,0.999),-0.999)); se_ld = 1/sqrt(N_ld-3);
        ci_ld = tanh([z_ld-1.96*se_ld, z_ld+1.96*se_ld]);
    else
        ci_ld = [NaN NaN];
    end

    if ~any(isnan(ci_sd))
        plot(ci_sd, [y_pos_sd(iV) y_pos_sd(iV)], '-', 'Color', col_sd, 'LineWidth', 2);
    end
    if ~any(isnan(ci_ld))
        plot(ci_ld, [y_pos_ld(iV) y_pos_ld(iV)], '-', 'Color', col_ld, 'LineWidth', 2);
    end

    if ~isnan(rho_sd)
        plot(rho_sd, y_pos_sd(iV), 'o', 'MarkerSize',8, 'MarkerFaceColor',col_sd, ...
            'MarkerEdgeColor','white', 'LineWidth',1);
    end
    if ~isnan(rho_ld)
        plot(rho_ld, y_pos_ld(iV), 's', 'MarkerSize',8, 'MarkerFaceColor',col_ld, ...
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
xlabel('Corrélation Spearman ρ (dérive ↔ ΔRPE)', 'FontSize',10);
title(sprintf('Spearman dérive ↔ ΔRPE — SD vs LD, Piano Normal%s (IC95%% Fisher z)', title_suffix), ...
    'FontSize',12, 'FontWeight','bold');

yline(3.5, ':', 'Color',[0.6 0.6 0.6]);
yline(8.5, ':', 'Color',[0.6 0.6 0.6]);

plot(NaN,NaN,'o','MarkerSize',8,'MarkerFaceColor',col_sd,'MarkerEdgeColor','white','LineWidth',1,'DisplayName','SD (petites mains)');
plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor',col_ld,'MarkerEdgeColor','white','LineWidth',1,'DisplayName','LD (grandes mains)');
legend('Location','southoutside','Orientation','horizontal','FontSize',9);

grid on; box on;
set(gca,'FontSize',9);

fname = fullfile(path_save, 'Fig_Forest_RPE_SD_LD.png');
saveas(fig, fname);
fprintf('Fig_Forest_RPE_SD_LD.png sauvegardée\n');
close(fig);

%% -----------------------------------------------------------------------
%  FIGURES SCATTER TRIÉES PAR DÉRIVE — POINTS RELIÉS PAR LIGNE FINE
%  -----------------------------------------------------------------------
fprintf('\n=== Génération des figures scatter triées par dérive ===\n');

col_sd = [0.20 0.40 0.70];
col_ld = [0.80 0.40 0.10];

for ib = 1:length(blocs)
    bloc     = blocs(ib);
    vars_fig = [bloc.vars, {['Index_' bloc.name]}];
    N_v      = length(vars_fig);

    fig = figure('Name', sprintf('RPE_SortedScatter_%s', bloc.name), ...
        'NumberTitle','off', 'Position',[30 30 380*N_v 700], 'Color','white');

    for iV = 1:N_v
        vn = vars_fig{iV};

        % --- SD ---
        x_sd=[]; y_sd=[]; lbl_sd={};
        for iS = 1:length(Subject_SD)
            subj = Subject_SD{iS};
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn) ...
                    && ~isnan(Derive_Results.(subj).(vn))
                x_sd(end+1) = Derive_Results.(subj).(vn);
                y_sd(end+1) = RPE.(subj).delta;
                lbl_sd{end+1} = subj;
            end
        end
        [x_sd_sorted, ord_sd] = sort(x_sd);
        y_sd_sorted = y_sd(ord_sd);
        lbl_sd_sorted = lbl_sd(ord_sd);

        % --- LD ---
        x_ld=[]; y_ld=[]; lbl_ld={};
        for iS = 1:length(Subject_LD)
            subj = Subject_LD{iS};
            if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn) ...
                    && ~isnan(Derive_Results.(subj).(vn))
                x_ld(end+1) = Derive_Results.(subj).(vn);
                y_ld(end+1) = RPE.(subj).delta;
                lbl_ld{end+1} = subj;
            end
        end
        [x_ld_sorted, ord_ld] = sort(x_ld);
        y_ld_sorted = y_ld(ord_ld);
        lbl_ld_sorted = lbl_ld(ord_ld);

        % --- Subplot SD (haut) ---
        subplot(2, N_v, iV); hold on;
        plot(1:length(x_sd_sorted), y_sd_sorted, '-', 'Color',[col_sd 0.5], 'LineWidth',1);
        scatter(1:length(x_sd_sorted), y_sd_sorted, 60, col_sd, 'o', 'filled');
        for k = 1:length(x_sd_sorted)
            text(k, y_sd_sorted(k), ['  ' lbl_sd_sorted{k}], 'FontSize',6, 'Rotation',45);
        end
        if isfield(RPE_Corr_Groups, vn)
            rho = RPE_Corr_Groups.(vn).rho_sd; p = RPE_Corr_Groups.(vn).p_sd;
            title(sprintf('%s\nSD ρ=%+.2f (%s)', strrep(vn,'_',' '), rho, p2star_full(p)), ...
                'FontSize',8,'FontWeight','bold');
        end
        xlabel('Sujets triés par Dérive ↑', 'FontSize',7);
        ylabel('ΔRPE', 'FontSize',7);
        xticks(1:length(x_sd_sorted));
        xticklabels(arrayfun(@(v) sprintf('%+.1f',v), x_sd_sorted, 'UniformOutput',false));
        xtickangle(45);
        set(gca,'FontSize',6);
        grid on; box on;

        % --- Subplot LD (bas) ---
        subplot(2, N_v, iV + N_v); hold on;
        plot(1:length(x_ld_sorted), y_ld_sorted, '-', 'Color',[col_ld 0.5], 'LineWidth',1);
        scatter(1:length(x_ld_sorted), y_ld_sorted, 60, col_ld, 's', 'filled');
        for k = 1:length(x_ld_sorted)
            text(k, y_ld_sorted(k), ['  ' lbl_ld_sorted{k}], 'FontSize',6, 'Rotation',45);
        end
        if isfield(RPE_Corr_Groups, vn)
            rho = RPE_Corr_Groups.(vn).rho_ld; p = RPE_Corr_Groups.(vn).p_ld;
            title(sprintf('LD ρ=%+.2f (%s)', rho, p2star_full(p)), 'FontSize',8,'FontWeight','bold');
        end
        xlabel('Sujets triés par Dérive ↑', 'FontSize',7);
        ylabel('ΔRPE', 'FontSize',7);
        xticks(1:length(x_ld_sorted));
        xticklabels(arrayfun(@(v) sprintf('%+.1f',v), x_ld_sorted, 'UniformOutput',false));
        xtickangle(45);
        set(gca,'FontSize',6);
        grid on; box on;
    end

    sgtitle(sprintf('Bloc %s — ΔRPE en fonction des sujets triés par dérive (SD haut, LD bas)%s', bloc.name, title_suffix), ...
        'FontWeight','bold','FontSize',11);

    fname = fullfile(path_save, sprintf('Fig_RPE_SortedScatter_%s.png', bloc.name));
    saveas(fig, fname);
    fprintf('  Fig_RPE_SortedScatter_%s.png sauvegardée\n', bloc.name);
    close(fig);
end

%% -----------------------------------------------------------------------
%  ÉTAPE 3 — SPEARMAN RPE_DÉBUT ↔ DÉRIVE, PAR GROUPE (SD / LD)
%  -----------------------------------------------------------------------
fprintf('\n=== ÉTAPE 3 : Spearman RPE_début ↔ Dérive, par groupe (SD / LD) ===\n\n');

RPE0_Corr_Groups = struct();

for iV = 1:length(vars_to_test)
    vn = vars_to_test{iV};

    X_sd=[]; Y_sd=[];
    for iS = 1:length(Subject_SD)
        subj = Subject_SD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn)
            drv = Derive_Results.(subj).(vn);
            if ~isnan(drv)
                X_sd(end+1) = RPE.(subj).debut;
                Y_sd(end+1) = drv;
            end
        end
    end

    X_ld=[]; Y_ld=[];
    for iS = 1:length(Subject_LD)
        subj = Subject_LD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn)
            drv = Derive_Results.(subj).(vn);
            if ~isnan(drv)
                X_ld(end+1) = RPE.(subj).debut;
                Y_ld(end+1) = drv;
            end
        end
    end

    rho_sd=NaN; p_sd=NaN; rho_ld=NaN; p_ld=NaN;
    if length(X_sd) >= 4, [rho_sd, p_sd] = corr(X_sd(:), Y_sd(:), 'Type','Spearman'); end
    if length(X_ld) >= 4, [rho_ld, p_ld] = corr(X_ld(:), Y_ld(:), 'Type','Spearman'); end

    RPE0_Corr_Groups.(vn).rho_sd = rho_sd; RPE0_Corr_Groups.(vn).p_sd = p_sd; RPE0_Corr_Groups.(vn).N_sd = length(X_sd);
    RPE0_Corr_Groups.(vn).rho_ld = rho_ld; RPE0_Corr_Groups.(vn).p_ld = p_ld; RPE0_Corr_Groups.(vn).N_ld = length(X_ld);

    fprintf('%-35s  SD (N=%d) ρ=%+.3f (%s)  |  LD (N=%d) ρ=%+.3f (%s)\n', ...
        vn, length(X_sd), rho_sd, p2star_full(p_sd), length(X_ld), rho_ld, p2star_full(p_ld));
end

%% -----------------------------------------------------------------------
%  FOREST PLOT — SPEARMAN RPE_DÉBUT ↔ DÉRIVE, PAR GROUPE
%  -----------------------------------------------------------------------
fprintf('\n=== Génération du forest plot RPE_début ↔ Dérive ===\n');

fig2 = figure('Name','Forest_RPE0_Derive_SD_LD','NumberTitle','off', ...
    'Position',[50 30 900 700], 'Color','white');
hold on;

for iV = 1:N_v
    vn = plot_order{iV};
    if ~isfield(RPE0_Corr_Groups, vn), continue; end

    rho_sd = RPE0_Corr_Groups.(vn).rho_sd; N_sd2 = RPE0_Corr_Groups.(vn).N_sd;
    rho_ld = RPE0_Corr_Groups.(vn).rho_ld; N_ld2 = RPE0_Corr_Groups.(vn).N_ld;

    if ~isnan(rho_sd) && N_sd2 > 3
        z_sd = atanh(max(min(rho_sd,0.999),-0.999)); se_sd = 1/sqrt(N_sd2-3);
        ci_sd = tanh([z_sd-1.96*se_sd, z_sd+1.96*se_sd]);
    else
        ci_sd = [NaN NaN];
    end
    if ~isnan(rho_ld) && N_ld2 > 3
        z_ld = atanh(max(min(rho_ld,0.999),-0.999)); se_ld = 1/sqrt(N_ld2-3);
        ci_ld = tanh([z_ld-1.96*se_ld, z_ld+1.96*se_ld]);
    else
        ci_ld = [NaN NaN];
    end

    if ~any(isnan(ci_sd))
        plot(ci_sd, [y_pos_sd(iV) y_pos_sd(iV)], '-', 'Color', col_sd, 'LineWidth', 2);
    end
    if ~any(isnan(ci_ld))
        plot(ci_ld, [y_pos_ld(iV) y_pos_ld(iV)], '-', 'Color', col_ld, 'LineWidth', 2);
    end

    if ~isnan(rho_sd)
        plot(rho_sd, y_pos_sd(iV), 'o', 'MarkerSize',8, 'MarkerFaceColor',col_sd, ...
            'MarkerEdgeColor','white', 'LineWidth',1);
    end
    if ~isnan(rho_ld)
        plot(rho_ld, y_pos_ld(iV), 's', 'MarkerSize',8, 'MarkerFaceColor',col_ld, ...
            'MarkerEdgeColor','white', 'LineWidth',1);
    end
end

xline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);

yticks(1:N_v);
yticklabels(ylabels);
ylim([0.5 N_v+0.5]);
xlim([-1 1]);
xlabel('Corrélation Spearman ρ (RPE_début ↔ Dérive)', 'FontSize',10);
title(sprintf('Spearman RPE_début ↔ Dérive — SD vs LD, Piano Normal%s (IC95%% Fisher z)', title_suffix), ...
    'FontSize',12, 'FontWeight','bold');

yline(3.5, ':', 'Color',[0.6 0.6 0.6]);
yline(8.5, ':', 'Color',[0.6 0.6 0.6]);

plot(NaN,NaN,'o','MarkerSize',8,'MarkerFaceColor',col_sd,'MarkerEdgeColor','white','LineWidth',1,'DisplayName','SD (petites mains)');
plot(NaN,NaN,'s','MarkerSize',8,'MarkerFaceColor',col_ld,'MarkerEdgeColor','white','LineWidth',1,'DisplayName','LD (grandes mains)');
legend('Location','southoutside','Orientation','horizontal','FontSize',9);

grid on; box on;
set(gca,'FontSize',9);

fname2 = fullfile(path_save, 'Fig_Forest_RPE0_Derive_SD_LD.png');
saveas(fig2, fname2);
fprintf('Fig_Forest_RPE0_Derive_SD_LD.png sauvegardée\n');
close(fig2);

%% -----------------------------------------------------------------------
%  VÉRIFICATION — RPE_DÉBUT ↔ ΔRPE, PAR GROUPE
%  -----------------------------------------------------------------------
fprintf('\n=== RPE_début ↔ ΔRPE, par groupe (SD / LD) ===\n\n');

rpe0_sd  = cellfun(@(s) RPE.(s).debut, Subject_SD)';
drpe_sd  = cellfun(@(s) RPE.(s).delta, Subject_SD)';
rpe0_ld  = cellfun(@(s) RPE.(s).debut, Subject_LD)';
drpe_ld  = cellfun(@(s) RPE.(s).delta, Subject_LD)';

[rho_rpe0_sd, p_rpe0_sd] = corr(rpe0_sd, drpe_sd, 'Type','Spearman');
[rho_rpe0_ld, p_rpe0_ld] = corr(rpe0_ld, drpe_ld, 'Type','Spearman');

fprintf('SD (N=%d) : RPE_début ↔ ΔRPE  ρ=%+.3f (%s)\n', length(rpe0_sd), rho_rpe0_sd, p2star_full(p_rpe0_sd));
fprintf('LD (N=%d) : RPE_début ↔ ΔRPE  ρ=%+.3f (%s)\n\n', length(rpe0_ld), rho_rpe0_ld, p2star_full(p_rpe0_ld));

%% -----------------------------------------------------------------------
%  ÉTAPE 4 — CORRÉLATION PARTIELLE : Dérive ↔ ΔRPE, CONTRÔLÉE PAR RPE_DÉBUT
%  Méthode : résidus de Dérive~RPE_début et résidus de ΔRPE~RPE_début
%  (basé sur les rangs, cohérent avec l'approche Spearman utilisée ailleurs)
%  -----------------------------------------------------------------------
fprintf('\n=== ÉTAPE 4 : Corrélation partielle Dérive ↔ ΔRPE | RPE_début ===\n\n');

Partial_Corr_Groups = struct();

for iV = 1:length(vars_to_test)
    vn = vars_to_test{iV};

    % --- SD ---
    X_sd=[]; Y_sd=[]; R0_sd=[];
    for iS = 1:length(Subject_SD)
        subj = Subject_SD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn)
            drv = Derive_Results.(subj).(vn);
            if ~isnan(drv)
                X_sd(end+1)  = drv;
                Y_sd(end+1)  = RPE.(subj).delta;
                R0_sd(end+1) = RPE.(subj).debut;
            end
        end
    end

    % --- LD ---
    X_ld=[]; Y_ld=[]; R0_ld=[];
    for iS = 1:length(Subject_LD)
        subj = Subject_LD{iS};
        if isfield(Derive_Results,subj) && isfield(Derive_Results.(subj),vn)
            drv = Derive_Results.(subj).(vn);
            if ~isnan(drv)
                X_ld(end+1)  = drv;
                Y_ld(end+1)  = RPE.(subj).delta;
                R0_ld(end+1) = RPE.(subj).debut;
            end
        end
    end

    rho_p_sd=NaN; p_p_sd=NaN; rho_p_ld=NaN; p_p_ld=NaN;
    if length(X_sd) >= 5
        [rho_p_sd, p_p_sd] = partialcorr(X_sd(:), Y_sd(:), R0_sd(:), 'Type','Spearman');
    end
    if length(X_ld) >= 5
        [rho_p_ld, p_p_ld] = partialcorr(X_ld(:), Y_ld(:), R0_ld(:), 'Type','Spearman');
    end

    Partial_Corr_Groups.(vn).rho_sd = rho_p_sd; Partial_Corr_Groups.(vn).p_sd = p_p_sd; Partial_Corr_Groups.(vn).N_sd = length(X_sd);
    Partial_Corr_Groups.(vn).rho_ld = rho_p_ld; Partial_Corr_Groups.(vn).p_ld = p_p_ld; Partial_Corr_Groups.(vn).N_ld = length(X_ld);

    fprintf('%-35s  SD (N=%d) ρ_partial=%+.3f (%s)  |  LD (N=%d) ρ_partial=%+.3f (%s)\n', ...
        vn, length(X_sd), rho_p_sd, p2star_full(p_p_sd), length(X_ld), rho_p_ld, p2star_full(p_p_ld));
end

%% -----------------------------------------------------------------------
%  COMPARAISON BRUTE vs PARTIELLE — TABLEAU RÉCAPITULATIF
%  -----------------------------------------------------------------------
fprintf('\n=== Comparaison ρ brut (Dérive↔ΔRPE) vs ρ partiel (contrôlé RPE_début) ===\n\n');
fprintf('%-35s  %-22s  %-22s\n', 'Variable', 'SD : brut → partiel', 'LD : brut → partiel');
fprintf('%s\n', repmat('-',1,90));
for iV = 1:length(vars_to_test)
    vn = vars_to_test{iV};
    if ~isfield(RPE_Corr_Groups, vn) || ~isfield(Partial_Corr_Groups, vn), continue; end
    rb_sd = RPE_Corr_Groups.(vn).rho_sd; pb_sd = RPE_Corr_Groups.(vn).p_sd;
    rp_sd = Partial_Corr_Groups.(vn).rho_sd; pp_sd = Partial_Corr_Groups.(vn).p_sd;
    rb_ld = RPE_Corr_Groups.(vn).rho_ld; pb_ld = RPE_Corr_Groups.(vn).p_ld;
    rp_ld = Partial_Corr_Groups.(vn).rho_ld; pp_ld = Partial_Corr_Groups.(vn).p_ld;
    fprintf('%-35s  %+.2f(%s) → %+.2f(%s)     %+.2f(%s) → %+.2f(%s)\n', ...
        vn, rb_sd, p2star(pb_sd), rp_sd, p2star(pp_sd), rb_ld, p2star(pb_ld), rp_ld, p2star(pp_ld));
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, 'RPE_SD_LD_Results.mat'), ...
    'RPE', 'Derive_Results', 'RPE_Corr_Groups', 'RPE0_Corr_Groups', 'Partial_Corr_Groups', 'Subject_SD', 'Subject_LD');
fprintf('\nRésultats sauvegardés : RPE_SD_LD_Results.mat\n');
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