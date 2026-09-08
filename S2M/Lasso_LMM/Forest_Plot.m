%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  FOREST PLOTS DE STABILITÉ — SCRIPT AUTONOME                      %%%%
%%%%  Charge Stability_Archive.mat et génère les forest plots           %%%%
%%%%  Sans avoir besoin de relancer les 100 seeds                      %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
% Chemin vers Stability_Archive.mat — adapter si nécessaire
path_archive = fileparts(mfilename('fullpath'));
if isempty(path_archive), path_archive = pwd; end

PathSave = path_archive;  % dossier de sauvegarde des PNG

THRESH_STAB = 50;         % seuil de fréquence de sélection (%)
bloc_names  = {'SegMod','SegAxis','Goubault','EMG'};
col_G1c     = [0.72 0.11 0.11];
col_G2c     = [0.15 0.32 0.60];

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
archive_file = fullfile(path_archive, 'Stability_Archive.mat');
if ~exist(archive_file, 'file')
    error('Fichier Stability_Archive.mat introuvable dans : %s', path_archive);
end

fprintf('Chargement de Stability_Archive.mat...\n');
tmp = load(archive_file);
stab_G1 = tmp.stability_archive.stab_G1;
stab_G2 = tmp.stability_archive.stab_G2;
N_stab  = tmp.stability_archive.N_stab;
fprintf('  N_stab = %d runs | Seuil = %d%%\n', N_stab, THRESH_STAB);

%% -----------------------------------------------------------------------
%  GÉNÉRATION DES FOREST PLOTS
%  -----------------------------------------------------------------------
fprintf('\n=== Génération des Forest Plots ===\n');

for ib = 1:4
    bn = bloc_names{ib};
    for iGrp = 1:2
        if iGrp==1, st=stab_G1.(bn); gname='G1'; col_g=col_G1c;
        else,        st=stab_G2.(bn); gname='G2'; col_g=col_G2c; end

        if isempty(fieldnames(st)), continue; end
        fns = fieldnames(st);

        % Collecter variables au-dessus du seuil
        vars_fp = {}; mu_fp = []; ci_lo_fp = []; ci_hi_fp = [];
        freq_fp = []; min_fp = []; max_fp = [];

        for iv = 1:length(fns)
            fn = fns{iv};
            s  = st.(fn);
            freq = s.n_sel / N_stab * 100;
            if freq < THRESH_STAB || isempty(s.beta_std), continue; end

            mu_b  = mean(s.beta_std, 'omitnan');
            sd_b  = std(s.beta_std, 0, 'omitnan');
            n_b   = sum(~isnan(s.beta_std));

            % Nettoyer le nom de variable
            vn_clean = strrep(strrep(fn, '___', ' — '), '_', ' ');

            vars_fp{end+1}  = vn_clean;
            mu_fp(end+1)    = mu_b;
            ci_lo_fp(end+1) = mu_b - 1.96*sd_b/sqrt(n_b);
            ci_hi_fp(end+1) = mu_b + 1.96*sd_b/sqrt(n_b);
            freq_fp(end+1)  = freq;
            min_fp(end+1)   = min(s.beta_std, [], 'omitnan');
            max_fp(end+1)   = max(s.beta_std, [], 'omitnan');
        end

        if isempty(vars_fp)
            fprintf('  [%s | %s] : aucune variable stable (f > %d%%)\n', bn, gname, THRESH_STAB);
            continue;
        end

        % Trier par beta_std décroissant
        [~, isort]  = sort(mu_fp, 'descend');
        vars_fp     = vars_fp(isort);
        mu_fp       = mu_fp(isort);
        ci_lo_fp    = ci_lo_fp(isort);
        ci_hi_fp    = ci_hi_fp(isort);
        freq_fp     = freq_fp(isort);
        min_fp      = min_fp(isort);
        max_fp      = max_fp(isort);
        N_fp        = length(vars_fp);

        % Figure
        fig_fp = figure('Name', sprintf('ForestPlot_%s_%s', bn, gname), ...
            'NumberTitle','off', ...
            'Position', [50 50 900 max(400, N_fp*80+150)], ...
            'Color','white');
        hold on;

        for iv = 1:N_fp
            y = N_fp - iv + 1;

            % Étendue min/max en ligne bleue fine
            plot([min_fp(iv) max_fp(iv)], [y y], '-', ...
                'Color', [0.15 0.32 0.60], 'LineWidth', 1.2);

            % IC 95% en ligne pleine épaisse ROUGE
            plot([ci_lo_fp(iv) ci_hi_fp(iv)], [y y], '-', ...
                'Color', [0.80 0.10 0.10], 'LineWidth', 4.0);

            % Point central blanc avec contour rouge
            scatter(mu_fp(iv), y, 80, 'white', 'filled', ...
                'MarkerEdgeColor', [0.80 0.10 0.10], 'LineWidth', 1.5);

            % Fréquence — proche du point central, en haut à droite
            text(mu_fp(iv) + 0.03, y + 0.03 , sprintf('%.0f%%', freq_fp(iv)), ...
                'FontSize', 9, 'Color', col_g, ...
                'VerticalAlignment','bottom', 'HorizontalAlignment','left',...
                'FontWeight','bold', 'FontName','Times New Roman');
        end

        % Ligne zéro
        xline(0, '--k', 'LineWidth', 0.9);

        % Axes
        yticks(1:N_fp);
        yticklabels(fliplr(vars_fp));
        xlabel('\beta_{std} (Moy ± IC 95%)', 'FontSize',11, 'FontName','Times New Roman');
        title(sprintf('Forest Plot — %s | %s', bn, gname), ...
            'FontWeight','bold', 'FontSize',12, 'FontName','Times New Roman');
        set(gca, 'FontName','Times New Roman', 'FontSize',9, ...
            'Box','off', 'TickDir','out', 'LineWidth',0.8, ...
            'YGrid','off', 'XGrid','on', ...
            'GridColor',[0.88 0.88 0.88], 'GridAlpha',1);

        % Légende dans coin supérieur gauche
        lx = 0.02; ly = 0.97;
        text(lx, ly, '— IC 95% (rouge)', ...
            'Units','normalized', 'FontSize',8, 'FontWeight','bold',...
            'Color',[0.80 0.10 0.10], 'FontName','Times New Roman',...
            'VerticalAlignment','top');
        text(lx, ly-0.06, '— Étendue min/max — 100 runs (bleu)', ...
            'Units','normalized', 'FontSize',8, ...
            'Color',[0.15 0.32 0.60], 'FontName','Times New Roman',...
            'VerticalAlignment','top');
        text(lx, ly-0.12, '○  β_{std} moyen', ...
            'Units','normalized', 'FontSize',8, ...
            'Color',[0.80 0.10 0.10], 'FontName','Times New Roman',...
            'VerticalAlignment','top');

        % xlim dynamique — basé sur min/max pour tout voir
        margin     = 0.15;
        txt_margin = 0.18;
        xlim([min(min_fp) - margin, max(ci_hi_fp) + margin + txt_margin]);

        % Sauvegarde
        fname = fullfile(PathSave, sprintf('ForestPlot_%s_%s.png', bn, gname));
        exportgraphics(fig_fp, fname, 'BackgroundColor','white', 'Resolution',200);
        close(fig_fp);
        fprintf('  Sauvegardé : ForestPlot_%s_%s.png  (%d variables)\n', bn, gname, N_fp);
    end
end

fprintf('\n=== Terminé ===\n');