%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  ANALYSE IM + GAM — FILTRE EN DEUX ÉTAPES                        %%%%
%%%%                                                                    %%%%
%%%%  Étape 1 (IM) : Information Mutuelle Kraskov (KSG-1) entre chaque %%%%
%%%%  dérive biomécanique et ΔRPE, par groupe (G1,G2,Normal,SD,LD).    %%%%
%%%%  Test de permutation (1000 tirages) + correction FDR (BH) sur     %%%%
%%%%  les 93×5 = 465 p-values. k adaptatif : k=3 si N<15, k=5 sinon.   %%%%
%%%%                                                                    %%%%
%%%%  Étape 2 (GAM) : pour les variables survivant au FDR (q<0.05)     %%%%
%%%%  dans au moins un groupe, ajustement d'un GAM (spline de lissage) %%%%
%%%%  ΔRPE ~ s(dérive_X) pour visualiser la forme de la relation.      %%%%
%%%%                                                                    %%%%
%%%%  Entrée  : Derive_DRPE_Results.mat (produit par                   %%%%
%%%%            Analyse_Derive_G1G2.m, doit contenir les champs raw_*) %%%%
%%%%  Sorties : Tableau_IM_FDR.png (résumé des scores IM + FDR)        %%%%
%%%%            Fig_GAM_<variable>_<groupe>.png (une par variable      %%%%
%%%%            retenue × groupe où elle est significative)            %%%%
%%%%            IM_GAM_Results.mat (sauvegarde complète)               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS & PARAMÈTRES
%  -----------------------------------------------------------------------
PathSave = fileparts(mfilename('fullpath'));
if isempty(PathSave), PathSave = pwd; end

N_PERM       = 1000;     % permutations pour le test contre l'IM aléatoire
FDR_ALPHA    = 0.05;     % seuil FDR (Benjamini-Hochberg), toujours calculé et affiché
SELECTION_MODE = 'raw';  % 'fdr' (corrigé, conservateur) ou 'raw' (p<0.05 brut, EXPLORATOIRE)
RAW_ALPHA    = 0.05;     % seuil utilisé si SELECTION_MODE = 'raw'
K_SMALL      = 3;        % k Kraskov si N < 15
K_LARGE      = 5;        % k Kraskov si N >= 15
N_THRESH_K   = 15;
MIN_N_GROUP  = 6;        % effectif minimum pour calculer l'IM sur un groupe
MAX_GAM_VARS = 10;       % nombre max de variables envoyées au GAM

groupes      = {'G1','G2','Normal','SD','LD'};        % labels d'affichage
groupes_fld  = {'G1','G2','Norm','SD','LD'};           % suffixes des champs raw_* (cohérent avec Analyse_Derive_G1G2.m)
N_groupes = length(groupes);

fprintf('=== Paramètres ===\n');
fprintf('  N_PERM=%d | FDR alpha=%.3f | k=%d (N<%d) ou k=%d (N>=%d)\n', ...
    N_PERM, FDR_ALPHA, K_SMALL, N_THRESH_K, K_LARGE, N_THRESH_K);
fprintf('  MODE DE SÉLECTION : %s', upper(SELECTION_MODE));
if strcmp(SELECTION_MODE, 'raw')
    fprintf(' (p<%.2f brut, SANS correction multiple — EXPLORATOIRE, ~%d faux positifs attendus sous H0 sur 465 tests)\n', ...
        RAW_ALPHA, round(465*RAW_ALPHA));
else
    fprintf(' (q<%.2f après FDR, conservateur)\n', FDR_ALPHA);
end

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
fprintf('\nChargement Derive_DRPE_Results.mat...\n');
load(fullfile(PathSave, 'Derive_DRPE_Results.mat'));   % -> Results_corr

if ~isfield(Results_corr, 'raw_dX_G1')
    error(['Derive_DRPE_Results.mat ne contient pas les champs raw_* (vecteurs ' ...
           'bruts dérive/ΔRPE/RPE0). Ce fichier provient d''une exécution ANTÉRIEURE ' ...
           'du script de calcul des dérives (le script qui produit Derive_DRPE_Results.mat, ' ...
           'quel que soit son nom local sur le disque), avant l''ajout de ces champs. ' ...
           'Relancer ce script une fois pour régénérer le fichier .mat avant de relancer Analyse_IM_GAM.m.']);
end

N_vars = length(Results_corr);
fprintf('  %d variables chargées.\n', N_vars);

%% -----------------------------------------------------------------------
%  ÉTAPE 1 : INFORMATION MUTUELLE (KRASKOV) + TEST DE PERMUTATION
%  -----------------------------------------------------------------------
fprintf('\n=== Étape 1 : Information Mutuelle (Kraskov KSG-1) ===\n');

IM_vals  = NaN(N_vars, N_groupes);
IM_pvals = NaN(N_vars, N_groupes);
IM_N     = NaN(N_vars, N_groupes);

rng(42);  % reproductibilité des permutations

t_start_global = tic;
for iv = 1:N_vars
    r = Results_corr(iv);
    if mod(iv, 10) == 0 || iv == 1
        fprintf('  Variable %d/%d (%s, %s)...\n', iv, N_vars, r.vname, r.source);
    end

    for ig = 1:N_groupes
        gn = groupes_fld{ig};
        dX   = r.(['raw_dX_'   gn]);
        dRPE = r.(['raw_dRPE_' gn]);

        ok = ~isnan(dX) & ~isnan(dRPE) & ~isinf(dX) & ~isinf(dRPE);
        N_ok = sum(ok);
        IM_N(iv, ig) = N_ok;
        if N_ok < MIN_N_GROUP, continue; end

        x_ok = dX(ok); y_ok = dRPE(ok);

        k_use = K_SMALL;
        if N_ok >= N_THRESH_K, k_use = K_LARGE; end
        % k doit rester < N pour que l'estimateur soit défini
        k_use = min(k_use, N_ok - 2);
        if k_use < 1, continue; end

        im_obs = kraskov_mi(x_ok, y_ok, k_use);
        IM_vals(iv, ig) = im_obs;

        % --- Test de permutation : mélange y, recalcule l'IM N_PERM fois ---
        im_perm = NaN(N_PERM, 1);
        for ip = 1:N_PERM
            y_perm = y_ok(randperm(N_ok));
            im_perm(ip) = kraskov_mi(x_ok, y_perm, k_use);
        end
        IM_pvals(iv, ig) = (sum(im_perm >= im_obs) + 1) / (N_PERM + 1);
    end
end
fprintf('  Temps total IM + permutations : %.1f s\n', toc(t_start_global));

%% -----------------------------------------------------------------------
%  CORRECTION FDR (Benjamini-Hochberg) sur les 93×5 = 465 p-values
%  -----------------------------------------------------------------------
fprintf('\n=== Correction FDR (Benjamini-Hochberg) ===\n');

p_flat = IM_pvals(:);
valid_idx = ~isnan(p_flat);
p_valid = p_flat(valid_idx);

[adj_p, ~] = bh_fdr_local(p_valid, FDR_ALPHA);

q_flat = NaN(size(p_flat));
q_flat(valid_idx) = adj_p;
IM_qvals = reshape(q_flat, size(IM_pvals));

n_raw_sig = sum(p_valid < 0.05);
n_fdr_sig = sum(adj_p < FDR_ALPHA);
fprintf('  Tests valides : %d / %d\n', sum(valid_idx), numel(p_flat));
fprintf('  Significatifs avant FDR (p<0.05) : %d\n', n_raw_sig);
fprintf('  Significatifs après  FDR (q<%.2f) : %d\n', FDR_ALPHA, n_fdr_sig);

%% -----------------------------------------------------------------------
%  SÉLECTION DES VARIABLES SURVIVANTES
%  Mode 'fdr' : au moins un groupe avec q < FDR_ALPHA (conservateur)
%  Mode 'raw' : au moins un groupe avec p < RAW_ALPHA (EXPLORATOIRE —
%               aucune correction pour les 465 tests multiples ; sur ce
%               nombre de tests, ~23 faux positifs sont attendus sous H0
%               avec un seuil de 0.05, donc les variables retenues ici ne
%               sont PAS des résultats confirmés, seulement des pistes à
%               visualiser et à valider sur un jeu de données indépendant)
%  -----------------------------------------------------------------------
if strcmp(SELECTION_MODE, 'fdr')
    SEL_CRIT = IM_qvals;
    SEL_ALPHA = FDR_ALPHA;
else
    SEL_CRIT = IM_pvals;
    SEL_ALPHA = RAW_ALPHA;
end

survives_any = any(SEL_CRIT < SEL_ALPHA, 2);
idx_survive  = find(survives_any);

fprintf('\n=== Variables survivant au filtre IM (mode %s, seuil %.3f) : %d ===\n', ...
    upper(SELECTION_MODE), SEL_ALPHA, length(idx_survive));
if strcmp(SELECTION_MODE, 'raw') && length(idx_survive) > 0
    fprintf('  ATTENTION : sélection EXPLORATOIRE sans correction multiple.\n');
end
for ii = 1:length(idx_survive)
    iv = idx_survive(ii);
    sig_groups = groupes(SEL_CRIT(iv,:) < SEL_ALPHA);
    fprintf('  %-40s (%-10s) — significatif dans : %s\n', ...
        Results_corr(iv).vname, Results_corr(iv).source, strjoin(sig_groups, ', '));
end

if length(idx_survive) > MAX_GAM_VARS
    % Si plus de MAX_GAM_VARS survivent, prioriser par critère minimal
    % (q-value en mode FDR, p-value brute en mode raw)
    [~, ord] = sort(min(SEL_CRIT(idx_survive,:), [], 2));
    idx_survive = idx_survive(ord(1:MAX_GAM_VARS));
    fprintf('  --> Limité aux %d variables avec le critère le plus faible pour le GAM.\n', MAX_GAM_VARS);
end

%% -----------------------------------------------------------------------
%  TABLEAU RÉCAPITULATIF IM + FDR (PNG)
%  -----------------------------------------------------------------------
fprintf('\n=== Génération du tableau récapitulatif ===\n');
drawIMTable(Results_corr, idx_survive, IM_vals, SEL_CRIT, SEL_ALPHA, SELECTION_MODE, IM_N, groupes, PathSave);

%% -----------------------------------------------------------------------
%  ÉTAPE 2 : GAM SUR LES VARIABLES RETENUES
%  -----------------------------------------------------------------------
fprintf('\n=== Étape 2 : Ajustement GAM sur les variables retenues ===\n');

if isempty(idx_survive)
    fprintf('  Aucune variable n''a survécu au filtre IM (mode %s) — pas de GAM à ajuster.\n', upper(SELECTION_MODE));
else
    if ~exist('fitrgam', 'file')
        fprintf(['  ATTENTION : fitrgam introuvable (nécessite Statistics and ' ...
                 'Machine Learning Toolbox, R2021a+). Le GAM ne peut pas être ajusté ' ...
                 'dans cet environnement — le filtre IM est néanmoins sauvegardé.\n']);
    else
        for ii = 1:length(idx_survive)
            iv = idx_survive(ii);
            r  = Results_corr(iv);
            sig_groups_idx = find(SEL_CRIT(iv,:) < SEL_ALPHA);

            for ig = sig_groups_idx
                gn = groupes_fld{ig};
                gn_lbl = groupes{ig};
                dX   = r.(['raw_dX_'   gn]);
                dRPE = r.(['raw_dRPE_' gn]);
                ok = ~isnan(dX) & ~isnan(dRPE) & ~isinf(dX) & ~isinf(dRPE);
                if sum(ok) < 8, continue; end  % GAM peu fiable sous N=8

                x_ok = dX(ok); y_ok = dRPE(ok);

                try
                    gam_mdl = fitrgam(x_ok, y_ok);

                    % Grille de prédiction pour visualiser la forme
                    x_grid = linspace(min(x_ok), max(x_ok), 100)';
                    y_pred = predict(gam_mdl, x_grid);

                    if strcmp(SELECTION_MODE, 'raw')
                        crit_display = IM_pvals(iv,ig);
                    else
                        crit_display = IM_qvals(iv,ig);
                    end
                    drawGAMFigure(x_ok, y_ok, x_grid, y_pred, r.vname, r.source, gn_lbl, ...
                        IM_vals(iv,ig), crit_display, SELECTION_MODE, PathSave);
                catch ME
                    fprintf('  GAM échoué pour %s [%s] : %s\n', r.vname, gn, ME.message);
                end
            end
        end
    end
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'IM_GAM_Results.mat'), ...
    'IM_vals','IM_pvals','IM_qvals','IM_N','groupes','idx_survive', ...
    'N_PERM','FDR_ALPHA','K_SMALL','K_LARGE','N_THRESH_K');
fprintf('\nSauvegardé : IM_GAM_Results.mat\n');
fprintf('\n=== PIPELINE TERMINÉ ===\n');

%% =========================================================================
%% FONCTIONS LOCALES
%% =========================================================================

function im = kraskov_mi(x, y, k)
% Estimateur Kraskov-Stögbauer-Grassberger (KSG), variante 1.
% Référence : Kraskov, Stögbauer & Grassberger (2004), Phys. Rev. E 69.
%
% Principe : pour chaque point i, trouve la distance epsilon_i au k-ième
% plus proche voisin dans l'espace joint (x,y) (norme infinie / Chebyshev).
% Compte ensuite, dans les espaces marginaux x et y séparément, le nombre
% de points à distance < epsilon_i (n_x(i), n_y(i)). L'IM s'estime par :
%   I(X;Y) = psi(k) - <psi(n_x+1) + psi(n_y+1)> + psi(N)
% où psi est la fonction digamma et <.> la moyenne sur tous les points i.

    x = x(:); y = y(:);
    N = length(x);
    if N <= k+1
        im = NaN; return;
    end

    % Normalisation (évite qu'une variable à grande échelle domine la
    % distance jointe) — sans incidence sur l'IM théorique car KSG utilise
    % les rangs de distance, mais améliore la stabilité numérique.
    x = (x - mean(x)) / (std(x) + eps);
    y = (y - mean(y)) / (std(y) + eps);

    n_x = zeros(N,1);
    n_y = zeros(N,1);

    for i = 1:N
        % Distance de Chebyshev (norme infinie) dans l'espace joint (x,y)
        dx = abs(x - x(i));
        dy = abs(y - y(i));
        d_joint = max(dx, dy);
        d_joint(i) = Inf;  % exclut le point lui-même

        d_sorted = sort(d_joint);
        eps_i = d_sorted(k);  % distance au k-ième plus proche voisin (strict <)

        % Comptage dans les marginales : points à distance < eps_i
        % (strictement, convention KSG-1, exclut le point i lui-même)
        n_x(i) = sum(dx < eps_i) - 1;
        n_y(i) = sum(dy < eps_i) - 1;
        n_x(i) = max(n_x(i), 0);
        n_y(i) = max(n_y(i), 0);
    end

    im = psi_local(k) - mean(psi_local(n_x+1) + psi_local(n_y+1)) + psi_local(N);
    im = max(im, 0);  % l'IM théorique est >=0 ; les estimateurs bruités
                       % peuvent légèrement dépasser sous zéro — on tronque
end

% -------------------------------------------------------------------------
function y = psi_local(x)
% Fonction digamma — wrapper autour de psi() si disponible (Symbolic Math
% Toolbox / Statistics), sinon approximation d'Euler-Maclaurin (précise
% pour x entiers positifs, le seul cas utilisé ici).
    x = x(:);
    y = zeros(size(x));
    for i = 1:length(x)
        xi = x(i);
        if xi <= 0
            y(i) = -Inf; continue;
        end
        % Récurrence digamma : psi(x) = psi(x+1) - 1/x, on remonte vers x>=6
        % pour que l'approximation asymptotique soit précise.
        val = 0;
        while xi < 6
            val = val - 1/xi;
            xi = xi + 1;
        end
        % Approximation asymptotique (Abramowitz & Stegun 6.3.18)
        val = val + log(xi) - 1/(2*xi) - 1/(12*xi^2) + 1/(120*xi^4) - 1/(252*xi^6);
        y(i) = val;
    end
end

% -------------------------------------------------------------------------
function [adj_p, h] = bh_fdr_local(p_vals, alpha)
% Correction FDR de Benjamini-Hochberg, implémentation locale (cf. les
% scripts SpearmanCorr_*_FDR.m du jeu 1, même algorithme step-up).
    n = numel(p_vals);
    [p_sorted, sort_idx] = sort(p_vals(:)');
    adj_sorted = p_sorted;
    adj_sorted(n) = p_sorted(n);
    for i = n-1:-1:1
        adj_sorted(i) = min(p_sorted(i) * n / i, adj_sorted(i+1));
    end
    adj_sorted = min(adj_sorted, 1);
    adj_p = NaN(1, n);
    adj_p(sort_idx) = adj_sorted;
    h = adj_p < alpha;
end

% -------------------------------------------------------------------------
function drawIMTable(Results_corr, idx_survive, IM_vals, SEL_CRIT, SEL_ALPHA, SELECTION_MODE, IM_N, groupes, PathSave)
% Tableau récapitulatif : une ligne par variable survivante, colonnes =
% score IM par groupe. SEL_CRIT est IM_pvals (mode 'raw') ou IM_qvals
% (mode 'fdr') — déjà choisi par l'appelant.

    N_v = max(length(idx_survive), 1);
    row_h_px = 20; title_h_px = 30; header_h_px = 24; legend_h_px = 18; margin_px = 6;
    fig_h = title_h_px + header_h_px + N_v*row_h_px + legend_h_px + 2*margin_px;
    fig_w = 1000;

    f_title_h  = title_h_px  / fig_h;
    f_header_h = header_h_px / fig_h;
    f_row_h    = row_h_px    / fig_h;
    f_legend_h = legend_h_px / fig_h;
    f_margin   = margin_px   / fig_h;

    y_title_top = 1 - f_margin;
    y_title_bot = y_title_top - f_title_h;
    y_header_top = y_title_bot;
    y_header_bot = y_header_top - f_header_h;
    y_table_top  = y_header_bot;

    fig = figure('Visible','off','Position',[50 50 fig_w fig_h],'Color','white');
    ax  = axes('Position',[0 0 1 1],'Visible','off');
    set(ax,'XLim',[0 1],'YLim',[0 1]); hold on;

    if strcmp(SELECTION_MODE, 'raw')
        title_str = 'Information Mutuelle (Kraskov) — sélection EXPLORATOIRE (p brut, non corrigé)';
        title_bg = [0.55 0.20 0.20];
    else
        title_str = 'Information Mutuelle (Kraskov) — variables survivant au FDR';
        title_bg = [0.12 0.22 0.42];
    end
    rectangle('Position',[0 y_title_bot 1 f_title_h],'FaceColor',title_bg,'EdgeColor','none');
    text(0.5, (y_title_top+y_title_bot)/2, title_str, ...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',10,'FontWeight','bold','Color','white');

    col_x = [0.01, 0.32, 0.44, 0.56, 0.68, 0.80];
    col_w = [0.30, 0.11, 0.11, 0.11, 0.11, 0.11];
    hdrs  = {'Variable','G1','G2','Normal','SD','LD'};

    rectangle('Position',[0 y_header_bot 1 f_header_h],'FaceColor',[0.78 0.80 0.88],'EdgeColor','none');
    for ic = 1:6
        text(col_x(ic)+col_w(ic)/2, (y_header_top+y_header_bot)/2, hdrs{ic}, ...
            'FontSize',8,'FontWeight','bold','Color',[0.1 0.1 0.3], ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    line([0 1],[y_header_bot y_header_bot],'Color',[0.5 0.5 0.6],'LineWidth',1);

    if ~isempty(idx_survive)
        for ii = 1:length(idx_survive)
            iv = idx_survive(ii);
            y_row = y_table_top - ii * f_row_h;
            y_txt = y_row + f_row_h/2;
            if mod(ii,2)==0
                rectangle('Position',[0 y_row 1 f_row_h],'FaceColor',[0.96 0.96 0.98],'EdgeColor','none');
            end
            vn_disp = strrep(Results_corr(iv).vname, '_', ' ');
            if length(vn_disp) > 42, vn_disp = [vn_disp(1:39) '...']; end
            text(col_x(1)+0.003, y_txt, vn_disp, ...
                'FontSize',6.5,'VerticalAlignment','middle','Interpreter','none','Color',[0.1 0.1 0.1]);

            for ig = 1:5
                im_v = IM_vals(iv,ig); crit_v = SEL_CRIT(iv,ig); n_v = IM_N(iv,ig);
                if isnan(im_v) || n_v < 6, continue; end
                sig = ''; fw = 'normal'; txt_col = [0.3 0.3 0.3];
                if ~isnan(crit_v) && crit_v < SEL_ALPHA
                    sig = '†'; fw = 'bold'; txt_col = [0.75 0.10 0.10];
                    bg = [0.95 0.85 0.85];
                    rectangle('Position',[col_x(ig+1) y_row col_w(ig+1) f_row_h],'FaceColor',bg,'EdgeColor','none');
                end
                txt = sprintf('%.3f%s', im_v, sig);
                text(col_x(ig+1)+col_w(ig+1)/2, y_txt, txt, ...
                    'FontSize',7,'FontWeight',fw,'Color',txt_col, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle');
            end
            line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
        end
    else
        if strcmp(SELECTION_MODE, 'raw')
            empty_msg = 'Aucune variable significative (p<0.05 brut)';
        else
            empty_msg = 'Aucune variable significative après correction FDR';
        end
        text(0.5, y_table_top - f_row_h, empty_msg, ...
            'FontSize',8,'Color',[0.5 0.5 0.5],'HorizontalAlignment','center');
    end

    line([0 1 1 0 0],[f_margin f_margin y_title_top y_title_top f_margin],'Color',[0.4 0.4 0.5],'LineWidth',1);
    if strcmp(SELECTION_MODE, 'raw')
        legend_txt = sprintf('† p<%.2f BRUT, non corrigé (EXPLORATOIRE) — valeur = IM Kraskov (nats)', SEL_ALPHA);
    else
        legend_txt = sprintf('† q<%.2f (FDR Benjamini-Hochberg) — valeur = IM Kraskov (nats)', SEL_ALPHA);
    end
    text(0.5, f_margin/2, legend_txt, ...
        'FontSize',6,'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);

    hold off;
    exportgraphics(fig, fullfile(PathSave,'Tableau_IM_FDR.png'), 'Resolution', 150);
    close(fig);
    fprintf('  Tableau sauvegardé : Tableau_IM_FDR.png\n');
end

% -------------------------------------------------------------------------
function drawGAMFigure(x_obs, y_obs, x_grid, y_pred, vname, source, gname, im_val, q_val, selection_mode, PathSave)
    fig = figure('Visible','off','Position',[50 50 700 500],'Color','white');
    hold on;
    scatter(x_obs, y_obs, 40, [0.2 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
    plot(x_grid, y_pred, '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 2.5);
    xlabel('Dérive (X_{fin} - X_{début})', 'FontSize', 11);
    ylabel('\Delta RPE', 'FontSize', 11);

    if strcmp(selection_mode, 'raw')
        subtitle_txt = sprintf('IM=%.3f nats | p brut=%.4f — EXPLORATOIRE, non corrigé (465 tests)', im_val, q_val);
        title_color = [0.75 0.15 0.15];
    else
        subtitle_txt = sprintf('IM=%.3f nats | FDR q=%.4f', im_val, q_val);
        title_color = [0 0 0];
    end

    title({sprintf('%s — %s', strrep(vname,'_',' '), gname), subtitle_txt}, ...
        'FontSize', 10, 'Interpreter', 'none', 'Color', title_color);
    grid on; box on; hold off;

    fname = sprintf('Fig_GAM_%s_%s_%s.png', source, vname, gname);
    fname = regexprep(fname, '[^a-zA-Z0-9_.]', '_');
    exportgraphics(fig, fullfile(PathSave, fname), 'Resolution', 150);
    close(fig);
    fprintf('  Figure GAM sauvegardée : %s\n', fname);
end