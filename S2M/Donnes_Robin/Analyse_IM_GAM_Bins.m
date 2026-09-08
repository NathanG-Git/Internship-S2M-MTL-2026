%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  ANALYSE GAMM PAR BINS — EFFET ALÉATOIRE SUJET                   %%%%
%%%%                                                                    %%%%
%%%%  Plutôt qu'un seul point dérive par sujet (N=8 à 26 selon groupe, %%%%
%%%%  cause du surajustement du GAM observé dans Analyse_IM_GAM.m),    %%%%
%%%%  on empile les différences consécutives bin-à-bin (lag-1) de TOUS %%%%
%%%%  les sujets d'un groupe (~N_sujets × 9 observations), avec un      %%%%
%%%%  effet aléatoire (1|Sujet) pour respecter la non-indépendance des %%%%
%%%%  observations intra-sujet :                                       %%%%
%%%%                                                                    %%%%
%%%%    ΔRPE(t) ~ base_spline(ΔX(t)) + (1|Sujet)      [fitglme]        %%%%
%%%%                                                                    %%%%
%%%%  Base de spline : B-splines cubiques, 5 nœuds internes (quantiles %%%%
%%%%  de ΔX), approximant un GAM avec flexibilité contrôlée — un seul  %%%%
%%%%  paramètre (nb de nœuds) plutôt que la flexibilité libre de       %%%%
%%%%  fitrgam, qui surajustait sur les 8-26 points dérive uniques.     %%%%
%%%%                                                                    %%%%
%%%%  RPE jeu 2 (Normal/SD/LD) : interpolé linéairement entre début et  %%%%
%%%%  fin (cf. extractBinsJeu2 dans Analyse_Derive_G1G2.m).             %%%%
%%%%                                                                    %%%%
%%%%  BRANCHEMENT G1/G2 vs Normal/SD/LD (cf. RCOND~1e-16/-18 et AIC     %%%%
%%%%  quasi constant entre variables observés juin 2026) :             %%%%
%%%%  le delta lag-1 d'une fonction interpolée LINÉAIREMENT est une     %%%%
%%%%  CONSTANTE — donc ΔRPE(t) est identique pour tous les bins t d'un  %%%%
%%%%  même sujet sur Normal/SD/LD. Le GAMM lag-1 + (1|Sujet) y est      %%%%
%%%%  alors structurellement dégénéré : l'effet aléatoire absorbe       %%%%
%%%%  exactement cette constante, rendant les coefficients spline non   %%%%
%%%%  identifiables (matrice quasi-singulière, p_global/AIC artefacts). %%%%
%%%%  Pour ces trois groupes, on revient donc à un point par sujet :    %%%%
%%%%  corrélation de Spearman entre dérive totale de X sur la session   %%%%
%%%%  (dernier bin valide - premier bin valide) et ΔRPE total (même     %%%%
%%%%  logique que Analyse_IM_GAM.m, N=8-26). G1/G2 (RPE mesuré          %%%%
%%%%  réellement bin par bin, pas d'interpolation) gardent le GAMM      %%%%
%%%%  lag-1 inchangé.                                                   %%%%
%%%%                                                                    %%%%
%%%%  Entrée  : Derive_DRPE_Results.mat (doit contenir les champs       %%%%
%%%%            bins_X_*/bins_RPE_* — relancer Analyse_Derive_G1G2.m   %%%%
%%%%            version à jour si absent)                              %%%%
%%%%  Sortie  : Fig_GAMM_<bloc>_<variable>_<groupe>.png  (G1/G2)        %%%%
%%%%            Fig_Drift_<bloc>_<variable>_<groupe>.png (Normal/SD/LD) %%%%
%%%%            GAMM_Bins_Results.mat                                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS & PARAMÈTRES
%  -----------------------------------------------------------------------
PathSave = fileparts(mfilename('fullpath'));
if isempty(PathSave), PathSave = pwd; end

N_SPLINE_KNOTS = 2;      % nœuds internes de la base de spline (cf. discussion) — réduit de 5 à 2 (q=9->6) pour le ratio q/G du test robuste en cluster (G=13-25 sujets)
MIN_N_SUBJ     = 6;      % nombre minimum de sujets pour ajuster un GAMM
MIN_N_OBS      = 20;     % nombre minimum d'observations (sujets x bins) après lag-1
VARS_TO_FIT    = {};     % si vide : utilise idx_survive d'Analyse_IM_GAM.m (s'il existe)
                          % sinon : liste manuelle de noms de variables à ajuster

groupes      = {'G1','G2','Normal','SD','LD'};
groupes_fld  = {'G1','G2','Norm','SD','LD'};   % suffixes des champs bins_* (cohérent avec Analyse_Derive_G1G2.m)
groupes_RPE_interp = {'Normal','SD','LD'};     % RPE interpolé linéairement -> branche corrélation drift (pas de GAMM lag-1, cf. en-tête)
N_groupes    = length(groupes);

fprintf('=== Paramètres ===\n');
fprintf('  N_SPLINE_KNOTS=%d | MIN_N_SUBJ=%d | MIN_N_OBS=%d\n', N_SPLINE_KNOTS, MIN_N_SUBJ, MIN_N_OBS);

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
fprintf('\nChargement Derive_DRPE_Results.mat...\n');
load(fullfile(PathSave, 'Derive_DRPE_Results.mat'));   % -> Results_corr

if ~isfield(Results_corr, 'bins_X_G1')
    error(['Derive_DRPE_Results.mat ne contient pas les champs bins_* (séries ' ...
           'complètes par bin). Ce fichier provient d''une exécution ANTÉRIEURE ' ...
           'du script de calcul des dérives, avant l''ajout de ces champs. ' ...
           'Relancer ce script (celui qui produit Derive_DRPE_Results.mat) une fois ' ...
           'avant de relancer Analyse_GAMM_Bins.m.']);
end

N_vars = length(Results_corr);
fprintf('  %d variables chargées.\n', N_vars);

% --- Sélection des variables à ajuster ---
if isempty(VARS_TO_FIT)
    fname_im = fullfile(PathSave, 'IM_GAM_Results.mat');
    if exist(fname_im, 'file')
        S = load(fname_im, 'idx_survive');
        idx_to_fit = S.idx_survive;
        fprintf('  %d variables reprises depuis IM_GAM_Results.mat (idx_survive).\n', length(idx_to_fit));
    else
        warning(['IM_GAM_Results.mat introuvable et VARS_TO_FIT vide — ' ...
                  'ajustement sur TOUTES les %d variables (peut être long).'], N_vars);
        idx_to_fit = 1:N_vars;
    end
else
    idx_to_fit = find(ismember({Results_corr.vname}, VARS_TO_FIT));
    fprintf('  %d variables sélectionnées manuellement via VARS_TO_FIT.\n', length(idx_to_fit));
end

if isempty(idx_to_fit)
    fprintf('Aucune variable à ajuster — fin du script.\n');
    return;
end

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE : GAMM par variable × groupe
%  -----------------------------------------------------------------------
fprintf('\n=== Ajustement GAMM (splines + effet aléatoire sujet) ===\n');

GAMM_Results = struct([]);
i_res = 0;

for ii = 1:length(idx_to_fit)
    iv = idx_to_fit(ii);
    r  = Results_corr(iv);
    fprintf('\nVariable %d/%d : %s (%s)\n', ii, length(idx_to_fit), r.vname, r.source);

    for ig = 1:N_groupes
        gn_fld = groupes_fld{ig};
        gn_lbl = groupes{ig};

        X_bins   = r.(['bins_X_'   gn_fld]);
        RPE_bins = r.(['bins_RPE_' gn_fld]);

        if isempty(X_bins) || isempty(RPE_bins), continue; end

        [N_subj_loc, N_bins_loc] = size(X_bins);
        if N_subj_loc < MIN_N_SUBJ || N_bins_loc < 3, continue; end

        if ismember(gn_lbl, groupes_RPE_interp)
            % --- Branche corrélation drift (RPE interpolé -> GAMM lag-1 dégénéré, cf. en-tête) ---
            try
                [rho, pval, N_subj_ok, dX_drift, dRPE_drift] = computeDriftCorrelation(X_bins, RPE_bins, MIN_N_SUBJ);
                if isnan(rho)
                    fprintf('  [%s] Corrélation drift ignorée (N_subj ou variance insuffisante)\n', gn_lbl);
                    continue;
                end

                i_res = i_res + 1;
                GAMM_Results(i_res).vname            = r.vname;
                GAMM_Results(i_res).source           = r.source;
                GAMM_Results(i_res).group            = gn_lbl;
                GAMM_Results(i_res).N_subj           = N_subj_ok;
                GAMM_Results(i_res).N_obs            = N_subj_ok;
                GAMM_Results(i_res).method           = 'Spearman_drift';
                GAMM_Results(i_res).p_global         = pval;
                GAMM_Results(i_res).AIC              = NaN;
                GAMM_Results(i_res).RPE_interpolated = true;
                GAMM_Results(i_res).Spearman_rho     = rho;
                GAMM_Results(i_res).resid_autocorr_r = NaN;   % non pertinent : pas de GAMM lag-1 ici
                GAMM_Results(i_res).resid_autocorr_p = NaN;
                GAMM_Results(i_res).resid_autocorr_N = NaN;
                GAMM_Results(i_res).p_global_notime = NaN;
                GAMM_Results(i_res).p_time = NaN;
                GAMM_Results(i_res).AIC_notime = NaN;
                GAMM_Results(i_res).resid_autocorr_r_notime = NaN;
                GAMM_Results(i_res).resid_autocorr_p_notime = NaN;
                GAMM_Results(i_res).resid_autocorr_N_notime = NaN;
                GAMM_Results(i_res).p_global_modelbased = NaN;
                GAMM_Results(i_res).F_robust   = NaN;
                GAMM_Results(i_res).df1_robust = NaN;
                GAMM_Results(i_res).df2_robust = NaN;

                fprintf('  [%s] N_subj=%d, rho_Spearman=%.3f, p=%.4f [point unique/sujet — RPE interpolé]\n', ...
                    gn_lbl, N_subj_ok, rho, pval);

                drawDriftFigure(dX_drift, dRPE_drift, rho, pval, r.vname, r.source, gn_lbl, N_subj_ok, PathSave);

            catch ME
                fprintf('  [%s] Corrélation drift échouée : %s\n', gn_lbl, ME.message);
            end
            continue;   % groupe traité (branche drift) -> on saute le GAMM lag-1 ci-dessous
        end

        % --- Construction des différences lag-1 empilées, avec ID sujet ---
        % dX(t) = X(t) - X(t-1), dRPE(t) = RPE(t) - RPE(t-1), pour t=2..N_bins
        N_obs_max = N_subj_loc * (N_bins_loc - 1);
        dX_all   = NaN(N_obs_max, 1);
        dRPE_all = NaN(N_obs_max, 1);
        subj_all = NaN(N_obs_max, 1);
        bin_all  = NaN(N_obs_max, 1);   % indice de bin t -> permet de vérifier la consécutivité réelle pour le diagnostic d'autocorrélation des résidus
        row = 0;
        for is = 1:N_subj_loc
            x_s   = X_bins(is, :);
            rpe_s = RPE_bins(is, :);
            if all(isnan(x_s)) || all(isnan(rpe_s)), continue; end
            for it = 2:N_bins_loc
                row = row + 1;
                dX_all(row)   = x_s(it) - x_s(it-1);
                dRPE_all(row) = rpe_s(it) - rpe_s(it-1);
                subj_all(row) = is;
                bin_all(row)  = it;
            end
        end
        dX_all   = dX_all(1:row);
        dRPE_all = dRPE_all(1:row);
        subj_all = subj_all(1:row);
        bin_all  = bin_all(1:row);

        ok = ~isnan(dX_all) & ~isnan(dRPE_all) & ~isinf(dX_all) & ~isinf(dRPE_all);
        dX_ok   = dX_all(ok);
        dRPE_ok = dRPE_all(ok);
        subj_ok = subj_all(ok);
        bin_ok  = bin_all(ok);

        N_obs_ok   = length(dX_ok);
        N_subj_ok  = length(unique(subj_ok));
        if N_obs_ok < MIN_N_OBS || N_subj_ok < MIN_N_SUBJ, continue; end
        if std(dX_ok) < 1e-12, continue; end  % variable constante, spline impossible

        try
            % --- Construction de la base de spline (B-splines cubiques) ---
            [B, knots_used] = buildSplineBasis(dX_ok, N_SPLINE_KNOTS);
            N_basis = size(B, 2);

            % --- Table pour fitglme : colonnes spline + Subject + dRPE ---
            varnames_B = arrayfun(@(k) sprintf('B%d', k), 1:N_basis, 'UniformOutput', false);
            Tbl = array2table(B, 'VariableNames', varnames_B);
            Tbl.dRPE = dRPE_ok;
            Tbl.Subject = categorical(subj_ok);
            Tbl.Time = zscore(bin_ok);   % indice de bin centré-réduit : covariable de contrôle pour une dérive temporelle commune à dX et dRPE (cf. autocorrélation résidus systématique observée)

            % Pas d'intercept explicite : la base de B-spline complète
            % vérifie sum(B,2)=1 partout (partition de l'unité), donc un
            % intercept séparé serait exactement colinéaire avec la somme
            % des colonnes spline -> matrice de design non de rang plein
            % (erreur "full column rank" rencontrée et corrigée ici).

            % --- Modèle A (référence, SANS contrôle du temps — celui d'avant) ---
            formula_raw = ['dRPE ~ -1 + ' strjoin(varnames_B, ' + ') ' + (1|Subject)'];
            gamm_mdl_raw = fitglme(Tbl, formula_raw, 'Distribution', 'normal');
            p_global_raw = waldTestSplineCoefs(gamm_mdl_raw);
            resid_raw = residuals(gamm_mdl_raw);
            [acr_raw, acp_raw, acN_raw] = checkResidAutocorr(resid_raw, subj_ok, bin_ok);

            % --- Modèle B (corrigé : + terme linéaire Time) ---
            % H0 testée pour p_global : tous les coefficients spline B# sont nuls,
            % EN CONTRÔLANT pour une dérive temporelle linéaire commune (Time).
            % Si le signal du modèle A disparaît ici, dX et dRPE partageaient une
            % dérive temporelle commune plutôt qu'un lien direct dX -> dRPE.
            formula_time = ['dRPE ~ -1 + Time + ' strjoin(varnames_B, ' + ') ' + (1|Subject)'];
            gamm_mdl_time = fitglme(Tbl, formula_time, 'Distribution', 'normal');
            p_global_modelbased = waldTestSplineCoefs(gamm_mdl_time);   % covariance NON robuste (suppose le modèle de covariance bien spécifié)
            resid_time = residuals(gamm_mdl_time);
            [resid_ac_r, resid_ac_p, resid_ac_N] = checkResidAutocorr(resid_time, subj_ok, bin_ok);

            time_coef_idx = find(strcmp(gamm_mdl_time.CoefficientNames, 'Time'));
            p_time = NaN;
            if ~isempty(time_coef_idx)
                p_time = gamm_mdl_time.Coefficients.pValue(time_coef_idx);
            end

            % --- Test de Wald robuste en cluster (sujet), F(q, G-1) (Cameron & Miller) ---
            % La covariance de coefTest ci-dessus suppose le modèle de covariance du
            % GAMM (intercept aléatoire + résidu homogène) correctement spécifié.
            % L'autocorrélation résiduelle qui persiste (surtout G2) indique que ce
            % n'est pas le cas -> covariance "sandwich" robuste en cluster (sujet),
            % test F(q,G-1), référence recommandée à G modéré plutôt qu'un chi².
            H_time = findSplineContrastH(gamm_mdl_time);
            [F_robust, p_global, df1_robust, df2_robust] = clusterRobustWaldTest(gamm_mdl_time, H_time, subj_ok, dRPE_ok);

            % --- Grille de prédiction pour visualiser la forme (modèle corrigé, Time fixé à sa moyenne = 0 car z-scoré) ---
            x_grid = linspace(min(dX_ok), max(dX_ok), 100)';
            B_grid = evalSplineBasis(x_grid, knots_used);
            Tbl_grid = array2table(B_grid, 'VariableNames', varnames_B);
            Tbl_grid.Time = zeros(100, 1);
            Tbl_grid.Subject = repmat(categorical(subj_ok(1)), 100, 1);  % valeur arbitraire, effet aléatoire neutralisé via predict 'Conditional',false
            y_grid = predict(gamm_mdl_time, Tbl_grid, 'Conditional', false);

            % Toujours faux ici : Normal/SD/LD sont désormais interceptés plus haut
            % (branche corrélation drift) avant d'atteindre ce bloc GAMM lag-1.
            is_interpolated_rpe = false;

            i_res = i_res + 1;
            GAMM_Results(i_res).vname    = r.vname;
            GAMM_Results(i_res).source   = r.source;
            GAMM_Results(i_res).group    = gn_lbl;
            GAMM_Results(i_res).N_subj   = N_subj_ok;
            GAMM_Results(i_res).N_obs    = N_obs_ok;
            GAMM_Results(i_res).p_global = p_global;                          % test ROBUSTE en cluster (F(q,G-1)) -> alimente le FDR
            GAMM_Results(i_res).p_global_modelbased = p_global_modelbased;    % test de Wald model-based (+Time, covariance non-robuste)
            GAMM_Results(i_res).F_robust   = F_robust;
            GAMM_Results(i_res).df1_robust = df1_robust;
            GAMM_Results(i_res).df2_robust = df2_robust;
            GAMM_Results(i_res).p_global_notime = p_global_raw;       % test brut (sans Time, sans robustesse), conservé pour comparaison
            GAMM_Results(i_res).p_time   = p_time;                    % significativité de la dérive temporelle elle-même
            GAMM_Results(i_res).AIC      = gamm_mdl_time.ModelCriterion.AIC;
            GAMM_Results(i_res).AIC_notime = gamm_mdl_raw.ModelCriterion.AIC;
            GAMM_Results(i_res).RPE_interpolated = is_interpolated_rpe;
            GAMM_Results(i_res).method   = 'GAMM_lag1_timeadj_robust';
            GAMM_Results(i_res).Spearman_rho = NaN;
            GAMM_Results(i_res).resid_autocorr_r = resid_ac_r;            % autocorr résidus, modèle CORRIGÉ
            GAMM_Results(i_res).resid_autocorr_p = resid_ac_p;
            GAMM_Results(i_res).resid_autocorr_N = resid_ac_N;
            GAMM_Results(i_res).resid_autocorr_r_notime = acr_raw;        % autocorr résidus, modèle brut (référence)
            GAMM_Results(i_res).resid_autocorr_p_notime = acp_raw;
            GAMM_Results(i_res).resid_autocorr_N_notime = acN_raw;

            fprintf('  [%s] N_subj=%d, N_obs=%d\n', gn_lbl, N_subj_ok, N_obs_ok);
            fprintf('       SANS contrôle temps                   : p_global=%.4f, AIC=%.1f, autocorr_resid(lag1) r=%.3f p=%.4f (N=%d)\n', ...
                p_global_raw, gamm_mdl_raw.ModelCriterion.AIC, acr_raw, acp_raw, acN_raw);
            fprintf('       AVEC temps (model-based, non-robuste) : p_global=%.4f, AIC=%.1f, p_Time=%.4f\n', ...
                p_global_modelbased, gamm_mdl_time.ModelCriterion.AIC, p_time);
            if isnan(resid_ac_r)
                fprintf('       AVEC temps, autocorr_resid(lag1)      : N insuffisant (%d)\n', resid_ac_N);
            else
                fprintf('       AVEC temps, autocorr_resid(lag1)      : r=%.3f p=%.4f (N=%d)\n', ...
                    resid_ac_r, resid_ac_p, resid_ac_N);
            end
            fprintf('       AVEC temps, ROBUSTE cluster F(%d,%d)   : p_global=%.4f (F=%.3f)\n', ...
                df1_robust, df2_robust, p_global, F_robust);

            drawGAMMFigure(dX_ok, dRPE_ok, subj_ok, x_grid, y_grid, ...
                r.vname, r.source, gn_lbl, N_subj_ok, p_global, is_interpolated_rpe, PathSave);

        catch ME
            fprintf('  [%s] GAMM échoué : %s\n', gn_lbl, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
%  CORRECTION FDR (Benjamini-Hochberg), séparément par groupe
%  -----------------------------------------------------------------------
% Une famille de tests par groupe (G1, G2, Normal, SD, LD), sur les variables
% effectivement ajustées dans cette exécution. Le champ p_global contient le
% test de Wald global (méthode GAMM_lag1, G1/G2) ou le p de Spearman (méthode
% Spearman_drift, Normal/SD/LD) — traités comme une seule famille de décision
% par groupe, puisqu'ils répondent à la même question ("ce candidat est-il
% associé au RPE dans ce groupe ?"), même si la méthode diffère.
% NB : implémentation BH locale (bh_fdr ci-dessous) — à remplacer par le
% bh_fdr.m existant si sa signature diffère.
fprintf('\n=== Correction FDR (Benjamini-Hochberg, par groupe) ===\n');
all_groups_res = {GAMM_Results.group};
for ig = 1:N_groupes
    gn_lbl = groupes{ig};
    idx_g = find(strcmp(all_groups_res, gn_lbl));
    if isempty(idx_g), continue; end

    p_vals_g = [GAMM_Results(idx_g).p_global];
    q_vals_g = bh_fdr(p_vals_g);
    for k = 1:length(idx_g)
        GAMM_Results(idx_g(k)).q_FDR = q_vals_g(k);
    end

    n_sig_raw = sum(p_vals_g < 0.05);
    n_sig_fdr = sum(q_vals_g < 0.05);
    fprintf('  [%s] %d tests | p<0.05 (brut) : %d | q<0.05 (FDR) : %d\n', ...
        gn_lbl, length(idx_g), n_sig_raw, n_sig_fdr);
    for k = 1:length(idx_g)
        sig_mark = '';
        if q_vals_g(k) < 0.05, sig_mark = '  *'; end
        fprintf('       %-45s p=%.4f -> q=%.4f%s\n', ...
            GAMM_Results(idx_g(k)).vname, p_vals_g(k), q_vals_g(k), sig_mark);
    end
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'GAMM_Bins_Results.mat'), 'GAMM_Results', ...
    'N_SPLINE_KNOTS', 'MIN_N_SUBJ', 'MIN_N_OBS');
fprintf('\nSauvegardé : GAMM_Bins_Results.mat (%d ajustements réussis)\n', length(GAMM_Results));
fprintf('\n=== PIPELINE TERMINÉ ===\n');

%% =========================================================================
%% FONCTIONS LOCALES
%% =========================================================================

% -------------------------------------------------------------------------
function H = findSplineContrastH(mdl)
% Matrice de contraste (une ligne par coefficient spline B#, 0 ailleurs) :
% partagée par le test de Wald model-based (waldTestSplineCoefs) et le test
% robuste en cluster (clusterRobustWaldTest) ci-dessous.
    coef_names = mdl.CoefficientNames;
    spline_idx = find(startsWith(coef_names, 'B') & ...
        ~cellfun(@isempty, regexp(coef_names, '^B\d+$', 'once')));
    n_coef = length(coef_names);
    H = zeros(length(spline_idx), n_coef);
    for hk = 1:length(spline_idx)
        H(hk, spline_idx(hk)) = 1;
    end
end

% -------------------------------------------------------------------------
function p = waldTestSplineCoefs(mdl)
% Test de Wald global model-based (H0 : tous les coefficients spline B#
% sont nuls), covariance NON robuste (celle de coefTest/CoefficientCovariance,
% qui suppose le modèle de covariance du GAMM correctement spécifié).
    H = findSplineContrastH(mdl);
    p = NaN;
    if isempty(H), return; end
    try
        [p, ~] = coefTest(mdl, H);
    catch ME_ct
        fprintf('  (coefTest a échoué : %s — p=NaN)\n', ME_ct.message);
    end
end

% -------------------------------------------------------------------------
function [F_stat, pval, df1, df2] = clusterRobustWaldTest(mdl, H, subj_ok, y)
% Test de Wald sur H*beta=0 avec covariance "sandwich" robuste en cluster
% (un cluster = un sujet), et test F(q, G-1) (Cameron & Miller, 2015) plutôt
% qu'un chi² asymptotique — recommandé quand le nombre de clusters G est
% modéré (ici G = N_sujets, 10 à 25), ce qui est le cas pour G1/G2.
%
% Sandwich : V_robuste = Bread * (somme_g X_g' V_g^-1 e_g e_g' V_g^-1 X_g) * Bread
%   - Bread = mdl.CoefficientCovariance = (X' V^-1 X)^-1 du modèle ajusté
%   - e_g   = résidus MARGINAUX (y - X*beta, SANS les effets aléatoires) du
%             cluster g (à ne pas confondre avec residuals(mdl), qui inclut
%             les BLUP et sert au diagnostic d'autocorrélation, pas ici)
%   - V_g   = sigma_b^2 * J_ng + sigma_e^2 * I_ng (structure compound-symmetry
%             d'un intercept aléatoire (1|Subject), mêmes sigma_b^2/sigma_e^2
%             que ceux estimés par le modèle (covarianceParameters)
%
% NB : pas de correction CR1 multiplicative supplémentaire ici (G/(G-1) etc.)
% — la référence F(q,G-1) porte déjà la correction petit-échantillon
% recommandée par Cameron & Miller pour G modéré.
    X     = designMatrix(mdl, 'Fixed');
    beta  = fixedEffects(mdl);
    Bread = mdl.CoefficientCovariance;
    e     = y - fitted(mdl, 'Conditional', false);   % résidus marginaux

    [psi, mse] = covarianceParameters(mdl);
    sigma_b2 = psi{1}(1,1);
    sigma_e2 = mse;

    uSubj = unique(subj_ok);
    G = length(uSubj);
    n_par = size(X, 2);
    meat = zeros(n_par, n_par);
    for ig = 1:G
        idx_g = (subj_ok == uSubj(ig));
        n_g   = sum(idx_g);
        X_g   = X(idx_g, :);
        e_g   = e(idx_g);
        V_g   = sigma_b2 * ones(n_g, n_g) + sigma_e2 * eye(n_g);
        score_g = X_g' * (V_g \ e_g);
        meat = meat + score_g * score_g';
    end

    V_robust = Bread * meat * Bread;

    beta_H = H * beta;
    M = H * V_robust * H';
    q = size(H, 1);
    W = beta_H' * (M \ beta_H);
    F_stat = W / q;
    df1 = q;
    df2 = G - 1;
    pval = 1 - fcdf(F_stat, df1, df2);
end

% -------------------------------------------------------------------------
function [r, pval, N_pairs] = checkResidAutocorr(resid, subj_ok, bin_ok)
% Autocorrélation lag-1 des résidus intra-sujet (diagnostic, pas un test
% formel) : on ne pairer que des bins réellement consécutifs (bin_ok(k+1) ==
% bin_ok(k)+1) au sein d'un même sujet, pour ne pas créer de fausse paire à
% travers un bin manquant (cf. filtrage NaN amont qui peut casser la suite).
    uSubj = unique(subj_ok);
    r0 = [];
    r1 = [];
    for is = 1:length(uSubj)
        idx_s = find(subj_ok == uSubj(is));
        [bin_sorted, order] = sort(bin_ok(idx_s));
        idx_s_sorted = idx_s(order);
        resid_s = resid(idx_s_sorted);
        for k = 1:length(bin_sorted)-1
            if bin_sorted(k+1) == bin_sorted(k) + 1
                r0(end+1,1) = resid_s(k);   %#ok<AGROW>
                r1(end+1,1) = resid_s(k+1); %#ok<AGROW>
            end
        end
    end
    N_pairs = length(r0);
    if N_pairs < 5
        r = NaN; pval = NaN;
        return;
    end
    [r, pval] = corr(r0, r1, 'Type', 'Pearson');
end

% -------------------------------------------------------------------------
function q = bh_fdr(p)
% Correction FDR de Benjamini-Hochberg (implémentation locale standard).
% q(i) = q-value (p-value ajustée) associée à p(i), même ordre que p en entrée.
    p = p(:);
    n = length(p);
    [p_sorted, idx_sort] = sort(p);
    q_sorted = p_sorted .* n ./ (1:n)';
    for i = n-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end
    q_sorted = min(q_sorted, 1);
    q = NaN(n, 1);
    q(idx_sort) = q_sorted;
end

% -------------------------------------------------------------------------
function [B, knots] = buildSplineBasis(x, n_knots)
% Construit une base de B-splines cubiques avec n_knots nœuds internes
% placés aux quantiles de x (répartition adaptative aux données plutôt
% qu'une grille régulière, plus robuste sur petits échantillons).
%
% Implémentation locale (algorithme de Cox-de Boor), sans dépendance à la
% Curve Fitting Toolbox.
    x = x(:);
    q = linspace(0, 1, n_knots+2);
    q = q(2:end-1);  % on exclut 0 et 1 (les bords sont gérés par les degrés du spline)
    inner_knots = quantile(x, q);
    inner_knots = unique(inner_knots);  % évite les nœuds dupliqués sur données discrètes

    x_min = min(x); x_max = max(x);
    degree = 3;  % cubique
    % Nœuds étendus : répétition des bords (degree+1 fois) + nœuds internes
    knots = [repmat(x_min, 1, degree+1), inner_knots, repmat(x_max, 1, degree+1)];

    % Clip x strictement à l'intérieur de [x_min, x_max] : la récursion de
    % Cox-de Boor utilise des intervalles semi-ouverts [knot_j, knot_j+1),
    % donc x_max lui-même tomberait hors de tout intervalle sans ce clip
    % (approche standard, plus robuste qu'un cas particulier sur le dernier
    % nœud, qui échoue avec des nœuds répétés en bord — bug détecté en test).
    eps_clip = (x_max - x_min) * 1e-10;
    if eps_clip == 0, eps_clip = 1e-10; end
    x_eval = min(max(x, x_min), x_max - eps_clip);

    n_basis = length(knots) - degree - 1;
    B = zeros(length(x), n_basis);
    for j = 1:n_basis
        B(:,j) = bsplineBasisFunc(x_eval, knots, j, degree);
    end
end

% -------------------------------------------------------------------------
function B = evalSplineBasis(x, knots)
% Évalue la même base de spline (mêmes nœuds) sur une nouvelle grille x —
% utilisé pour la prédiction sur x_grid avec les noeuds appris sur x_ok.
    degree = 3;
    n_basis = length(knots) - degree - 1;
    x = x(:);

    x_min = knots(1); x_max = knots(end);
    eps_clip = (x_max - x_min) * 1e-10;
    if eps_clip == 0, eps_clip = 1e-10; end
    x_eval = min(max(x, x_min), x_max - eps_clip);

    B = zeros(length(x), n_basis);
    for j = 1:n_basis
        B(:,j) = bsplineBasisFunc(x_eval, knots, j, degree);
    end
end

% -------------------------------------------------------------------------
function y = bsplineBasisFunc(x, knots, j, degree)
% Fonction de base B-spline j de degré 'degree', récurrence de Cox-de Boor.
% j est 1-based (MATLAB) ; converti en 0-based pour la récurrence classique.
    j0 = j - 1;
    y = bsplineRecursive(x, knots, j0, degree);
end

function y = bsplineRecursive(x, knots, j, p)
    x = x(:);
    if p == 0
        y = double(x >= knots(j+1) & x < knots(j+2));
        return;
    end

    denom1 = knots(j+p+1) - knots(j+1);
    denom2 = knots(j+p+2) - knots(j+2);

    if denom1 > 0
        term1 = (x - knots(j+1)) / denom1 .* bsplineRecursive(x, knots, j, p-1);
    else
        term1 = zeros(size(x));
    end

    if denom2 > 0
        term2 = (knots(j+p+2) - x) / denom2 .* bsplineRecursive(x, knots, j+1, p-1);
    else
        term2 = zeros(size(x));
    end

    y = term1 + term2;
end

% -------------------------------------------------------------------------
function [rho, pval, N_subj_ok, dX_drift_ok, dRPE_drift_ok] = computeDriftCorrelation(X_bins, RPE_bins, MIN_N_SUBJ)
% Corrélation de Spearman à un point par sujet : dérive totale de X sur la
% session (dernier bin valide - premier bin valide) vs dérive totale de RPE
% (même logique). Remplace le GAMM lag-1 pour les groupes où le RPE est
% interpolé linéairement (cf. en-tête du script) : dans ce cas, le delta
% lag-1 de RPE est une constante par sujet et le GAMM lag-1 + (1|Sujet) y
% est structurellement dégénéré.
    N_subj_loc = size(X_bins, 1);
    dX_drift   = NaN(N_subj_loc, 1);
    dRPE_drift = NaN(N_subj_loc, 1);
    for is = 1:N_subj_loc
        x_start   = firstValid(X_bins(is, :));
        x_end     = lastValid(X_bins(is, :));
        rpe_start = firstValid(RPE_bins(is, :));
        rpe_end   = lastValid(RPE_bins(is, :));
        dX_drift(is)   = x_end - x_start;
        dRPE_drift(is) = rpe_end - rpe_start;
    end

    ok = ~isnan(dX_drift) & ~isnan(dRPE_drift);
    dX_drift_ok   = dX_drift(ok);
    dRPE_drift_ok = dRPE_drift(ok);
    N_subj_ok     = length(dX_drift_ok);

    if N_subj_ok < MIN_N_SUBJ || std(dX_drift_ok) < 1e-12
        rho = NaN; pval = NaN;
        return;
    end

    [rho, pval] = corr(dX_drift_ok, dRPE_drift_ok, 'Type', 'Spearman');
end

% -------------------------------------------------------------------------
function v = firstValid(row)
% Première valeur non-NaN d'un vecteur ligne (NaN si aucune).
    idx = find(~isnan(row), 1, 'first');
    if isempty(idx), v = NaN; else, v = row(idx); end
end

% -------------------------------------------------------------------------
function v = lastValid(row)
% Dernière valeur non-NaN d'un vecteur ligne (NaN si aucune).
    idx = find(~isnan(row), 1, 'last');
    if isempty(idx), v = NaN; else, v = row(idx); end
end

% -------------------------------------------------------------------------
function drawDriftFigure(dX_drift, dRPE_drift, rho, pval, vname, source, gname, n_subj, PathSave)
% Figure pour la branche corrélation drift (1 point par sujet) : nuage de
% points + droite de régression linéaire en repère visuel uniquement — le
% test statistique reporté (titre) est la corrélation de Spearman, robuste
% aux valeurs extrêmes contrairement à cette droite.
    fig = figure('Visible','off','Position',[50 50 750 550],'Color','white');
    hold on;

    scatter(dX_drift, dRPE_drift, 50, [0.2 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.7);

    p_fit  = polyfit(dX_drift, dRPE_drift, 1);
    x_line = linspace(min(dX_drift), max(dX_drift), 100);
    y_line = polyval(p_fit, x_line);
    plot(x_line, y_line, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2);

    xlabel('Dérive totale de X (dernier bin valide - premier bin valide)', 'FontSize', 11);
    ylabel('\Delta RPE total (session)', 'FontSize', 11);

    title({sprintf('%s — %s (corrélation drift, 1 point/sujet)', strrep(vname,'_',' '), gname), ...
           sprintf('N_{sujets}=%d | \\rho_{Spearman}=%.3f | p=%.4f | RPE INTERPOLÉ (approximation)', n_subj, rho, pval)}, ...
        'FontSize', 10, 'Interpreter', 'tex', 'Color', [0.75 0.15 0.15]);
    grid on; box on; hold off;

    fname = sprintf('Fig_Drift_%s_%s_%s.png', source, vname, gname);
    fname = regexprep(fname, '[^a-zA-Z0-9_.]', '_');
    exportgraphics(fig, fullfile(PathSave, fname), 'Resolution', 150);
    close(fig);
    fprintf('  Figure drift sauvegardée : %s\n', fname);
end

% -------------------------------------------------------------------------
function drawGAMMFigure(x_obs, y_obs, subj_obs, x_grid, y_pred, vname, source, gname, n_subj, p_global, is_interpolated_rpe, PathSave)
    fig = figure('Visible','off','Position',[50 50 750 550],'Color','white');
    hold on;

    % Couleur par sujet pour visualiser le clustering intra-sujet
    uSubj = unique(subj_obs);
    cmap = lines(length(uSubj));
    for is = 1:length(uSubj)
        idx_s = subj_obs == uSubj(is);
        scatter(x_obs(idx_s), y_obs(idx_s), 30, cmap(is,:), 'filled', 'MarkerFaceAlpha', 0.5);
    end

    plot(x_grid, y_pred, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.5);
    xlabel('\Delta X (lag-1, bin t - bin t-1)', 'FontSize', 11);
    ylabel('\Delta RPE (lag-1)', 'FontSize', 11);

    if isnan(p_global)
        sig_str = 'p_{global}=NaN';
    else
        sig_str = sprintf('p_{global}=%.4f', p_global);
    end

    title_color = [0 0 0];
    subtitle2 = sprintf('N_{sujets}=%d | %s | points colorés par sujet', n_subj, sig_str);
    if is_interpolated_rpe
        subtitle2 = [subtitle2 ' — RPE INTERPOLÉ (approximation)'];
        title_color = [0.75 0.15 0.15];
    end

    title({sprintf('%s — %s (GAMM, effet aléatoire sujet)', strrep(vname,'_',' '), gname), ...
           subtitle2}, ...
        'FontSize', 10, 'Interpreter', 'none', 'Color', title_color);
    grid on; box on; hold off;

    fname = sprintf('Fig_GAMM_%s_%s_%s.png', source, vname, gname);
    fname = regexprep(fname, '[^a-zA-Z0-9_.]', '_');
    exportgraphics(fig, fullfile(PathSave, fname), 'Resolution', 150);
    close(fig);
    fprintf('  Figure GAMM sauvegardée : %s\n', fname);
end