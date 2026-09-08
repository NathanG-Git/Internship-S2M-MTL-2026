%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  PLSR + VIP ITERATIF — Dataset Robin                              %%%%
%%%%  4 modèles : M1=SH/Std  M2=SH/Ergo  M3=LH/Std  M4=Pool/Std      %%%%
%%%%  1280 variables (doublons IMU/Goubault et DeltMed supprimés)        %%%%
%%%%  VIP + signe B calculés sur jeu TRAIN de chaque fold             %%%%
%%%%  B moyenné sur folds où la feature est sélectionnée uniquement   %%%%
%%%%  Grouped k-fold CV (par sujet) × 50 répétitions                  %%%%
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
path_plsr  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\PLSR\';
F_spectral = load(fullfile(path_plsr, 'SpectralFeatures_Ergo3.mat'));
F_emg      = load(fullfile(path_data, 'EMG_Features_NewDataset_25vars.mat'));
fprintf('  OK\n\n');

%% -----------------------------------------------------------------------
%  SUJETS
%  -----------------------------------------------------------------------
Subject_Both = {'P02','P03','P04','P05','P08','P10','P12',...
                'P15','P17','P18','P20','P21','P24'};
Subject_LD   = {'P01','P06','P09','P11','P13','P14',...
                'P16','P19','P22','P23'};

%% -----------------------------------------------------------------------
%  RPE
%  -----------------------------------------------------------------------
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
    s = RPE_data{i,1};
    RPE.(s).Norm.debut = RPE_data{i,2};
    RPE.(s).Norm.fin   = RPE_data{i,3};
    if ~isnan(RPE_data{i,4})
        RPE.(s).CS60.debut = RPE_data{i,4};
        RPE.(s).CS60.fin   = RPE_data{i,5};
    end
end

%% -----------------------------------------------------------------------
%  FORME RPE NORMALISÉE G1+G2
%  -----------------------------------------------------------------------
RPE_shape = [0, 0.2757, 0.4889, 0.6171, 0.7147, 0.7717, 0.8557, 0.9065, 0.9555, 1.0];

%% -----------------------------------------------------------------------
%  PARAMÈTRES PLSR
%  -----------------------------------------------------------------------
N_bins   = 10;
N_rep    = 50;    % répétitions Passe 2
N_rep_p1 = 50;   % répétitions Passe 1 (moyennage AbsErr par niveau, comme Goubault)
K_fold   = 5;
VIP_init = 1.0;
VIP_step = 0.1;

%% -----------------------------------------------------------------------
%  LISTE DES 1280 VARIABLES
%  1260 features spectrales (SpectralFeatures_Ergo) + 20 EMG
%  Variables cinématiques (IMU/Goubault/SegMod/SegAxis) supprimées :
%  remplacées par le pipeline Goubault dans SpectralFeatures_Ergo
%  -----------------------------------------------------------------------

% ── Features spectrales (1260) ───────────────────────────────────────────
seg_names_sp = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
sig_types_sp = {'Accel','AngVel','Jerk'};
comp_sp      = {'X','Y','Z','Mod'};
feat_sp      = {'Mean','Std','Median','Max','p10','p25','p75','p90',...
                'MedianFreq','SpectralEntropy','Power_below4Hz',...
                'Power_above4Hz','TotalPower','PeakPower','PeakPowerFreq'};

feat_spectral = {};
for iSig = 1:numel(sig_types_sp)
    for iCo = 1:numel(comp_sp)
        for iFt = 1:numel(feat_sp)
            for iSeg = 1:numel(seg_names_sp)
                feat_spectral{end+1} = sprintf('%s_%s_%s_%s', ...
                    sig_types_sp{iSig}, comp_sp{iCo}, ...
                    feat_sp{iFt}, seg_names_sp{iSeg});
            end
        end
    end
end
src_spectral = repmat({'Spectral'}, 1, numel(feat_spectral));

% ── Features EMG (20) — DeltMed supprimé (doublon confirmé de DeltAnt) ──
chans     = {'Biceps','Triceps','DeltAnt','SupTrap'};
emg_types = {'Activity','Mobility','SampleEntropy','TFR_MedianFreq','TFR_SpectralEntropy'};
feat_emg  = {};
for et = emg_types
    for ch = chans
        feat_emg{end+1} = sprintf('%s_%s', et{1}, ch{1});
    end
end
src_emg = repmat({'EMG'}, 1, numel(feat_emg));

% ── Concaténation ────────────────────────────────────────────────────────
all_feat   = [feat_spectral, feat_emg];
all_source = [src_spectral,  src_emg];
N_feat     = numel(all_feat);

fprintf('DEBUG : spectral=%d  emg=%d  total=%d\n', numel(feat_spectral), numel(feat_emg), N_feat);

% ── Renommage algorithmique des features spectrales ─────────────────────
% Règle : <Signal>_<Comp>_<Feature>_<Segment> → <Segment>_<Feature>_<Signal><Comp>
% Mod → Magnitude | X/Y/Z conservés en majuscule
% EMG : noms inchangés
rename_map = containers.Map('KeyType','char','ValueType','char');
for iF = 1:numel(feat_spectral)
    old_name = feat_spectral{iF};
    % Parser : Signal_Comp_Feature_Segment
    parts = strsplit(old_name, '_');
    % parts{1}=Signal, parts{2}=Comp, parts{3..end-1}=Feature, parts{end}=Segment
    sig_p  = parts{1};                          % Accel / AngVel / Jerk
    comp_p = parts{2};                          % X Y Z Mod
    seg_p  = parts{end};                        % L5 T8 Head ...
    feat_p = strjoin(parts(3:end-1), '_');      % Mean / MedianFreq / Power_below4Hz ...

    % Substitutions composante
    comp_p = strrep(comp_p, 'Mod', 'Magnitude');

    % Nouveau nom : Segment_Feature_SignalComp
    new_name = sprintf('%s_%s_%s%s', seg_p, feat_p, sig_p, comp_p);
    rename_map(old_name) = new_name;
end

% Appliquer le renommage sur all_feat
for iF = 1:numel(all_feat)
    if rename_map.isKey(all_feat{iF})
        all_feat{iF} = rename_map(all_feat{iF});
    end
end

% Reverse map : nouveau_nom → ancien_nom (pour l'extraction depuis Results_Spectral)
reverse_map = containers.Map('KeyType','char','ValueType','char');
old_keys = keys(rename_map);
for iK = 1:numel(old_keys)
    reverse_map(rename_map(old_keys{iK})) = old_keys{iK};
end

fprintf('Renommage spectral appliqué. Exemples :\n');
ex_old = {'Accel_Y_MedianFreq_Hand','AngVel_Mod_SpectralEntropy_Forearm',...
          'Jerk_X_PeakPower_Arm','Accel_Mod_Mean_T8'};
for i = 1:numel(ex_old)
    if rename_map.isKey(ex_old{i})
        fprintf('  %s → %s\n', ex_old{i}, rename_map(ex_old{i}));
    end
end
fprintf('\n');

fprintf('Variables : %d  (Spectral=%d  EMG=%d)\n\n', ...
    N_feat, numel(feat_spectral), numel(feat_emg));

%% -----------------------------------------------------------------------
%  MODÈLES
%  -----------------------------------------------------------------------
models = {
    'M1_SH_Std',  Subject_Both,               'Norm';
    'M2_SH_Ergo', Subject_Both,               'CS60';
    'M3_LH_Std',  Subject_LD,                 'Norm';
    'M4_Pool_Std',[Subject_Both, Subject_LD],  'Norm';
};

Results_PLSR = struct();

%% -----------------------------------------------------------------------

%%  BOUCLE PRINCIPALE PAR MODÈLE — Méthode Goubault
%  BOUCLE PRINCIPALE PAR MODÈLE — Méthode Goubault
%  Passe 1 : modèle complet → VIP → niveaux de réduction successifs
%  Pour chaque niveau : erreur CV → tableau comparatif → meilleur modèle
%  Passe 2 : 50 répétitions sur le meilleur niveau → VIP stables + B signés
%  -----------------------------------------------------------------------
for iM = 1:size(models,1)
    model_name = models{iM,1};
    subjects   = models{iM,2};
    cond       = models{iM,3};
    N_subj     = numel(subjects);

    fprintf('=== %s (N=%d, cond=%s) ===\n', model_name, N_subj, cond);

    %% --- Construction X, Y, G ------------------------------------------
    X_raw = NaN(N_subj * N_bins, N_feat);
    Y_raw = NaN(N_subj * N_bins, 1);
    G_raw = NaN(N_subj * N_bins, 1);

    for iS = 1:N_subj
        subj = subjects{iS};
        rows = (iS-1)*N_bins + (1:N_bins);
        if ~isfield(RPE,subj) || ~isfield(RPE.(subj),cond)
            fprintf('  [SKIP] %s — pas de RPE\n', subj);
            continue;
        end
        rpe0 = RPE.(subj).(cond).debut;
        rpe1 = RPE.(subj).(cond).fin;
        Y_raw(rows) = rpe0 + (rpe1-rpe0) .* RPE_shape(:);
        G_raw(rows) = iS;

        for iF = 1:N_feat
            fname = all_feat{iF};
            src   = all_source{iF};
            val   = NaN(N_bins,1);
            try
                switch src
                    case 'Spectral'
                        % fname est le nom renommé → traduire vers l'ancien nom
                        % pour accéder à Results_Spectral qui garde les anciens noms
                        if reverse_map.isKey(fname)
                            fname_orig = reverse_map(fname);
                        else
                            fname_orig = fname;
                        end
                        val = F_spectral.Results_Spectral.(subj).(cond).(fname_orig)(:);
                    case 'EMG'
                        val = F_emg.Results_EMG.(subj).(cond).(fname)(:);
                end
            catch, end
            X_raw(rows,iF) = val;
        end
    end

    %% --- Nettoyage -----------------------------------------------------
    valid     = ~isnan(Y_raw) & ~isnan(G_raw);
    X         = X_raw(valid,:);
    Y         = Y_raw(valid);
    G         = G_raw(valid);

    keep      = mean(isnan(X),1) < 0.2;
    X         = X(:,keep);
    feat_used = all_feat(keep);
    src_used  = all_source(keep);
    N_fu      = sum(keep);
    fprintf('  Features utilisées : %d / %d\n', N_fu, N_feat);

    for iF = 1:N_fu
        col = X(:,iF);
        col(isnan(col)) = nanmedian(col);
        X(:,iF) = col;
    end

    mu_X = nanmean(X);   sd_X = nanstd(X);   sd_X(sd_X==0) = 1;
    mu_Y = nanmean(Y);   sd_Y = nanstd(Y);
    X_z  = (X - mu_X) ./ sd_X;
    Y_z  = (Y - mu_Y) ./ sd_Y;

    N_obs    = size(X_z,1);
    subj_ids = unique(G);
    N_s      = numel(subj_ids);

    %% -------------------------------------------------------------------
    %  PASSE 1 — Modèle complet + tableau de réduction (méthode Goubault)
    %  N_latent = 20 fixe pour cohérence avec Goubault
    %  Pour chaque seuil VIP : évaluer N_input, N_latent_opt, %RPE, AbsErr
    %  -------------------------------------------------------------------
    N_latent_goub = min(20, N_fu-1);   % N_latent fixe comme Goubault

    fprintf('  Passe 1 : calcul VIP modèle complet...\n');

    % VIP sur modèle complet (une seule passe, toutes observations)
    [XL_full, YL_full, XS_full, ~, ~, PCTVAR_full, ~, stats_full] = ...
        plsregress(X_z, Y_z, N_latent_goub);
    W0_full  = stats_full.W ./ sqrt(sum(stats_full.W.^2,1));
    SC_full  = sum(XS_full.^2,1) .* (YL_full.^2);
    VIP_full = sqrt(N_fu * ((W0_full.^2) * SC_full') / sum(SC_full));

    % Variance Y expliquée sur train (modèle complet)
    pct_train_full = sum(PCTVAR_full(2,:)) * 100;

    % Niveaux de seuil à tester
    thresh_levels = VIP_init : VIP_step : max(VIP_full) + VIP_step;
    N_levels      = numel(thresh_levels);

    % Tableau résultats par niveau — accumulé sur N_rep_p1 répétitions
    tbl_thresh   = thresh_levels;
    tbl_ninput   = zeros(1, N_levels);
    tbl_nlatent  = zeros(1, N_levels);
    tbl_pct_train= zeros(1, N_levels);
    tbl_mae_all  = nan(N_rep_p1, N_levels);   % AbsErr par run × niveau
    tbl_r2       = nan(1, N_levels);

    % Niveau 0 = modèle complet (N_input fixe)
    tbl_ninput(1)    = N_fu;
    tbl_nlatent(1)   = N_latent_goub;
    tbl_pct_train(1) = pct_train_full;

    fprintf('  Passe 1 : %d répétitions par niveau...\n', N_rep_p1);

    for iRep_p1 = 1:N_rep_p1

        % Partition aléatoire indépendante à chaque répétition
        subj_order_p1  = subj_ids(randperm(N_s));
        fold_assign_p1 = mod((1:N_s)'-1, K_fold) + 1;

        % ── Niveau 0 : modèle complet ────────────────────────────────
        Y_pred_full = NaN(N_obs,1);
        for iFold = 1:K_fold
            te_s = subj_order_p1(fold_assign_p1 == iFold);
            te   = ismember(G, te_s);   tr = ~te;
            if sum(tr)<10 || sum(te)<1, continue; end
            n_c = min(N_latent_goub, numel(unique(G(tr)))-1);
            if n_c < 1, continue; end
            [~,~,~,~,B_c] = plsregress(X_z(tr,:), Y_z(tr), n_c);
            Y_pred_full(te) = [ones(sum(te),1), X_z(te,:)] * B_c;
        end
        Yp = Y_pred_full * sd_Y + mu_Y;
        ok = ~isnan(Yp);
        if sum(ok)>2
            tbl_mae_all(iRep_p1,1) = mean(abs(Yp(ok)-Y(ok)));
        end

        % ── Niveaux de réduction ─────────────────────────────────────
        for iLvl = 2:N_levels
            thresh = thresh_levels(iLvl);
            sel    = VIP_full >= thresh;
            n_sel  = sum(sel);
            if n_sel < 2
                break;
            end

            X_sel = X_z(:, sel);

            % N_latent optimal (cherché une fois par répétition)
            n_lat_max    = min(N_latent_goub, min(n_sel-1, N_s-2));
            best_mae_lat = Inf;
            best_nlat    = 1;

            for n_lat = 1:n_lat_max
                Y_pred_lat = NaN(N_obs,1);
                for iFold = 1:K_fold
                    te_s = subj_order_p1(fold_assign_p1 == iFold);
                    te   = ismember(G, te_s);   tr = ~te;
                    if sum(tr)<10 || sum(te)<1, continue; end
                    n_c = min(n_lat, numel(unique(G(tr)))-1);
                    if n_c < 1, continue; end
                    [~,~,~,~,B_c] = plsregress(X_sel(tr,:), Y_z(tr), n_c);
                    Y_pred_lat(te) = [ones(sum(te),1), X_sel(te,:)] * B_c;
                end
                Yp_lat = Y_pred_lat * sd_Y + mu_Y;
                ok_lat = ~isnan(Yp_lat);
                if sum(ok_lat) > 2
                    mae_lat = mean(abs(Yp_lat(ok_lat) - Y(ok_lat)));
                    if mae_lat < best_mae_lat
                        best_mae_lat = mae_lat;
                        best_nlat    = n_lat;
                    end
                end
            end

            % AbsErr finale avec N_latent optimal
            Y_pred_lv = NaN(N_obs,1);
            for iFold = 1:K_fold
                te_s = subj_order_p1(fold_assign_p1 == iFold);
                te   = ismember(G, te_s);   tr = ~te;
                if sum(tr)<10 || sum(te)<1, continue; end
                n_c = min(best_nlat, numel(unique(G(tr)))-1);
                if n_c < 1, continue; end
                [~,~,~,~,B_c] = plsregress(X_sel(tr,:), Y_z(tr), n_c);
                Y_pred_lv(te) = [ones(sum(te),1), X_sel(te,:)] * B_c;
            end
            Yp_lv = Y_pred_lv * sd_Y + mu_Y;
            ok_lv = ~isnan(Yp_lv);
            if sum(ok_lv) > 2
                tbl_mae_all(iRep_p1,iLvl) = mean(abs(Yp_lv(ok_lv) - Y(ok_lv)));
            end

            % Stocker N_input, N_latent, %RPE train (identiques quelle que soit la répétition)
            if iRep_p1 == 1
                [~,~,~,~,~,PV_sel] = plsregress(X_sel, Y_z, best_nlat);
                tbl_ninput(iLvl)    = n_sel;
                tbl_nlatent(iLvl)   = best_nlat;
                tbl_pct_train(iLvl) = sum(PV_sel(2,:)) * 100;
            end
        end

        if mod(iRep_p1,10)==0
            fprintf('    Passe 1 : %d/%d répétitions\n', iRep_p1, N_rep_p1);
        end
    end   % fin boucle N_rep_p1

    % Moyenner les AbsErr sur les N_rep_p1 répétitions
    tbl_mae = nanmean(tbl_mae_all, 1);
    tbl_r2  = nan(1, N_levels);   % R² non moyenné (non utilisé pour le choix)

    %% --- Sélection du meilleur niveau — méthode Goubault -------------------
    % 1. ANOVA à un facteur sur les AbsErr des N_rep_p1 runs × N_niveaux
    % 2. Post-hoc Tukey : niveaux non significativement différents du minimum
    % 3. Parmi ces niveaux : choisir celui avec le %RPE train le plus élevé

    % Niveaux valides (au moins une valeur non-NaN)
    valid_lvls = find(~all(isnan(tbl_mae_all), 1));
    N_valid    = numel(valid_lvls);

    if N_valid < 2
        % Pas assez de niveaux pour ANOVA → fallback sur min AbsErr
        [~, best_idx_v] = min(tbl_mae(valid_lvls));
        best_lvl        = valid_lvls(best_idx_v);
        fprintf('\n  → [Fallback] Meilleur modèle : seuil=%.1f  N=%d  Nlat=%d  AbsErr=%.2f\n',...
            tbl_thresh(best_lvl), tbl_ninput(best_lvl), tbl_nlatent(best_lvl), tbl_mae(best_lvl));
    else
        % Construire la matrice ANOVA : N_rep_p1 × N_valid_levels
        X_anova = tbl_mae_all(:, valid_lvls);   % (N_rep_p1 × N_valid)

        % ANOVA à un facteur (silencieuse)
        grp_labels = arrayfun(@(k) sprintf('L%d',k), 1:N_valid, 'UniformOutput',false);
        [~, anova_tbl, anova_stats] = anova1(X_anova, grp_labels, 'off');

        % Post-hoc Tukey
        [mc_results, ~] = multcompare(anova_stats, 'CType','tukey-kramer', 'Display','off');
        % mc_results : [grp_i, grp_j, lower, diff, upper, p_value]

        % Niveau avec AbsErr minimale (dans les niveaux valides)
        [~, best_idx_v] = min(tbl_mae(valid_lvls));
        best_idx_global = best_idx_v;   % index dans valid_lvls

        % Trouver les niveaux non significativement différents du meilleur (p > 0.05)
        equiv_idx = best_idx_global;   % commence avec le meilleur seul
        for iRow = 1:size(mc_results,1)
            gi = mc_results(iRow,1);
            gj = mc_results(iRow,2);
            pv = mc_results(iRow,6);
            % Si une comparaison implique best_idx_global et p > 0.05 → ajouter l'autre
            if gi == best_idx_global && pv > 0.05
                equiv_idx(end+1) = gj;
            elseif gj == best_idx_global && pv > 0.05
                equiv_idx(end+1) = gi;
            end
        end
        equiv_idx = unique(equiv_idx);

        % Parmi les niveaux équivalents → choisir celui avec %RPE train le plus élevé
        pct_equiv = tbl_pct_train(valid_lvls(equiv_idx));
        [~, best_pct_idx] = max(pct_equiv);
        best_idx_v   = equiv_idx(best_pct_idx);
        best_lvl     = valid_lvls(best_idx_v);

        fprintf('\n  ANOVA p=%.4f | %d niveaux équivalents au minimum AbsErr (Tukey p>0.05)\n',...
            anova_tbl{2,6}, numel(equiv_idx));
        fprintf('  → Meilleur modèle (min AbsErr + max %%RPE parmi équivalents) :\n');
        fprintf('     seuil=%.1f  N=%d  Nlat=%d  AbsErr=%.2f  %%RPE=%.1f%% (moy. %d runs)\n',...
            tbl_thresh(best_lvl), tbl_ninput(best_lvl), tbl_nlatent(best_lvl), ...
            tbl_mae(best_lvl), tbl_pct_train(best_lvl), N_rep_p1);
    end

    best_thresh = tbl_thresh(best_lvl);
    best_sel    = VIP_full >= best_thresh;

    % Stocker tableau réduction
    Results_PLSR.(model_name).reduction.thresh    = tbl_thresh;
    Results_PLSR.(model_name).reduction.n_input   = tbl_ninput;
    Results_PLSR.(model_name).reduction.n_latent  = tbl_nlatent;
    Results_PLSR.(model_name).reduction.pct_train = tbl_pct_train;
    Results_PLSR.(model_name).reduction.mae       = tbl_mae;        % AbsErr moyennée
    Results_PLSR.(model_name).reduction.mae_all   = tbl_mae_all;    % AbsErr par run
    Results_PLSR.(model_name).reduction.mae_std   = nanstd(tbl_mae_all, 0, 1); % ±SD
    Results_PLSR.(model_name).reduction.r2        = tbl_r2;
    Results_PLSR.(model_name).reduction.best_lvl  = best_lvl;
    Results_PLSR.(model_name).reduction.best_thresh = best_thresh;

    %% -------------------------------------------------------------------
    %  PASSE 2 — 50 répétitions sur le meilleur sous-ensemble
    %  VIP calculé sur jeu TRAIN de chaque fold → variance réelle
    %  B moyenné sur folds où feature sélectionnée uniquement
    %  -------------------------------------------------------------------
    fprintf('  Passe 2 : 50 répétitions sur meilleur sous-ensemble...\n');

    % Sous-ensemble optimal
    feat_opt = feat_used(best_sel);
    src_opt  = src_used(best_sel);
    X_opt    = X_z(:, best_sel);
    N_opt    = sum(best_sel);
    n_lat_opt= tbl_nlatent(best_lvl);

    VIP_all  = NaN(N_rep, N_opt);
    B_all    = NaN(N_rep, N_opt);
    MAE_all  = NaN(N_rep, 1);
    RMSE_all = NaN(N_rep, 1);
    R2_all   = NaN(N_rep, 1);

    for iRep = 1:N_rep
        subj_order  = subj_ids(randperm(N_s));
        fold_assign = mod((1:N_s)'-1, K_fold) + 1;

        VIP_folds = NaN(K_fold, N_opt);
        B_folds   = NaN(K_fold, N_opt);
        Y_pred_cv = NaN(N_obs, 1);

        for iFold = 1:K_fold
            te_subjs = subj_order(fold_assign == iFold);
            te       = ismember(G, te_subjs);
            tr       = ~te;
            if sum(tr) < 10 || sum(te) < 1, continue; end

            X_tr = X_opt(tr,:);   Y_tr = Y_z(tr);
            X_te = X_opt(te,:);

            % N_latent fixe (optimal déterminé en Passe 1)
            n_c = min(n_lat_opt, numel(unique(G(tr)))-1);
            if n_c < 1, continue; end

            % VIP sur train
            [XL,YL,XS,~,~,~,~,stats] = plsregress(X_tr, Y_tr, n_c);
            W0  = stats.W ./ sqrt(sum(stats.W.^2,1));
            SC  = sum(XS.^2,1) .* (YL.^2);
            vip = sqrt(N_opt * ((W0.^2) * SC') / sum(SC));

            vip_fold = vip;   % pas de nouvelle réduction en Passe 2
            VIP_folds(iFold,:) = vip_fold';

            % Coefficients B
            [~,~,~,~,B_cv] = plsregress(X_tr, Y_tr, n_c);
            B_folds(iFold,:) = B_cv(2:end)';

            % Prédiction test
            Y_pred_cv(te) = [ones(sum(te),1), X_te] * B_cv;
        end

        VIP_all(iRep,:) = nanmean(VIP_folds, 1);
        B_all(iRep,:)   = nanmean(B_folds,   1);

        Yp = Y_pred_cv * sd_Y + mu_Y;
        ok = ~isnan(Yp);
        if sum(ok) > 2
            MAE_all(iRep)  = mean(abs(Yp(ok) - Y(ok)));
            RMSE_all(iRep) = sqrt(mean((Yp(ok)-Y(ok)).^2));
            R2_all(iRep)   = 1 - sum((Y(ok)-Yp(ok)).^2) / ...
                                 sum((Y(ok)-mean(Y(ok))).^2);
        end

        if mod(iRep,10)==0
            fprintf('  Rep %d/%d  AbsErr=%.2f  R²=%.3f\n', iRep, N_rep, ...
                nanmean(MAE_all(1:iRep)), nanmean(R2_all(1:iRep)));
        end
    end

    %% --- Agrégation ----------------------------------------------------
    VIP_mean = nanmean(VIP_all)';
    VIP_std  = nanstd(VIP_all)';
    VIP_p025 = prctile(VIP_all, 2.5)';
    VIP_p975 = prctile(VIP_all, 97.5)';
    sel_freq = mean(~isnan(VIP_all))';
    B_mean   = nanmean(B_all)';
    B_std    = nanstd(B_all)';
    B_p025   = prctile(B_all, 2.5)';
    B_p975   = prctile(B_all, 97.5)';

    fprintf('\n--- %s ---\n', model_name);
    fprintf('  AbsErr=%.2f±%.2f  RMSE=%.2f±%.2f  R²=%.3f±%.3f\n', ...
        nanmean(MAE_all), nanstd(MAE_all), ...
        nanmean(RMSE_all),nanstd(RMSE_all), ...
        nanmean(R2_all),  nanstd(R2_all));

    [~,ord] = sort(VIP_mean,'descend','MissingPlacement','last');
    fprintf('\n  Top 15 (VIP | B_mean | signe) :\n');
    for k = 1:min(15,N_opt)
        fi = ord(k);
        if B_mean(fi) >= 0, sgn = '+'; else, sgn = '-'; end
        fprintf('    %-42s  VIP=%.3f±%.3f  B=%+.4f±%.4f  [%s]\n', ...
            feat_opt{fi}, VIP_mean(fi), VIP_std(fi), B_mean(fi), B_std(fi), sgn);
    end

    % VIP_full : VIP de la Passe 1 sur le modèle complet (N_fu variables)
    % Utilisé pour les figures comparatives (équivalent Table III Goubault)
    Results_PLSR.(model_name).feat_full = feat_used;   % toutes les N_fu variables
    Results_PLSR.(model_name).VIP_full  = VIP_full;    % VIP modèle complet Passe 1

    Results_PLSR.(model_name).feat_used = feat_opt;
    Results_PLSR.(model_name).src_used  = src_opt;
    Results_PLSR.(model_name).VIP_mean  = VIP_mean;
    Results_PLSR.(model_name).VIP_std   = VIP_std;
    Results_PLSR.(model_name).VIP_p025  = VIP_p025;
    Results_PLSR.(model_name).VIP_p975  = VIP_p975;
    Results_PLSR.(model_name).sel_freq  = sel_freq;
    Results_PLSR.(model_name).B_mean    = B_mean;
    Results_PLSR.(model_name).B_std     = B_std;
    Results_PLSR.(model_name).B_p025    = B_p025;
    Results_PLSR.(model_name).B_p975    = B_p975;
    Results_PLSR.(model_name).AbsErr    = nanmean(MAE_all);
    Results_PLSR.(model_name).RMSE      = nanmean(RMSE_all);
    Results_PLSR.(model_name).R2        = nanmean(R2_all);
    Results_PLSR.(model_name).VIP_all   = VIP_all;
    Results_PLSR.(model_name).B_all     = B_all;
    Results_PLSR.(model_name).MAE_all  = MAE_all;
    Results_PLSR.(model_name).RMSE_all = RMSE_all;
    Results_PLSR.(model_name).R2_all   = R2_all;
    fprintf('\n');
end


%%%  -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save,'PLSR_Results.mat'), 'Results_PLSR');
fprintf('Sauvegardé : PLSR_Results.mat\n\n');

% Propager le renommage aux feat_used stockés dans Results_PLSR
% Doit être fait APRÈS la boucle PLSR, avant les figures
mn_list = fieldnames(Results_PLSR);
for iMn = 1:numel(mn_list)
    mn = mn_list{iMn};
    if ~isfield(Results_PLSR.(mn),'feat_used'), continue; end
    fu = Results_PLSR.(mn).feat_used;
    for iF = 1:numel(fu)
        if rename_map.isKey(fu{iF})
            fu{iF} = rename_map(fu{iF});
        end
    end
    Results_PLSR.(mn).feat_used = fu;
    if isfield(Results_PLSR.(mn),'feat_full')
        ff = Results_PLSR.(mn).feat_full;
        for iF = 1:numel(ff)
            if rename_map.isKey(ff{iF})
                ff{iF} = rename_map(ff{iF});
            end
        end
        Results_PLSR.(mn).feat_full = ff;
    end
end
fprintf('Renommage propagé aux feat_used de Results_PLSR\n\n');


%% -----------------------------------------------------------------------
%  FIGURE TABLEAU RÉDUCTION — Style Goubault Table II
%  Pour chaque modèle : N_input, N_latent, %RPE train, AbsErr test
%  -----------------------------------------------------------------------
fprintf('Génération Fig_Tableau_Reduction...\n');

model_names_r  = {'M1_SH_Std','M2_SH_Ergo','M3_LH_Std','M4_Pool_Std'};
model_labels_r = {'M1 SH Std','M2 SH Ergo','M3 LH Std','M4 Pool Std'};
colors_r       = {[0.2 0.4 0.8],[0.15 0.6 0.35],[0.85 0.65 0.0],[0.8 0.3 0.1]};

fig_tbl = figure('Name','Tableau_Reduction','NumberTitle','off',...
    'Position',[30 30 1400 900],'Color','w');

for iM = 1:4
    mn  = model_names_r{iM};
    if ~isfield(Results_PLSR,mn) || ~isfield(Results_PLSR.(mn),'reduction')
        continue;
    end
    rd = Results_PLSR.(mn).reduction;

    % Nettoyer les NaN
    valid_lvl = ~isnan(rd.n_input) & ~isnan(rd.mae);
    thresh_v  = rd.thresh(valid_lvl);
    ninput_v  = rd.n_input(valid_lvl);
    nlatent_v = rd.n_latent(valid_lvl);
    pct_v     = rd.pct_train(valid_lvl);
    mae_v     = rd.mae(valid_lvl);
    r2_v      = rd.r2(valid_lvl);
    N_lvl     = sum(valid_lvl);
    best_idx  = rd.best_lvl;
    if best_idx > N_lvl, best_idx = N_lvl; end

    ax = subplot(2,2,iM);
    hold(ax,'on');

    % ── Courbe Absolute Error (axe gauche) ───────────────────────────
    yyaxis(ax,'left');
    mae_sd = rd.mae_std(valid_lvl);
    % Enveloppe ±SD
    fill(ax, [1:N_lvl, fliplr(1:N_lvl)], ...
        [max(mae_v-mae_sd,0), fliplr(mae_v+mae_sd)], ...
        colors_r{iM}, 'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
    plot(ax, 1:N_lvl, mae_v, 'o-', ...
        'Color', colors_r{iM}, 'LineWidth', 2.0, ...
        'MarkerSize', 5, 'MarkerFaceColor', colors_r{iM}, ...
        'HandleVisibility','off');
    % Meilleur modèle — étoile
    plot(ax, best_idx, mae_v(best_idx), 'p', ...
        'MarkerSize', 14, 'MarkerFaceColor', [1 0.8 0], ...
        'MarkerEdgeColor', [0.6 0.4 0], 'LineWidth', 1.5, ...
        'HandleVisibility','off');
    ylabel(ax, 'Absolute Error (unités RPE)', 'FontSize', 8);
    ax.YColor = colors_r{iM};

    % ── Courbe %RPE train (axe droit) ────────────────────────────────
    yyaxis(ax,'right');
    plot(ax, 1:N_lvl, pct_v, 's--', ...
        'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.7 0.7 0.7], ...
        'HandleVisibility','off');
    ylabel(ax, '% RPE expliqué (train)', 'FontSize', 8);
    ax.YColor = [0.3 0.3 0.3];
    ylim(ax, [0 105]);

    % ── Annotations N_input ──────────────────────────────────────────
    yyaxis(ax,'left');
    for k = 1:N_lvl
        text(ax, k, min(mae_v)*0.92, sprintf('%d', ninput_v(k)), ...
            'HorizontalAlignment','center','FontSize',6,...
            'Color',[0.3 0.3 0.3],'VerticalAlignment','top');
    end
    text(ax, N_lvl/2+0.5, min(mae_v)*0.86, 'N variables →', ...
        'HorizontalAlignment','center','FontSize',6,...
        'Color',[0.5 0.5 0.5],'FontAngle','italic');

    % ── Axes ─────────────────────────────────────────────────────────
    set(ax,'XTick',1:N_lvl,...
        'XTickLabel',arrayfun(@(t) sprintf('%.1f',t), thresh_v,...
        'UniformOutput',false),...
        'XTickLabelRotation',45,'FontSize',7,'XGrid','on','Box','on');
    xlabel(ax, 'Seuil VIP', 'FontSize', 8);
    xlim(ax, [0.5, N_lvl+0.5]);

    title(ax, sprintf('%s  |  Meilleur : seuil=%.1f  N=%d  Nlat=%d  AbsErr=%.2f',...
        model_labels_r{iM}, thresh_v(best_idx), ninput_v(best_idx),...
        nlatent_v(best_idx), mae_v(best_idx)),...
        'FontSize', 8, 'FontWeight','bold', 'Color', colors_r{iM});

    % ── Légende manuelle (entrées propres uniquement) ─────────────────
    yyaxis(ax,'left');
    h1 = plot(ax,NaN,NaN,'o-','Color',colors_r{iM},'LineWidth',2);
    h2 = plot(ax,NaN,NaN,'p','MarkerSize',10,'MarkerFaceColor',[1 0.8 0],...
        'MarkerEdgeColor',[0.6 0.4 0]);
    yyaxis(ax,'right');
    h3 = plot(ax,NaN,NaN,'s--','Color',[0.4 0.4 0.4],'LineWidth',1.2,...
        'MarkerFaceColor',[0.7 0.7 0.7]);
    legend(ax,[h1 h2 h3],{'Absolute Error (CV)','Meilleur modèle','% RPE train'},...
        'Location','northeast','FontSize',7);
end

sgtitle('Tableau de réduction itérative — Absolute Error (test CV) et %RPE (train) par seuil VIP',...
    'FontSize',11,'FontWeight','bold');

saveas(fig_tbl, fullfile(path_save,'Fig_Tableau_Reduction.png')); close(fig_tbl);
fprintf('Figure sauvegardée : Fig_Tableau_Reduction.png\n\n');

%% -----------------------------------------------------------------------
%  PRÉPARATION COMMUNES AUX FIGURES
%  -----------------------------------------------------------------------
model_names  = {'M1_SH_Std','M2_SH_Ergo','M3_LH_Std','M4_Pool_Std'};
model_labels = {'M1 SH Std','M2 SH Ergo','M3 LH Std','M4 Pool Std'};
colors       = {[0.2 0.4 0.8],[0.15 0.6 0.35],[0.85 0.65 0.0],[0.8 0.3 0.1]};
N_top        = 10;

% ── Passe 2 : VIP modèle optimal (50 répétitions) ──────────────────────
% Utilisé pour TOUTES les figures sauf Fig_Tableau_Reduction
% B signé depuis Passe 2 également
all_vip      = NaN(N_feat, 4);   % VIP_mean Passe 2
all_vip_std  = NaN(N_feat, 4);
all_vip_p025 = NaN(N_feat, 4);
all_vip_p975 = NaN(N_feat, 4);
all_B        = NaN(N_feat, 4);
all_B_std    = NaN(N_feat, 4);

for iM = 1:4
    mn = model_names{iM};
    if ~isfield(Results_PLSR,mn), continue; end
    fu = Results_PLSR.(mn).feat_used;
    for iF = 1:numel(fu)
        idx = find(strcmp(all_feat,fu{iF}),1);
        if isempty(idx), continue; end
        all_vip(idx,iM)      = Results_PLSR.(mn).VIP_mean(iF);
        all_vip_std(idx,iM)  = Results_PLSR.(mn).VIP_std(iF);
        all_vip_p025(idx,iM) = Results_PLSR.(mn).VIP_p025(iF);
        all_vip_p975(idx,iM) = Results_PLSR.(mn).VIP_p975(iF);
        all_B(idx,iM)        = Results_PLSR.(mn).B_mean(iF);
        all_B_std(idx,iM)    = Results_PLSR.(mn).B_std(iF);
    end
end

% ── Variables avec VIP ≥ 1 dans le modèle optimal Passe 2 ───────────────
% Critère : VIP ≥ 1 (significativement explicatives dans le meilleur modèle)
% Remplace le top 10 arbitraire — le nombre varie selon le modèle
vip1_idx    = cell(1, 4);   % indices dans all_feat
vip1_labels = cell(1, 4);   % labels correspondants

for iM = 1:4
    % Sélectionner les variables avec VIP ≥ 1, triées par VIP décroissant
    sel_vip1 = find(~isnan(all_vip(:,iM)) & all_vip(:,iM) >= 1.0);
    [~, ord_v1] = sort(all_vip(sel_vip1,iM), 'descend');
    sel_sorted = sel_vip1(ord_v1);
    vip1_idx{iM}    = sel_sorted;
    vip1_labels{iM} = strrep(all_feat(sel_sorted),'_',' ')';
    fprintf('  %s : %d variables VIP≥1 dans modèle optimal\n', ...
        model_names{iM}, numel(sel_sorted));
end

% Conserver top10 pour Fig_VIP_IC95 (top 10 du modèle optimal pour lisibilité)
top10_idx    = NaN(N_top, 4);
top10_labels = cell(N_top, 4);
for iM = 1:4
    n_avail = numel(vip1_idx{iM});
    if n_avail >= N_top
        top10_idx(:,iM)    = vip1_idx{iM}(1:N_top);
        top10_labels(:,iM) = vip1_labels{iM}(1:N_top);
    else
        % Moins de 10 variables VIP≥1 : compléter avec les meilleures VIP<1
        [~,ord_all] = sort(all_vip(:,iM),'descend','MissingPlacement','last');
        top10_idx(:,iM)    = ord_all(1:N_top);
        top10_labels(:,iM) = strrep(all_feat(ord_all(1:N_top)),'_',' ')';
    end
end
% Alias pour Fig_VIP_IC95
top10_idx_p2    = top10_idx;
top10_labels_p2 = top10_labels;
all_vip_p2      = all_vip;

%% -----------------------------------------------------------------------
%  FIGURE 2 — IC 95% VIP par modèle
%  -----------------------------------------------------------------------
fig2 = figure('Name','VIP_IC95','NumberTitle','off',...
    'Position',[50 50 1400 720],'Color','w');
for iM = 1:4
    subplot(2,2,iM); hold on;
    ord = top10_idx_p2(:,iM);           % ← Passe 2
    vm  = all_vip_p2(ord,iM);           % ← Passe 2
    vlo = all_vip_p025(ord,iM);
    vhi = all_vip_p975(ord,iM);
    for k = 1:N_top
        pos = N_top+1-k;
        xl  = max(vlo(k),0);  xh = vhi(k);
        if ~isnan(xl) && ~isnan(xh)
            fill([xl xh xh xl],[pos-0.40 pos-0.40 pos+0.40 pos+0.40],...
                colors{iM},'FaceAlpha',0.15,'EdgeColor','none');
        end
    end
    barh(N_top:-1:1, vm, 0.28, 'FaceColor',colors{iM},'EdgeColor','none');
    for k = 1:N_top
        pos = N_top+1-k;
        xl  = max(vlo(k),0);  xh = vhi(k);
        if ~isnan(xl) && ~isnan(xh)
            plot([xl xh],[pos pos],'k-','LineWidth',2.0);
            plot([xl xl],[pos-0.30 pos+0.30],'k-','LineWidth',1.4);
            plot([xh xh],[pos-0.30 pos+0.30],'k-','LineWidth',1.4);
        end
    end
    xline(1,'--k','LineWidth',1.0);
    set(gca,'YTick',1:N_top,'YTickLabel',flip(top10_labels_p2(:,iM)),...
        'FontSize',8,'XGrid','on','TickDir','out');
    xlabel('VIP moyen  [IC 95%]','FontSize',8);
    title(model_labels{iM},'FontSize',10,'FontWeight','bold','Color',colors{iM});
    xlim([0, max(nanmax(vhi)*1.08, 1.8)]);
end
sgtitle('Stabilité VIP — IC 95% sur 50 répétitions (top 10 par modèle)',...
    'FontSize',11,'FontWeight','bold');
saveas(fig2,fullfile(path_save,'Fig_VIP_IC95.png')); close(fig2);









fprintf('Figures sauvegardées :\n');
fprintf('  Fig_VIP_IC95.png\n');
fprintf('  Fig_VIP_Intersection_Groupes.png\n');
fprintf('  Fig_Heatmap_Convergence.png\n');
%% -----------------------------------------------------------------------
%  TABLEAU 1 — Top 10 PLSR Robin (4 modèles)
%  -----------------------------------------------------------------------
model_names_t  = {'M1_SH_Std','M2_SH_Ergo','M3_LH_Std','M4_Pool_Std'};
model_labels_t = {'M1 SH Std','M2 SH Ergo','M3 LH Std','M4 Pool Std'};
colors_t = {[0.2 0.4 0.8],[0.15 0.6 0.35],[0.85 0.65 0.0],[0.8 0.3 0.1]};
N_top_t  = 10;

fig_t1 = figure('Name','Tableau_Top10_PLSR','NumberTitle','off',...
    'Position',[30 30 1400 420],'Color','w');

sgtitle('Top 10 — PLSR Robin  (↑ B>0 | ↓ B<0  |  VIP Passe 2)',...
    'FontSize',11,'FontWeight','bold');

for iM = 1:4
    mn = model_names_t{iM};
    if ~isfield(Results_PLSR,mn), continue; end
    fu = Results_PLSR.(mn).feat_used;
    vm = Results_PLSR.(mn).VIP_mean;
    bm = Results_PLSR.(mn).B_mean;
    [~,ord] = sort(vm,'descend','MissingPlacement','last');
    n_show  = min(N_top_t, numel(fu));

    ax = subplot(1,4,iM);
    axis(ax,'off');
    set(ax,'XLim',[0 1],'YLim',[0 1]);

    % Titre colonne
    text(ax,0.5,0.97, model_labels_t{iM},...
        'FontSize',9,'FontWeight','bold','Color',colors_t{iM},...
        'HorizontalAlignment','center','VerticalAlignment','top');

    % Ligne séparatrice
    plot(ax,[0 1],[0.93 0.93],'Color',colors_t{iM},'LineWidth',1.2);
    hold(ax,'on');

    % En-têtes
    row_h = 0.082;
    y0    = 0.91;
    text(ax,0.04,y0,'#',   'FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    text(ax,0.14,y0,'Variable','FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    text(ax,0.80,y0,'VIP', 'FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    text(ax,0.93,y0,'Dir', 'FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    plot(ax,[0 1],[y0-0.02 y0-0.02],'Color',[0.8 0.8 0.8],'LineWidth',0.5);

    for k = 1:n_show
        fi  = ord(k);
        vip = vm(fi);
        b   = bm(fi);
        nm  = strrep(fu{fi},'_',' ');
        if length(nm)>26, nm=[nm(1:23) '...']; end
        if b>=0, dir_s='↑'; dc=[0.7 0.1 0.1]; else, dir_s='↓'; dc=[0.1 0.1 0.7]; end

        yr = y0 - k*row_h;
        if mod(k,2)==0
            fill(ax,[0 1 1 0],[yr-0.01 yr-0.01 yr+row_h-0.01 yr+row_h-0.01],...
                [0.96 0.96 0.96],'EdgeColor','none');
        end
        text(ax,0.04,yr,sprintf('%d',k),    'FontSize',7,'Color',[0.5 0.5 0.5]);
        text(ax,0.14,yr,nm,                  'FontSize',6.8,'Color',[0.1 0.1 0.1]);
        text(ax,0.80,yr,sprintf('%.3f',vip),'FontSize',7,'Color',[0.3 0.3 0.3]);
        text(ax,0.93,yr,dir_s,'FontSize',9,'FontWeight','bold','Color',dc);
    end
end

saveas(fig_t1,fullfile(path_save,'Tableau_Top10_PLSR.png')); close(fig_t1);
fprintf('Tableau sauvegardé : Tableau_Top10_PLSR.png\n');

%% -----------------------------------------------------------------------
%  TABLEAU 2 — Référence : LMM Robin + Goubault Chord + Digital
%  -----------------------------------------------------------------------
lmm_rows = {
    'Hand MedianFreq AccelY',              '↑','β=+0.96';
    'Head Mean AccelMagnitude',            '↑','β=+1.09';
    'Forearm SpEnt AngVelMagnitude',       '↑','β=+0.90';
    'Head PeakPower AngVelX',             '↑','β=+0.56';
    'Hand Mean AccelMagnitude',            '↓','β=−1.12';
    'Hand PeakPower AccelMagnitude',       '↓','β=−0.98';
    'SampleEntropy Biceps',               '↓','β=−0.88';
    'TFR SpEnt DeltAnt',                  '↓','β=−0.59';
    'TFR MedianFreq Triceps',             '↓','β=−0.58';
};

chord_rows = {
    'Median freq / Acc.y / Hand',          '↑','';
    '90th pct / Acc.x / Hand',             '↓','';
    'Peak-Power / AngVel.x / Head',        '↑','';
    'Mean / Acc.y / Forearm',              '↓','';
    'SpEnt / AngVel.mag / Forearm',        '↑','';
    'Mean / AngVel.x / Shoulder',          '↑','';
    'PeakPowerFreq / AngVel.mag / Forearm','↓','';
    'Peak-Power / Acc.mag / Hand',         '↓','';
    'Mean / AngVel.y / Shoulder',          '↑','';
    'SpEnt / Acc.y / Hand',                '↑','';
};

digital_rows = {
    'SpEnt / Acc.x / Head',                '↑','';
    'PeakPowerFreq / AngVel.mag / Forearm','↓','';
    'PeakPowerFreq / Acc.mag / Shoulder',  '↓','';
    '75th pct / Jerk.y / Hand',            '↑','';
    'PeakPowerFreq / AngVel.mag / Arm',    '↓','';
    'Peak-Power / AngVel.x / Arm',         '↓','';
    'Peak-Power / Acc.z / Trunk',          '↑','';
    'SpEnt / Jerk.z / Shoulder',           '↑','';
    'SpEnt / Acc.mag / Arm',               '↓','';
    'Median / AngVel.y / Head',            '~','';
};

ref_data   = {lmm_rows, chord_rows, digital_rows};
ref_titles = {'LMM Robin  (9 variables)','Goubault — Chord','Goubault — Digital'};
ref_cols   = {[0.2 0.5 0.2],[0.7 0.35 0.1],[0.2 0.2 0.6]};

fig_t2 = figure('Name','Tableau_Reference','NumberTitle','off',...
    'Position',[30 30 1300 440],'Color','w');

sgtitle('Référence — LMM Robin et PLSR Goubault  (↑ augmente | ↓ diminue avec RPE)',...
    'FontSize',11,'FontWeight','bold');

for iR = 1:3
    rows = ref_data{iR};
    rc   = ref_cols{iR};
    N_r  = size(rows,1);
    has_beta = iR==1;

    ax = subplot(1,3,iR);
    axis(ax,'off');
    set(ax,'XLim',[0 1],'YLim',[0 1]);
    hold(ax,'on');

    text(ax,0.5,0.97,ref_titles{iR},...
        'FontSize',9,'FontWeight','bold','Color',rc,...
        'HorizontalAlignment','center','VerticalAlignment','top');
    plot(ax,[0 1],[0.93 0.93],'Color',rc,'LineWidth',1.2);

    row_h = 0.082;
    y0    = 0.91;
    text(ax,0.04,y0,'#',       'FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    text(ax,0.12,y0,'Variable','FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    text(ax,0.80,y0,'Dir',     'FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    if has_beta
        text(ax,0.88,y0,'|β|','FontSize',7,'FontWeight','bold','Color',[0.4 0.4 0.4]);
    end
    plot(ax,[0 1],[y0-0.02 y0-0.02],'Color',[0.8 0.8 0.8],'LineWidth',0.5);

    for k = 1:N_r
        nm  = rows{k,1};
        dir_s = rows{k,2};
        val   = rows{k,3};
        if length(nm)>32, nm=[nm(1:29) '...']; end
        if strcmp(dir_s,'↑'), dc=[0.7 0.1 0.1];
        elseif strcmp(dir_s,'↓'), dc=[0.1 0.1 0.7];
        else, dc=[0.5 0.5 0.5]; end

        yr = y0 - k*row_h;
        if mod(k,2)==0
            fill(ax,[0 1 1 0],[yr-0.01 yr-0.01 yr+row_h-0.01 yr+row_h-0.01],...
                [0.96 0.96 0.96],'EdgeColor','none');
        end
        text(ax,0.04,yr,sprintf('%d',k),   'FontSize',7,'Color',[0.5 0.5 0.5]);
        text(ax,0.12,yr,nm,                 'FontSize',6.8,'Color',[0.1 0.1 0.1]);
        text(ax,0.80,yr,dir_s,'FontSize',9,'FontWeight','bold','Color',dc);
        if has_beta && ~isempty(val)
            text(ax,0.88,yr,val,'FontSize',7,'Color',[0.3 0.3 0.3]);
        end
    end
end

saveas(fig_t2,fullfile(path_save,'Tableau_Reference.png')); close(fig_t2);
fprintf('Tableau sauvegardé : Tableau_Reference.png\n');

fprintf('\n=== TERMINÉ ===\n');