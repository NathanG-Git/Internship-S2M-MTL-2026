%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LASSO + LMM - Version Finale Corrigée et Sécurisée %%%%
%%%% Version restructurée - Group LASSO FISTA + GroupKFold %%%%
%%%% %%%%
%%%% Pipeline : %%%%
%%%% 1) Chargement de toutes les variables (sans filtrage Spearman) %%%%
%%%% 2) Construction des tableaux longs par bloc %%%%
%%%% 3b) Group LASSO FISTA + GroupKFold par sujet + Mundlak HAC %%%%
%%%% 4) LMM par bloc sur les variables sélectionnées par 3b %%%%
%%%%    - Correction autocorrélation résiduelle : Cochrane-Orcutt %%%%
%%%%      itératif (Prais-Winsten, convergence phi < 0.001) %%%%
%%%%    - Vérification homoscédasticité : Breusch-Pagan simplifié %%%%
%%%%      (corrélation Spearman résidus² ~ temps) %%%%
%%%% 5) Figures LMM %%%%
%%%% 6) Boucle de stabilité Group LASSO (10 graines) %%%%
%%%% 7) Sauvegarde %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
rng(42);

%% -----------------------------------------------------------------------
% CHEMINS ET CONFIGURATION
% -----------------------------------------------------------------------
PathSave = fileparts(mfilename('fullpath'));
if isempty(PathSave), PathSave = pwd; end

% Modifie ces chemins si nécessaire pour pointer vers tes répertoires locaux
PathG1_IMU = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';
PathG2_IMU = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_long_duration\';
PathG1_EMG = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\1_EMG\';
PathG2_EMG = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_long_duration\2_EMG\';

%% -----------------------------------------------------------------------
% CHARGEMENT DES DONNÉES
% -----------------------------------------------------------------------
fprintf('Chargement des données en cours...\n');

try
    d1 = load(fullfile(PathG1_IMU, 'Workload_Li_TimeNormalised.mat'));
    d2 = load(fullfile(PathG2_IMU, 'Workload_Li_TimeNormalised_G2.mat'));
    dG1 = load(fullfile(PathG1_IMU, 'Workload_Li_TimeNormalised_Goubault.mat'));
    dG2 = load(fullfile(PathG2_IMU, 'Workload_Li_TimeNormalised_Goubault_G2.mat'));
    eG1 = load(fullfile(PathG1_EMG, 'EMG_Li_TimeNormalised_G1.mat'));
    eG2 = load(fullfile(PathG2_EMG, 'EMG_Li_TimeNormalised_G2.mat'));
    fprintf('--> Chargement des fichiers réussi.\n');
catch ME
    warning('Erreur de chargement des chemins OneDrive. Utilisation du dossier courant.');
    % Fallback automatique sur le dossier courant si OneDrive est inaccessible
    d1 = load('Workload_Li_TimeNormalised.mat');
    d2 = load('Workload_Li_TimeNormalised_G2.mat');
    dG1 = load('Workload_Li_TimeNormalised_Goubault.mat');
    dG2 = load('Workload_Li_TimeNormalised_Goubault_G2.mat');
    eG1 = load('EMG_Li_TimeNormalised_G1.mat');
    eG2 = load('EMG_Li_TimeNormalised_G2.mat');
end

Vi1_26 = find(d1.valid_26);
Vi2_26 = find(d2.valid_26);
Vi1 = find(d1.valid_subj);
Vi2 = find(d2.valid_subj);
Vi1g = find(dG1.valid_subj);
Vi2g = find(dG2.valid_subj);
Vi1e = find(eG1.valid_subj);
Vi2e = find(eG2.valid_subj);

N_bins = 10;
Segments = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Axes_lbl = {'X','Y','Z'};
N_seg = 7;
Channels = eG1.Channels;
VarNames_emg = eG1.VarNames;
N_ch = length(Channels);
MIN_BINS = 8;

%% -----------------------------------------------------------------------
%  CORRECTION INF/NAN DANS MAT_EMG (SampleEntropy)
%  Etape 1 : Remplacer Inf par NaN puis imputer par moyenne intra-individuelle
%  Etape 2 : Exclure les sujets avec 100% de valeurs invalides sur un canal
%            (signal EMG absent sur toute la session — pas imputable)
%            Ces sujets sont exclus uniquement du bloc EMG,
%            ils restent dans SegMod et Goubault.
%  -----------------------------------------------------------------------
fprintf('\n=== Correction Inf/NaN SampleEntropy dans Mat_EMG ===\n');
for iG = 1:2
    if iG==1, Mat=eG1.Mat_EMG; VN=eG1.VarNames; gname='G1';
    else,     Mat=eG2.Mat_EMG; VN=eG2.VarNames; gname='G2'; end

    % Etape 1 : Inf → NaN puis imputation intra-individuelle
    n_inf = sum(isinf(Mat(:)));
    Mat(isinf(Mat)) = NaN;
    [NS,NT,NC,NV] = size(Mat);
    for iP = 1:NS
        for iC = 1:NC
            for iV = 1:NV
                x = squeeze(Mat(iP,:,iC,iV));
                bad = isnan(x);
                if any(bad) && ~all(bad)
                    x(bad) = mean(x(~bad));
                    Mat(iP,:,iC,iV) = x;
                end
            end
        end
    end

    % Etape 2 : identifier les sujets avec 100% invalide sur au moins un canal
    iV_se = find(strcmp(VN, 'SampleEntropy'));
    excl = [];
    for iP = 1:NS
        for iC = 1:NC
            x = squeeze(Mat(iP,:,iC,iV_se));
            if all(isnan(x))
                excl(end+1) = iP;
                break
            end
        end
    end
    excl = unique(excl);
    fprintf('  %s : %d Inf corrigés | %d sujets exclus du bloc EMG : [%s]\n',...
        gname, n_inf, length(excl), num2str(excl));

    if iG==1
        eG1.Mat_EMG = Mat;
        % Exclure de Vi1e les sujets invalides
        % excl contient les indices dans Mat_EMG (1..N_subj_emg)
        % Vi1e contient les indices absolus dans la matrice originale
        Vi1e_valid = Vi1e;
        Vi1e_valid(excl) = [];
        Vi1e = Vi1e_valid;
    else
        eG2.Mat_EMG = Mat;
        Vi2e_valid = Vi2e;
        Vi2e_valid(excl) = [];
        Vi2e = Vi2e_valid;
    end
end

fprintf('  Vi1e après exclusion : %d sujets EMG G1\n', length(Vi1e));
fprintf('  Vi2e après exclusion : %d sujets EMG G2\n', length(Vi2e));

fprintf('G1: Features3 n=%d | Goubault n=%d | EMG n=%d\n', length(Vi1_26), length(Vi1g), length(Vi1e));
fprintf('G2: Features3 n=%d | Goubault n=%d | EMG n=%d\n', length(Vi2_26), length(Vi2g), length(Vi2e));

%% -----------------------------------------------------------------------
% ÉTAPE 1 : EXTRACTION DE TOUTES LES VARIABLES SANS FILTRAGE SPEARMAN
% -----------------------------------------------------------------------
fprintf('\n=== Étape 1 : Indexation des variables brutes ===\n');
all_vnames_test = {};
all_pG1 = []; all_pG2 = [];
all_rhoG1 = []; all_rhoG2 = [];
all_source = {};
all_Vi1 = {}; all_Vi2 = {};
all_Mat1 = {}; all_Mat2 = {};
all_RPE1 = {}; all_RPE2 = {};

% --- BySegMod (Features3) ---
for iSig = 1:2
    if iSig==1
        sn='Accel'; m1=d1.Mat_Seg_Accel_26; m2=d2.Mat_Seg_Accel_26;
    else
        sn='Jerk'; m1=d1.Mat_Seg_Jerk_26; m2=d2.Mat_Seg_Jerk_26;
    end
    for iSeg = 1:N_seg
        [r1,p1] = fisherTtest(squeeze(m1(:,:,iSeg)), d1.Mat_RPE_26, Vi1_26, MIN_BINS);
        [r2,p2] = fisherTtest(squeeze(m2(:,:,iSeg)), d2.Mat_RPE_26, Vi2_26, MIN_BINS);
        vn = [sn '_' Segments{iSeg}];
        all_vnames_test{end+1} = vn;
        all_pG1(end+1) = p1; all_pG2(end+1) = p2;
        all_rhoG1(end+1) = r1; all_rhoG2(end+1) = r2;
        all_source{end+1} = 'SegMod';
        all_Vi1{end+1} = Vi1_26; all_Vi2{end+1} = Vi2_26;
        all_Mat1{end+1} = squeeze(m1(:,:,iSeg)); all_Mat2{end+1} = squeeze(m2(:,:,iSeg));
        all_RPE1{end+1} = d1.Mat_RPE_26; all_RPE2{end+1} = d2.Mat_RPE_26;
    end
end

% --- Total (Features3) ---
for iSig = 1:2
    if iSig==1
        sn='Total_Accel'; m1=d1.Mat_Total_Accel_26; m2=d2.Mat_Total_Accel_26;
    else
        sn='Total_Jerk'; m1=d1.Mat_Total_Jerk_26; m2=d2.Mat_Total_Jerk_26;
    end
    [r1,p1] = fisherTtest(m1, d1.Mat_RPE_26, Vi1_26, MIN_BINS);
    [r2,p2] = fisherTtest(m2, d2.Mat_RPE_26, Vi2_26, MIN_BINS);
    all_vnames_test{end+1} = sn;
    all_pG1(end+1)=p1; all_pG2(end+1)=p2;
    all_rhoG1(end+1)=r1; all_rhoG2(end+1)=r2;
    all_source{end+1}='SegMod';
    all_Vi1{end+1}=Vi1_26; all_Vi2{end+1}=Vi2_26;
    all_Mat1{end+1}=m1; all_Mat2{end+1}=m2;
    all_RPE1{end+1}=d1.Mat_RPE_26; all_RPE2{end+1}=d2.Mat_RPE_26;
end

% --- BySegAxis ---
for iSig = 1:2
    if iSig==1
        sn='Accel'; m1=d1.Mat_SegAxis_Accel; m2=d2.Mat_SegAxis_Accel;
    else
        sn='Jerk'; m1=d1.Mat_SegAxis_Jerk; m2=d2.Mat_SegAxis_Jerk;
    end
    for iSeg = 1:N_seg
        for iAx = 1:3
            [r1,p1] = fisherTtest(squeeze(m1(:,:,iSeg,iAx)), d1.Mat_RPE, Vi1, MIN_BINS);
            [r2,p2] = fisherTtest(squeeze(m2(:,:,iSeg,iAx)), d2.Mat_RPE, Vi2, MIN_BINS);
            vn = [sn '_Axe' Axes_lbl{iAx} '_' Segments{iSeg}];
            all_vnames_test{end+1} = vn;
            all_pG1(end+1)=p1; all_pG2(end+1)=p2;
            all_rhoG1(end+1)=r1; all_rhoG2(end+1)=r2;
            all_source{end+1}='SegAxis';
            all_Vi1{end+1}=Vi1; all_Vi2{end+1}=Vi2;
            all_Mat1{end+1}=squeeze(m1(:,:,iSeg,iAx));
            all_Mat2{end+1}=squeeze(m2(:,:,iSeg,iAx));
            all_RPE1{end+1}=d1.Mat_RPE; all_RPE2{end+1}=d2.Mat_RPE;
        end
    end
end

% --- Goubault ---
for iV = 1:size(dG1.Vars,1)
    [r1,p1] = fisherTtest(squeeze(dG1.Mat_Vars(:,:,iV)), dG1.Mat_RPE, Vi1g, MIN_BINS);
    [r2,p2] = fisherTtest(squeeze(dG2.Mat_Vars(:,:,iV)), dG2.Mat_RPE, Vi2g, MIN_BINS);
    vn = strrep(dG1.Vars{iV,1},' ','_');
    vn = strrep(vn,'-','_');
    vn = strrep(vn,'—','_');
    vn = ['Goub_' vn];
    all_vnames_test{end+1} = vn;
    all_pG1(end+1)=p1; all_pG2(end+1)=p2;
    all_rhoG1(end+1)=r1; all_rhoG2(end+1)=r2;
    all_source{end+1}='Goubault';
    all_Vi1{end+1}=Vi1g; all_Vi2{end+1}=Vi2g;
    all_Mat1{end+1}=squeeze(dG1.Mat_Vars(:,:,iV));
    all_Mat2{end+1}=squeeze(dG2.Mat_Vars(:,:,iV));
    all_RPE1{end+1}=dG1.Mat_RPE; all_RPE2{end+1}=dG2.Mat_RPE;
end

% --- EMG sans Amplitude ---
for iV = 1:length(VarNames_emg)
    if strcmp(VarNames_emg{iV},'Amplitude'), continue; end
    for iC = 1:N_ch
        [r1,p1] = fisherTtest(squeeze(eG1.Mat_EMG(:,:,iC,iV)), eG1.Mat_RPE, Vi1e, MIN_BINS);
        [r2,p2] = fisherTtest(squeeze(eG2.Mat_EMG(:,:,iC,iV)), eG2.Mat_RPE, Vi2e, MIN_BINS);
        vn = [VarNames_emg{iV} '_' Channels{iC}];
        all_vnames_test{end+1} = vn;
        all_pG1(end+1)=p1; all_pG2(end+1)=p2;
        all_rhoG1(end+1)=r1; all_rhoG2(end+1)=r2;
        all_source{end+1}='EMG';
        all_Vi1{end+1}=Vi1e; all_Vi2{end+1}=Vi2e;
        all_Mat1{end+1}=squeeze(eG1.Mat_EMG(:,:,iC,iV));
        all_Mat2{end+1}=squeeze(eG2.Mat_EMG(:,:,iC,iV));
        all_RPE1{end+1}=eG1.Mat_RPE; all_RPE2{end+1}=eG2.Mat_RPE;
    end
end

N_sig = length(all_vnames_test);
fprintf('--> %d variables brutes transmises au Group LASSO (Aucun filtrage Spearman préalable).\n', N_sig);

sel_vnames = all_vnames_test;
sel_source = all_source;
sel_Vi1 = all_Vi1;
sel_Vi2 = all_Vi2;
sel_Mat1 = all_Mat1;
sel_Mat2 = all_Mat2;
sel_RPE1 = all_RPE1;
sel_RPE2 = all_RPE2;

%% -----------------------------------------------------------------------
% BOUCLE EXTERNE SUR LES 4 VARIANTES DE FENÊTRE DE DIFFÉRENCE
% lag1 : DeltaRPE(t)=RPE(t)-RPE(t-1)  | dynamique fine, bruitée, 9 obs/sujet
% lag2 : DeltaRPE(t)=RPE(t)-RPE(t-2)  | lissage modéré, 8 obs/sujet
% lag3 : DeltaRPE(t)=RPE(t)-RPE(t-3)  | lissage fort,   7 obs/sujet
% global : DeltaRPE = RPE(10)-RPE(1)  | cohérent avec validation externe,
%          1 obs/sujet, perd la dynamique intra-session
% -----------------------------------------------------------------------
lag_modes      = {'lag1','lag2','lag3','global'};
lag_mode_lbls  = {'Lag-1','Lag-2','Lag-3','Fenêtre globale (10-1)'};

Results_AllLagModes = struct();

for iLag = 1:length(lag_modes)
lag_mode = lag_modes{iLag};
lag_lbl  = lag_mode_lbls{iLag};

fprintf('\n\n');
fprintf('#########################################################################\n');
fprintf('###  VARIANTE FENÊTRE : %-20s                                 ###\n', lag_lbl);
fprintf('#########################################################################\n');

%% -----------------------------------------------------------------------
% ÉTAPE 2 : CONSTRUCTION DES TABLEAUX LONGS PAR BLOC (DeltaX_within + X_between)
% -----------------------------------------------------------------------
fprintf('\n=== Étape 2 : Construction des tableaux longs par bloc [%s] ===\n', lag_lbl);
sources_uniq = unique(sel_source);
bloc_names = {'SegMod','SegAxis','Goubault','EMG'};
bloc_labels = {'BySegMod+Total','BySegAxis','Goubault','EMG'};

vnames_bloc = struct('SegMod',{{}},'SegAxis',{{}},'Goubault',{{}},'EMG',{{}});
% Xw = DeltaX_within (dynamique), Xb = X_between (niveau moyen sujet, constant)
Xw_G1_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
Xw_G2_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
Xb_G1_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
Xb_G2_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
RPE_G1_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
RPE_G2_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
subj_G1_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
subj_G2_bloc = struct('SegMod',[],'SegAxis',[],'Goubault',[],'EMG',[]);
first_bloc = struct('SegMod',true,'SegAxis',true,'Goubault',true,'EMG',true);

for iv = 1:N_sig
    src = sel_source{iv};
    vn = sel_vnames{iv};
    Vi_1 = sel_Vi1{iv}; Vi_2 = sel_Vi2{iv};
    M1 = sel_Mat1{iv}; M2 = sel_Mat2{iv};
    R1 = sel_RPE1{iv}; R2 = sel_RPE2{iv};

    [xw1,xb1,r1,s1] = buildBlock(Vi_1, M1, R1, N_bins, lag_mode);
    [xw2,xb2,r2,s2] = buildBlock(Vi_2, M2, R2, N_bins, lag_mode);

    vnames_bloc.(src){end+1} = vn;
    Xw_G1_bloc.(src) = [Xw_G1_bloc.(src) xw1];
    Xw_G2_bloc.(src) = [Xw_G2_bloc.(src) xw2];
    Xb_G1_bloc.(src) = [Xb_G1_bloc.(src) xb1];
    Xb_G2_bloc.(src) = [Xb_G2_bloc.(src) xb2];

    if first_bloc.(src)
        RPE_G1_bloc.(src) = r1; RPE_G2_bloc.(src) = r2;
        subj_G1_bloc.(src) = s1; subj_G2_bloc.(src) = s2;
        first_bloc.(src) = false;
    end
end

for ib = 1:4
    bn = bloc_names{ib};
    fprintf('  Bloc %-10s : %3d variables parentes | G1 n=%4d obs | G2 n=%4d obs\n', ...
        bn, length(vnames_bloc.(bn)), size(Xw_G1_bloc.(bn),1), size(Xw_G2_bloc.(bn),1));
end

%% -----------------------------------------------------------------------
% ÉTAPE 3b : GROUP LASSO FISTA + GROUP K-FOLD PAR SUJET + MUNDLAK HAC
% -----------------------------------------------------------------------
fprintf('\n=== Étape 3b : Group LASSO Mundlak & HAC Newey-West ===\n');
K_FOLDS   = 5;
N_LAMBDA  = 100;
% Note : la sélection du lambda utilise lambda_min_mse (pas la règle 1SE)
% ce qui est plus adapté quand le signal biomécanique est faible.
HAC_BW    = 2;
MAX_ITER = 1000;
TOL_FISTA = 1e-6;

mundlak_G1 = struct();
mundlak_G2 = struct();

for ib = 1:4
    bn = bloc_names{ib};
    vnames = vnames_bloc.(bn);
    if isempty(vnames), continue; end
    N_vars = length(vnames);

    for iGrp = 1:2
        if iGrp == 1
            Xw = Xw_G1_bloc.(bn); Xb = Xb_G1_bloc.(bn);
            Y = RPE_G1_bloc.(bn); S = subj_G1_bloc.(bn); gname = 'G1';
        else
            Xw = Xw_G2_bloc.(bn); Xb = Xb_G2_bloc.(bn);
            Y = RPE_G2_bloc.(bn); S = subj_G2_bloc.(bn); gname = 'G2';
        end

        ok_all = ~any(isnan(Xw)|isinf(Xw),2) & ~any(isnan(Xb)|isinf(Xb),2) & ~isnan(Y) & ~isinf(Y);
        Xw = Xw(ok_all,:); Xb = Xb(ok_all,:); Yok = Y(ok_all); Sok = S(ok_all);
        N_obs = length(Yok);

        % Diagnostic effectif attendu — dépend du lag_mode (géré par buildBlock)
        uSubj_tmp = unique(Sok);
        N_subj_tmp = length(uSubj_tmp);
        if N_subj_tmp <= 15
            K_FOLDS_loc = N_subj_tmp;   % LOSO
            fprintf('\n  [%s | %s | %s] : %d obs, %d vars, LOSO (%d folds)\n', bn, gname, lag_lbl, N_obs, N_vars, K_FOLDS_loc);
        elseif N_subj_tmp <= 20
            K_FOLDS_loc = 3;
            fprintf('\n  [%s | %s | %s] : %d obs, %d vars, K=%d folds\n', bn, gname, lag_lbl, N_obs, N_vars, K_FOLDS_loc);
        else
            K_FOLDS_loc = K_FOLDS;      % K=5 par défaut
            fprintf('\n  [%s | %s | %s] : %d obs, %d vars, K=%d folds\n', bn, gname, lag_lbl, N_obs, N_vars, K_FOLDS_loc);
        end

        % --- Mundlak déjà décomposé en amont par buildBlock ---
        % Xw = DeltaX_within(t) — dynamique de changement, varie dans le temps
        % Xb = X_between        — niveau moyen du sujet sur la session, constant
        uSubj = unique(Sok);
        N_subj = length(uSubj);
        X_within  = Xw;
        X_between = Xb;

        % --- Matrice Mundlak ordonnée [X_w1, X_b1, X_w2, X_b2, ...] ---
        N_cols = 2 * N_vars;
        X_mundlak = NaN(N_obs, N_cols);
        group_idx = NaN(1, N_cols);
        col_names = cell(1, N_cols);
        col_names_raw = cellfun(@(x) regexprep(x,'[^a-zA-Z0-9]','_'), ...
            vnames, 'UniformOutput', false);

        for iv = 1:N_vars
            cw = 2*iv-1; cb = 2*iv;
            X_mundlak(:,cw) = X_within(:,iv);
            X_mundlak(:,cb) = X_between(:,iv);
            group_idx(cw) = iv; group_idx(cb) = iv;
            col_names{cw} = [col_names_raw{iv} '_w'];
            col_names{cb} = [col_names_raw{iv} '_b'];
        end

        ok2 = ~any(isnan(X_mundlak)|isinf(X_mundlak), 2);
        X_mundlak = X_mundlak(ok2,:);
        Y_cv = Yok(ok2); S_cv = Sok(ok2);
        N_final = sum(ok2);

        if N_final < 20
            fprintf('  [%s | %s] : Pas assez de données (%d obs), skip\n', bn, gname, N_final);
            continue;
        end

        % --- Partition Group K-Fold par sujet (zéro leakage) ---
        fold_id = groupKFoldBySubject(S_cv, uSubj, K_FOLDS_loc);

        % --- Calcul analytique de Lambda Max ---
        mu_all = mean(X_mundlak,1);
        sd_all = std(X_mundlak,0,1); sd_all(sd_all==0) = 1;
        Xs_all = (X_mundlak - mu_all) ./ sd_all;
        Yc_all = Y_cv - mean(Y_cv);
        lambda_max = computeLambdaMax(Xs_all, Yc_all, group_idx, N_vars);
        lambda_min = 1e-4 * lambda_max;
        lambdas = exp(linspace(log(lambda_max), log(lambda_min), N_LAMBDA));

        % --- Calcul de la constante de Lipschitz pour l'ajustement final ---
        L_f = norm(Xs_all' * Xs_all, 'fro') / N_final;
        if L_f < 1e-10, L_f = 1; end

        % --- Cross-Validation Group K-Fold ---
        mse_folds = NaN(K_FOLDS_loc, N_LAMBDA);
        for k = 1:K_FOLDS_loc
            idx_tr = fold_id ~= k;
            idx_val = fold_id == k;
            X_tr = X_mundlak(idx_tr,:); Y_tr = Y_cv(idx_tr);
            X_val= X_mundlak(idx_val,:); Y_val= Y_cv(idx_val);

            % Standardisation interne locale pour chaque pli (sans contamination)
            mu_tr = mean(X_tr,1); sd_tr = std(X_tr,0,1); sd_tr(sd_tr==0)=1;
            Xs_tr = (X_tr - mu_tr) ./ sd_tr;
            Xs_val = (X_val - mu_tr) ./ sd_tr;
            Yc_tr = Y_tr - mean(Y_tr);
            Y_mean_tr = mean(Y_tr);

            % Calcul de la constante de Lipschitz pré-boucle pour accélérer FISTA
            L_cv = norm(Xs_tr' * Xs_tr, 'fro') / sum(idx_tr);
            if L_cv < 1e-10, L_cv = 1; end

            beta_warm = zeros(N_cols,1);
            for il = 1:N_LAMBDA
                beta_warm = fitGroupLassoFISTA(Xs_tr, Yc_tr, group_idx, N_vars, ...
                    lambdas(il), MAX_ITER, TOL_FISTA, beta_warm, L_cv);
                Y_pred = Xs_val * beta_warm + Y_mean_tr;
                mse_folds(k,il) = mean((Y_val - Y_pred).^2);
            end
        end

        % --- Sélection de Lambda par règle SE adaptative
        % 1SE standard pour K>=5 (SE stable)
        % 2SE pour K<=3 (petit échantillon, SE instable et gonflé)
        mse_mean = mean(mse_folds, 1, 'omitnan');
        mse_se   = std(mse_folds, 0, 1) / sqrt(K_FOLDS_loc);
        [mse_min_val, idx_min] = min(mse_mean);

        if K_FOLDS_loc <= 3
            se_mult = 2.0;   % 2SE pour petits échantillons
        else
            se_mult = 1.0;   % 1SE standard
        end

        threshold_1se = mse_min_val + se_mult * mse_se(idx_min);
        idx_1se = find(mse_mean <= threshold_1se, 1, 'first');
        if isempty(idx_1se), idx_1se = idx_min; end
        lambda_opt = lambdas(idx_1se);
        fprintf('    CV => lambda_min=%.4f (idx=%d) | lambda_%dSE=%.4f (idx=%d)\n',...
            lambdas(idx_min), idx_min, se_mult, lambda_opt, idx_1se);

        % --- Ajustement final ---
        mu_f = mean(X_mundlak,1); sd_f = std(X_mundlak,0,1); sd_f(sd_f==0)=1;
        Xs_f = (X_mundlak - mu_f) ./ sd_f;
        Yc_f = Y_cv - mean(Y_cv);
        beta_final = fitGroupLassoFISTA(Xs_f, Yc_f, group_idx, N_vars, ...
            lambda_opt, MAX_ITER, TOL_FISTA, zeros(N_cols,1), L_f);

        groups_selected = false(1, N_vars);
        for iv = 1:N_vars
            cols_iv = find(group_idx == iv);
            if any(abs(beta_final(cols_iv)) > 1e-10)
                groups_selected(iv) = true;
            end
        end
        n_sel = sum(groups_selected);
        fprintf('  [%s | %s] : %d/%d variables sélectionnées par Group LASSO (lambda = %.4f)\n', ...
            bn, gname, n_sel, N_vars, lambda_opt);

        if n_sel == 0, continue; end

        % --- Ajustement Final fitlm + HAC Newey-West ---
        sel_idx_vars = find(groups_selected);
        sel_col_w = cell(1,n_sel); sel_col_b = cell(1,n_sel);
        X_sel = NaN(N_final, 2*n_sel);
        for jj = 1:n_sel
            iv = sel_idx_vars(jj);
            cols_iv = find(group_idx == iv);
            X_sel(:,2*jj-1) = X_mundlak(:,cols_iv(1));
            X_sel(:,2*jj)   = X_mundlak(:,cols_iv(2));
            sel_col_w{jj} = matlab.lang.makeValidName(col_names{cols_iv(1)});
            sel_col_b{jj} = matlab.lang.makeValidName(col_names{cols_iv(2)});
        end
        all_pred_names = [sel_col_w sel_col_b];

        ok3 = ~any(isnan(X_sel),2) & ~isnan(Y_cv);
        X_sel = X_sel(ok3,:); Y_hac = Y_cv(ok3); N_hac = sum(ok3);

        % --- Option A : Filtrage de colinéarité avant fitlm ---
        % Supprimer les variables avec |r| > 0.90 (garder la plus corrélée au RPE)
        if size(X_sel,2) > 2
            R_mat = corr(X_sel, 'Rows','complete');
            R_mat = abs(R_mat);
            keep_mask = true(1, size(X_sel,2));
            for ic = 1:size(R_mat,2)-1
                if ~keep_mask(ic), continue; end
                for jc = ic+1:size(R_mat,2)
                    if keep_mask(jc) && R_mat(ic,jc) > 0.90
                        % Garder la variable la plus corrélée au RPE
                        r_ic = abs(corr(X_sel(:,ic), Y_hac, 'Rows','complete'));
                        r_jc = abs(corr(X_sel(:,jc), Y_hac, 'Rows','complete'));
                        if r_ic >= r_jc
                            keep_mask(jc) = false;
                        else
                            keep_mask(ic) = false;
                            break
                        end
                    end
                end
            end
            n_removed = sum(~keep_mask);
            if n_removed > 0
                fprintf('    Filtrage colinéarité : %d colonnes supprimées (|r|>0.90)\n', n_removed);
                X_sel = X_sel(:, keep_mask);
                all_pred_names = all_pred_names(keep_mask);
            end
        end

        tbl_hac = array2table([Y_hac X_sel], 'VariableNames', ['RPE' all_pred_names]);
        formula_hac = ['RPE ~ ' strjoin(all_pred_names, ' + ')];

        try
            mdl = fitlm(tbl_hac, formula_hac);
            try
                [~, SE_hac, ~] = hac(mdl, 'Bandwidth', HAC_BW, 'Weights','BT','Display','off');
            catch ME_hac
                fprintf('    HAC échoué (%s) — SE classiques utilisés\n', ME_hac.message);
                SE_hac = mdl.Coefficients.SE;
            end
            coef_vals = mdl.Coefficients.Estimate;
            coef_names = mdl.CoefficientNames';
            df = N_hac - length(coef_vals);
            t_hac = coef_vals ./ SE_hac;
            p_hac = 2 * tcdf(-abs(t_hac), df);

            res = struct();
            res.var_names_raw   = col_names_raw(sel_idx_vars);
            res.col_w           = sel_col_w;
            res.col_b           = sel_col_b;
            res.coef_names      = coef_names;
            res.coef            = coef_vals;
            res.SE_hac          = SE_hac;
            res.t_hac           = t_hac;
            res.p_hac           = p_hac;
            res.R2_adj          = mdl.Rsquared.Adjusted;
            res.lambda_opt      = lambda_opt;
            res.lambda_max      = lambda_max;
            res.mse_cv_mean     = mse_mean;
            res.mse_cv_se       = mse_se;
            res.lambdas         = lambdas;
            res.idx_1se         = idx_1se;
            res.idx_min         = idx_min;
            res.groups_selected = groups_selected;
            res.n_vars_parent   = N_vars;
            res.n_sel           = n_sel;

            if iGrp==1, mundlak_G1.(bn) = res;
            else, mundlak_G2.(bn) = res; end
        catch ME
            fprintf('  [ERREUR] %s %s : %s\n', bn, gname, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
% FIGURES GROUP LASSO MUNDLAK
% -----------------------------------------------------------------------
fprintf('\n=== Étape 3c : Tracé des figures Group LASSO Mundlak ===\n');
col_G1c = [0.85 0.15 0.15];
col_G2c = [0.15 0.45 0.85];

for ib = 1:4
    bn = bloc_names{ib};
    has_G1 = isfield(mundlak_G1, bn);
    has_G2 = isfield(mundlak_G2, bn);
    if ~has_G1 && ~has_G2, continue; end

    % 1. Courbes d'optimisation de Lambda (MSE de validation)
    fig_cv = figure('Name',sprintf('CV_Path_%s',bn),'NumberTitle','off',...
        'Position',[50 50 900 400],'Color','white');
    hold on;
    if has_G1
        r = mundlak_G1.(bn);
        shadedErrorBar_simple(log(r.lambdas), r.mse_cv_mean, r.mse_cv_se, col_G1c);
        xline(log(r.lambdas(r.idx_1se)),'--','Color',col_G1c,'LineWidth',1.5,...
            'Label',sprintf('G1 1SE \\lambda=%.4f',r.lambda_opt));
    end
    if has_G2
        r = mundlak_G2.(bn);
        shadedErrorBar_simple(log(r.lambdas), r.mse_cv_mean, r.mse_cv_se, col_G2c);
        xline(log(r.lambdas(r.idx_1se)),'--','Color',col_G2c,'LineWidth',1.5,...
            'Label',sprintf('G2 1SE \\lambda=%.4f',r.lambda_opt));
    end
    xlabel('log(\lambda)'); ylabel('MSE moyen (Group K-Fold)');
    title(sprintf('Chemin de régularisation Group LASSO - %s',bn),'FontWeight','bold');
    h1=patch(NaN,NaN,col_G1c,'FaceAlpha',0.4,'EdgeColor','none');
    h2=patch(NaN,NaN,col_G2c,'FaceAlpha',0.4,'EdgeColor','none');
    legend([h1 h2],{'G1 Short Duration','G2 Long Duration'},'Location','best');
    grid on; box on; hold off;
    saveas(fig_cv, fullfile(PathSave, sprintf('%s_GroupLASSO_CVPath_%s.png',lag_mode,bn)));
    close(fig_cv);

    % 2. Histogrammes comparatifs des bêtas (Within vs Between)
    all_vn_raw = {};
    if has_G1, all_vn_raw = [all_vn_raw; mundlak_G1.(bn).var_names_raw(:)]; end
    if has_G2, all_vn_raw = [all_vn_raw; mundlak_G2.(bn).var_names_raw(:)]; end
    all_vn_raw = unique(all_vn_raw,'stable');
    N_v = length(all_vn_raw);
    if N_v == 0, continue; end

    fig_beta = figure('Name',sprintf('Beta_Mundlak_%s',bn),'NumberTitle','off',...
        'Position',[50 50 max(900,N_v*100+200) 700],'Color','white');
    tiledlayout(2,1,'Padding','compact','TileSpacing','compact');
    subtitles = {'Effet intra-individuel (Within \beta_w) — Dynamique de fatigue',...
                 'Effet inter-individuel (Between \beta_b) — Niveau technique moyen'};
    suffixes = {'_w','_b'};

    for iSub = 1:2
        nexttile; hold on; offset = 0.2;
        for i = 1:N_v
            vn_raw = all_vn_raw{i};
            cn_suf = [matlab.lang.makeValidName(vn_raw) suffixes{iSub}];
            for iGrp = 1:2
                if iGrp==1 && ~has_G1, continue; end
                if iGrp==2 && ~has_G2, continue; end

                if iGrp==1
                    res=mundlak_G1.(bn); col_g=col_G1c; off=-offset;
                else
                    res=mundlak_G2.(bn); col_g=col_G2c; off= offset;
                end

                idx = find(strcmp(res.coef_names, cn_suf));
                if isempty(idx), continue; end
                b_=res.coef(idx); se_=res.SE_hac(idx); p_=res.p_hac(idx);

                bar(i+off, b_, 0.35, 'FaceColor', col_g, 'EdgeColor', 'none', 'FaceAlpha',0.8);
                errorbar(i+off, b_, se_, se_, 'k', 'LineWidth', 1.1, 'CapSize', 4);

                sig = p2star(p_);
                if ~isempty(sig)
                    text(i+off, b_+sign(b_)*(se_+0.05), sig, 'HorizontalAlignment','center',...
                        'FontSize',8, 'Color',col_g, 'FontWeight','bold');
                end
            end
        end
        yline(0,'-k','LineWidth',0.8);
        xticks(1:N_v); xticklabels(cellfun(@normLabel, all_vn_raw, 'UniformOutput', false));
        xtickangle(35); ylabel('\beta (± SE HAC NW)');
        title(subtitles{iSub},'FontWeight','bold');
        grid on; box on; hold off; xlim([0.5 N_v+0.5]);
    end
    saveas(fig_beta, fullfile(PathSave, sprintf('%s_GroupLASSO_Mundlak_Beta_%s.png',lag_mode,bn)));
    close(fig_beta);
end

%% -----------------------------------------------------------------------
% ÉTAPE 4 : MODÈLES MIXTES LINÉAIRES (LMM) PAR BLOC
% -----------------------------------------------------------------------
fprintf('\n=== Étape 4 : Ajustement des modèles mixtes (LMM) ===\n');
lmm_G1 = struct(); lmm_G2 = struct();

for ib = 1:4
    bn = bloc_names{ib};
    vnames = vnames_bloc.(bn);
    if isempty(vnames), continue; end

    for iGrp = 1:2
        if iGrp==1
            Xw=Xw_G1_bloc.(bn); Xb=Xb_G1_bloc.(bn); Y=RPE_G1_bloc.(bn); S=subj_G1_bloc.(bn); gname='G1';
            has_sel = isfield(mundlak_G1, bn);
            if has_sel, sel_v = mundlak_G1.(bn).var_names_raw; else, sel_v = {}; end
        else
            Xw=Xw_G2_bloc.(bn); Xb=Xb_G2_bloc.(bn); Y=RPE_G2_bloc.(bn); S=subj_G2_bloc.(bn); gname='G2';
            has_sel = isfield(mundlak_G2, bn);
            if has_sel, sel_v = mundlak_G2.(bn).var_names_raw; else, sel_v = {}; end
        end

        if isempty(sel_v)
            fprintf('  LMM %s %s : Aucune variable sélectionnée par Group LASSO\n', bn, gname);
            continue;
        end

        % Remapping des indices dans vnames_bloc par nom nettoyé
        vnames_clean = cellfun(@(x) regexprep(x,'[^a-zA-Z0-9]','_'), vnames, 'UniformOutput',false);
        sel_idx = find(ismember(vnames_clean, sel_v));

        if isempty(sel_idx)
            fprintf('  LMM %s %s : correspondance noms introuvable\n', bn, gname);
            continue;
        end

        % Reconstruire les colonnes Within (DeltaX) ET Between (X niveau
        % moyen) pour chaque variable sélectionnée par le groupe — même
        % logique que l'étape 3b (sel_col_w + sel_col_b), pas le X brut.
        Xw_sel = Xw(:, sel_idx);
        Xb_sel = Xb(:, sel_idx);
        Xsel   = [Xw_sel Xb_sel];
        ok = ~any(isnan(Xsel)|isinf(Xsel),2) & ~isnan(Y) & ~isinf(Y);
        Xsel = Xsel(ok,:); Yok = Y(ok); Sok = S(ok);
        mu_s = mean(Xsel,1); sd_s = std(Xsel,0,1); sd_s(sd_s==0) = 1;
        Xstd = (Xsel - mu_s) ./ sd_s;

        vnames_sel_clean = cellfun(@(x) regexprep(x,'[^a-zA-Z0-9]','_'), ...
            vnames(sel_idx), 'UniformOutput', false);
        col_names_lmm = [cellfun(@(x) [x '_w'], vnames_sel_clean, 'UniformOutput',false), ...
                          cellfun(@(x) [x '_b'], vnames_sel_clean, 'UniformOutput',false)];
        % DeltaX_within(t) capture la dynamique de changement récente ;
        % X_between capture le niveau moyen du sujet sur la session
        % (confondant potentiel de type RPE_début) — les deux entrent
        % comme prédicteurs séparés du modèle de changement.

        T = array2table([Yok Xstd double(Sok)], ...
            'VariableNames', ['RPE' col_names_lmm {'Participant'}]);
        T.Participant = categorical(T.Participant);

        % Vecteur temps par sujet pour la structure AR(1)
        uSok = unique(Sok);
        temps = zeros(length(Yok),1);
        for iS = 1:length(uSok)
            idx_s = find(Sok == uSok(iS));
            temps(idx_s) = 1:length(idx_s);
        end
        T.temps = temps;

        try
            % ---------------------------------------------------------------
            % LMM STANDARD via fitlme
            % Structure : RPE ~ variables + (1|Participant)
            % Autocorrélation résiduelle documentée via phi AR(1)
            % mais non corrigée structurellement (limite documentée)
            % Homoscédasticité vérifiée via Breusch-Pagan simplifié
            % ---------------------------------------------------------------
            formula = ['RPE ~ ' strjoin(col_names_lmm,' + ') ' + (1|Participant)'];
            lme_final = fitlme(T, formula);

            [~,~,FEStats] = fixedEffects(lme_final, 'DFMethod','satterthwaite');
            [R2m_full, R2c_full] = computeNakagawaR2(lme_final, T);

            % --- Phi AR(1) résiduel (diagnostic) ---
            resid_final = residuals(lme_final, 'ResidualType', 'Raw');
            phi_post_vec = NaN(length(uSok), 1);
            for iS = 1:length(uSok)
                idx_s = find(Sok == uSok(iS));
                r_s   = resid_final(idx_s);
                if numel(r_s) >= 2
                    phi_post_vec(iS) = corr(r_s(1:end-1), r_s(2:end), 'Rows','complete');
                end
            end
            phi_ar1 = mean(phi_post_vec, 'omitnan');

            % --- Test Breusch-Pagan simplifié (homoscédasticité) ---
            resid_sq = resid_final .^ 2;
            try
                [bp_rho, bp_p] = corr(resid_sq, temps, 'Type', 'Spearman');
            catch
                bp_rho = NaN; bp_p = NaN;
            end

            % --- Normalité des résidus ---
            try
                [sw_h, sw_p] = swtest(resid_final);
            catch
                try
                    [sw_h, sw_p] = lillietest(resid_final);
                catch
                    sw_h = NaN; sw_p = NaN;
                end
            end

            % --- Coefficients standardisés ---
            sd_Y     = std(Yok, 0, 'omitnan');
            beta_std = FEStats.Estimate(2:end) / sd_Y;

            res = struct();
            res.coef                = FEStats.Estimate(2:end);
            res.coef_std            = beta_std;
            res.ci_lo               = FEStats.Lower(2:end);
            res.ci_hi               = FEStats.Upper(2:end);
            res.pval                = FEStats.pValue(2:end);
            res.names               = col_names_lmm;
            res.sd_X                = ones(length(col_names_lmm),1);
            res.sd_Y                = sd_Y;
            res.R2_adj              = lme_final.Rsquared.Adjusted;
            res.R2_marginal_full    = R2m_full;
            res.R2_conditional_full = R2c_full;
            res.lme                 = lme_final;
            res.phi_ar1             = phi_ar1;
            res.bp_rho              = bp_rho;
            res.bp_p                = bp_p;
            res.sw_p                = sw_p;
            res.sw_h                = sw_h;

            if numel(resid_final) > 1
                res.r_lag1 = corr(resid_final(1:end-1), resid_final(2:end), 'Rows','complete');
            else
                res.r_lag1 = NaN;
            end

            % Contribution unique Drop-One dR2
            dR2_drop1 = NaN(length(col_names_lmm),1);
            for iVar = 1:length(col_names_lmm)
                col_red = col_names_lmm([1:iVar-1, iVar+1:end]);
                try
                    if isempty(col_red)
                        lme_red = fitlme(T, 'RPE ~ 1 + (1|Participant)');
                    else
                        lme_red = fitlme(T, ['RPE ~ ' strjoin(col_red,' + ') ' + (1|Participant)']);
                    end
                    dR2_drop1(iVar) = res.R2_adj - lme_red.Rsquared.Adjusted;
                catch; end
            end
            res.dR2_drop1 = dR2_drop1;

            if iGrp==1, lmm_G1.(bn)=res; else, lmm_G2.(bn)=res; end

            sw_str = 'N/A';
            if ~isnan(sw_p), sw_str = sprintf('p=%.3f%s', sw_p, p2star(sw_p)); end
            bp_str = 'N/A';
            if ~isnan(bp_p), bp_str = sprintf('p=%.3f%s (rho=%.2f)', bp_p, p2star(bp_p), bp_rho); end
            fprintf('  LMM [%s | %s] : R2_adj=%.3f | R2m=%.3f | R2c=%.3f | phi_AR1=%.2f | BP %s | SW %s\n', ...
                bn, gname, res.R2_adj, res.R2_marginal_full, res.R2_conditional_full, ...
                phi_ar1, bp_str, sw_str);
            fprintf('    %-35s  beta_raw   beta_std   p\n', 'Variable');
            for iv = 1:length(col_names_lmm)
                fprintf('    %-35s  %8.4f   %8.4f   %s\n', ...
                    col_names_lmm{iv}, res.coef(iv), res.coef_std(iv), p2star(res.pval(iv)));
            end

        catch ME
            fprintf('  [LMM ERREUR] %s %s : %s\n', bn, gname, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
% ÉTAPE 5 : FIGURES COMPARATIVES DES MODÈLES MIXTES (LMM)
% -----------------------------------------------------------------------
fprintf('\n=== Étape 5 : Tracé des graphiques comparatifs LMM ===\n');

for ib = 1:4
    bn = bloc_names{ib};
    has_G1 = isfield(lmm_G1,bn); has_G2 = isfield(lmm_G2,bn);
    if ~has_G1 && ~has_G2, continue; end

    all_names = {};
    if has_G1, all_names = [all_names; lmm_G1.(bn).names(:)]; end
    if has_G2, all_names = [all_names; lmm_G2.(bn).names(:)]; end
    all_names = unique(all_names,'stable'); N_v = length(all_names);

    fig_name = sprintf('%s_Indiv_LMM_Bloc_%s_G1vsG2', lag_mode, bn);
    figure('Name',fig_name,'NumberTitle','off','Position',[50 50 max(800,N_v*80+200) 550],'Color','white');
    hold on; offset=0.2;
    for i = 1:N_v
        vn=all_names{i};
        if has_G1
            idx1=find(strcmp(lmm_G1.(bn).names,vn));
            if ~isempty(idx1)
                b1=lmm_G1.(bn).coef(idx1); lo1=lmm_G1.(bn).ci_lo(idx1);
                hi1=lmm_G1.(bn).ci_hi(idx1); p1=lmm_G1.(bn).pval(idx1);
                bar(i-offset,b1,0.35,'FaceColor',col_G1c,'EdgeColor','none','FaceAlpha',0.8);
                errorbar(i-offset,b1,b1-lo1,hi1-b1,'k','LineWidth',1.1,'CapSize',4);
                sig=p2star(p1); if ~isempty(sig), text(i-offset,hi1+0.05,sig,'HorizontalAlignment','center','FontSize',8,'Color',col_G1c,'FontWeight','bold'); end
            end
        end
        if has_G2
            idx2=find(strcmp(lmm_G2.(bn).names,vn));
            if ~isempty(idx2)
                b2=lmm_G2.(bn).coef(idx2); lo2=lmm_G2.(bn).ci_lo(idx2);
                hi2=lmm_G2.(bn).ci_hi(idx2); p2=lmm_G2.(bn).pval(idx2);
                bar(i+offset,b2,0.35,'FaceColor',col_G2c,'EdgeColor','none','FaceAlpha',0.8);
                errorbar(i+offset,b2,b2-lo2,hi2-b2,'k','LineWidth',1.1,'CapSize',4);
                sig=p2star(p2); if ~isempty(sig), text(i+offset,hi2+0.05,sig,'HorizontalAlignment','center','FontSize',8,'Color',col_G2c,'FontWeight','bold'); end
            end
        end
    end
    yline(0,'-k','LineWidth',0.8);
    xticks(1:N_v); xticklabels(cellfun(@normLabel,all_names,'UniformOutput',false)); xtickangle(35);
    ylabel('Coefficient \beta standardisé (± IC 95%)');
    title(sprintf('Coefficients des effets fixes LMM - Bloc %s',bn),'FontWeight','bold');
    h1=patch(NaN,NaN,col_G1c,'FaceAlpha',0.8,'EdgeColor','none');
    h2=patch(NaN,NaN,col_G2c,'FaceAlpha',0.8,'EdgeColor','none');
    legend([h1 h2],{'G1 Short Duration','G2 Long Duration'},'Location','best');
    grid on; box on; hold off; xlim([0.5 N_v+0.5]);
    saveas(gcf, fullfile(PathSave,[fig_name '.png'])); close(gcf);
end

% --- Figure comparative du R² ajusté par Bloc ---
figure('Name','Indiv_R2Adj_ParBloc','NumberTitle','off','Position',[50 50 700 450],'Color','white');
hold on;
R2_G1_vec=NaN(4,1); R2_G2_vec=NaN(4,1);
for ib=1:4
    bn=bloc_names{ib};
    if isfield(lmm_G1,bn), R2_G1_vec(ib)=lmm_G1.(bn).R2_adj; end
    if isfield(lmm_G2,bn), R2_G2_vec(ib)=lmm_G2.(bn).R2_adj; end
end
x_coords = 1:4;
b1=bar(x_coords-0.2,R2_G1_vec,0.35,'FaceColor',col_G1c,'EdgeColor','none');
b2=bar(x_coords+0.2,R2_G2_vec,0.35,'FaceColor',col_G2c,'EdgeColor','none');
for ib=1:4
    if ~isnan(R2_G1_vec(ib)), text(ib-0.2,R2_G1_vec(ib)+0.01,sprintf('%.2f',R2_G1_vec(ib)),'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',col_G1c); end
    if ~isnan(R2_G2_vec(ib)), text(ib+0.2,R2_G2_vec(ib)+0.01,sprintf('%.2f',R2_G2_vec(ib)),'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',col_G2c); end
end
yline(0.65,'--k','Goubault max (0.65)','LineWidth',1.2);
xticks(1:4); xticklabels(bloc_labels);
ylabel('R² ajusté (Modèle complet)'); ylim([0 1]);
title('R² ajusté global par bloc de variables','FontWeight','bold');
legend([b1 b2],{'G1 Short Duration','G2 Long Duration'},'Location','northeast');
grid on; box on; hold off;
saveas(gcf, fullfile(PathSave,sprintf('%s_Indiv_LMM_R2Adj_ParBloc.png',lag_mode))); close(gcf);

%% -----------------------------------------------------------------------
% ÉTAPE 6 : BOUCLE DE STABILITÉ GROUP LASSO (10 GRAINES — réduit depuis 100
%           pour ce premier passage exploratoire à 4 variantes de fenêtre)
% -----------------------------------------------------------------------
fprintf('\n=== Étape 6 : Analyse de stabilité sur 10 graines [%s] ===\n', lag_lbl);
seeds_stab = 1:10;
N_stab = length(seeds_stab);
thresh_pct = 50;

stab_G1 = struct('SegMod',struct(),'SegAxis',struct(),'Goubault',struct(),'EMG',struct());
stab_G2 = struct('SegMod',struct(),'SegAxis',struct(),'Goubault',struct(),'EMG',struct());

% --- sd_Y global par bloc et groupe (Point 3 : normalisation constante) ---
% Calculé une seule fois avant la boucle pour que les beta_std soient
% comparables entre tous les runs — indépendant de la composition des folds.
sd_Y_global = struct();
for ib = 1:4
    bn = bloc_names{ib};
    sd_Y_global.(bn).G1 = std(RPE_G1_bloc.(bn)(~isnan(RPE_G1_bloc.(bn))), 0, 'omitnan');
    sd_Y_global.(bn).G2 = std(RPE_G2_bloc.(bn)(~isnan(RPE_G2_bloc.(bn))), 0, 'omitnan');
end

for iSeed = 1:N_stab
    seed_val = seeds_stab(iSeed);
    fprintf('  -> Exécution de la simulation : Run %d/%d (Graine active : %d) [%s]...\n', iSeed, N_stab, seed_val, lag_lbl);

    for ib = 1:4
        bn = bloc_names{ib};
        vnames = vnames_bloc.(bn);
        if isempty(vnames), continue; end
        N_vars = length(vnames);

        for iGrp = 1:2
            if iGrp==1
                Xw=Xw_G1_bloc.(bn); Xb=Xb_G1_bloc.(bn); Y=RPE_G1_bloc.(bn); S=subj_G1_bloc.(bn);
            else
                Xw=Xw_G2_bloc.(bn); Xb=Xb_G2_bloc.(bn); Y=RPE_G2_bloc.(bn); S=subj_G2_bloc.(bn);
            end

            ok = ~any(isnan(Xw)|isinf(Xw),2) & ~any(isnan(Xb)|isinf(Xb),2) & ~isnan(Y) & ~isinf(Y);
            if sum(ok) < 20, continue; end
            Xw_ok=Xw(ok,:); Xb_ok=Xb(ok,:); Yok=Y(ok); Sok=S(ok); N_obs=length(Yok);

            % Mundlak déjà décomposé en amont (DeltaX_within + X_between)
            uSubj = unique(Sok); N_subj = length(uSubj);
            X_within  = Xw_ok;
            X_between = Xb_ok;

            N_cols_s = 2*N_vars;
            X_mnd = NaN(N_obs, N_cols_s);
            grp_s = NaN(1, N_cols_s);
            for iv = 1:N_vars
                X_mnd(:, 2*iv-1) = X_within(:, iv);
                X_mnd(:, 2*iv)   = X_between(:, iv);
                grp_s(2*iv-1) = iv; grp_s(2*iv) = iv;
            end

            ok2 = ~any(isnan(X_mnd)|isinf(X_mnd),2);
            X_mnd=X_mnd(ok2,:); Y_mnd=Yok(ok2); S_mnd=Sok(ok2);
            if length(Y_mnd) < 20, continue; end

            % K adaptatif dans la boucle de stabilité
            if N_subj <= 15
                K_FOLDS_stab = N_subj;
            elseif N_subj <= 20
                K_FOLDS_stab = 3;
            else
                K_FOLDS_stab = K_FOLDS;
            end

            % CORRECTION 2 : rng placé AVANT groupKFoldBySubject
            % pour garantir la reproductibilité de la permutation par graine
            rng(seed_val);
            fold_id_s = groupKFoldBySubject(S_mnd, uSubj, K_FOLDS_stab);

            % Grille et lambdas locaux
            mu_s2=mean(X_mnd,1); sd_s2=std(X_mnd,0,1); sd_s2(sd_s2==0)=1;
            Xs_s2=(X_mnd-mu_s2)./sd_s2; Yc_s2=Y_mnd-mean(Y_mnd);
            lmax_s = computeLambdaMax(Xs_s2, Yc_s2, grp_s, N_vars);
            lmin_s = 1e-4 * lmax_s;
            lams_s = exp(linspace(log(lmax_s), log(lmin_s), N_LAMBDA));

            % Constante de Lipschitz finale globale pré-calculée pour cette graine
            L_glob_s = norm(Xs_s2' * Xs_s2, 'fro') / length(Y_mnd);
            if L_glob_s < 1e-10, L_glob_s = 1; end

            try
                % CV Group K-Fold
                mse_s = NaN(K_FOLDS_stab, N_LAMBDA);
                for k = 1:K_FOLDS_stab
                    idx_tr=fold_id_s~=k; idx_val=fold_id_s==k;
                    X_tr=X_mnd(idx_tr,:); Y_tr=Y_mnd(idx_tr);
                    X_vl=X_mnd(idx_val,:); Y_vl=Y_mnd(idx_val);

                    mu_t=mean(X_tr,1); sd_t=std(X_tr,0,1); sd_t(sd_t==0)=1;
                    Xs_t=(X_tr-mu_t)./sd_t; Xs_v=(X_vl-mu_t)./sd_t;
                    Yc_t=Y_tr-mean(Y_tr); Ym_t=mean(Y_tr);

                    % Constante de Lipschitz locale pré-calculée pour ce pli
                    L_cv_s = norm(Xs_t' * Xs_t, 'fro') / sum(idx_tr);
                    if L_cv_s < 1e-10, L_cv_s = 1; end

                    bw=zeros(N_cols_s,1);
                    for il=1:N_LAMBDA
                        bw=fitGroupLassoFISTA(Xs_t,Yc_t,grp_s,N_vars,lams_s(il),200,1e-5,bw,L_cv_s);
                        mse_s(k,il)=mean((Y_vl - Xs_v*bw - Ym_t).^2);
                    end
                end

                mse_m=mean(mse_s,1,'omitnan'); mse_e=std(mse_s,0,1)/sqrt(K_FOLDS_stab);
                [mm,im]=min(mse_m);
                % Règle SE adaptative : 2SE si K<=3, 1SE sinon
                se_mult_stab = 1.0;
                if K_FOLDS_stab <= 3, se_mult_stab = 2.0; end
                thr=mm + se_mult_stab * mse_e(im);
                i1se=find(mse_m<=thr,1,'first');
                if isempty(i1se), i1se=im; end
                lopt_s=lams_s(i1se);

                % Ajustement final de la graine
                bf=fitGroupLassoFISTA(Xs_s2,Yc_s2,grp_s,N_vars,lopt_s,300,1e-5,zeros(N_cols_s,1),L_glob_s);

                groups_sel_s = false(1,N_vars);
                for iv=1:N_vars
                    cols_iv=find(grp_s==iv);
                    if any(abs(bf(cols_iv))>1e-10), groups_sel_s(iv)=true; end
                end
                idx_sel_s = find(groups_sel_s);
                if isempty(idx_sel_s), continue; end

                % Calcul LMM unifié pour stocker les métriques associées à la sélection
                % Variable cible : DeltaRPE (différence lag-k), prédicteurs DeltaX_within + X_between
                vnames_clean_s = cellfun(@(x) regexprep(x,'[^a-zA-Z0-9]','_'),vnames,'UniformOutput',false);
                % Reconstruire Xsel depuis Xw_ok (DeltaX_within) et Xb_ok (X_between)
                % pour chaque variable parente sélectionnée par le groupe — colonnes _w et _b
                Xw_sel_s = Xw_ok(:, idx_sel_s);
                Xb_sel_s = Xb_ok(:, idx_sel_s);
                Xsel_s   = [Xw_sel_s Xb_sel_s];
                ok_lmm = ~any(isnan(Xsel_s)|isinf(Xsel_s),2) & ~isnan(Yok) & ~isinf(Yok);
                Xsel_s=Xsel_s(ok_lmm,:); Ylmm=Yok(ok_lmm); Slmm=Sok(ok_lmm);

                mu_l=mean(Xsel_s,1); sd_l=std(Xsel_s,0,1); sd_l(sd_l==0)=1;
                Xs_l=(Xsel_s-mu_l)./sd_l;
                % Noms : [var1_w, var2_w, ..., var1_b, var2_b, ...]
                col_run = [cellfun(@(x) [x '_w'], vnames_clean_s(idx_sel_s), 'UniformOutput',false), ...
                           cellfun(@(x) [x '_b'], vnames_clean_s(idx_sel_s), 'UniformOutput',false)];
                % Normalisation constante : sd_Y calculé sur l'ensemble du bloc
                % (pas sur le sous-ensemble du run) pour comparabilité entre runs
                if iGrp==1, sd_Y_cst = sd_Y_global.(bn).G1;
                else,        sd_Y_cst = sd_Y_global.(bn).G2; end

                T_run=array2table([Ylmm Xs_l double(Slmm)],...
                    'VariableNames',['RPE' col_run {'Participant'}]);
                T_run.Participant=categorical(T_run.Participant);
                lme_run=fitlme(T_run,['RPE ~ ' strjoin(col_run,' + ') ' + (1|Participant)']);
                [R2m_run,R2c_run]=computeNakagawaR2(lme_run,T_run);
                R2full_run=lme_run.Rsquared.Adjusted;

                % Extraire coefficients et p-values du LMM
                [~,~,FE_run] = fixedEffects(lme_run,'DFMethod','satterthwaite');
                beta_run_raw = FE_run.Estimate(2:end);   % excl. intercept
                pval_run     = FE_run.pValue(2:end);
                beta_run_std = beta_run_raw / sd_Y_cst;  % standardisés (sd_Y constant)

                dR2_run=NaN(length(idx_sel_s),1);
                for iVar=1:length(idx_sel_s)
                    col_red=col_run([1:iVar-1,iVar+1:end]);
                    try
                        if isempty(col_red)
                            lme_red=fitlme(T_run,'RPE ~ 1 + (1|Participant)');
                        else
                            lme_red=fitlme(T_run,['RPE ~ ' strjoin(col_red,' + ') ' + (1|Participant)']);
                        end
                        dR2_run(iVar)=R2full_run-lme_red.Rsquared.Adjusted;
                    catch; end
                end

                % Stockage unifié dans la structure de stabilité par variable parente
                for jj=1:length(idx_sel_s)
                    iv = idx_sel_s(jj);
                    vn = vnames{iv};
                    fn = matlab.lang.makeValidName(vn);
                    cols_iv = find(grp_s==iv);

                    % Trouver l'indice de la colonne _w de cette variable dans col_run
                    idx_col_run = find(strcmp(col_run, [vnames_clean_s{iv} '_w']));

                    if iGrp == 1
                        if ~isfield(stab_G1.(bn), fn)
                            stab_G1.(bn).(fn) = struct('coefs_w',[],'coefs_b',[],...
                                'beta_std',[],'pval',[],'dR2',[],'R2m_z',[],'R2c_z',[],'n_sel',0);
                        end
                        stab_G1.(bn).(fn).coefs_w(end+1) = bf(cols_iv(1));
                        stab_G1.(bn).(fn).coefs_b(end+1) = bf(cols_iv(2));
                        if ~isempty(idx_col_run)
                            stab_G1.(bn).(fn).beta_std(end+1) = beta_run_std(idx_col_run);
                            stab_G1.(bn).(fn).pval(end+1)     = pval_run(idx_col_run);
                        end
                        stab_G1.(bn).(fn).dR2(end+1) = dR2_run(jj);
                        if R2m_run>0 && R2m_run<1, stab_G1.(bn).(fn).R2m_z(end+1)=atanh(sqrt(R2m_run)); end
                        if R2c_run>0 && R2c_run<1, stab_G1.(bn).(fn).R2c_z(end+1)=atanh(sqrt(R2c_run)); end
                        stab_G1.(bn).(fn).n_sel = stab_G1.(bn).(fn).n_sel + 1;
                    else
                        if ~isfield(stab_G2.(bn), fn)
                            stab_G2.(bn).(fn) = struct('coefs_w',[],'coefs_b',[],...
                                'beta_std',[],'pval',[],'dR2',[],'R2m_z',[],'R2c_z',[],'n_sel',0);
                        end
                        stab_G2.(bn).(fn).coefs_w(end+1) = bf(cols_iv(1));
                        stab_G2.(bn).(fn).coefs_b(end+1) = bf(cols_iv(2));
                        if ~isempty(idx_col_run)
                            stab_G2.(bn).(fn).beta_std(end+1) = beta_run_std(idx_col_run);
                            stab_G2.(bn).(fn).pval(end+1)     = pval_run(idx_col_run);
                        end
                        stab_G2.(bn).(fn).dR2(end+1) = dR2_run(jj);
                        if R2m_run>0 && R2m_run<1, stab_G2.(bn).(fn).R2m_z(end+1)=atanh(sqrt(R2m_run)); end
                        if R2c_run>0 && R2c_run<1, stab_G2.(bn).(fn).R2c_z(end+1)=atanh(sqrt(R2c_run)); end
                        stab_G2.(bn).(fn).n_sel = stab_G2.(bn).(fn).n_sel + 1;
                    end
                end
            catch
                % Sécurité si une graine ou FISTA diverge localement
            end
        end
    end
end

% --- Génération des figures synthétiques de stabilité sous format PNG ---
% CORRECTION 1 : blocs_labels_stab manquait dans le document 4
blocs_labels_stab = {'BySegMod+Total','BySegAxis','Goubault','EMG'};

%% -----------------------------------------------------------------------
%  TABLEAU MAÎTRE DE STABILITÉ (seuil f > 50%)
%  -----------------------------------------------------------------------
THRESH_STAB = 50;  % seuil de fréquence de sélection (%)
fprintf('\n=== TABLEAU MAÎTRE DE STABILITÉ (seuil f > %d%%) ===\n', THRESH_STAB);

for ib = 1:4
    bn = bloc_names{ib};
    for iGrp = 1:2
        if iGrp==1, st=stab_G1.(bn); gname='G1';
        else,        st=stab_G2.(bn); gname='G2'; end
        fns = fieldnames(st);
        if isempty(fns), continue; end

        % Calculer métriques pour chaque variable
        rows = {};
        for iv = 1:length(fns)
            fn = fns{iv};
            s  = st.(fn);
            freq = s.n_sel / N_stab * 100;
            if freq < THRESH_STAB, continue; end

            % Beta_std : moyenne et IC95% sur les runs où sélectionnée
            if ~isempty(s.beta_std)
                mu_b  = mean(s.beta_std, 'omitnan');
                sd_b  = std(s.beta_std,  0, 'omitnan');
                n_b   = sum(~isnan(s.beta_std));
                ci_lo = mu_b - 1.96*sd_b/sqrt(n_b);
                ci_hi = mu_b + 1.96*sd_b/sqrt(n_b);
            else
                mu_b=NaN; ci_lo=NaN; ci_hi=NaN;
            end

            % p-value médiane
            if ~isempty(s.pval)
                pmed = median(s.pval, 'omitnan');
            else
                pmed = NaN;
            end

            rows{end+1} = {normLabel(fn), freq, mu_b, ci_lo, ci_hi, pmed};
        end

        if isempty(rows), continue; end

        % Trier par fréquence décroissante
        freqs = cellfun(@(r) r{2}, rows);
        [~, isort] = sort(freqs, 'descend');
        rows = rows(isort);

        fprintf('\n  [%s | %s] — %d variables stables\n', bn, gname, length(rows));
        fprintf('  %-40s  %8s  %9s  %18s  %10s\n', ...
            'Variable', 'Stab(%)', 'βstd(moy)', 'IC95%', 'p(med)');
        fprintf('  %s\n', repmat('-',1,92));
        for ir = 1:length(rows)
            r = rows{ir};
            p_str = '';
            if ~isnan(r{6})
                if r{6}<0.001, p_str='< 0.001';
                elseif r{6}<0.01, p_str=sprintf('%.3f**',r{6});
                elseif r{6}<0.05, p_str=sprintf('%.3f*',r{6});
                else, p_str=sprintf('%.3f',r{6}); end
            end
            fprintf('  %-40s  %7.1f%%  %9.4f  [%6.3f; %6.3f]  %10s\n', ...
                r{1}, r{2}, r{3}, r{4}, r{5}, p_str);
        end
    end
end

%% -----------------------------------------------------------------------
%  FOREST PLOT DE STABILITÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Génération des Forest Plots de stabilité ===\n');

% Variables nécessaires — redéfinies ici pour pouvoir lancer cette section seule
bloc_names        = {'SegMod','SegAxis','Goubault','EMG'};
blocs_labels_stab = {'BySegMod+Total','BySegAxis','Goubault','EMG'};
col_G1c           = [0.72 0.11 0.11];
col_G2c           = [0.15 0.32 0.60];
thresh_pct        = THRESH_STAB;
PathSave_fp = fileparts(mfilename('fullpath'));
if isempty(PathSave_fp), PathSave_fp = pwd; end
PathSave = PathSave_fp;

for ib = 1:4
    bn = bloc_names{ib};
    for iGrp = 1:2
        if iGrp==1, st=stab_G1.(bn); gname='G1'; col_g=col_G1c;
        else,        st=stab_G2.(bn); gname='G2'; col_g=col_G2c; end
        fns = fieldnames(st);
        if isempty(fns), continue; end

        % Collecter variables au-dessus du seuil
        vars_fp = {}; mu_fp = []; ci_lo_fp = []; ci_hi_fp = [];
        freq_fp = []; min_fp = []; max_fp = [];
        for iv = 1:length(fns)
            fn = fns{iv};
            s  = st.(fn);
            freq = s.n_sel / N_stab * 100;
            if freq < THRESH_STAB || isempty(s.beta_std), continue; end
            mu_b  = mean(s.beta_std,'omitnan');
            sd_b  = std(s.beta_std,0,'omitnan');
            n_b   = sum(~isnan(s.beta_std));
            vars_fp{end+1}  = normLabel(fn);
            mu_fp(end+1)    = mu_b;
            ci_lo_fp(end+1) = mu_b - 1.96*sd_b/sqrt(n_b);
            ci_hi_fp(end+1) = mu_b + 1.96*sd_b/sqrt(n_b);
            freq_fp(end+1)  = freq;
            min_fp(end+1)   = min(s.beta_std,[],'omitnan');
            max_fp(end+1)   = max(s.beta_std,[],'omitnan');
        end

        if isempty(vars_fp), continue; end

        % Trier par beta_std
        [~, isort] = sort(mu_fp, 'descend');
        vars_fp=vars_fp(isort); mu_fp=mu_fp(isort);
        ci_lo_fp=ci_lo_fp(isort); ci_hi_fp=ci_hi_fp(isort);
        freq_fp=freq_fp(isort);
        min_fp=min_fp(isort); max_fp=max_fp(isort);
        N_fp = length(vars_fp);

        fig_fp = figure('Name', sprintf('ForestPlot_%s_%s',bn,gname), ...
            'NumberTitle','off','Position',[50 50 700 max(300,N_fp*45+100)], ...
            'Color','white');
        hold on;

        for iv = 1:N_fp
            y = N_fp - iv + 1;
            alpha_val = 0.4 + 0.6*(freq_fp(iv)-THRESH_STAB)/(100-THRESH_STAB);
            col_alpha = col_g * alpha_val + (1-alpha_val)*[1 1 1];

            % Étendue min/max en pointillés fins
            plot([min_fp(iv) max_fp(iv)], [y y], ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0);
            % IC 95% en trait plein épais — couleur bloc, bien visible
            plot([ci_lo_fp(iv) ci_hi_fp(iv)], [y y], '-', 'Color', col_g, 'LineWidth', 3.5);
            % Marqueurs aux extrémités des IC
            plot([ci_lo_fp(iv) ci_hi_fp(iv)], [y y], '|', 'Color', col_g, 'LineWidth', 2.5, 'MarkerSize', 10);
            % Point central — blanc pour ressortir sur la barre
            scatter(mu_fp(iv), y, 80, col_g, 'filled');
            scatter(mu_fp(iv), y, 30, 'white', 'filled');
            % Fréquence en texte
            x_txt = max(ci_hi_fp(iv), max_fp(iv)) + 0.03;
            text(x_txt, y, sprintf('%.0f%%', freq_fp(iv)), ...
                'FontSize', 9, 'Color', col_g, 'VerticalAlignment','middle', 'FontWeight','bold');
        end

        xline(0, '--k', 'LineWidth', 0.8);
        yticks(1:N_fp);
        yticklabels(fliplr(vars_fp));
        xlabel('\beta_{std} (Moy ± IC 95%)','FontSize',10,'FontName','Times New Roman');
        title(sprintf('Forest Plot — %s | %s', bn, gname), 'FontWeight','bold','FontName','Times New Roman');
        set(gca,'FontName','Times New Roman','FontSize',9,'Box','off','TickDir','out','LineWidth',0.8);
        grid on;
        % Élargir xlim pour que les IC très étroits restent visibles
        ic_range = max(ci_hi_fp) - min(ci_lo_fp);
        margin = max(0.15, ic_range * 0.2);
        xlim([min(ci_lo_fp) - margin, max(ci_hi_fp) + margin + 0.12]);
        saveas(fig_fp, fullfile(PathSave, sprintf('%s_ForestPlot_%s_%s.png',lag_mode,bn,gname)));
        close(fig_fp);
        fprintf('  Forest plot sauvegardé : ForestPlot_%s_%s.png\n', bn, gname);
    end
end

for ib = 1:4
    bn=bloc_names{ib};
    st1=stab_G1.(bn); st2=stab_G2.(bn);
    fn1=fieldnames(st1); fn2=fieldnames(st2);
    all_fn=unique([fn1;fn2]); if isempty(all_fn), continue; end
    N_v=length(all_fn);

    freq1=zeros(N_v,1); freq2=zeros(N_v,1);
    med_w1=NaN(N_v,1); med_w2=NaN(N_v,1);
    med_b1=NaN(N_v,1); med_b2=NaN(N_v,1);
    iqr_w1=NaN(N_v,1); iqr_w2=NaN(N_v,1);
    dR2_1=NaN(N_v,1); dR2_2=NaN(N_v,1);
    R2m_1=NaN(N_v,1); R2m_2=NaN(N_v,1);
    R2c_1=NaN(N_v,1); R2c_2=NaN(N_v,1);
    vn_disp=cell(N_v,1);

    for iv=1:N_v
        fn=all_fn{iv}; vn_disp{iv}=normLabel(fn);
        if isfield(st1,fn)
            s=st1.(fn); freq1(iv)=s.n_sel/N_stab*100;
            if ~isempty(s.coefs_w), med_w1(iv)=median(s.coefs_w); iqr_w1(iv)=iqr(s.coefs_w); end
            if ~isempty(s.coefs_b), med_b1(iv)=median(s.coefs_b); end
            if ~isempty(s.dR2), dR2_1(iv)=mean(s.dR2,'omitnan'); end
            if ~isempty(s.R2m_z), R2m_1(iv)=tanh(mean(s.R2m_z,'omitnan'))^2; end
            if ~isempty(s.R2c_z), R2c_1(iv)=tanh(mean(s.R2c_z,'omitnan'))^2; end
        end
        if isfield(st2,fn)
            s=st2.(fn); freq2(iv)=s.n_sel/N_stab*100;
            if ~isempty(s.coefs_w), med_w2(iv)=median(s.coefs_w); iqr_w2(iv)=iqr(s.coefs_w); end
            if ~isempty(s.coefs_b), med_b2(iv)=median(s.coefs_b); end
            if ~isempty(s.dR2), dR2_2(iv)=mean(s.dR2,'omitnan'); end
            if ~isempty(s.R2m_z), R2m_2(iv)=tanh(mean(s.R2m_z,'omitnan'))^2; end
            if ~isempty(s.R2c_z), R2c_2(iv)=tanh(mean(s.R2c_z,'omitnan'))^2; end
        end
    end

    [~,si]=sort(max(freq1,freq2),'descend');
    freq1=freq1(si); freq2=freq2(si);
    med_w1=med_w1(si); med_w2=med_w2(si);
    med_b1=med_b1(si); med_b2=med_b2(si);
    iqr_w1=iqr_w1(si); iqr_w2=iqr_w2(si);
    dR2_1=dR2_1(si); dR2_2=dR2_2(si);
    R2m_1=R2m_1(si); R2m_2=R2m_2(si);
    R2c_1=R2c_1(si); R2c_2=R2c_2(si);
    vn_disp=vn_disp(si);

    % --- Génération de l'image de synthèse en tableau de bord ---
    col_hdrs={'Variable',...
        'Freq G1','betaW G1','betaB G1','IQRw G1','dR2 G1','R2m G1','R2c G1',...
        'Freq G2','betaW G2','betaB G2','IQRw G2','dR2 G2','R2m G2','R2c G2'};
    col_x=[0.01,0.17,0.24,0.30,0.36,0.42,0.48,0.54,...
        0.61,0.68,0.74,0.80,0.86,0.91,0.96];

    row_h=max(0.010,min(0.025,0.80/max(N_v,1)));
    fig_h=max(400,100+N_v*22); fig_w=1600;
    fig=figure('Visible','off','Position',[50 50 fig_w fig_h],'Color','white');
    ax=axes('Position',[0 0 1 1],'Visible','off'); ax.XLim=[0 1]; ax.YLim=[0 1]; hold on;

    rectangle('Position',[0 0.96 1 0.04],'FaceColor',[0.15 0.25 0.45],'EdgeColor','none');
    text(0.5,0.98,sprintf('Stabilité de sélection des biomarqueurs parentiels - %s (%d runs)',...
        blocs_labels_stab{ib},N_stab),...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10,'FontWeight','bold','Color','white');

    rectangle('Position',[col_x(2) 0.92 col_x(9)-col_x(2) 0.04],'FaceColor',[0.95 0.88 0.88],'EdgeColor','none');
    text(mean([col_x(2) col_x(8)+0.05]),0.94,'G1 - Short Duration (Sujets fatigables)','HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',[0.7 0.1 0.1]);

    rectangle('Position',[col_x(9) 0.92 1-col_x(9) 0.04],'FaceColor',[0.88 0.91 0.97],'EdgeColor','none');
    text(mean([col_x(9) 1.0]),0.94,'G2 - Long Duration (Sujets résistants)','HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',[0.1 0.2 0.7]);

    head_y=0.89;
    rectangle('Position',[0 head_y-0.005 1 0.030],'FaceColor',[0.75 0.78 0.85],'EdgeColor','none');
    for ic=1:15, text(col_x(ic),head_y+0.010,col_hdrs{ic},'FontSize',6.5,'FontWeight','bold','Color',[0.1 0.1 0.3],'VerticalAlignment','middle'); end
    line([0 1],[head_y-0.005 head_y-0.005],'Color',[0.5 0.5 0.6],'LineWidth',1);
    line([col_x(9) col_x(9)],[0 1],'Color',[0.5 0.5 0.6],'LineWidth',0.8);

    col_r=[0.7 0.1 0.1]; col_b=[0.1 0.2 0.7]; col_gr=[0.2 0.5 0.2]; col_gy=[0.4 0.4 0.4];
    for iv=1:N_v
        y_row=head_y-0.005-iv*row_h; y_txt=y_row+row_h*0.5;
        if mod(iv,2)==0, rectangle('Position',[0 y_row 1 row_h],'FaceColor',[0.97 0.97 0.99],'EdgeColor','none'); end
        if freq1(iv)>=thresh_pct, rectangle('Position',[col_x(2) y_row col_x(9)-col_x(2) row_h],'FaceColor',[1.0 0.93 0.88],'EdgeColor','none'); end
        if freq2(iv)>=thresh_pct, rectangle('Position',[col_x(9) y_row 1-col_x(9) row_h],'FaceColor',[0.88 0.93 1.0],'EdgeColor','none'); end
        text(col_x(1)+0.002,y_txt,vn_disp{iv},'FontSize',6,'VerticalAlignment','middle','Interpreter','none','Color',[0.1 0.1 0.1]);

        % Rendu G1
        if freq1(iv)>0
            cf=col_gy; if freq1(iv)>=thresh_pct,cf=col_r; end
            text(col_x(2)+0.03,y_txt,sprintf('%.0f%%',freq1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',cf,'FontWeight','bold');
        end
        if ~isnan(med_w1(iv)), text(col_x(3)+0.03,y_txt,sprintf('%.3f',med_w1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_r,'FontWeight','bold'); end
        if ~isnan(med_b1(iv)), text(col_x(4)+0.03,y_txt,sprintf('%.3f',med_b1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_r); end
        if ~isnan(iqr_w1(iv)), text(col_x(5)+0.03,y_txt,sprintf('%.3f',iqr_w1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_gy); end
        if ~isnan(dR2_1(iv))
            cd=col_gy; if dR2_1(iv)>0.03,cd=col_gr; end
            text(col_x(6)+0.03,y_txt,sprintf('%.4f',dR2_1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',cd,'FontWeight','bold');
        end
        if ~isnan(R2m_1(iv)), text(col_x(7)+0.03,y_txt,sprintf('%.3f',R2m_1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_r); end
        if ~isnan(R2c_1(iv)), text(col_x(8)+0.03,y_txt,sprintf('%.3f',R2c_1(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_r); end

        % Rendu G2
        if freq2(iv)>0
            cf=col_gy; if freq2(iv)>=thresh_pct,cf=col_b; end
            text(col_x(9)+0.03,y_txt,sprintf('%.0f%%',freq2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',cf,'FontWeight','bold');
        end
        if ~isnan(med_w2(iv)), text(col_x(10)+0.03,y_txt,sprintf('%.3f',med_w2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_b,'FontWeight','bold'); end
        if ~isnan(med_b2(iv)), text(col_x(11)+0.03,y_txt,sprintf('%.3f',med_b2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_b); end
        if ~isnan(iqr_w2(iv)), text(col_x(12)+0.03,y_txt,sprintf('%.3f',iqr_w2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_gy); end
        if ~isnan(dR2_2(iv))
            cd=col_gy; if dR2_2(iv)>0.03,cd=col_gr; end
            text(col_x(13)+0.03,y_txt,sprintf('%.4f',dR2_2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',cd,'FontWeight','bold');
        end
        if ~isnan(R2m_2(iv)), text(col_x(14)+0.03,y_txt,sprintf('%.3f',R2m_2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_b); end
        if ~isnan(R2c_2(iv)), text(col_x(15)+0.03,y_txt,sprintf('%.3f',R2c_2(iv)),'FontSize',6.5,'VerticalAlignment','middle','HorizontalAlignment','center','Color',col_b); end

        line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
    end
    line([0 1 1 0 0],[0 0 1 1 0],'Color',[0.4 0.4 0.5],'LineWidth',1); hold off;
    fig_name=sprintf('%s_Indiv_GroupLASSO_Stabilite_%s.png',lag_mode,bn);
    exportgraphics(fig,fullfile(PathSave,fig_name),'Resolution',150); close(fig);
end

%% -----------------------------------------------------------------------
% ÉTAPE 7 : EXPORTATION ET SAUVEGARDE (variante courante : lag_mode/lag_lbl)
% -----------------------------------------------------------------------

% --- Archivage des résultats de stabilité (fichier séparé) ---
% Permet de régénérer les figures sans relancer les 10 runs
stability_archive = struct();
stability_archive.stab_G1      = stab_G1;
stability_archive.stab_G2      = stab_G2;
stability_archive.N_stab       = N_stab;
stability_archive.THRESH_STAB  = THRESH_STAB;
stability_archive.sd_Y_global  = sd_Y_global;
stability_archive.bloc_names   = bloc_names;
stability_archive.vnames_bloc  = vnames_bloc;
stability_archive.seeds_stab   = seeds_stab;
stability_archive.lag_mode     = lag_mode;
stability_archive.date_run     = datestr(now, 'yyyy-mm-dd HH:MM');

save(fullfile(PathSave,sprintf('%s_Stability_Archive.mat',lag_mode)), 'stability_archive', '-v7.3');
fprintf('Stabilité archivée : %s_Stability_Archive.mat\n', lag_mode);

% --- Sauvegarde principale pour cette variante ---
save(fullfile(PathSave,sprintf('%s_Indiv_LASSO_LMM_Results.mat',lag_mode)),...
    'sel_vnames','sel_source','mundlak_G1','mundlak_G2','lmm_G1','lmm_G2',...
    'stab_G1','stab_G2','sd_Y_global');
fprintf('\n========================================================================\n');
fprintf('  VARIANTE [%s] EXÉCUTÉE AVEC SUCCÈS | RÉSULTATS SAUVEGARDÉS\n', lag_lbl);
fprintf('========================================================================\n');

% --- Stockage du résumé de cette variante pour la comparaison finale ---
% Recense, par bloc et par groupe, les variables sélectionnées (LASSO),
% leur phi AR(1) résiduel post-LMM, et le R² ajusté — pour comparer les
% 4 fenêtres de différence en fin de script.
Results_AllLagModes.(lag_mode).lag_lbl   = lag_lbl;
Results_AllLagModes.(lag_mode).mundlak_G1 = mundlak_G1;
Results_AllLagModes.(lag_mode).mundlak_G2 = mundlak_G2;
Results_AllLagModes.(lag_mode).lmm_G1     = lmm_G1;
Results_AllLagModes.(lag_mode).lmm_G2     = lmm_G2;
Results_AllLagModes.(lag_mode).stab_G1    = stab_G1;
Results_AllLagModes.(lag_mode).stab_G2    = stab_G2;

end % for iLag — fin de la boucle externe sur les 4 variantes de fenêtre

%% -----------------------------------------------------------------------
% ÉTAPE 8 : COMPARAISON DES 4 VARIANTES DE FENÊTRE
% -----------------------------------------------------------------------
fprintf('\n\n');
fprintf('#########################################################################\n');
fprintf('###  ÉTAPE 8 : TABLEAU COMPARATIF DES 4 VARIANTES DE FENÊTRE         ###\n');
fprintf('#########################################################################\n\n');

fprintf('%-10s  %-10s  %-6s  %-35s  %-8s  %-8s\n', ...
    'Fenêtre', 'Bloc', 'Grp', 'Variable', 'phi_AR1', 'R2_adj');
fprintf('%s\n', repmat('-', 1, 90));

for iLagCmp = 1:length(lag_modes)
    lag_mode = lag_modes{iLagCmp};
    if ~isfield(Results_AllLagModes, lag_mode), continue; end
    R = Results_AllLagModes.(lag_mode);

    for ib = 1:4
        bn = bloc_names{ib};
        for iGrp = 1:2
            if iGrp==1, gname='G1'; lmm_s = R.lmm_G1; else, gname='G2'; lmm_s = R.lmm_G2; end
            if ~isfield(lmm_s, bn), continue; end
            res = lmm_s.(bn);
            if ~isfield(res,'names') || isempty(res.names), continue; end
            for iv = 1:length(res.names)
                phi_v = NaN;
                if isfield(res,'phi_ar1'), phi_v = res.phi_ar1; end
                fprintf('%-10s  %-10s  %-6s  %-35s  %-8.3f  %-8.3f\n', ...
                    lag_mode, bn, gname, res.names{iv}, phi_v, res.R2_adj);
            end
        end
    end
end

save(fullfile(PathSave,'Comparaison_4_Fenetres.mat'), 'Results_AllLagModes', 'lag_modes', 'lag_mode_lbls');
fprintf('\nComparaison sauvegardée : Comparaison_4_Fenetres.mat\n');
fprintf('\n========================================================================\n');
fprintf('  PIPELINE COMPLET (4 VARIANTES) EXÉCUTÉ AVEC SUCCÈS\n');
fprintf('========================================================================\n');

%% =======================================================================
%% FONCTIONS LOCALES
%% =======================================================================

function beta = fitGroupLassoFISTA(X, Y, group_idx, N_groups, lambda, max_iter, tol, beta_init, L_step)
    % Optimiseur autonome Group LASSO via FISTA (Lipschitz pré-calculée)
    N = length(Y);
    step_size = 1/L_step;
    beta = beta_init;
    z = beta;
    t_prev = 1;
    for iter = 1:max_iter
        beta_prev = beta;
        grad = -(X'*(Y-X*z))/N;
        u = z - step_size*grad;
        beta = groupSoftThreshold(u, step_size*lambda, group_idx, N_groups, size(X,2));
        t_curr = (1 + sqrt(1 + 4 * t_prev^2)) / 2;
        z = beta + ((t_prev - 1) / t_curr) * (beta - beta_prev);
        t_prev = t_curr;
        if max(abs(beta-beta_prev)) < tol, break; end
    end
end

function beta_prox = groupSoftThreshold(u, gamma, group_idx, N_groups, P)
    % Opérateur de projection proximal Group Soft-Thresholding
    beta_prox = zeros(P, 1);
    for g = 1:N_groups
        idx_g = (group_idx == g);
        u_g = u(idx_g);
        norm_g = norm(u_g, 2);
        if norm_g > gamma
            beta_prox(idx_g) = (1 - gamma/norm_g) * u_g;
        end
    end
end

function lambda_max = computeLambdaMax(X, Y, group_idx, N_groups)
    % Calcule analytiquement le point exact d'annulation de tous les coefficients
    N = length(Y);
    lambda_max = 0;
    for g = 1:N_groups
        idx_g = (group_idx == g);
        ng = norm(X(:,idx_g)'*Y, 2)/N;
        if ng > lambda_max, lambda_max = ng; end
    end
    lambda_max = lambda_max * 1.01;
end

function fold_id = groupKFoldBySubject(S, uSubj, K)
    % Stratification stricte par sujet (zéro leakage)
    N_subj = length(uSubj);
    fold_id = zeros(length(S), 1);
    perm = randperm(N_subj);
    subj_fold = mod((0:N_subj-1), K) + 1;
    subj_fold = subj_fold(perm);
    for iS = 1:N_subj
        idx_s = (S == uSubj(iS));
        fold_id(idx_s) = subj_fold(iS);
    end
end

function shadedErrorBar_simple(x, y, err, col)
    x = x(:)'; y = y(:)'; err = err(:)';
    x_patch = [x, fliplr(x)];
    y_patch = [y + err, fliplr(y - err)];
    patch(x_patch, y_patch, col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(x, y, 'Color', col, 'LineWidth', 1.8);
end

function [rho_moy, p_ttest] = fisherTtest(Mat_var, Mat_RPE, Vi, MIN_BINS)
    N_subj = length(Vi);
    z_vals = NaN(N_subj, 1);
    for iG = 1:N_subj
        iP = Vi(iG);
        x = Mat_var(iP,:)';
        rpe = Mat_RPE(iP,:)';
        ok = ~isnan(x) & ~isnan(rpe);
        if sum(ok) < MIN_BINS, continue; end
        [rho, ~] = corr(x(ok), rpe(ok), 'Type', 'Spearman');
        rho = max(min(rho, 0.9999), -0.9999);
        z_vals(iG) = atanh(rho);
    end
    z_ok = z_vals(~isnan(z_vals));
    if length(z_ok) < 3
        rho_moy = NaN; p_ttest = NaN; return;
    end
    [~, p_ttest] = ttest(z_ok);
    rho_moy = tanh(mean(z_ok));
end

function [Xw_long, Xb_long, RPE_long, subj_long] = buildBlock(Vi, Mat_data, Mat_RPE, N_bins, lag_mode)
% Construit le tableau long pour la formulation changement<->changement.
%
% lag_mode : 'lag1' | 'lag2' | 'lag3' | 'global'
%   lag1/lag2/lag3 : DeltaRPE(t) = RPE(t) - RPE(t-lag), pour t = lag+1..N_bins
%                    DeltaX(t)   = X(t)   - X(t-lag)
%   global         : DeltaRPE = RPE(N_bins) - RPE(1) (une seule obs/sujet)
%                    DeltaX   = X(N_bins) - X(1)
%
% Décomposition Mundlak appliquée à DeltaX :
%   Xb_long (Between) = moyenne sur la session de X(t) brut (niveau moyen
%                       du sujet, constant dans le temps, calculé sur tout
%                       l'intervalle 1..N_bins avant différenciation)
%   Xw_long (Within)  = DeltaX_within(t) = DeltaX(t) lui-même, car
%                       soustraire deux fois la même moyenne_sujet à t et
%                       t-lag l'annule mathématiquement — DeltaX brut EST
%                       déjà la composante Within de la différence.
%
% Note : subj_long stocke iG (indice local 1..N_subj), pas iP (indice
% absolu), pour garantir la cohérence entre blocs avec des Vi différents.

    N_subj = length(Vi);

    switch lag_mode
        case 'lag1', lag = 1;
        case 'lag2', lag = 2;
        case 'lag3', lag = 3;
        case 'global', lag = NaN;  % traité à part
        otherwise, error('lag_mode inconnu : %s', lag_mode);
    end

    if strcmp(lag_mode, 'global')
        % --- Une seule observation par sujet : intervalle N_bins - intervalle 1 ---
        N_obs = N_subj;
        Xw_long   = NaN(N_obs, 1);
        Xb_long   = NaN(N_obs, 1);
        RPE_long  = NaN(N_obs, 1);
        subj_long = NaN(N_obs, 1);

        for iG = 1:N_subj
            iP       = Vi(iG);
            x_subj   = Mat_data(iP, :);
            rpe_subj = Mat_RPE(iP, :);

            Xw_long(iG)   = x_subj(N_bins) - x_subj(1);     % DeltaX global
            Xb_long(iG)   = mean(x_subj, 'omitnan');        % niveau moyen sujet
            RPE_long(iG)  = rpe_subj(N_bins) - rpe_subj(1); % DeltaRPE global
            subj_long(iG) = iG;
        end
    else
        % --- Lag-k : N_bins - lag observations par sujet ---
        N_obs   = N_subj * (N_bins - lag);
        Xw_long   = NaN(N_obs, 1);
        Xb_long   = NaN(N_obs, 1);
        RPE_long  = NaN(N_obs, 1);
        subj_long = NaN(N_obs, 1);
        row_idx = 0;

        for iG = 1:N_subj
            iP       = Vi(iG);
            x_subj   = Mat_data(iP, :);
            rpe_subj = Mat_RPE(iP, :);
            mu_subj  = mean(x_subj, 'omitnan');   % X_between : niveau moyen du sujet

            for k = (lag+1):N_bins
                row_idx = row_idx + 1;
                Xw_long(row_idx)   = x_subj(k) - x_subj(k-lag);       % DeltaX_within(t)
                Xb_long(row_idx)   = mu_subj;                          % X_between (constant/sujet)
                RPE_long(row_idx)  = rpe_subj(k) - rpe_subj(k-lag);   % DeltaRPE(t)
                subj_long(row_idx) = iG;
            end
        end
    end
end

function RPE_lag = makeLagBySubject(Y, S)
    RPE_lag = NaN(size(Y));
    uS = unique(S(:))';
    for iSubj = uS
        idx_s = find(S == iSubj);
        if numel(idx_s) >= 2
            RPE_lag(idx_s(2:end)) = Y(idx_s(1:end-1));
        end
    end
end

function [R2m, R2c] = computeNakagawaR2(lme, T)
    yhat_fix = fitted(lme, 'Conditional', false);
    varF = var(yhat_fix, 1, 'omitnan');
    try
        Z = designMatrix(lme, 'Random');
        b = randomEffects(lme);
        yhat_rand = Z * b;
        varR = var(yhat_rand, 1, 'omitnan');
    catch
        varR = NaN;
    end
    sigma2 = lme.MSE;
    if isnan(varR)
        R2m = varF / (varF + sigma2);
        R2c = NaN;
    else
        denom = varF + varR + sigma2;
        if denom <= 0 || isnan(denom)
            R2m = NaN; R2c = NaN;
        else
            R2m = varF / denom;
            R2c = (varF + varR) / denom;
        end
    end
end

function sig = p2star(p)
    if isnan(p), sig = '';
    elseif p < 0.001, sig = '***';
    elseif p < 0.01, sig = '**';
    elseif p < 0.05, sig = '*';
    else, sig = ''; end
end

function lbl = normLabel(lbl)
    lbl = strrep(lbl, '_', ' ');
    lbl = strrep(lbl, ' — ', ' - ');
    lbl = strrep(lbl, '—', ' - ');
    lbl = strrep(lbl, 'TFR ', '');
    lbl = strrep(lbl, 'Goub ', '');
    lbl = strtrim(lbl);
end

function Sigma = buildAR1blockDiag(phi, N_subj, T_i)
% Construit la matrice de covariance AR(1) bloc-diagonale
% pour N_subj sujets avec T_i observations chacun.
% Chaque bloc T_i×T_i : Sigma_ij = phi^|i-j|
    N_total = N_subj * T_i;
    Sigma   = zeros(N_total, N_total);
    exp_mat = abs((0:T_i-1)' - (0:T_i-1));
    bloc    = phi .^ exp_mat;
    for iS = 1:N_subj
        r0 = (iS-1)*T_i + 1;
        r1 = iS*T_i;
        Sigma(r0:r1, r0:r1) = bloc;
    end
    Sigma = Sigma + 1e-8 * eye(N_total);
end

function nll = ar1NegLogLik(phi, Y, Xfixed, Zrand, N_subj, T_i)
% Log-vraisemblance REML négative calculée manuellement.
%
% Modèle : Y = X*β + Z*b + ε
%   b ~ N(0, σ²_b * I_q)
%   ε ~ N(0, σ²_e * R)   où R = AR(1) bloc-diagonale
%
% Vraisemblance marginale (intégrée sur b) :
%   Y ~ N(X*β, V)  avec V = σ²_e * R + Z * (σ²_b * I) * Z'
%
% Estimation profilée :
%   1) Pour phi fixé, estimer β par GLS : β = (X' V^{-1} X)^{-1} X' V^{-1} Y
%   2) Estimer σ²_e et σ²_b par grille 1D sur ratio τ = σ²_b / σ²_e
%   3) Calculer log|V| + résidus'V^{-1}résidus et retourner -REML

    if phi <= 0 || phi >= 0.99
        nll = Inf; return;
    end

    N   = length(Y);
    p   = size(Xfixed, 2);  % nb paramètres fixes (incl. intercept)
    R   = buildAR1blockDiag(phi, N_subj, T_i);

    % --- Grille sur tau = sigma2_b / sigma2_e ---
    tau_grid = [0, logspace(-3, 3, 30)];
    best_nll = Inf;

    for tau = tau_grid
        try
            % V = sigma2_e * (R + tau * Z*Z')
            % On travaille sur V_norm = R + tau*Z*Z' (facteur sigma2_e profilé)
            V_norm = R + tau * (Zrand * Zrand');
            V_norm = (V_norm + V_norm') / 2;  % forcer symétrie

            % Décomposition de Cholesky de V_norm
            [L, flag] = chol(V_norm, 'lower');
            if flag ~= 0, continue; end

            % Résolution par substitution
            LinvX = L \ Xfixed;   % L^{-1} X
            LinvY = L \ Y(:);     % L^{-1} Y

            % GLS : β = (X' V^{-1} X)^{-1} X' V^{-1} Y
            XtVinvX = LinvX' * LinvX;
            XtVinvY = LinvX' * LinvY;
            [~, rflag] = chol(XtVinvX);
            if rflag ~= 0, continue; end
            beta_gls = XtVinvX \ XtVinvY;

            % Résidus GLS
            resid_gls = LinvY - LinvX * beta_gls;

            % Estimation profilée de sigma2_e
            % REML : sigma2_e = RSS / (N - p)
            RSS      = resid_gls' * resid_gls;
            sigma2_e = RSS / (N - p);
            if sigma2_e <= 0, continue; end

            % Log-vraisemblance REML
            % logL_REML = -0.5 * [(N-p)*log(sigma2_e) + log|V_norm|
            %              + log|X'V^{-1}X| + RSS/sigma2_e]
            logdetV   = 2 * sum(log(diag(L)));
            logdetXVX = 2 * sum(log(diag(chol(XtVinvX))));
            nll_cur   = 0.5 * ((N-p)*log(sigma2_e) + logdetV ...
                         + logdetXVX + (N-p));

            if nll_cur < best_nll
                best_nll = nll_cur;
            end
        catch
            continue;
        end
    end

    nll = best_nll;
end

function [beta_out, ci_lo, ci_hi, pval_out, sigma2_e, sigma2_b] = ...
        fitLMM_AR1(phi, Y, Xfixed, Zrand, N_subj, T_i, col_names)
% Ajustement final du LMM AR(1) avec phi fixé.
% Retourne β, IC 95% (Wald), p-valeurs, σ²_e, σ²_b.
%
% Variance de β : Var(β) = σ²_e * (X' V^{-1} X)^{-1}
% Test de Wald  : t = β / se(β), ddl ≈ N - p (Satterthwaite simplifié)

    N   = length(Y);
    p   = size(Xfixed, 2);
    R   = buildAR1blockDiag(phi, N_subj, T_i);

    % Optimiser tau par grille fine
    tau_grid = [0, logspace(-3, 3, 50)];
    best_nll = Inf; best_tau = 0;

    for tau = tau_grid
        try
            V_norm = R + tau * (Zrand * Zrand');
            V_norm = (V_norm + V_norm') / 2;
            [~, flag] = chol(V_norm, 'lower');
            if flag ~= 0, continue; end
            nll_cur = ar1NegLogLik_tau(phi, tau, Y, Xfixed, Zrand, N_subj, T_i, R);
            if nll_cur < best_nll
                best_nll = nll_cur; best_tau = tau;
            end
        catch; continue; end
    end

    % Ajustement final avec best_tau
    V_norm = R + best_tau * (Zrand * Zrand');
    V_norm = (V_norm + V_norm') / 2;
    L      = chol(V_norm, 'lower');
    LinvX  = L \ Xfixed;
    LinvY  = L \ Y(:);
    XtVinvX = LinvX' * LinvX;
    XtVinvY = LinvX' * LinvY;
    beta_gls = XtVinvX \ XtVinvY;
    resid_gls_raw = Y(:) - Xfixed * beta_gls;
    sigma2_e = (resid_gls_raw' * (V_norm \ resid_gls_raw)) / (N - p);
    sigma2_b = best_tau * sigma2_e;

    % Variance de β et inférence
    Var_beta = sigma2_e * inv(XtVinvX);
    se_beta  = sqrt(diag(Var_beta));
    df       = max(N - p, 1);
    t_vals   = beta_gls ./ se_beta;
    pval_out = 2 * (1 - tcdf(abs(t_vals), df));
    t_crit   = tinv(0.975, df);
    ci_lo    = beta_gls - t_crit * se_beta;
    ci_hi    = beta_gls + t_crit * se_beta;

    % Retourner uniquement les effets fixes (sans intercept)
    beta_out = beta_gls(2:end);
    ci_lo    = ci_lo(2:end);
    ci_hi    = ci_hi(2:end);
    pval_out = pval_out(2:end);
end

function nll = ar1NegLogLik_tau(phi, tau, Y, Xfixed, Zrand, N_subj, T_i, R)
% Vraisemblance pour un (phi, tau) donné — fonction auxiliaire
    N = length(Y); p = size(Xfixed,2);
    V_norm = R + tau * (Zrand * Zrand');
    V_norm = (V_norm + V_norm') / 2;
    [L, flag] = chol(V_norm, 'lower');
    if flag ~= 0, nll = Inf; return; end
    LinvX = L \ Xfixed; LinvY = L \ Y(:);
    XtVinvX = LinvX' * LinvX;
    beta_gls = XtVinvX \ (LinvX' * LinvY);
    resid = LinvY - LinvX * beta_gls;
    RSS = resid' * resid;
    sigma2_e = RSS / (N - p);
    if sigma2_e <= 0, nll = Inf; return; end
    logdetV = 2 * sum(log(diag(L)));
    logdetXVX = 2 * sum(log(diag(chol(XtVinvX))));
    nll = 0.5 * ((N-p)*log(sigma2_e) + logdetV + logdetXVX + (N-p));
end