%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  ANALYSE DÉRIVE BIOMARQUEURS ↔ ΔRPE                              %%%%
%%%%  Métrique : X(bin10) − X(bin1) vs RPE(bin10) − RPE(bin1)        %%%%
%%%%  Corrélation Spearman pour 5 groupes :                           %%%%
%%%%    G1 (Short Duration, jeu 1)                                    %%%%
%%%%    G2 (Long Duration, jeu 1)                                     %%%%
%%%%    Normal (SD+LD confondus, jeu 2, condition Normal)             %%%%
%%%%    SD (petites mains, jeu 2, condition Normal)                   %%%%
%%%%    LD (grandes mains, jeu 2, condition Normal)                   %%%%
%%%%  Variables : 93 variables (SegMod, SegAxis, Goubault, EMG)      %%%%
%%%%  Sorties   : tableaux PNG par bloc                               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
% CHEMINS
% -----------------------------------------------------------------------
PathSave = fileparts(mfilename('fullpath'));
if isempty(PathSave), PathSave = pwd; end

PathG1_IMU = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';
PathG2_IMU = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_long_duration\';
PathG1_EMG = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\1_EMG\';
PathG2_EMG = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_long_duration\2_EMG\';
PathJeu2   = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';

%% -----------------------------------------------------------------------
% CHARGEMENT JEU 1 (Goubault, G1/G2)
% -----------------------------------------------------------------------
fprintf('Chargement jeu 1 (Goubault G1/G2)...\n');
try
    d1  = load(fullfile(PathG1_IMU, 'Workload_Li_TimeNormalised.mat'));
    d2  = load(fullfile(PathG2_IMU, 'Workload_Li_TimeNormalised_G2.mat'));
    dG1 = load(fullfile(PathG1_IMU, 'Workload_Li_TimeNormalised_Goubault.mat'));
    dG2 = load(fullfile(PathG2_IMU, 'Workload_Li_TimeNormalised_Goubault_G2.mat'));
    eG1 = load(fullfile(PathG1_EMG, 'EMG_Li_TimeNormalised_G1.mat'));
    eG2 = load(fullfile(PathG2_EMG, 'EMG_Li_TimeNormalised_G2.mat'));
    fprintf('--> Jeu 1 chargé.\n');
catch
    warning('Chemin OneDrive inaccessible, tentative dossier courant.');
    d1  = load('Workload_Li_TimeNormalised.mat');
    d2  = load('Workload_Li_TimeNormalised_G2.mat');
    dG1 = load('Workload_Li_TimeNormalised_Goubault.mat');
    dG2 = load('Workload_Li_TimeNormalised_Goubault_G2.mat');
    eG1 = load('EMG_Li_TimeNormalised_G1.mat');
    eG2 = load('EMG_Li_TimeNormalised_G2.mat');
end

Vi1_26 = find(d1.valid_26);
Vi2_26 = find(d2.valid_26);
Vi1    = find(d1.valid_subj);
Vi2    = find(d2.valid_subj);
Vi1g   = find(dG1.valid_subj);
Vi2g   = find(dG2.valid_subj);
Vi1e   = find(eG1.valid_subj);
Vi2e   = find(eG2.valid_subj);

N_bins   = 10;
Segments = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Axes_lbl = {'X','Y','Z'};
N_seg    = 7;
Channels     = eG1.Channels;
VarNames_emg = eG1.VarNames;
N_ch     = length(Channels);
MIN_BINS = 8;

% --- Correction Inf/NaN SampleEntropy EMG ---
fprintf('\n=== Correction Inf/NaN SampleEntropy Mat_EMG ===\n');
for iG = 1:2
    if iG==1, Mat=eG1.Mat_EMG; VN=eG1.VarNames; gname='G1';
    else,     Mat=eG2.Mat_EMG; VN=eG2.VarNames; gname='G2'; end
    n_inf = sum(isinf(Mat(:))); Mat(isinf(Mat)) = NaN;
    [NS,NT,NC,NV] = size(Mat);
    for iP=1:NS
        for iC=1:NC
            for iV=1:NV
                x=squeeze(Mat(iP,:,iC,iV)); bad=isnan(x);
                if any(bad)&&~all(bad), x(bad)=mean(x(~bad)); Mat(iP,:,iC,iV)=x; end
            end
        end
    end
    iV_se = find(strcmp(VN,'SampleEntropy')); excl=[];
    for iP=1:NS
        for iC=1:NC
            if all(isnan(squeeze(Mat(iP,:,iC,iV_se)))), excl(end+1)=iP; break; end
        end
    end
    excl=unique(excl);
    fprintf('  %s : %d Inf corrigés | %d sujets exclus EMG : [%s]\n',gname,n_inf,length(excl),num2str(excl));
    if iG==1, eG1.Mat_EMG=Mat; Vi1e_v=Vi1e; Vi1e_v(excl)=[]; Vi1e=Vi1e_v;
    else,     eG2.Mat_EMG=Mat; Vi2e_v=Vi2e; Vi2e_v(excl)=[]; Vi2e=Vi2e_v; end
end

%% -----------------------------------------------------------------------
% CHARGEMENT JEU 2 (Robin, Normal/CS60, SD/LD)
% -----------------------------------------------------------------------
fprintf('\nChargement jeu 2 (Robin SD/LD)...\n');
try
    load(fullfile(PathJeu2, 'IMU_Features_NewDataset.mat'));   % Results_IMU
    load(fullfile(PathJeu2, 'EMG_Features_NewDataset_25vars.mat'));   % Results_EMG (25 vars)
    load(fullfile(PathJeu2, 'SegMod_SegAxis_Features_NewDataset.mat'));  % Results_SegMod, Results_SegAxis
    load(fullfile(PathJeu2, 'Goubault_Features_NewDataset.mat'));  % Results_Goubault (10 vars)
    fprintf('--> Jeu 2 chargé.\n');
catch
    load('IMU_Features_NewDataset.mat');
    load('EMG_Features_NewDataset_25vars.mat');
    load('SegMod_SegAxis_Features_NewDataset.mat');
    load('Goubault_Features_NewDataset.mat');
end

% Exclusion P07
if isfield(Results_IMU,'P07'), Results_IMU = rmfield(Results_IMU,'P07'); end
if exist('Results_EMG','var') && isfield(Results_EMG,'P07'), Results_EMG = rmfield(Results_EMG,'P07'); end
if exist('Results_Goubault','var') && isfield(Results_Goubault,'P07'), Results_Goubault = rmfield(Results_Goubault,'P07'); end

% Groupes SD / LD (condition Normal uniquement)
SD_subjects = {'P02','P03','P04','P05','P08','P10','P12','P15','P17','P18','P20','P21','P24'};  % N=13
LD_subjects = {'P01','P06','P09','P11','P13','P14','P16','P19','P22','P23'};                      % N=10

% --- RPE début/fin (condition Normal) — pas de champ RPE dans Results_IMU,
% valeurs en dur identiques à Analyse_RPE_SD_LD.m ---
RPE_data_j2 = {
    'P02', 30, 75; 'P03', 25, 85; 'P04', 15, 23; 'P05', 35, 75; 'P08',  7, 18;
    'P10',  7, 30; 'P12', 13, 95; 'P15',  5, 70; 'P17', 13, 60; 'P18', 13, 60;
    'P20', 15, 55; 'P21', 10, 50; 'P24', 10, 25;
    'P01',  0, 30; 'P06', 10, 55; 'P09',  5, 60; 'P11',  7, 60; 'P13', 17, 55;
    'P14', 22, 28; 'P16', 20, 60; 'P19',  0, 25; 'P22',  1, 20; 'P23',  3, 35;
};
RPE_j2 = struct();
for iRpe = 1:size(RPE_data_j2,1)
    s = RPE_data_j2{iRpe,1};
    RPE_j2.(s).debut = RPE_data_j2{iRpe,2};
    RPE_j2.(s).fin   = RPE_data_j2{iRpe,3};
end
% Filtrer les sujets présents dans Results_IMU
all_subj_j2 = fieldnames(Results_IMU);
SD_subjects = SD_subjects(ismember(SD_subjects, all_subj_j2));
LD_subjects = LD_subjects(ismember(LD_subjects, all_subj_j2));
Norm_subjects = [SD_subjects, LD_subjects];

fprintf('Jeu 2 — Normal: N=%d | SD: N=%d | LD: N=%d\n', ...
    length(Norm_subjects), length(SD_subjects), length(LD_subjects));

%% -----------------------------------------------------------------------
% FONCTION LOCALE : CALCUL DÉRIVE + ΔRPE (jeu 1)
% dérive_X = X(bin10) - X(bin1) par sujet
% dérive_RPE = RPE(bin10) - RPE(bin1) par sujet
% -----------------------------------------------------------------------
% (utilisée ci-dessous dans la boucle de calcul)

%% -----------------------------------------------------------------------
% BOUCLE PRINCIPALE : calcul dérive pour les 93 variables
% -----------------------------------------------------------------------
fprintf('\n=== Calcul des dérives et corrélations Spearman ===\n');

% Structures de stockage
% Pré-allocation à taille fixe (93 variables) avec ordre de champs explicite
% — plus robuste que struct([]) face aux variations de forme entre appels.
% IMPORTANT : doit contenir EXACTEMENT les mêmes champs que ceux retournés
% par calcCorr, sinon "incompatible fields in struct assignment".
N_VARS_TOTAL = 93;
Results_corr(N_VARS_TOTAL) = struct( ...
    'vname','', 'source','', ...
    'rho_G1',NaN, 'p_G1',NaN, 'N_G1',0, ...
    'rho_G2',NaN, 'p_G2',NaN, 'N_G2',0, ...
    'rho_Norm',NaN, 'p_Norm',NaN, 'N_Norm',0, ...
    'rho_SD',NaN, 'p_SD',NaN, 'N_SD',0, ...
    'rho_LD',NaN, 'p_LD',NaN, 'N_LD',0, ...
    'rho_partial_G1',NaN, 'p_partial_G1',NaN, ...
    'rho_partial_G2',NaN, 'p_partial_G2',NaN, ...
    'rho_partial_Norm',NaN, 'p_partial_Norm',NaN, ...
    'rho_partial_SD',NaN, 'p_partial_SD',NaN, ...
    'rho_partial_LD',NaN, 'p_partial_LD',NaN, ...
    'raw_dX_G1',  [], 'raw_dRPE_G1',  [], 'raw_RPE0_G1',  [], ...
    'raw_dX_G2',  [], 'raw_dRPE_G2',  [], 'raw_RPE0_G2',  [], ...
    'raw_dX_Norm',[], 'raw_dRPE_Norm',[], 'raw_RPE0_Norm',[], ...
    'raw_dX_SD',  [], 'raw_dRPE_SD',  [], 'raw_RPE0_SD',  [], ...
    'raw_dX_LD',  [], 'raw_dRPE_LD',  [], 'raw_RPE0_LD',  [], ...
    'bins_X_G1',  [], 'bins_RPE_G1',  [], ...
    'bins_X_G2',  [], 'bins_RPE_G2',  [], ...
    'bins_X_Norm',[], 'bins_RPE_Norm',[], ...
    'bins_X_SD',  [], 'bins_RPE_SD',  [], ...
    'bins_X_LD',  [], 'bins_RPE_LD',  []);
% Champs : vname, source, rho_G1, p_G1, N_G1, rho_G2, p_G2, N_G2,
%          rho_Norm, p_Norm, N_Norm, rho_SD, p_SD, N_SD, rho_LD, p_LD, N_LD,
%          rho_partial_{G1,G2,Norm,SD,LD}, p_partial_{G1,G2,Norm,SD,LD},
%          raw_{dX,dRPE,RPE0}_{G1,G2,Norm,SD,LD} (vecteurs bruts, NaN inclus,
%          réutilisés par Analyse_IM_GAM.m), bins_{X,RPE}_{G1,G2,Norm,SD,LD}
%          (matrices [N_sujets x N_bins], réutilisées par Analyse_GAMM_Bins.m)

idx_var = 0;

% =========================================================
% BLOC 1 : SegMod (Features3 — par segment, signal modulé)
% =========================================================
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_Seg_Accel_26; m2=d2.Mat_Seg_Accel_26; RPE1=d1.Mat_RPE_26; RPE2=d2.Mat_RPE_26; Vi1_loc=Vi1_26; Vi2_loc=Vi2_26;
    else,        sn='Jerk';  m1=d1.Mat_Seg_Jerk_26;  m2=d2.Mat_Seg_Jerk_26;  RPE1=d1.Mat_RPE_26; RPE2=d2.Mat_RPE_26; Vi1_loc=Vi1_26; Vi2_loc=Vi2_26; end
    for iSeg = 1:N_seg
        vn = [sn '_' Segments{iSeg}];
        M1 = squeeze(m1(:,:,iSeg)); M2 = squeeze(m2(:,:,iSeg));
        idx_var = idx_var + 1;
        Results_corr(idx_var) = calcCorr(vn, 'SegMod', M1, RPE1, Vi1_loc, M2, RPE2, Vi2_loc, ...
            Results_SegMod, [], SD_subjects, LD_subjects, RPE_j2, vn, 'IMU');
    end
end

% Total (SegMod)
for iSig = 1:2
    if iSig==1, sn='Total_Accel'; m1=d1.Mat_Total_Accel_26; m2=d2.Mat_Total_Accel_26; RPE1=d1.Mat_RPE_26; RPE2=d2.Mat_RPE_26; Vi1_loc=Vi1_26; Vi2_loc=Vi2_26;
    else,        sn='Total_Jerk';  m1=d1.Mat_Total_Jerk_26;  m2=d2.Mat_Total_Jerk_26;  RPE1=d1.Mat_RPE_26; RPE2=d2.Mat_RPE_26; Vi1_loc=Vi1_26; Vi2_loc=Vi2_26; end
    idx_var = idx_var + 1;
    Results_corr(idx_var) = calcCorr(sn, 'SegMod', m1, RPE1, Vi1_loc, m2, RPE2, Vi2_loc, ...
        Results_SegMod, [], SD_subjects, LD_subjects, RPE_j2, sn, 'IMU');
end

% =========================================================
% BLOC 2 : SegAxis (par segment × axe)
% =========================================================
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_SegAxis_Accel; m2=d2.Mat_SegAxis_Accel; RPE1=d1.Mat_RPE; RPE2=d2.Mat_RPE; Vi1_loc=Vi1; Vi2_loc=Vi2;
    else,        sn='Jerk';  m1=d1.Mat_SegAxis_Jerk;  m2=d2.Mat_SegAxis_Jerk;  RPE1=d1.Mat_RPE; RPE2=d2.Mat_RPE; Vi1_loc=Vi1; Vi2_loc=Vi2; end
    for iSeg = 1:N_seg
        for iAx = 1:3
            vn = [sn '_Axe' Axes_lbl{iAx} '_' Segments{iSeg}];
            M1 = squeeze(m1(:,:,iSeg,iAx)); M2 = squeeze(m2(:,:,iSeg,iAx));
            idx_var = idx_var + 1;
            Results_corr(idx_var) = calcCorr(vn, 'SegAxis', M1, RPE1, Vi1_loc, M2, RPE2, Vi2_loc, ...
                Results_SegAxis, [], SD_subjects, LD_subjects, RPE_j2, vn, 'IMU');
        end
    end
end

% =========================================================
% BLOC 3 : Goubault
% =========================================================
% Mapping explicite par index — les noms de champs produits par
% Compute_Goubault_NewDataset.m ne suivent pas la même convention que
% dG1.Vars{iV,1} (nom d'affichage avec tirets), donc on mappe directement
% par position (même ordre confirmé dans les deux scripts).
vn_goubault_j2 = {
    'Hand_MedianFreq_AccelY', 'Hand_Prc90_AccelX', 'Hand_PeakPower_AccelModule', ...
    'Hand_SpectralEntropy_AccelY', 'Head_PeakPower_AngVelX', 'Forearm_Mean_AccelY', ...
    'Forearm_SpectralEntropy_AngVelMod', 'Forearm_PeakPowerFreq_AngVelMod', ...
    'Shoulder_Mean_AngVelX', 'Shoulder_Mean_AngVelY'
};

for iV = 1:size(dG1.Vars,1)
    vn_raw = strrep(strrep(strrep(dG1.Vars{iV,1},' ','_'),'-','_'),'—','_');
    vn = ['Goub_' vn_raw];
    M1 = squeeze(dG1.Mat_Vars(:,:,iV)); M2 = squeeze(dG2.Mat_Vars(:,:,iV));
    RPE1 = dG1.Mat_RPE; RPE2 = dG2.Mat_RPE;
    idx_var = idx_var + 1;
    vn_j2 = vn_goubault_j2{iV};
    Results_corr(idx_var) = calcCorr(vn, 'Goubault', M1, RPE1, Vi1g, M2, RPE2, Vi2g, ...
        Results_Goubault, [], SD_subjects, LD_subjects, RPE_j2, vn_j2, 'IMU');
end

% =========================================================
% BLOC 4 : EMG
% =========================================================
for iV = 1:length(VarNames_emg)
    if strcmp(VarNames_emg{iV},'Amplitude'), continue; end
    for iC = 1:N_ch
        vn = [VarNames_emg{iV} '_' Channels{iC}];
        M1 = squeeze(eG1.Mat_EMG(:,:,iC,iV)); M2 = squeeze(eG2.Mat_EMG(:,:,iC,iV));
        RPE1 = eG1.Mat_RPE; RPE2 = eG2.Mat_RPE;
        idx_var = idx_var + 1;
        Results_corr(idx_var) = calcCorr(vn, 'EMG', M1, RPE1, Vi1e, M2, RPE2, Vi2e, ...
            [], Results_EMG, SD_subjects, LD_subjects, RPE_j2, vn, 'EMG');
    end
end

fprintf('--> %d variables calculées.\n', idx_var);

%% -----------------------------------------------------------------------
% GÉNÉRATION DES TABLEAUX PNG PAR BLOC
% -----------------------------------------------------------------------
fprintf('\n=== Génération des tableaux PNG ===\n');

blocs      = {'SegMod','SegAxis','Goubault','EMG'};
blocs_lbl  = {'BySegMod','BySegAxis','Goubault','EMG'};
groupes    = {'G1','G2','Normal','SD','LD'};
col_grp    = {[0.80 0.15 0.15],[0.15 0.40 0.85],[0.20 0.55 0.20],[0.55 0.20 0.70],[0.90 0.50 0.05]};

% --- Tableaux corrélation SIMPLE (déjà existants) ---
for ib = 1:length(blocs)
    bn = blocs{ib};
    idx_bloc = find(strcmp({Results_corr.source}, bn));
    if isempty(idx_bloc), continue; end
    drawCorrTable(Results_corr, idx_bloc, bn, blocs_lbl{ib}, groupes, col_grp, ...
        'rho_G1','p_G1','N_G1','rho_G2','p_G2','N_G2', ...
        'rho_Norm','p_Norm','N_Norm','rho_SD','p_SD','N_SD','rho_LD','p_LD','N_LD', ...
        'Corrélations Spearman Dérive↔ΔRPE', ...
        sprintf('Derive_DRPE_Spearman_%s.png', bn), PathSave);
end

% --- Tableaux corrélation PARTIELLE (contrôlant RPE0) ---
fprintf('\n=== Génération des tableaux PNG — corrélation partielle (contrôle RPE0) ===\n');
for ib = 1:length(blocs)
    bn = blocs{ib};
    idx_bloc = find(strcmp({Results_corr.source}, bn));
    if isempty(idx_bloc), continue; end
    drawCorrTable(Results_corr, idx_bloc, bn, blocs_lbl{ib}, groupes, col_grp, ...
        'rho_partial_G1','p_partial_G1','N_G1', ...
        'rho_partial_G2','p_partial_G2','N_G2', ...
        'rho_partial_Norm','p_partial_Norm','N_Norm', ...
        'rho_partial_SD','p_partial_SD','N_SD', ...
        'rho_partial_LD','p_partial_LD','N_LD', ...
        'Corrélation partielle Spearman Dérive↔ΔRPE (contrôle RPE0)', ...
        sprintf('Derive_DRPE_SpearmanPartial_%s.png', bn), PathSave);
end

% Sauvegarde des résultats numériques
save(fullfile(PathSave,'Derive_DRPE_Results.mat'),'Results_corr');
fprintf('\n=== PIPELINE TERMINÉ ===\n');

%% =======================================================================
%% FONCTIONS LOCALES
%% =======================================================================

function drawCorrTable(Results_corr, idx_bloc, bn, bn_lbl, groupes, col_grp, ...
    f_rho1,f_p1,f_N1, f_rho2,f_p2,f_N2, f_rho3,f_p3,f_N3, f_rho4,f_p4,f_N4, f_rho5,f_p5,f_N5, ...
    title_txt, fname, PathSave)
% Dessine un tableau de corrélations 5-groupes pour un bloc donné.
% Les noms de champs (f_rho*, f_p*, f_N*) permettent de réutiliser cette
% fonction pour la corrélation simple ET la corrélation partielle, sans
% dupliquer 130 lignes de code de mise en page.

    N_v = length(idx_bloc);

    % --- Dimensionnement en pixels absolus, converti en fractions ensuite ---
    row_h_px    = 20;
    title_h_px  = 30;
    header_h_px = 24;
    legend_h_px = 26;
    margin_px   = 6;

    fig_h = title_h_px + header_h_px + N_v*row_h_px + legend_h_px + 2*margin_px;
    fig_w = 1100;

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
    y_table_bot  = f_legend_h + f_margin;

    fig = figure('Visible','off','Position',[50 50 fig_w fig_h],'Color','white');
    ax  = axes('Position',[0 0 1 1],'Visible','off');
    set(ax,'XLim',[0 1],'YLim',[0 1]); hold on;

    rectangle('Position',[0 y_title_bot 1 f_title_h],'FaceColor',[0.12 0.22 0.42],'EdgeColor','none');
    text(0.5, (y_title_top+y_title_bot)/2, sprintf('%s — %s', title_txt, bn_lbl), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',10,'FontWeight','bold','Color','white');

    col_x   = [0.01, 0.28, 0.40, 0.52, 0.64, 0.76];
    col_w   = [0.27, 0.11, 0.11, 0.11, 0.11, 0.11];
    hdrs    = {'Variable','G1','G2','Normal','SD','LD'};

    rectangle('Position',[0 y_header_bot 1 f_header_h],'FaceColor',[0.78 0.80 0.88],'EdgeColor','none');
    for ic = 1:6
        text(col_x(ic)+col_w(ic)/2, (y_header_top+y_header_bot)/2, hdrs{ic}, ...
            'FontSize', 8, 'FontWeight','bold','Color',[0.1 0.1 0.3],...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    line([0 1],[y_header_bot y_header_bot],'Color',[0.5 0.5 0.6],'LineWidth',1);

    field_rho = {f_rho1,f_rho2,f_rho3,f_rho4,f_rho5};
    field_p   = {f_p1,f_p2,f_p3,f_p4,f_p5};
    field_N   = {f_N1,f_N2,f_N3,f_N4,f_N5};

    for iv = 1:N_v
        r  = Results_corr(idx_bloc(iv));
        y_row = y_table_top - iv * f_row_h;
        y_txt = y_row + f_row_h/2;

        if mod(iv,2)==0
            rectangle('Position',[0 y_row 1 f_row_h],'FaceColor',[0.96 0.96 0.98],'EdgeColor','none');
        end

        vn_disp = strrep(r.vname, '_', ' ');
        if length(vn_disp) > 38, vn_disp = [vn_disp(1:35) '...']; end
        text(col_x(1)+0.003, y_txt, vn_disp, ...
            'FontSize',6.5,'VerticalAlignment','middle','Interpreter','none','Color',[0.1 0.1 0.1]);

        for ig = 1:5
            rho = r.(field_rho{ig}); p = r.(field_p{ig}); N = r.(field_N{ig});
            if isnan(rho) || N < 6, continue; end

            sig = pstar(p);
            txt = sprintf('%.2f%s', rho, sig);
            col_txt = col_grp{ig};
            fw = 'normal';
            if p < 0.05, fw = 'bold'; end

            if p < 0.05
                alpha_bg = min(0.3, abs(rho) * 0.4);
                bg_col = col_txt * alpha_bg + [1 1 1] * (1 - alpha_bg);
                rectangle('Position',[col_x(ig+1) y_row col_w(ig+1) f_row_h],...
                    'FaceColor',bg_col,'EdgeColor','none');
            end

            text(col_x(ig+1)+col_w(ig+1)/2, y_txt, txt, ...
                'FontSize',7,'FontWeight',fw,'Color',col_txt,...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end

        line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
    end

    line([0 1 1 0 0],[y_table_bot y_table_bot y_title_top y_title_top y_table_bot],...
        'Color',[0.4 0.4 0.5],'LineWidth',1);

    Ns = [Results_corr(idx_bloc(1)).N_G1, Results_corr(idx_bloc(1)).N_G2, ...
          Results_corr(idx_bloc(1)).N_Norm, Results_corr(idx_bloc(1)).N_SD, Results_corr(idx_bloc(1)).N_LD];
    leg_str = '';
    for ig = 1:5
        leg_str = [leg_str sprintf('%s N=%d  ', groupes{ig}, Ns(ig))];
    end
    text(0.5, f_margin + f_legend_h/2, ['* p<.05  ** p<.01  *** p<.001    |    ' leg_str], ...
        'FontSize',6,'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);

    hold off;
    exportgraphics(fig, fullfile(PathSave, fname), 'Resolution', 150);
    close(fig);
    fprintf('  Tableau sauvegardé : %s\n', fname);
end

function R = calcCorr(vname, source, M1, RPE1, Vi1, M2, RPE2, Vi2, ...
    Results_IMU, Results_EMG, SD_subjects, LD_subjects, RPE_j2, vn_j2, type_j2)
% Calcule les corrélations Spearman dérive↔ΔRPE pour les 5 groupes,
% en simple ET en partielle (contrôlant RPE0, le niveau de fatigue initial).
% Jeu 1 : matrices [N_sujets x N_bins], indices Vi1/Vi2
% Jeu 2 : structures Results_IMU ou Results_EMG, condition 'Norm' ;
%         RPE_j2.(sujet).debut/.fin — table en dur (pas de champ RPE
%         dans Results_IMU/Results_EMG pour ce dataset).
%
% IMPORTANT : tous les champs sont initialisés explicitement à des types
% fixes (char pour vname/source, double scalaire pour le reste) avant tout
% calcul, pour garantir une structure identique à chaque appel — un struct
% array MATLAB/Octave exige la même "forme" (mêmes champs, mêmes types)
% sur tous ses éléments, sinon "dissimilar structures" lève une erreur.

    R = struct( ...
        'vname',   char(vname), ...
        'source',  char(source), ...
        'rho_G1',  NaN, 'p_G1',  NaN, 'N_G1',  0, ...
        'rho_G2',  NaN, 'p_G2',  NaN, 'N_G2',  0, ...
        'rho_Norm',NaN, 'p_Norm',NaN, 'N_Norm',0, ...
        'rho_SD',  NaN, 'p_SD',  NaN, 'N_SD',  0, ...
        'rho_LD',  NaN, 'p_LD',  NaN, 'N_LD',  0, ...
        'rho_partial_G1',  NaN, 'p_partial_G1',  NaN, ...
        'rho_partial_G2',  NaN, 'p_partial_G2',  NaN, ...
        'rho_partial_Norm',NaN, 'p_partial_Norm',NaN, ...
        'rho_partial_SD',  NaN, 'p_partial_SD',  NaN, ...
        'rho_partial_LD',  NaN, 'p_partial_LD',  NaN, ...
        'raw_dX_G1',  [], 'raw_dRPE_G1',  [], 'raw_RPE0_G1',  [], ...
        'raw_dX_G2',  [], 'raw_dRPE_G2',  [], 'raw_RPE0_G2',  [], ...
        'raw_dX_Norm',[], 'raw_dRPE_Norm',[], 'raw_RPE0_Norm',[], ...
        'raw_dX_SD',  [], 'raw_dRPE_SD',  [], 'raw_RPE0_SD',  [], ...
        'raw_dX_LD',  [], 'raw_dRPE_LD',  [], 'raw_RPE0_LD',  [], ...
        'bins_X_G1',  [], 'bins_RPE_G1',  [], ...
        'bins_X_G2',  [], 'bins_RPE_G2',  [], ...
        'bins_X_Norm',[], 'bins_RPE_Norm',[], ...
        'bins_X_SD',  [], 'bins_RPE_SD',  [], ...
        'bins_X_LD',  [], 'bins_RPE_LD',  [] ...
    );
    % Les champs raw_* stockent les vecteurs bruts dérive/ΔRPE/RPE0 par
    % groupe (avec NaN inclus), pour réutilisation par Analyse_IM_GAM.m.
    % Les champs bins_* stockent les séries complètes [N_sujets x N_bins]
    % (X et RPE, RPE interpolé linéairement pour le jeu 2 — cf.
    % extractBinsJeu2), pour réutilisation par Analyse_GAMM_Bins.m.

    % --- Jeu 1 : G1 ---
    [d_X1, d_RPE1, RPE0_1] = deriveJeu1(M1, RPE1, Vi1);
    [rho, p] = spearman_safe(d_X1, d_RPE1);
    R.rho_G1 = rho; R.p_G1 = p;
    R.N_G1 = sum(~isnan(d_X1) & ~isnan(d_RPE1));
    [rho_p, p_p] = partial_spearman_safe(d_X1, d_RPE1, RPE0_1);
    R.rho_partial_G1 = rho_p; R.p_partial_G1 = p_p;
    R.raw_dX_G1 = d_X1; R.raw_dRPE_G1 = d_RPE1; R.raw_RPE0_G1 = RPE0_1;
    [Xb1, RPEb1] = extractBinsJeu1(M1, RPE1, Vi1);
    R.bins_X_G1 = Xb1; R.bins_RPE_G1 = RPEb1;

    % --- Jeu 1 : G2 ---
    [d_X2, d_RPE2, RPE0_2] = deriveJeu1(M2, RPE2, Vi2);
    [rho, p] = spearman_safe(d_X2, d_RPE2);
    R.rho_G2 = rho; R.p_G2 = p;
    R.N_G2 = sum(~isnan(d_X2) & ~isnan(d_RPE2));
    [rho_p, p_p] = partial_spearman_safe(d_X2, d_RPE2, RPE0_2);
    R.rho_partial_G2 = rho_p; R.p_partial_G2 = p_p;
    R.raw_dX_G2 = d_X2; R.raw_dRPE_G2 = d_RPE2; R.raw_RPE0_G2 = RPE0_2;
    [Xb2, RPEb2] = extractBinsJeu1(M2, RPE2, Vi2);
    R.bins_X_G2 = Xb2; R.bins_RPE_G2 = RPEb2;

    % --- Jeu 2 : Normal (SD+LD) / SD / LD ---
    Norm_subjects = [SD_subjects, LD_subjects];
    [d_Xn, d_RPEn, RPE0_n] = deriveJeu2(Norm_subjects, Results_IMU, Results_EMG, RPE_j2, vn_j2, type_j2);
    [rho, p] = spearman_safe(d_Xn, d_RPEn);
    R.rho_Norm = rho; R.p_Norm = p;
    R.N_Norm = sum(~isnan(d_Xn) & ~isnan(d_RPEn));
    [rho_p, p_p] = partial_spearman_safe(d_Xn, d_RPEn, RPE0_n);
    R.rho_partial_Norm = rho_p; R.p_partial_Norm = p_p;
    R.raw_dX_Norm = d_Xn; R.raw_dRPE_Norm = d_RPEn; R.raw_RPE0_Norm = RPE0_n;
    [Xbn, RPEbn] = extractBinsJeu2(Norm_subjects, Results_IMU, Results_EMG, RPE_j2, vn_j2, type_j2);
    R.bins_X_Norm = Xbn; R.bins_RPE_Norm = RPEbn;

    [d_Xsd, d_RPEsd, RPE0_sd] = deriveJeu2(SD_subjects, Results_IMU, Results_EMG, RPE_j2, vn_j2, type_j2);
    [rho, p] = spearman_safe(d_Xsd, d_RPEsd);
    R.rho_SD = rho; R.p_SD = p;
    R.N_SD = sum(~isnan(d_Xsd) & ~isnan(d_RPEsd));
    [rho_p, p_p] = partial_spearman_safe(d_Xsd, d_RPEsd, RPE0_sd);
    R.rho_partial_SD = rho_p; R.p_partial_SD = p_p;
    R.raw_dX_SD = d_Xsd; R.raw_dRPE_SD = d_RPEsd; R.raw_RPE0_SD = RPE0_sd;
    [Xbsd, RPEbsd] = extractBinsJeu2(SD_subjects, Results_IMU, Results_EMG, RPE_j2, vn_j2, type_j2);
    R.bins_X_SD = Xbsd; R.bins_RPE_SD = RPEbsd;

    [d_Xld, d_RPEld, RPE0_ld] = deriveJeu2(LD_subjects, Results_IMU, Results_EMG, RPE_j2, vn_j2, type_j2);
    [rho, p] = spearman_safe(d_Xld, d_RPEld);
    R.rho_LD = rho; R.p_LD = p;
    R.N_LD = sum(~isnan(d_Xld) & ~isnan(d_RPEld));
    [rho_p, p_p] = partial_spearman_safe(d_Xld, d_RPEld, RPE0_ld);
    R.rho_partial_LD = rho_p; R.p_partial_LD = p_p;
    R.raw_dX_LD = d_Xld; R.raw_dRPE_LD = d_RPEld; R.raw_RPE0_LD = RPE0_ld;
    [Xbld, RPEbld] = extractBinsJeu2(LD_subjects, Results_IMU, Results_EMG, RPE_j2, vn_j2, type_j2);
    R.bins_X_LD = Xbld; R.bins_RPE_LD = RPEbld;
end

% -------------------------------------------------------------------------
function [rho, p] = partial_spearman_safe(x, y, z)
% Corrélation partielle de Spearman entre x et y, contrôlant z (RPE0).
% Équivalent à partialcorr(x,y,z,'Type','Spearman') mais codé localement
% pour ne dépendre d'aucune toolbox : on range-transforme x, y, z, puis
% on calcule la corrélation partielle de Pearson sur les rangs.
    ok = ~isnan(x) & ~isnan(y) & ~isnan(z) & ~isinf(x) & ~isinf(y) & ~isinf(z);
    if sum(ok) < 6
        rho = NaN; p = NaN; return;
    end
    xo = x(ok); yo = y(ok); zo = z(ok);
    n = length(xo);

    % Transformation en rangs (équivalent Spearman)
    rx = tiedrank_local(xo);
    ry = tiedrank_local(yo);
    rz = tiedrank_local(zo);

    % Corrélation partielle de Pearson sur les rangs :
    % r_xy.z = (r_xy - r_xz*r_yz) / sqrt((1-r_xz^2)(1-r_yz^2))
    r_xy = corrcoef_local(rx, ry);
    r_xz = corrcoef_local(rx, rz);
    r_yz = corrcoef_local(ry, rz);

    denom = sqrt((1 - r_xz^2) * (1 - r_yz^2));
    if denom < 1e-10
        rho = NaN; p = NaN; return;
    end
    rho = (r_xy - r_xz*r_yz) / denom;
    rho = max(-1, min(1, rho));  % sécurité numérique

    % Test de significativité : df = n - 3 (un degré de liberté de moins
    % que la corrélation simple, car on contrôle une 3e variable)
    df = n - 3;
    if df < 1 || abs(rho) >= 1
        p = NaN;
    else
        t_stat = rho * sqrt(df / (1 - rho^2));
        p = 2 * (1 - tcdf_local(abs(t_stat), df));
    end
end

% -------------------------------------------------------------------------
function r = tiedrank_local(x)
% Rang moyen en cas d'égalités (équivalent tiedrank de MATLAB)
    [sorted_x, idx] = sort(x(:));
    n = length(x);
    ranks = zeros(n,1);
    i = 1;
    while i <= n
        j = i;
        while j < n && sorted_x(j+1) == sorted_x(i)
            j = j + 1;
        end
        avg_rank = (i + j) / 2;
        ranks(i:j) = avg_rank;
        i = j + 1;
    end
    r = zeros(n,1);
    r(idx) = ranks;
end

% -------------------------------------------------------------------------
function r = corrcoef_local(x, y)
    x = x(:) - mean(x); y = y(:) - mean(y);
    denom = sqrt(sum(x.^2) * sum(y.^2));
    if denom < 1e-10
        r = 0;
    else
        r = sum(x.*y) / denom;
    end
end

% -------------------------------------------------------------------------
function p = tcdf_local(t, df)
% CDF de Student approximée via la fonction beta incomplète régularisée.
% Évite toute dépendance à la Statistics Toolbox (betainc est en Core MATLAB).
    x = df / (df + t^2);
    p_two_tail = betainc(x, df/2, 0.5);
    p = 1 - 0.5 * p_two_tail;
end
function [dX, dRPE, RPE0] = deriveJeu1(Mat, Mat_RPE, Vi)
% Calcule dérive X, ΔRPE, et RPE0 (niveau initial) pour les sujets Vi
    N = length(Vi);
    dX   = NaN(N,1);
    dRPE = NaN(N,1);
    RPE0 = NaN(N,1);
    for i = 1:N
        iP = Vi(i);
        x   = Mat(iP,:);
        rpe = Mat_RPE(iP,:);
        if sum(~isnan(x)) >= 2 && sum(~isnan(rpe)) >= 2
            dX(i)   = x(end)   - x(1);
            dRPE(i) = rpe(end) - rpe(1);
            RPE0(i) = rpe(1);
        end
    end
end

% -------------------------------------------------------------------------
function [X_bins, RPE_bins] = extractBinsJeu1(Mat, Mat_RPE, Vi)
% Extrait les séries complètes par bin (pas seulement dérive bin10-bin1),
% pour réutilisation par Analyse_GAMM_Bins.m (GAMM avec effet aléatoire
% sujet sur les 10 bins, au lieu d'un seul point dérive par sujet).
% X_bins, RPE_bins : [N_sujets x N_bins], NaN si sujet/bin manquant.
    N = length(Vi);
    N_bins_loc = size(Mat, 2);
    X_bins   = NaN(N, N_bins_loc);
    RPE_bins = NaN(N, N_bins_loc);
    for i = 1:N
        iP = Vi(i);
        X_bins(i,:)   = Mat(iP,:);
        RPE_bins(i,:) = Mat_RPE(iP,:);
    end
end

% -------------------------------------------------------------------------
function [dX, dRPE, RPE0] = deriveJeu2(subjects, Results_IMU, Results_EMG, RPE_j2, vn, type_j2)
% Calcule dérive X, ΔRPE, et RPE0 (RPE_début) pour les sujets du jeu 2,
% condition Normal. RPE_j2.(sujet).debut/.fin : table en dur (pas de champ
% RPE dans Results_IMU/Results_EMG pour ce dataset).
    N    = length(subjects);
    dX   = NaN(N,1);
    dRPE = NaN(N,1);
    RPE0 = NaN(N,1);
    for i = 1:N
        subj = subjects{i};
        try
            % --- ΔRPE et RPE0 : toujours depuis la table en dur RPE_j2 ---
            if ~isfield(RPE_j2, subj), continue; end
            d_rpe = RPE_j2.(subj).fin - RPE_j2.(subj).debut;
            rpe0_val = RPE_j2.(subj).debut;

            % --- Variable X : condition Normal ---
            if strcmp(type_j2,'IMU')
                if ~isfield(Results_IMU, subj), continue; end
                kbds = fieldnames(Results_IMU.(subj));
                cond = 'Norm';
                if ~isfield(Results_IMU.(subj), cond)
                    cond_candidates = kbds(contains(kbds,'Norm','IgnoreCase',true));
                    if isempty(cond_candidates), continue; end
                    cond = cond_candidates{1};
                end
                sub = Results_IMU.(subj).(cond);
                if ~isfield(sub, vn), continue; end
                x = sub.(vn)(:)';
            else % EMG
                if isempty(Results_EMG) || ~isfield(Results_EMG, subj), continue; end
                cond = 'Norm';
                if ~isfield(Results_EMG.(subj), cond)
                    kbds2 = fieldnames(Results_EMG.(subj));
                    cond_candidates = kbds2(contains(kbds2,'Norm','IgnoreCase',true));
                    if isempty(cond_candidates), continue; end
                    cond = cond_candidates{1};
                end
                sub = Results_EMG.(subj).(cond);
                if ~isfield(sub, vn), continue; end
                x = sub.(vn)(:)';
            end

            if length(x) >= 2
                dX(i)   = x(end) - x(1);
                dRPE(i) = d_rpe;
                RPE0(i) = rpe0_val;
            elseif length(x)==1
                % X est un scalaire (pas un vecteur de bins) : x(end)-x(1)
                % vaudrait 0 silencieusement — on laisse NaN.
                continue;
            end
        catch
            % Sujet ou variable absent — on laisse NaN
        end
    end
end

% -------------------------------------------------------------------------
function [X_bins, RPE_bins, subj_used] = extractBinsJeu2(subjects, Results_IMU, Results_EMG, RPE_j2, vn, type_j2)
% Extrait les séries complètes par bin pour le jeu 2 (Robin), pour
% réutilisation par Analyse_GAMM_Bins.m.
%
% APPROXIMATION IMPORTANTE : le jeu 2 ne fournit que RPE_début/RPE_fin
% (table en dur RPE_j2), pas de RPE intermédiaire mesuré par bin. Le RPE
% par bin est donc approximé par INTERPOLATION LINÉAIRE entre RPE_début
% et RPE_fin sur les N_bins de la variable X — ce qui suppose une montée
% de fatigue parfaitement constante dans le temps, hypothèse non vérifiée
% et probablement fausse (la fatigue s'accélère souvent en fin de tâche).
% À traiter comme une approximation, pas une mesure, dans l'interprétation.
%
% X_bins, RPE_bins : [N_sujets x N_bins], NaN si sujet/variable absent.
% subj_used : liste des sujets effectivement résolus (même ordre que les lignes).
    N = length(subjects);
    N_bins_loc = NaN;
    X_tmp = cell(N,1);
    RPE0_tmp = NaN(N,1);
    RPE1_tmp = NaN(N,1);
    ok_subj = false(N,1);

    for i = 1:N
        subj = subjects{i};
        try
            if ~isfield(RPE_j2, subj), continue; end
            rpe_debut = RPE_j2.(subj).debut;
            rpe_fin   = RPE_j2.(subj).fin;

            if strcmp(type_j2,'IMU')
                if ~isfield(Results_IMU, subj), continue; end
                kbds = fieldnames(Results_IMU.(subj));
                cond = 'Norm';
                if ~isfield(Results_IMU.(subj), cond)
                    cond_candidates = kbds(contains(kbds,'Norm','IgnoreCase',true));
                    if isempty(cond_candidates), continue; end
                    cond = cond_candidates{1};
                end
                sub = Results_IMU.(subj).(cond);
                if ~isfield(sub, vn), continue; end
                x = sub.(vn)(:)';
            else % EMG
                if isempty(Results_EMG) || ~isfield(Results_EMG, subj), continue; end
                cond = 'Norm';
                if ~isfield(Results_EMG.(subj), cond)
                    kbds2 = fieldnames(Results_EMG.(subj));
                    cond_candidates = kbds2(contains(kbds2,'Norm','IgnoreCase',true));
                    if isempty(cond_candidates), continue; end
                    cond = cond_candidates{1};
                end
                sub = Results_EMG.(subj).(cond);
                if ~isfield(sub, vn), continue; end
                x = sub.(vn)(:)';
            end

            if length(x) < 2, continue; end
            X_tmp{i} = x;
            RPE0_tmp(i) = rpe_debut;
            RPE1_tmp(i) = rpe_fin;
            ok_subj(i) = true;
            if isnan(N_bins_loc), N_bins_loc = length(x); end
        catch
            % Sujet ou variable absent — on laisse NaN
        end
    end

    if isnan(N_bins_loc), N_bins_loc = 10; end  % défaut si aucun sujet résolu

    X_bins   = NaN(N, N_bins_loc);
    RPE_bins = NaN(N, N_bins_loc);
    for i = 1:N
        if ~ok_subj(i), continue; end
        x = X_tmp{i};
        if length(x) ~= N_bins_loc, continue; end  % tailles incohérentes entre sujets
        X_bins(i,:) = x;
        % RPE interpolé linéairement entre début (bin 1) et fin (bin N_bins_loc)
        RPE_bins(i,:) = linspace(RPE0_tmp(i), RPE1_tmp(i), N_bins_loc);
    end
    subj_used = subjects;
end

% -------------------------------------------------------------------------
function [rho, p] = spearman_safe(x, y)
    ok = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
    if sum(ok) < 5
        rho = NaN; p = NaN; return;
    end
    [rho, p] = corr(x(ok), y(ok), 'Type','Spearman');
end

% -------------------------------------------------------------------------
function s = pstar(p)
    if isnan(p),    s = '';    return; end
    if p < 0.001,   s = '***';
    elseif p < 0.01, s = '**';
    elseif p < 0.05, s = '*';
    else,            s = '';
    end
end