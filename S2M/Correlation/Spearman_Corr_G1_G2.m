%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Correlations Spearman - Toutes variables vs RPE                  %%%%
%%%%  Sorties : figures PNG avec tableaux visuels                       %%%%
%%%%  Significativite : * p<0.05 | ** p<0.01 | *** p<0.001             %%%%
%%%%  Valeurs significatives en GRAS dans le tableau                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS
%  -----------------------------------------------------------------------
PathSave  = fileparts(mfilename('fullpath'));

PathG1_IMU = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';
PathG2_IMU = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_long_duration\';
PathG1_EMG = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\1_EMG\';
PathG2_EMG = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_long_duration\2_EMG\';

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
fprintf('Chargement des donnees...\n');

d1  = load(fullfile(PathG1_IMU, 'Workload_Li_TimeNormalised.mat'));
d2  = load(fullfile(PathG2_IMU, 'Workload_Li_TimeNormalised_G2.mat'));
dG1 = load(fullfile(PathG1_IMU, 'Workload_Li_TimeNormalised_Goubault.mat'));
dG2 = load(fullfile(PathG2_IMU, 'Workload_Li_TimeNormalised_Goubault_G2.mat'));
eG1 = load(fullfile(PathG1_EMG, 'EMG_Li_TimeNormalised_G1.mat'));
eG2 = load(fullfile(PathG2_EMG, 'EMG_Li_TimeNormalised_G2.mat'));

Vi1 = find(d1.valid_subj);   RPE_G1 = nanmean(d1.Mat_RPE(Vi1,:), 1);
Vi2 = find(d2.valid_subj);   RPE_G2 = nanmean(d2.Mat_RPE(Vi2,:), 1);
Vi1g = find(dG1.valid_subj); RPE_G1g = nanmean(dG1.Mat_RPE(Vi1g,:), 1);
Vi2g = find(dG2.valid_subj); RPE_G2g = nanmean(dG2.Mat_RPE(Vi2g,:), 1);
Vi1e = find(eG1.valid_subj); RPE_G1e = nanmean(eG1.Mat_RPE(Vi1e,:), 1);
Vi2e = find(eG2.valid_subj); RPE_G2e = nanmean(eG2.Mat_RPE(Vi2e,:), 1);

Segments = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Axes     = {'X','Y','Z'};
N_seg    = 7;

fprintf('Chargement OK.\n\n');

%% -----------------------------------------------------------------------
%  TABLEAU 1 : Total + BySegMod + ByAxis
%  -----------------------------------------------------------------------
labels1 = {};
data1   = [];
npart1  = [];

% --- Vi26 : participants valides pour Features3 ---
Vi1_26 = find(d1.valid_26);
Vi2_26 = find(d2.valid_26);
RPE_G1_26 = nanmean(d1.Mat_RPE_26(Vi1_26,:), 1);
RPE_G2_26 = nanmean(d2.Mat_RPE_26(Vi2_26,:), 1);

% Total (somme des normes de segments, Features3, n=26/25)
for iSig = 1:2
    if iSig==1, sn='Accel Total'; m1=d1.Mat_Total_Accel_26; m2=d2.Mat_Total_Accel_26; r1_rpe=RPE_G1_26; r2_rpe=RPE_G2_26;
    else,        sn='Jerk Total';  m1=d1.Mat_Total_Jerk_26;  m2=d2.Mat_Total_Jerk_26;  r1_rpe=RPE_G1_26; r2_rpe=RPE_G2_26; end
    n1_v = sum(~all(isnan(m1(Vi1_26,:)),2));
    n2_v = sum(~all(isnan(m2(Vi2_26,:)),2));
    [r1,p1]=spear(nanmean(m1(Vi1_26,:),1), r1_rpe);
    [r2,p2]=spear(nanmean(m2(Vi2_26,:),1), r2_rpe);
    labels1{end+1} = sn;
    data1(end+1,:)  = [r1 p1 r2 p2];
    npart1(end+1,:) = [n1_v n2_v];
end

% BySegMod (Features3, n=26/25)
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_Seg_Accel_26; m2=d2.Mat_Seg_Accel_26; r1_rpe=RPE_G1_26; r2_rpe=RPE_G2_26;
    else,        sn='Jerk';  m1=d1.Mat_Seg_Jerk_26;  m2=d2.Mat_Seg_Jerk_26;  r1_rpe=RPE_G1_26; r2_rpe=RPE_G2_26; end
    for iSeg=1:N_seg
        n1_v = sum(~all(isnan(squeeze(m1(Vi1_26,:,iSeg))),2));
        n2_v = sum(~all(isnan(squeeze(m2(Vi2_26,:,iSeg))),2));
        [r1,p1]=spear(nanmean(squeeze(m1(Vi1_26,:,iSeg)),1), r1_rpe);
        [r2,p2]=spear(nanmean(squeeze(m2(Vi2_26,:,iSeg)),1), r2_rpe);
        labels1{end+1} = [sn ' - ' Segments{iSeg}];
        data1(end+1,:)  = [r1 p1 r2 p2];
        npart1(end+1,:) = [n1_v n2_v];
    end
end

% ByAxis
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_Axis_Accel; m2=d2.Mat_Axis_Accel;
    else,        sn='Jerk';  m1=d1.Mat_Axis_Jerk;  m2=d2.Mat_Axis_Jerk; end
    for iAx=1:3
        n1_v = sum(~all(isnan(squeeze(m1(Vi1,:,iAx))),2));
        n2_v = sum(~all(isnan(squeeze(m2(Vi2,:,iAx))),2));
        [r1,p1]=spear(nanmean(squeeze(m1(Vi1,:,iAx)),1), RPE_G1);
        [r2,p2]=spear(nanmean(squeeze(m2(Vi2,:,iAx)),1), RPE_G2);
        labels1{end+1} = [sn ' - Axe ' Axes{iAx}];
        data1(end+1,:)  = [r1 p1 r2 p2];
        npart1(end+1,:) = [n1_v n2_v];
    end
end
[fdr1_pass, qvals1] = computeFDR(data1);
plotTable(labels1, data1, fdr1_pass, qvals1, npart1, ...
    'Spearman - Cinematique : Total, BySegMod, ByAxis', ...
    fullfile(PathSave, 'Spearman_Tab1_Cinematique.png'));

%% -----------------------------------------------------------------------
%  TABLEAU 2 : BySegAxis - une figure par axe
%  -----------------------------------------------------------------------
labels2_X = {}; data2_X = []; npart2_X = [];
labels2_Y = {}; data2_Y = []; npart2_Y = [];
labels2_Z = {}; data2_Z = []; npart2_Z = [];

for iAx = 1:3
    labels2 = {};
    data2   = [];
    npart2  = [];
    for iSig = 1:2
        if iSig==1, sn='Accel'; m1=d1.Mat_SegAxis_Accel; m2=d2.Mat_SegAxis_Accel;
        else,        sn='Jerk';  m1=d1.Mat_SegAxis_Jerk;  m2=d2.Mat_SegAxis_Jerk; end
        for iSeg=1:N_seg
            n1_v=sum(~all(isnan(squeeze(m1(Vi1,:,iSeg,iAx))),2));
            n2_v=sum(~all(isnan(squeeze(m2(Vi2,:,iSeg,iAx))),2));
            [r1,p1]=spear(nanmean(squeeze(m1(Vi1,:,iSeg,iAx)),1), RPE_G1);
            [r2,p2]=spear(nanmean(squeeze(m2(Vi2,:,iSeg,iAx)),1), RPE_G2);
            labels2{end+1} = [sn ' Axe' Axes{iAx} ' - ' Segments{iSeg}];
            data2(end+1,:)  = [r1 p1 r2 p2];
            npart2(end+1,:) = [n1_v n2_v];
        end
    end
    [fdr2_pass, qvals2] = computeFDR(data2);
    plotTable(labels2, data2, fdr2_pass, qvals2, npart2, ...
        ['Spearman - Cinematique BySegAxis - Axe ' Axes{iAx}], ...
        fullfile(PathSave, ['Spearman_Tab2_BySegAxis_' Axes{iAx} '.png']));
    % Stocker pour la synthese
    if iAx==1, labels2_X=labels2; data2_X=data2; fdr2_X_pass=fdr2_pass; qvals2_X=qvals2; npart2_X=npart2;
    elseif iAx==2, labels2_Y=labels2; data2_Y=data2; fdr2_Y_pass=fdr2_pass; qvals2_Y=qvals2; npart2_Y=npart2;
    else,          labels2_Z=labels2; data2_Z=data2; fdr2_Z_pass=fdr2_pass; qvals2_Z=qvals2; npart2_Z=npart2;
    end
end

%% -----------------------------------------------------------------------
%  TABLEAU 3 : EMG
%  -----------------------------------------------------------------------
VarNames_emg = eG1.VarNames;
Channels_emg = eG1.Channels;
N_ve = length(VarNames_emg);
N_ch = length(Channels_emg);

labels3 = {};
data3   = [];
npart3  = [];
for iV = 1:N_ve
    % Exclure la variable Amplitude
    if strcmp(VarNames_emg{iV}, 'Amplitude'), continue; end
    for iC = 1:N_ch
        x1 = squeeze(eG1.Mat_EMG(Vi1e,:,iC,iV));
        x2 = squeeze(eG2.Mat_EMG(Vi2e,:,iC,iV));
        n1_v = sum(~all(isnan(x1),2));
        n2_v = sum(~all(isnan(x2),2));
        [r1,p1] = spear(nanmean(x1,1), RPE_G1e);
        [r2,p2] = spear(nanmean(x2,1), RPE_G2e);
        labels3{end+1} = [VarNames_emg{iV} ' - ' Channels_emg{iC}];
        data3(end+1,:)  = [r1 p1 r2 p2];
        npart3(end+1,:) = [n1_v n2_v];
    end
end
[fdr3_pass, qvals3] = computeFDR(data3);
plotTable(labels3, data3, fdr3_pass, qvals3, npart3, ...
    'Spearman - Variables EMG', ...
    fullfile(PathSave, 'Spearman_Tab3_EMG.png'));

%% -----------------------------------------------------------------------
%  TABLEAU 4 : Goubault
%  -----------------------------------------------------------------------
N_vg = size(dG1.Vars, 1);
labels4 = {};
data4   = [];
npart4  = [];
for iV = 1:N_vg
    x1 = squeeze(dG1.Mat_Vars(Vi1g,:,iV));
    x2 = squeeze(dG2.Mat_Vars(Vi2g,:,iV));
    n1_v = sum(~all(isnan(x1),2));
    n2_v = sum(~all(isnan(x2),2));
    [r1,p1] = spear(nanmean(x1,1), RPE_G1g);
    [r2,p2] = spear(nanmean(x2,1), RPE_G2g);
    labels4{end+1} = dG1.Vars{iV,1};
    data4(end+1,:)  = [r1 p1 r2 p2];
    npart4(end+1,:) = [n1_v n2_v];
end
[fdr4_pass, qvals4] = computeFDR(data4);
plotTable(labels4, data4, fdr4_pass, qvals4, npart4, ...
    'Spearman - Variables Goubault', ...
    fullfile(PathSave, 'Spearman_Tab4_Goubault.png'));

fprintf('\nTermine. Figures PNG sauvegardees dans :\n%s\n', PathSave);

%% -----------------------------------------------------------------------
%  FIGURE SYNTHESE : variables avec q<0.01 ET |rho|>=0.9
%  Si G1 OU G2 remplit la condition, on affiche les deux groupes
%  -----------------------------------------------------------------------

% Collecter toutes les donnees de tous les tableaux
all_labels = [labels1(:); labels2_X(:); labels2_Y(:); labels2_Z(:); labels3(:); labels4(:)];
all_data   = [data1;   data2_X;   data2_Y;   data2_Z;   data3;   data4];
all_fdr    = [fdr1_pass; fdr2_X_pass; fdr2_Y_pass; fdr2_Z_pass; fdr3_pass; fdr4_pass];
all_qvals  = [qvals1;  qvals2_X;  qvals2_Y;  qvals2_Z;  qvals3;  qvals4];
all_npart  = [npart1;  npart2_X;  npart2_Y;  npart2_Z;  npart3;  npart4];

% Critere : q < 0.01 (2 etoiles) ET |rho| >= 0.9 pour G1 OU G2
% Garder toutes les variables significatives (q<0.05) pour G1 OU G2
crit_G1 = all_qvals(:,1) < 0.05;
crit_G2 = all_qvals(:,2) < 0.05;
keep    = crit_G1 | crit_G2;

sel_labels = all_labels(keep);
sel_data   = all_data(keep,:);
sel_fdr    = all_fdr(keep,:);
sel_qvals  = all_qvals(keep,:);
sel_npart  = all_npart(keep,:);

fprintf('\n%d variables significatives (q<0.05 pour G1 ou G2)\n', sum(keep));

if sum(keep) > 0
    % Tri par segment dans l'ordre anatomique proximal -> distal
    seg_order = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
    seg_rank  = zeros(length(sel_labels), 1);
    for iL = 1:length(sel_labels)
        lbl = sel_labels{iL};
        ranked = false;
        for iS = 1:length(seg_order)
            if contains(lbl, seg_order{iS})
                seg_rank(iL) = iS;
                ranked = true;
                break;
            end
        end
        if ~ranked
            seg_rank(iL) = length(seg_order) + 1;  % EMG/Goubault sans segment -> fin
        end
    end
    [~, sort_idx] = sort(seg_rank, 'ascend');
    sel_labels = sel_labels(sort_idx);
    sel_data   = sel_data(sort_idx, :);
    sel_fdr    = sel_fdr(sort_idx, :);
    sel_qvals  = sel_qvals(sort_idx, :);
    sel_npart  = sel_npart(sort_idx, :);

    plotTable(sel_labels, sel_data, sel_fdr, sel_qvals, sel_npart, ...
        'Synthese - Variables significatives (q<0.05)', ...
        fullfile(PathSave, 'Spearman_Synthese.png'));
    synthese_path = fullfile(PathSave, 'Spearman_Synthese.png');
else
    fprintf('Aucune variable ne remplit le critere.\n');
end

%% Ouvrir les figures automatiquement%% Ouvrir les figures automatiquement
figs_to_open = {
    fullfile(PathSave, 'Spearman_Tab1_Cinematique.png');
    fullfile(PathSave, 'Spearman_Tab2_BySegAxis_X.png');
    fullfile(PathSave, 'Spearman_Tab2_BySegAxis_Y.png');
    fullfile(PathSave, 'Spearman_Tab2_BySegAxis_Z.png');
    fullfile(PathSave, 'Spearman_Tab3_EMG.png');
    fullfile(PathSave, 'Spearman_Tab4_Goubault.png');
};
if exist('synthese_path', 'var') && exist(synthese_path, 'file')
    figs_to_open{end+1} = synthese_path;
end
for i = 1:length(figs_to_open)
    open(figs_to_open{i});
end

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------

function [rho, pval] = spear(x, y)
    x = x(:); y = y(:);
    ok = ~isnan(x) & ~isnan(y);
    if sum(ok) < 4
        rho = NaN; pval = NaN; return;
    end
    [rho, pval] = corr(x(ok), y(ok), 'Type', 'Spearman');
end

function stars = getStars(p)
    if isnan(p),      stars = '';
    elseif p < 0.001, stars = '***';
    elseif p < 0.01,  stars = '**';
    elseif p < 0.05,  stars = '*';
    else,             stars = '';
    end
end


function [fdr_pass, qvals] = computeFDR(data, alpha)
% Correction Benjamini-Hochberg (FDR) appliquee separement sur G1 et G2
% data     : N x 4 = [rho_G1, p_G1, rho_G2, p_G2]
% fdr_pass : N x 2 logical [G1_pass, G2_pass]
% qvals    : N x 2 q-values ajustees [q_G1, q_G2]
    if nargin < 2, alpha = 0.05; end
    N = size(data, 1);
    fdr_pass = false(N, 2);
    qvals    = NaN(N, 2);
    for iGrp = 1:2
        col_p = iGrp * 2;  % col 2 = p_G1, col 4 = p_G2
        pvals = data(:, col_p);
        ok    = ~isnan(pvals);
        if sum(ok) < 2, continue; end
        p_ok  = pvals(ok);
        n_ok  = length(p_ok);
        [p_sorted, sort_idx] = sort(p_ok);
        % Q-values BH : q(i) = p(i) * n / rang(i), en partant du rang le plus eleve
        q_sorted = p_sorted * n_ok ./ (1:n_ok)';
        % Assurer monotonie decroissante (cummin depuis la fin)
        for k = n_ok-1:-1:1
            q_sorted(k) = min(q_sorted(k), q_sorted(k+1));
        end
        q_sorted = min(q_sorted, 1);  % borner a 1
        % Repasser dans l'ordre original
        q_full = NaN(n_ok, 1);
        q_full(sort_idx) = q_sorted;
        idx_ok = find(ok);
        for j = 1:length(idx_ok)
            qvals(idx_ok(j), iGrp) = q_full(j);
        end
        % fdr_pass = q < alpha
        fdr_pass(:, iGrp) = qvals(:, iGrp) < alpha & ok;
    end
end

function plotTable(labels, data, fdr_pass, qvals, n_part, title_str, fpath, n_max_g1, n_max_g2)
% data : N x 4 = [rho_G1, p_G1, rho_G2, p_G2]
% Produit une figure PNG avec tableau visuel

    if nargin < 8, n_max_g1 = 26; end
    if nargin < 9, n_max_g2 = 23; end
    N = length(labels);
    col_headers = {'Variable', 'rho G1', 'p G1', 'q G1', 'sig G1', 'n G1', 'sc G1', 'rho G2', 'p G2', 'q G2', 'sig G2', 'n G2', 'sc G2'};
    N_cols = length(col_headers);

    % Hauteur adaptee au nombre de lignes
    row_h  = max(0.012, min(0.022, 0.82 / max(N,1)));  % adaptatif selon N
    fig_h  = max(400, 80 + N * 18);
    fig_w  = 1050;

    fig = figure('Visible','off','Position',[50 50 fig_w fig_h],...
        'Color','white');
    ax = axes('Position',[0 0 1 1],'Visible','off');
    ax.XLim = [0 1]; ax.YLim = [0 1];
    hold on;

    % Couleurs colonnes
    col_x = [0.02, 0.19, 0.26, 0.33, 0.40, 0.45, 0.50, 0.56, 0.63, 0.70, 0.76, 0.81, 0.87];
    % col_x: 1=Var 2=rhoG1 3=pG1 4=qG1 5=sigG1 6=nG1 7=scG1 8=rhoG2 9=pG2 10=qG2 11=sigG2 12=nG2 13=scG2

    top_y   = 0.97;
    head_y  = top_y - 0.03;

    % Fond titre
    rectangle('Position',[0 head_y+0.02 1 0.05],...
        'FaceColor',[0.15 0.25 0.45],'EdgeColor','none');
    text(0.5, head_y+0.045, title_str, ...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',11,'FontWeight','bold','Color','white');

    % Fond en-tete colonnes
    rectangle('Position',[0 head_y-0.005 1 0.03],...
        'FaceColor',[0.85 0.88 0.93],'EdgeColor','none');
    col_colors_head = [0.85 0.88 0.93];
    for iC = 1:N_cols
        text(col_x(iC), head_y+0.010, col_headers{iC}, ...
            'FontSize',8,'FontWeight','bold','Color',[0.1 0.1 0.3],...
            'VerticalAlignment','middle');
    end

    % Ligne separatrice
    line([0 1],[head_y-0.005 head_y-0.005],'Color',[0.5 0.5 0.6],'LineWidth',1);

    % Lignes de donnees
    for iR = 1:N
        y_row = head_y - 0.005 - iR * row_h;

        % Fond alterne
        if mod(iR,2)==0
            rectangle('Position',[0 y_row 1 row_h],...
                'FaceColor',[0.96 0.97 0.99],'EdgeColor','none');
        end

        rho1 = data(iR,1); p1 = data(iR,2);
        rho2 = data(iR,3); p2 = data(iR,4);
        q1   = qvals(iR,1); q2 = qvals(iR,2);
        % Etoiles et gras bases sur la q-value
        s1   = getStars(q1);
        s2   = getStars(q2);
        fdr1 = fdr_pass(iR, 1);  % q1 < 0.05
        fdr2 = fdr_pass(iR, 2);  % q2 < 0.05

        % Fond colore si q significatif (base q-value)
        if fdr1
            c1 = getSigColor(q1);
            rectangle('Position',[col_x(2) y_row col_x(8)-col_x(2) row_h],...
                'FaceColor',c1,'EdgeColor','none');
        end
        if fdr2
            c2 = getSigColor(q2);
            rectangle('Position',[col_x(8) y_row 1-col_x(8) row_h],...
                'FaceColor',c2,'EdgeColor','none');
        end

        y_txt = y_row + row_h*0.5;
        % Gras si q < 0.05 (passe FDR)
        fw1 = 'normal'; if fdr1, fw1='bold'; end
        fw2 = 'normal'; if fdr2, fw2='bold'; end

        % Colonne Variable
        text(col_x(1)+0.005, y_txt, labels{iR}, ...
            'FontSize',7.5,'VerticalAlignment','middle',...
            'FontWeight','normal','Color',[0.1 0.1 0.1],...
            'Interpreter','none');

        % rho G1
        if ~isnan(rho1)
            text(col_x(2)+0.013, y_txt, sprintf('%.3f', rho1), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight',fw1,'HorizontalAlignment','center');
        end
        % p G1
        if ~isnan(p1)
            text(col_x(3)+0.013, y_txt, sprintf('%.3f', p1), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight','normal','HorizontalAlignment','center');
        end
        % q G1
        if ~isnan(q1)
            text(col_x(4)+0.013, y_txt, sprintf('%.3f', q1), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight',fw1,'HorizontalAlignment','center');
        end
        % sig G1
        text(col_x(5)+0.008, y_txt, s1, ...
            'FontSize',7.5,'VerticalAlignment','middle',...
            'FontWeight','bold','HorizontalAlignment','center',...
            'Color',[0.7 0.1 0.1]);
        % n G1
        n1 = n_part(iR, 1);
        if ~isnan(n1)
            text(col_x(6)+0.010, y_txt, num2str(n1), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight','normal','HorizontalAlignment','center',...
                'Color',[0.3 0.3 0.3]);
        end
        % score G1 = |rho| x sqrt(n/n_max)
        if ~isnan(rho1) && ~isnan(n1) && n1 > 0
            sc1 = abs(rho1) * sqrt(n1 / n_max_g1);
            col_sc = [0.6 0.6 0.6];
            if sc1 >= 0.8, col_sc = [0.1 0.5 0.1]; end
            text(col_x(7)+0.010, y_txt, sprintf('%.2f', sc1), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight',fw1,'HorizontalAlignment','center',...
                'Color', col_sc);
        end

        % rho G2
        if ~isnan(rho2)
            text(col_x(8)+0.013, y_txt, sprintf('%.3f', rho2), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight',fw2,'HorizontalAlignment','center');
        end
        % p G2
        if ~isnan(p2)
            text(col_x(9)+0.013, y_txt, sprintf('%.3f', p2), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight','normal','HorizontalAlignment','center');
        end
        % q G2
        if ~isnan(q2)
            text(col_x(10)+0.013, y_txt, sprintf('%.3f', q2), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight',fw2,'HorizontalAlignment','center');
        end
        % sig G2
        text(col_x(11)+0.008, y_txt, s2, ...
            'FontSize',7.5,'VerticalAlignment','middle',...
            'FontWeight','bold','HorizontalAlignment','center',...
            'Color',[0.1 0.2 0.7]);
        % n G2
        n2 = n_part(iR, 2);
        if ~isnan(n2)
            text(col_x(12)+0.010, y_txt, num2str(n2), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight','normal','HorizontalAlignment','center',...
                'Color',[0.3 0.3 0.3]);
        end
        % score G2 = |rho| x sqrt(n/n_max)
        if ~isnan(rho2) && ~isnan(n2) && n2 > 0
            sc2 = abs(rho2) * sqrt(n2 / n_max_g2);
            col_sc2 = [0.6 0.6 0.6];
            if sc2 >= 0.8, col_sc2 = [0.1 0.1 0.6]; end
            text(col_x(13)+0.010, y_txt, sprintf('%.2f', sc2), ...
                'FontSize',6.5,'VerticalAlignment','middle',...
                'FontWeight',fw2,'HorizontalAlignment','center',...
                'Color', col_sc2);
        end

        % Ligne separatrice legere
        line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
    end

    % Bordure exterieure
    line([0 1 1 0 0],[0 0 1 1 0],'Color',[0.4 0.4 0.5],'LineWidth',1);

    % Legende significativite
    y_leg = head_y - 0.005 - (N+0.8)*row_h;
    text(0.02, y_leg, '* q<0.05  ** q<0.01  *** q<0.001  |  Gras = q<0.05  |  sc = |rho|*sqrt(n/n_max)  vert>=0.8  |  G1=rouge G2=bleu', ...
        'FontSize',7,'Color',[0.3 0.3 0.3],'FontAngle','italic');

    hold off;
    exportgraphics(fig, fpath, 'Resolution', 150);
    close(fig);
    fprintf('  Figure : %s\n', fpath);
end

function c = getSigColor(p)
    if p < 0.001,     c = [1.0 0.88 0.88];   % rouge pale fort
    elseif p < 0.01,  c = [1.0 0.93 0.88];   % orange pale
    else,             c = [1.0 0.98 0.88];   % jaune pale
    end
end