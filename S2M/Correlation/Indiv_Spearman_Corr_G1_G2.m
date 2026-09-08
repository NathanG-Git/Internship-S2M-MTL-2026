%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Correlations Spearman INDIVIDUELLES - Toutes variables vs RPE    %%%%
%%%%  Methode : correlation par participant puis moyenne via Fisher     %%%%
%%%%  FDR appliquee sur les p-values du test t sur les z de Fisher     %%%%
%%%%  Sorties : figures PNG Indiv_Spearman_Tab*.png                    %%%%
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

Vi1   = find(d1.valid_subj);
Vi2   = find(d2.valid_subj);
Vi1_26= find(d1.valid_26);
Vi2_26= find(d2.valid_26);
Vi1g  = find(dG1.valid_subj);
Vi2g  = find(dG2.valid_subj);
Vi1e  = find(eG1.valid_subj);
Vi2e  = find(eG2.valid_subj);

Segments = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Axes     = {'X','Y','Z'};
N_seg    = 7;

fprintf('Chargement OK.\n\n');

%% -----------------------------------------------------------------------
%  TABLEAU 1 : Total + BySegMod + ByAxis
%  -----------------------------------------------------------------------
labels1 = {}; data1 = []; npart1 = [];

% Total (Features3, Vi1_26/Vi2_26)
for iSig = 1:2
    if iSig==1, sn='Accel Total'; m1=d1.Mat_Total_Accel_26; m2=d2.Mat_Total_Accel_26;
                rpe1=d1.Mat_RPE_26; rpe2=d2.Mat_RPE_26; Vi_1=Vi1_26; Vi_2=Vi2_26;
    else,        sn='Jerk Total';  m1=d1.Mat_Total_Jerk_26;  m2=d2.Mat_Total_Jerk_26;
                rpe1=d1.Mat_RPE_26; rpe2=d2.Mat_RPE_26; Vi_1=Vi1_26; Vi_2=Vi2_26; end
    [r1,p1,n1_v] = spearIndiv(m1, rpe1, Vi_1);
    [r2,p2,n2_v] = spearIndiv(m2, rpe2, Vi_2);
    labels1{end+1} = sn;
    data1(end+1,:)  = [r1 p1 r2 p2];
    npart1(end+1,:) = [n1_v n2_v];
end

% BySegMod (Features3)
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_Seg_Accel_26; m2=d2.Mat_Seg_Accel_26;
                rpe1=d1.Mat_RPE_26; rpe2=d2.Mat_RPE_26; Vi_1=Vi1_26; Vi_2=Vi2_26;
    else,        sn='Jerk';  m1=d1.Mat_Seg_Jerk_26;  m2=d2.Mat_Seg_Jerk_26;
                rpe1=d1.Mat_RPE_26; rpe2=d2.Mat_RPE_26; Vi_1=Vi1_26; Vi_2=Vi2_26; end
    for iSeg = 1:N_seg
        [r1,p1,n1_v] = spearIndiv(squeeze(m1(:,:,iSeg)), rpe1, Vi_1);
        [r2,p2,n2_v] = spearIndiv(squeeze(m2(:,:,iSeg)), rpe2, Vi_2);
        labels1{end+1} = [sn ' - ' Segments{iSeg}];
        data1(end+1,:)  = [r1 p1 r2 p2];
        npart1(end+1,:) = [n1_v n2_v];
    end
end

% ByAxis (donnees brutes)
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_Axis_Accel; m2=d2.Mat_Axis_Accel;
    else,        sn='Jerk';  m1=d1.Mat_Axis_Jerk;  m2=d2.Mat_Axis_Jerk; end
    for iAx = 1:3
        [r1,p1,n1_v] = spearIndiv(squeeze(m1(:,:,iAx)), d1.Mat_RPE, Vi1);
        [r2,p2,n2_v] = spearIndiv(squeeze(m2(:,:,iAx)), d2.Mat_RPE, Vi2);
        labels1{end+1} = [sn ' - Axe ' Axes{iAx}];
        data1(end+1,:)  = [r1 p1 r2 p2];
        npart1(end+1,:) = [n1_v n2_v];
    end
end
[fdr1_pass, qvals1] = computeFDR(data1);
plotTable(labels1, data1, fdr1_pass, qvals1, npart1, ...
    'Indiv Spearman - Cinematique : Total, BySegMod, ByAxis', ...
    fullfile(PathSave, 'Indiv_Spearman_Tab1_Cinematique.png'));

%% -----------------------------------------------------------------------
%  TABLEAU 2 : BySegAxis
%  -----------------------------------------------------------------------
labels2_X = {}; data2_X = []; npart2_X = [];
labels2_Y = {}; data2_Y = []; npart2_Y = [];
labels2_Z = {}; data2_Z = []; npart2_Z = [];

for iAx = 1:3
    labels2 = {}; data2 = []; npart2 = [];
    for iSig = 1:2
        if iSig==1, sn='Accel'; m1=d1.Mat_SegAxis_Accel; m2=d2.Mat_SegAxis_Accel;
        else,        sn='Jerk';  m1=d1.Mat_SegAxis_Jerk;  m2=d2.Mat_SegAxis_Jerk; end
        for iSeg = 1:N_seg
            [r1,p1,n1_v] = spearIndiv(squeeze(m1(:,:,iSeg,iAx)), d1.Mat_RPE, Vi1);
            [r2,p2,n2_v] = spearIndiv(squeeze(m2(:,:,iSeg,iAx)), d2.Mat_RPE, Vi2);
            labels2{end+1} = [sn ' Axe' Axes{iAx} ' - ' Segments{iSeg}];
            data2(end+1,:)  = [r1 p1 r2 p2];
            npart2(end+1,:) = [n1_v n2_v];
        end
    end
    [fdr2_pass, qvals2] = computeFDR(data2);
    plotTable(labels2, data2, fdr2_pass, qvals2, npart2, ...
        ['Indiv Spearman - Cinematique BySegAxis - Axe ' Axes{iAx}], ...
        fullfile(PathSave, ['Indiv_Spearman_Tab2_BySegAxis_' Axes{iAx} '.png']));
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

labels3 = {}; data3 = []; npart3 = [];
for iV = 1:N_ve
    if strcmp(VarNames_emg{iV}, 'Amplitude'), continue; end
    for iC = 1:N_ch
        m1 = squeeze(eG1.Mat_EMG(:,:,iC,iV));
        m2 = squeeze(eG2.Mat_EMG(:,:,iC,iV));
        [r1,p1,n1_v] = spearIndiv(m1, eG1.Mat_RPE, Vi1e);
        [r2,p2,n2_v] = spearIndiv(m2, eG2.Mat_RPE, Vi2e);
        labels3{end+1} = [VarNames_emg{iV} ' - ' Channels_emg{iC}];
        data3(end+1,:)  = [r1 p1 r2 p2];
        npart3(end+1,:) = [n1_v n2_v];
    end
end
[fdr3_pass, qvals3] = computeFDR(data3);
plotTable(labels3, data3, fdr3_pass, qvals3, npart3, ...
    'Indiv Spearman - Variables EMG', ...
    fullfile(PathSave, 'Indiv_Spearman_Tab3_EMG.png'));

%% -----------------------------------------------------------------------
%  TABLEAU 4 : Goubault
%  -----------------------------------------------------------------------
N_vg = size(dG1.Vars, 1);
labels4 = {}; data4 = []; npart4 = [];
for iV = 1:N_vg
    [r1,p1,n1_v] = spearIndiv(squeeze(dG1.Mat_Vars(:,:,iV)), dG1.Mat_RPE, Vi1g);
    [r2,p2,n2_v] = spearIndiv(squeeze(dG2.Mat_Vars(:,:,iV)), dG2.Mat_RPE, Vi2g);
    labels4{end+1} = dG1.Vars{iV,1};
    data4(end+1,:)  = [r1 p1 r2 p2];
    npart4(end+1,:) = [n1_v n2_v];
end
[fdr4_pass, qvals4] = computeFDR(data4);
plotTable(labels4, data4, fdr4_pass, qvals4, npart4, ...
    'Indiv Spearman - Variables Goubault', ...
    fullfile(PathSave, 'Indiv_Spearman_Tab4_Goubault.png'));

fprintf('\nTermine. Figures PNG sauvegardees dans :\n%s\n', PathSave);

%% -----------------------------------------------------------------------
%  SYNTHESE
%  -----------------------------------------------------------------------
all_labels = [labels1(:); labels2_X(:); labels2_Y(:); labels2_Z(:); labels3(:); labels4(:)];
all_data   = [data1;   data2_X;   data2_Y;   data2_Z;   data3;   data4];
all_fdr    = [fdr1_pass; fdr2_X_pass; fdr2_Y_pass; fdr2_Z_pass; fdr3_pass; fdr4_pass];
all_qvals  = [qvals1;  qvals2_X;  qvals2_Y;  qvals2_Z;  qvals3;  qvals4];
all_npart  = [npart1;  npart2_X;  npart2_Y;  npart2_Z;  npart3;  npart4];

crit_G1 = all_qvals(:,1) < 0.05;
crit_G2 = all_qvals(:,2) < 0.05;
keep    = crit_G1 | crit_G2;
fprintf('\n%d variables significatives (q<0.05 pour G1 ou G2)\n', sum(keep));

%% Sauvegarder les variables individuellement significatives
Indiv_Sig = struct();
Indiv_Sig.labels   = all_labels(keep);
Indiv_Sig.data     = all_data(keep,:);
Indiv_Sig.qvals    = all_qvals(keep,:);
Indiv_Sig.npart    = all_npart(keep,:);
Indiv_Sig.fdr_pass = all_fdr(keep,:);
save(fullfile(PathSave, 'Indiv_Spearman_Results.mat'), 'Indiv_Sig');
fprintf('Variables sauvegardees : Indiv_Spearman_Results.mat\n');

if sum(keep) > 0
    sel_labels = all_labels(keep);
    sel_data   = all_data(keep,:);
    sel_fdr    = all_fdr(keep,:);
    sel_qvals  = all_qvals(keep,:);
    sel_npart  = all_npart(keep,:);

    seg_order = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
    seg_rank  = zeros(length(sel_labels), 1);
    for iL = 1:length(sel_labels)
        lbl = sel_labels{iL};
        ranked = false;
        for iS = 1:length(seg_order)
            if contains(lbl, seg_order{iS})
                seg_rank(iL) = iS; ranked = true; break;
            end
        end
        if ~ranked, seg_rank(iL) = length(seg_order) + 1; end
    end
    [~, sort_idx] = sort(seg_rank, 'ascend');
    sel_labels = sel_labels(sort_idx);
    sel_data   = sel_data(sort_idx,:);
    sel_fdr    = sel_fdr(sort_idx,:);
    sel_qvals  = sel_qvals(sort_idx,:);
    sel_npart  = sel_npart(sort_idx,:);

    plotTable(sel_labels, sel_data, sel_fdr, sel_qvals, sel_npart, ...
        'Indiv Synthese - Variables significatives (q<0.05)', ...
        fullfile(PathSave, 'Indiv_Spearman_Synthese.png'));
    synthese_path = fullfile(PathSave, 'Indiv_Spearman_Synthese.png');
end

%% Ouvrir les figures
figs_to_open = {
    fullfile(PathSave, 'Indiv_Spearman_Tab1_Cinematique.png');
    fullfile(PathSave, 'Indiv_Spearman_Tab2_BySegAxis_X.png');
    fullfile(PathSave, 'Indiv_Spearman_Tab2_BySegAxis_Y.png');
    fullfile(PathSave, 'Indiv_Spearman_Tab2_BySegAxis_Z.png');
    fullfile(PathSave, 'Indiv_Spearman_Tab3_EMG.png');
    fullfile(PathSave, 'Indiv_Spearman_Tab4_Goubault.png');
};
if exist('synthese_path', 'var') && exist(synthese_path, 'file')
    figs_to_open{end+1} = synthese_path;
end
for i = 1:length(figs_to_open)
    try; open(figs_to_open{i}); catch; end
end

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------

function [rho_moy, p_ttest, n_valid] = spearIndiv(Mat_var, Mat_RPE, Vi)
% Calcule la correlation de Spearman individuelle pour chaque participant
% puis moyenne via transformee de Fisher et teste via t-test
%
% Mat_var : N_all x N_bins
% Mat_RPE : N_all x N_bins
% Vi      : indices des participants valides
%
% rho_moy  : rho moyen retransforme depuis Fisher
% p_ttest  : p-value du t-test sur les z de Fisher (H0 : z_moy = 0)
% n_valid  : nombre de participants ayant un rho calculable

    MIN_BINS = 8;  % minimum de bins valides par participant

    N_subj  = length(Vi);
    z_vals  = NaN(N_subj, 1);
    n_valid = 0;

    for iG = 1:N_subj
        iP  = Vi(iG);
        x   = Mat_var(iP, :)';   % N_bins x 1
        rpe = Mat_RPE(iP, :)';

        % Verifier les NaN
        ok  = ~isnan(x) & ~isnan(rpe);
        n_ok = sum(ok);

        % Alerte si trop de NaN
        n_nan = sum(~ok);
        if n_nan > 2
            fprintf('  ATTENTION : participant %d a %d NaN sur %d bins pour cette variable\n', ...
                iP, n_nan, length(x));
        end

        if n_ok < MIN_BINS
            % Pas assez de bins valides — participant exclu
            continue;
        end

        [rho, ~] = corr(x(ok), rpe(ok), 'Type', 'Spearman');

        % Borner rho pour eviter Inf dans Fisher (|rho| = 1 exactement)
        rho = max(min(rho, 0.9999), -0.9999);

        % Transformee de Fisher : z = atanh(rho)
        z_vals(iG) = atanh(rho);
        n_valid = n_valid + 1;
    end

    % Garder uniquement les z valides
    z_ok = z_vals(~isnan(z_vals));
    n_valid = length(z_ok);

    if n_valid < 3
        rho_moy = NaN; p_ttest = NaN; return;
    end

    % Test t a un echantillon : H0 : moyenne(z) = 0
    [~, p_ttest, ~, stats] = ttest(z_ok);

    % Retransformer la moyenne des z en rho moyen
    z_mean  = mean(z_ok);
    rho_moy = tanh(z_mean);
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
    if nargin < 2, alpha = 0.05; end
    N = size(data, 1);
    fdr_pass = false(N, 2);
    qvals    = NaN(N, 2);
    for iGrp = 1:2
        col_p = iGrp * 2;
        pvals = data(:, col_p);
        ok    = ~isnan(pvals);
        if sum(ok) < 2, continue; end
        p_ok  = pvals(ok);
        n_ok  = length(p_ok);
        [p_sorted, sort_idx] = sort(p_ok);
        q_sorted = p_sorted * n_ok ./ (1:n_ok)';
        for k = n_ok-1:-1:1
            q_sorted(k) = min(q_sorted(k), q_sorted(k+1));
        end
        q_sorted = min(q_sorted, 1);
        q_full = NaN(n_ok, 1);
        q_full(sort_idx) = q_sorted;
        idx_ok = find(ok);
        for j = 1:length(idx_ok)
            qvals(idx_ok(j), iGrp) = q_full(j);
        end
        fdr_pass(:, iGrp) = qvals(:, iGrp) < alpha & ok;
    end
end

function plotTable(labels, data, fdr_pass, qvals, n_part, title_str, fpath, n_max_g1, n_max_g2)
    if nargin < 8, n_max_g1 = 26; end
    if nargin < 9, n_max_g2 = 23; end
    N = length(labels);
    col_headers = {'Variable', 'rho G1', 'p G1', 'q G1', 'sig G1', 'n G1', 'sc G1', 'rho G2', 'p G2', 'q G2', 'sig G2', 'n G2', 'sc G2'};

    row_h = max(0.012, min(0.022, 0.82 / max(N,1)));
    fig_h = max(400, 80 + N * 18);
    fig_w = 1300;

    fig = figure('Visible','off','Position',[50 50 fig_w fig_h],'Color','white');
    ax  = axes('Position',[0 0 1 1],'Visible','off');
    ax.XLim = [0 1]; ax.YLim = [0 1];
    hold on;

    col_x = [0.02, 0.19, 0.26, 0.33, 0.40, 0.45, 0.50, 0.56, 0.63, 0.70, 0.76, 0.81, 0.87];

    top_y  = 0.97;
    head_y = top_y - 0.03;

    rectangle('Position',[0 head_y+0.02 1 0.05],'FaceColor',[0.15 0.25 0.45],'EdgeColor','none');
    text(0.5, head_y+0.045, title_str, 'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',11,'FontWeight','bold','Color','white');

    rectangle('Position',[0 head_y-0.005 1 0.03],'FaceColor',[0.85 0.88 0.93],'EdgeColor','none');
    for iC = 1:length(col_headers)
        text(col_x(iC), head_y+0.010, col_headers{iC},'FontSize',8,'FontWeight','bold',...
            'Color',[0.1 0.1 0.3],'VerticalAlignment','middle');
    end
    line([0 1],[head_y-0.005 head_y-0.005],'Color',[0.5 0.5 0.6],'LineWidth',1);

    for iR = 1:N
        y_row = head_y - 0.005 - iR * row_h;
        if mod(iR,2)==0
            rectangle('Position',[0 y_row 1 row_h],'FaceColor',[0.96 0.97 0.99],'EdgeColor','none');
        end

        rho1=data(iR,1); p1=data(iR,2); rho2=data(iR,3); p2=data(iR,4);
        q1=qvals(iR,1);  q2=qvals(iR,2);
        s1=getStars(q1); s2=getStars(q2);
        fdr1=fdr_pass(iR,1); fdr2=fdr_pass(iR,2);

        if fdr1
            rectangle('Position',[col_x(2) y_row col_x(8)-col_x(2) row_h],...
                'FaceColor',getSigColor(q1),'EdgeColor','none');
        end
        if fdr2
            rectangle('Position',[col_x(8) y_row 1-col_x(8) row_h],...
                'FaceColor',getSigColor(q2),'EdgeColor','none');
        end

        y_txt = y_row + row_h*0.5;
        fw1='normal'; if fdr1, fw1='bold'; end
        fw2='normal'; if fdr2, fw2='bold'; end

        % Normaliser les separateurs : tirets longs, underscores -> tiret court ' - '
        lbl_display = labels{iR};
        lbl_display = strrep(lbl_display, '_', ' ');
        lbl_display = strrep(lbl_display, ' — ', ' - ');
        lbl_display = strrep(lbl_display, '—', ' - ');
        lbl_display = strrep(lbl_display, '  ', ' ');
        % Supprimer le prefixe TFR_ -> garder juste le nom de la feature
        lbl_display = strrep(lbl_display, 'TFR ', '');
        text(col_x(1)+0.005, y_txt, lbl_display,'FontSize',7.5,'VerticalAlignment','middle',...
            'FontWeight','normal','Color',[0.1 0.1 0.1],'Interpreter','none');

        if ~isnan(rho1)
            text(col_x(2)+0.013, y_txt, sprintf('%.3f',rho1),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight',fw1,'HorizontalAlignment','center');
        end
        if ~isnan(p1)
            text(col_x(3)+0.013, y_txt, sprintf('%.3f',p1),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight','normal','HorizontalAlignment','center');
        end
        if ~isnan(q1)
            text(col_x(4)+0.013, y_txt, sprintf('%.3f',q1),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight',fw1,'HorizontalAlignment','center');
        end
        text(col_x(5)+0.008, y_txt, s1,'FontSize',7.5,'VerticalAlignment','middle',...
            'FontWeight','bold','HorizontalAlignment','center','Color',[0.7 0.1 0.1]);

        n1_val = n_part(iR,1);
        if ~isnan(n1_val)
            text(col_x(6)+0.010, y_txt, num2str(n1_val),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight','normal','HorizontalAlignment','center',...
                'Color',[0.3 0.3 0.3]);
        end
        if ~isnan(rho1) && ~isnan(n1_val) && n1_val > 0
            sc1 = abs(rho1) * sqrt(n1_val / n_max_g1);
            col_sc = [0.6 0.6 0.6]; if sc1 >= 0.8, col_sc = [0.1 0.5 0.1]; end
            text(col_x(7)+0.010, y_txt, sprintf('%.2f',sc1),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight',fw1,'HorizontalAlignment','center',...
                'Color',col_sc);
        end

        if ~isnan(rho2)
            text(col_x(8)+0.013, y_txt, sprintf('%.3f',rho2),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight',fw2,'HorizontalAlignment','center');
        end
        if ~isnan(p2)
            text(col_x(9)+0.013, y_txt, sprintf('%.3f',p2),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight','normal','HorizontalAlignment','center');
        end
        if ~isnan(q2)
            text(col_x(10)+0.013, y_txt, sprintf('%.3f',q2),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight',fw2,'HorizontalAlignment','center');
        end
        text(col_x(11)+0.008, y_txt, s2,'FontSize',7.5,'VerticalAlignment','middle',...
            'FontWeight','bold','HorizontalAlignment','center','Color',[0.1 0.2 0.7]);

        n2_val = n_part(iR,2);
        if ~isnan(n2_val)
            text(col_x(12)+0.010, y_txt, num2str(n2_val),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight','normal','HorizontalAlignment','center',...
                'Color',[0.3 0.3 0.3]);
        end
        if ~isnan(rho2) && ~isnan(n2_val) && n2_val > 0
            sc2 = abs(rho2) * sqrt(n2_val / n_max_g2);
            col_sc2 = [0.6 0.6 0.6]; if sc2 >= 0.8, col_sc2 = [0.1 0.1 0.6]; end
            text(col_x(13)+0.010, y_txt, sprintf('%.2f',sc2),'FontSize',6.5,...
                'VerticalAlignment','middle','FontWeight',fw2,'HorizontalAlignment','center',...
                'Color',col_sc2);
        end

        line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
    end

    line([0 1 1 0 0],[0 0 1 1 0],'Color',[0.4 0.4 0.5],'LineWidth',1);
    y_leg = head_y - 0.005 - (N+0.8)*row_h;
    text(0.02, y_leg, ['* q<0.05  ** q<0.01  *** q<0.001  |  p = t-test sur z Fisher  |  ' ...
        'Gras = q<0.05  |  sc = |rho_moy|*sqrt(n/n_max)  vert>=0.8  |  G1=rouge G2=bleu'], ...
        'FontSize',7,'Color',[0.3 0.3 0.3],'FontAngle','italic');

    hold off;
    exportgraphics(fig, fpath, 'Resolution', 150);
    close(fig);
    fprintf('  Figure : %s\n', fpath);
end

function c = getSigColor(p)
    if p < 0.001,     c = [1.0 0.88 0.88];
    elseif p < 0.01,  c = [1.0 0.93 0.88];
    else,             c = [1.0 0.98 0.88];
    end
end