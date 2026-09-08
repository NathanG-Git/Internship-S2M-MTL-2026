%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        Spearman Correlation — Workload Scores vs RPE              %%%%
%%%%                WITH FDR CORRECTION (Benjamini-Hochberg)           %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For each variable (Total, ByAxis, BySegMod, BySegAxis),
% for each signal (Acceleration, Jerk),
% for each score (Mean, Median):
%   1. Compute Spearman rho between RPE and score for each participant
%   2. Compute group median rho
%   3. Test if rho significantly != 0 (Wilcoxon signed-rank test)
%
% FDR correction applied per score×signal combination (32 p-values each)
% AND globally across all 128 p-values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  LOAD
%  -----------------------------------------------------------------------
PathSave = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';
load(fullfile(PathSave, 'Workload_Li_G1.mat'));

N_participants = length(Workload);
Scores   = {'Mean', 'Median'};
Signals  = {'Acceleration', 'Jerk'};
Segments = Workload_Labels.BySegMod_cols;
Axes_lbl = Workload_Labels.ByAxis_cols;

validPart = false(N_participants, 1);
for iG = 1:N_participants
    if isfield(Workload(iG), 'SubjectID') && ~isempty(Workload(iG).SubjectID)
        validPart(iG) = true;
    end
end
validIdx = find(validPart);
N_valid  = length(validIdx);
fprintf('%d valid participants\n', N_valid);

%% -----------------------------------------------------------------------
%  HELPERS
%  -----------------------------------------------------------------------
function [rho_all, med_rho, p_val] = computeCorr(Workload, validIdx, ...
                                                   ScoreName, SigName, field, col)
    rho_all = NaN(length(validIdx), 1);
    for i = 1:length(validIdx)
        iG = validIdx(i);
        if ~isfield(Workload(iG), ScoreName), continue; end
        RPE   = Workload(iG).RPE_values;
        Score = Workload(iG).(ScoreName).(SigName).(field)(:, col);
        if length(RPE) < 3, continue; end
        rho_all(i) = corr(RPE(:), Score(:), 'Type', 'Spearman', 'Rows', 'complete');
    end
    rho_clean = rho_all(~isnan(rho_all));
    if length(rho_clean) < 3
        med_rho = NaN; p_val = NaN;
    else
        med_rho = median(rho_clean);
        p_val   = signrank(rho_clean);
    end
end

function sig = sigMarker(p)
    if isnan(p),      sig = 'NaN';
    elseif p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,             sig = 'ns';
    end
end

%% -----------------------------------------------------------------------
%  COMPUTE CORRELATIONS (same as before)
%  -----------------------------------------------------------------------
Results = struct();
VarNames         = {};
MedRho_Mean_Accel = []; MedRho_Mean_Jerk = [];
MedRho_Med_Accel  = []; MedRho_Med_Jerk  = [];
Pval_Mean_Accel   = []; Pval_Mean_Jerk   = [];
Pval_Med_Accel    = []; Pval_Med_Jerk    = [];

for iScore = 1:2
    ScoreName = Scores{iScore};
    for iSig = 1:2
        SigName = Signals{iSig};

        % Total
        [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'Total', 1);
        Results.(ScoreName).(SigName).Total = struct('rho',rho,'median',med,'pval',p);
        if iScore==1 && iSig==1, VarNames{end+1}='Total'; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
        if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
        if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
        if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end

        % By Axis (3)
        for iAx = 1:3
            [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'ByAxis', iAx);
            Results.(ScoreName).(SigName).ByAxis(iAx) = struct('rho',rho,'median',med,'pval',p);
            varname = ['Axis_' Axes_lbl{iAx}];
            if iScore==1 && iSig==1, VarNames{end+1}=varname; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
            if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
            if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
            if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end
        end

        % By Segment Module (7)
        for iSeg = 1:7
            [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'BySegMod', iSeg);
            Results.(ScoreName).(SigName).BySegMod(iSeg) = struct('rho',rho,'median',med,'pval',p);
            varname = ['Mod_' Segments{iSeg}];
            if iScore==1 && iSig==1, VarNames{end+1}=varname; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
            if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
            if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
            if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end
        end

        % By Segment + Axis (21)
        BySegAxis_cols = Workload_Labels.BySegAxis_cols;
        for iCol = 1:21
            [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'BySegAxis', iCol);
            Results.(ScoreName).(SigName).BySegAxis(iCol) = struct('rho',rho,'median',med,'pval',p);
            varname = BySegAxis_cols{iCol};
            if iScore==1 && iSig==1, VarNames{end+1}=varname; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
            if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
            if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
            if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end
        end
    end
end

%% -----------------------------------------------------------------------
%  FDR CORRECTION — applied globally across all 4×32 = 128 p-values
%  -----------------------------------------------------------------------
% Implémentation BH sans toolbox externe (aucune dépendance requise).
all_pvals = [Pval_Mean_Accel, Pval_Mean_Jerk, Pval_Med_Accel, Pval_Med_Jerk];
n_vars = length(Pval_Mean_Accel);   % = 32

valid_idx = ~isnan(all_pvals);

[adj_tmp, h_tmp] = bh_fdr(all_pvals(valid_idx), 0.05);
fprintf('FDR method: Benjamini-Hochberg (implémentation locale, sans toolbox)\n');

all_padj = NaN(size(all_pvals));
all_h    = false(size(all_pvals));
all_padj(valid_idx) = adj_tmp;
all_h(valid_idx)    = h_tmp;

% Split back into 4 groups of 32
Padj_Mean_Accel = all_padj(1          : n_vars);
Padj_Mean_Jerk  = all_padj(n_vars+1   : 2*n_vars);
Padj_Med_Accel  = all_padj(2*n_vars+1 : 3*n_vars);
Padj_Med_Jerk   = all_padj(3*n_vars+1 : 4*n_vars);

H_Mean_Accel = all_h(1          : n_vars);
H_Mean_Jerk  = all_h(n_vars+1   : 2*n_vars);
H_Med_Accel  = all_h(2*n_vars+1 : 3*n_vars);
H_Med_Jerk   = all_h(3*n_vars+1 : 4*n_vars);

fprintf('\n--- FDR summary (128 p-values total) ---\n');
fprintf('  Raw p<0.05 : %d / 128\n', sum(all_pvals(valid_idx) < 0.05));
fprintf('  FDR q<0.05 : %d / 128\n', sum(all_h));

%% -----------------------------------------------------------------------
%  PRINT SUMMARY TABLE
%  -----------------------------------------------------------------------
fprintf('\n%-20s | Mean_Accel         | Mean_Jerk          | Med_Accel          | Med_Jerk\n', 'Variable');
fprintf('%s\n', repmat('-', 1, 105));
for iV = 1:length(VarNames)
    fprintf('%-20s | %+.3f %s/%s       | %+.3f %s/%s      | %+.3f %s/%s      | %+.3f %s/%s\n', ...
        VarNames{iV}, ...
        MedRho_Mean_Accel(iV), sigMarker(Pval_Mean_Accel(iV)), sigMarker(Padj_Mean_Accel(iV)), ...
        MedRho_Mean_Jerk(iV),  sigMarker(Pval_Mean_Jerk(iV)),  sigMarker(Padj_Mean_Jerk(iV)), ...
        MedRho_Med_Accel(iV),  sigMarker(Pval_Med_Accel(iV)),  sigMarker(Padj_Med_Accel(iV)), ...
        MedRho_Med_Jerk(iV),   sigMarker(Pval_Med_Jerk(iV)),   sigMarker(Padj_Med_Jerk(iV)));
end
fprintf('(format: raw/FDR)\n');

%% -----------------------------------------------------------------------
%  SAVE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'Correlation_Workload_RPE_FDR.mat'), ...
    'Results', 'VarNames', ...
    'MedRho_Mean_Accel', 'MedRho_Mean_Jerk', 'MedRho_Med_Accel', 'MedRho_Med_Jerk', ...
    'Pval_Mean_Accel',   'Pval_Mean_Jerk',   'Pval_Med_Accel',   'Pval_Med_Jerk',   ...
    'Padj_Mean_Accel',   'Padj_Mean_Jerk',   'Padj_Med_Accel',   'Padj_Med_Jerk',   ...
    'H_Mean_Accel',      'H_Mean_Jerk',      'H_Med_Accel',      'H_Med_Jerk');

T = table(VarNames', ...
    MedRho_Mean_Accel', Pval_Mean_Accel', Padj_Mean_Accel', H_Mean_Accel', ...
    MedRho_Mean_Jerk',  Pval_Mean_Jerk',  Padj_Mean_Jerk',  H_Mean_Jerk', ...
    MedRho_Med_Accel',  Pval_Med_Accel',  Padj_Med_Accel',  H_Med_Accel', ...
    MedRho_Med_Jerk',   Pval_Med_Jerk',   Padj_Med_Jerk',   H_Med_Jerk', ...
    'VariableNames', {'Variable', ...
        'MA_rho','MA_p_raw','MA_p_adj','MA_FDR', ...
        'MJ_rho','MJ_p_raw','MJ_p_adj','MJ_FDR', ...
        'MdA_rho','MdA_p_raw','MdA_p_adj','MdA_FDR', ...
        'MdJ_rho','MdJ_p_raw','MdJ_p_adj','MdJ_FDR'});
writetable(T, fullfile(PathSave, 'Correlation_Workload_RPE_FDR.xlsx'));

%% -----------------------------------------------------------------------
%  HEATMAP — 4 colonnes (score×signal), 32 lignes (variables)
%  Légende des marqueurs :
%    † gras  = survit à la correction FDR (q<0.05)
%    * gris   = raw p<0.05 uniquement (ne survit pas au FDR)
%    rien     = non significatif
%  -----------------------------------------------------------------------
data_heatmap = [MedRho_Mean_Accel; MedRho_Mean_Jerk; ...
                MedRho_Med_Accel;  MedRho_Med_Jerk]';

praw_heatmap = [Pval_Mean_Accel; Pval_Mean_Jerk; ...
                Pval_Med_Accel;  Pval_Med_Jerk]';

padj_heatmap = [Padj_Mean_Accel; Padj_Mean_Jerk; ...
                Padj_Med_Accel;  Padj_Med_Jerk]';

figure('Name', 'Correlation Heatmap - Workload (FDR)', ...
    'NumberTitle', 'off', 'Position', [50 50 1400 700]);

imagesc(data_heatmap);
colormap(redblue_colormap());
cb = colorbar;
cb.Label.String = 'Spearman \rho';
cb.Label.FontSize = 10;
caxis([-1 1]);

[rows, cols] = size(data_heatmap);
for r = 1:rows
    for c = 1:cols
        if isnan(data_heatmap(r,c)), continue; end
        val  = data_heatmap(r,c);
        praw = praw_heatmap(r,c);
        padj = padj_heatmap(r,c);

        if ~isnan(padj) && padj < 0.05
            fw = 'bold'; fc = 'k';
            marker = '†';
        elseif ~isnan(praw) && praw < 0.05
            fw = 'bold'; fc = [0.3 0.3 0.3];
            marker = '*';
        else
            fw = 'normal'; fc = [0.5 0.5 0.5];
            marker = '';
        end

        text(c, r-0.2, sprintf('%.2f', val), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'FontWeight', fw, 'Color', fc);

        if ~isempty(marker)
            text(c, r+0.25, marker, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', fc);
        end
    end
end

set(gca, ...
    'XTick', 1:4, ...
    'XTickLabel', {'Mean Accel', 'Mean Jerk', 'Median Accel', 'Median Jerk'}, ...
    'YTick', 1:length(VarNames), 'YTickLabel', VarNames, ...
    'FontSize', 9, 'TickLabelInterpreter', 'none');

title({'Spearman \rho — Workload scores vs RPE', ...
       '\bf†\rm = FDR q<0.05   |   \bf*\rm = raw p<0.05 only'}, ...
    'FontSize', 12, 'FontWeight', 'bold');
xlabel('Score type', 'FontSize', 11);
ylabel('Variable', 'FontSize', 11);

saveas(gcf, fullfile(PathSave, 'Fig8_Correlation_Heatmap_FDR.png'));
fprintf('\nHeatmap saved: Fig8_Correlation_Heatmap_FDR.png\n');
fprintf('All results saved to %s\n', PathSave);

%% -----------------------------------------------------------------------
%  LOCAL FUNCTION
%  -----------------------------------------------------------------------
function cmap = redblue_colormap()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    cmap = [r', g', b'];
end

function [adj_p, h] = bh_fdr(p_vals, alpha)
% BH_FDR  Correction FDR de Benjamini-Hochberg — aucun toolbox requis.
%
%   [adj_p, h] = bh_fdr(p_vals, alpha)
%
%   p_vals : vecteur de p-values (pas de NaN — filtrer avant l'appel)
%   alpha  : seuil FDR souhaité (ex: 0.05)
%   adj_p  : q-values (p-values ajustées), même ordre que p_vals
%   h      : logical, 1 = significatif après correction FDR
%
%   Algorithme : tri croissant, ajustement step-up de Hochberg,
%                identique à fdr_bh(...,'pdep') et à scipy fdr_bh.

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