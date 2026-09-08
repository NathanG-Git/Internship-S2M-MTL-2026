%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Spearman Correlation — Goubault Features vs RPE (G1 Liszt)    %%%%
%%%%                  WITH FDR CORRECTION (Benjamini-Hochberg)         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes Spearman rho between each IMU feature and the RPE
% for each participant (intra-individual), then tests if the
% group median rho is significantly different from 0 (signrank test).
% FDR correction applied across all 128 p-values (Benjamini-Hochberg).
%
% Data source: datpart1_[Feature]_[Signal].mat
%   Structure: dat_part1 (N_obs x 7)
%     Col 1 : Participant ID
%     Col 2 : RPE (interpolated 0->max_RPE over 100 time points)
%     Col 3 : Timestamp
%     Col 4-7: 4 segments (Shoulder, Arm, Forearm, Hand)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS
%  -----------------------------------------------------------------------
PathData = 'J:\Piano_Fatigue\Data_Exported\Tables_XSENS_LMM\';
PathSave = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';

%% -----------------------------------------------------------------------
%  DEFINITIONS
%  -----------------------------------------------------------------------
Features = {'MedianFreq', 'PeakPower', 'PeakPower_Freq', 'Peak', ...
            'PowerLF', 'PowerHF', 'PowerTot', 'SpectralEntropy'};

Signals  = {'Acceleration', 'Angular_Velocity', ...
            'Module_Acceleration', 'Module_Angular_Velocity'};

Segments = {'Shoulder', 'Arm', 'Forearm', 'Hand'};

G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];

%% -----------------------------------------------------------------------
%  FIRST PASS: collect all rho and p-values (128 variables)
%  -----------------------------------------------------------------------
VarNames    = {};
SegNames    = {};
Rho_median  = [];
Pval        = [];
N_valid_arr = [];

for iF = 1:length(Features)
    feat = Features{iF};
    for iS = 1:length(Signals)
        sig = Signals{iS};

        fname = ['datpart1_' feat '_' sig '.mat'];
        fpath = fullfile(PathData, fname);

        if ~exist(fpath, 'file')
            fprintf('  MISSING: %s\n', fname);
            % Fill with NaN to preserve indexing
            for iSeg = 1:4
                VarNames{end+1}   = [feat '_' sig '_' Segments{iSeg}];
                SegNames{end+1}   = Segments{iSeg};
                Rho_median(end+1) = NaN;
                Pval(end+1)       = NaN;
                N_valid_arr(end+1) = 0;
            end
            continue;
        end

        load(fpath);  % -> dat_part1

        for iSeg = 1:4
            col = iSeg + 3;

            rho_all = NaN(length(G1_Liszt), 1);
            for iG = 1:length(G1_Liszt)
                iP  = G1_Liszt(iG);
                idx = dat_part1(:,1) == iP;
                if sum(idx) < 3, continue; end
                RPE   = dat_part1(idx, 2);
                Score = dat_part1(idx, col);
                valid = ~isnan(RPE) & ~isnan(Score);
                if sum(valid) < 3, continue; end
                rho_all(iG) = corr(RPE(valid), Score(valid), 'Type', 'Spearman');
            end

            rho_clean = rho_all(~isnan(rho_all));
            if length(rho_clean) < 3
                med_rho = NaN; p_val = NaN; n = 0;
            else
                med_rho = median(rho_clean);
                p_val   = signrank(rho_clean);
                n       = length(rho_clean);
            end

            VarNames{end+1}    = [feat '_' sig '_' Segments{iSeg}];
            SegNames{end+1}    = Segments{iSeg};
            Rho_median(end+1)  = med_rho;
            Pval(end+1)        = p_val;
            N_valid_arr(end+1) = n;
        end
    end
end

%% -----------------------------------------------------------------------
%  FDR CORRECTION (Benjamini-Hochberg) across all 128 p-values
%  -----------------------------------------------------------------------
% Implémentation BH sans toolbox externe (aucune dépendance requise).

% Indices de p-values non-NaN uniquement
valid_pval_idx = find(~isnan(Pval));
p_for_fdr = Pval(valid_pval_idx);

[adj_p, h_fdr] = bh_fdr(p_for_fdr, 0.05);
fprintf('FDR method: Benjamini-Hochberg (implémentation locale, sans toolbox)\n');

% Map back to full 128-length arrays
Pval_adj    = NaN(size(Pval));
H_fdr       = false(size(Pval));
Pval_adj(valid_pval_idx)  = adj_p;
H_fdr(valid_pval_idx)     = h_fdr;

n_raw = sum(Pval(valid_pval_idx) < 0.05);
n_fdr = sum(H_fdr);
fprintf('\n--- FDR summary ---\n');
fprintf('  Significant before FDR (p<0.05) : %d / %d\n', n_raw, length(valid_pval_idx));
fprintf('  Significant after  FDR (q<0.05) : %d / %d\n', n_fdr, length(valid_pval_idx));

%% -----------------------------------------------------------------------
%  PRINT RESULTS
%  -----------------------------------------------------------------------
fprintf('\n%-50s | %+8s | %8s | %8s | %4s | Raw | FDR\n', ...
    'Variable', 'Rho', 'p_raw', 'p_adj', 'N');
fprintf('%s\n', repmat('-', 1, 100));

for i = 1:length(VarNames)
    if isnan(Pval(i)), continue; end
    raw_sig = sigMarker(Pval(i));
    fdr_sig = sigMarker(Pval_adj(i));
    if ~strcmp(raw_sig, 'ns') || ~strcmp(fdr_sig, 'ns')
        fprintf('%-50s | %+8.3f | %8.4f | %8.4f | %3d | %-3s | %s\n', ...
            VarNames{i}, Rho_median(i), Pval(i), Pval_adj(i), ...
            N_valid_arr(i), raw_sig, fdr_sig);
    end
end

%% -----------------------------------------------------------------------
%  SAVE RESULTS
%  -----------------------------------------------------------------------
T = table(VarNames', Rho_median', Pval', Pval_adj', H_fdr', N_valid_arr', ...
    'VariableNames', {'Variable', 'Median_rho', 'p_raw', 'p_adj_FDR', 'Sig_FDR', 'N'});

[~, sortIdx] = sort(abs(Rho_median), 'descend');
T_sorted = T(sortIdx, :);

writetable(T_sorted, fullfile(PathSave, 'Correlation_Goubault_Features_RPE_FDR.xlsx'));
save(fullfile(PathSave, 'Correlation_Goubault_Features_RPE_FDR.mat'), ...
    'T_sorted', 'VarNames', 'Rho_median', 'Pval', 'Pval_adj', 'H_fdr', 'N_valid_arr');
fprintf('\nResults saved.\n');

%% -----------------------------------------------------------------------
%  HEATMAP — now with 3 significance levels:
%    bold + border : survives FDR  (q < 0.05)
%    bold only     : raw p < 0.05, does NOT survive FDR
%    normal        : not significant
%  -----------------------------------------------------------------------
N_feat = length(Features);
N_sig  = length(Signals);
N_seg  = 4;

RhoMatrix  = NaN(N_feat * N_sig, N_seg);
PvalMatrix = NaN(N_feat * N_sig, N_seg);   % raw
AdjMatrix  = NaN(N_feat * N_sig, N_seg);   % FDR-adjusted
RowNames   = {};

iRow = 0;
for iF = 1:N_feat
    for iS = 1:N_sig
        iRow = iRow + 1;
        RowNames{iRow} = [Features{iF} ' — ' Signals{iS}];
        for iSeg = 1:N_seg
            varname = [Features{iF} '_' Signals{iS} '_' Segments{iSeg}];
            idx = strcmp(VarNames, varname);
            if any(idx)
                RhoMatrix(iRow, iSeg)  = Rho_median(idx);
                PvalMatrix(iRow, iSeg) = Pval(idx);
                AdjMatrix(iRow, iSeg)  = Pval_adj(idx);
            end
        end
    end
end

validRows = any(~isnan(RhoMatrix), 2);
RhoMatrix  = RhoMatrix(validRows, :);
PvalMatrix = PvalMatrix(validRows, :);
AdjMatrix  = AdjMatrix(validRows, :);
RowNames   = RowNames(validRows);

figure('Name', 'Correlation Heatmap - Goubault Features (FDR)', ...
    'NumberTitle', 'off', 'Position', [50 50 980 1050]);

imagesc(RhoMatrix);
colormap(redblue_cmap());
cb = colorbar;
cb.Label.String = 'Spearman \rho';
cb.Label.FontSize = 10;
caxis([-1 1]);

[nr, nc] = size(RhoMatrix);
for r = 1:nr
    for c = 1:nc
        if isnan(RhoMatrix(r,c)), continue; end
        val  = RhoMatrix(r,c);
        praw = PvalMatrix(r,c);
        padj = AdjMatrix(r,c);

        % Determine significance level
        if ~isnan(padj) && padj < 0.05
            % Survives FDR: bold value + star marker
            str = sprintf('\\bf%.2f^{†}', val);
            fw  = 'bold';
            fc  = 'k';
        elseif ~isnan(praw) && praw < 0.05
            % Raw only: bold value + different marker
            str = sprintf('\\bf%.2f*', val);
            fw  = 'bold';
            fc  = [0.3 0.3 0.3];
        else
            str = sprintf('%.2f', val);
            fw  = 'normal';
            fc  = [0.5 0.5 0.5];
        end

        text(c, r, sprintf('%.2f', val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 8, 'FontWeight', fw, 'Color', fc);

        % Add marker below value for significance
        if ~isnan(padj) && padj < 0.05
            text(c, r+0.3, '†', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', 'k');
        elseif ~isnan(praw) && praw < 0.05
            text(c, r+0.3, '*', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);
        end
    end
end

set(gca, 'XTick', 1:4, 'XTickLabel', Segments, 'FontSize', 9, ...
    'YTick', 1:nr, 'YTickLabel', RowNames, 'TickLabelInterpreter', 'none');
xtickangle(0);

title({'Spearman \rho — Goubault Features vs RPE', ...
       '\bf†\rm = FDR q<0.05   |   \bf*\rm = raw p<0.05 only (does not survive FDR)'}, ...
    'FontSize', 11, 'FontWeight', 'bold');
xlabel('Segment', 'FontSize', 11);
ylabel('Feature — Signal', 'FontSize', 11);

% Draw grid lines between cells
hold on;
for r = 0.5:1:nr+0.5
    plot([0.5 nc+0.5], [r r], 'k-', 'LineWidth', 0.3);
end
for c = 0.5:1:nc+0.5
    plot([c c], [0.5 nr+0.5], 'k-', 'LineWidth', 0.3);
end
hold off;

saveas(gcf, fullfile(PathSave, 'Fig_Correlation_Goubault_Heatmap_FDR.png'));
fprintf('Heatmap saved: Fig_Correlation_Goubault_Heatmap_FDR.png\n');

%% -----------------------------------------------------------------------
%  LOCAL FUNCTIONS
%  -----------------------------------------------------------------------
function sig = sigMarker(p)
    if isnan(p),      sig = 'NaN';
    elseif p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,             sig = 'ns';
    end
end

function cmap = redblue_cmap()
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

    % Calcul des q-values step-up (de la plus grande à la plus petite)
    adj_sorted = p_sorted;
    adj_sorted(n) = p_sorted(n);
    for i = n-1:-1:1
        adj_sorted(i) = min(p_sorted(i) * n / i, adj_sorted(i+1));
    end
    adj_sorted = min(adj_sorted, 1);  % plafonner à 1

    % Remettre dans l'ordre original
    adj_p = NaN(1, n);
    adj_p(sort_idx) = adj_sorted;

    h = adj_p < alpha;
end