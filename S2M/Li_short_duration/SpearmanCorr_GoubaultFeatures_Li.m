%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Spearman Correlation — Goubault Features vs RPE (G1 Liszt)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes Spearman rho between each IMU feature and the RPE
% for each participant (intra-individual), then tests if the
% group median rho is significantly different from 0 (signrank test).
%
% Data source: datpart1_[Feature]_[Signal].mat
%   Structure: dat_part1 (N_obs x 7)
%     Col 1 : Participant ID
%     Col 2 : RPE (interpolated 0→max_RPE over 100 time points)
%     Col 3 : Timestamp
%     Col 4-7: 4 segments (Shoulder, Arm, Forearm, Hand)
%
% Features available:
%   MedianFreq, PeakPower, PeakPower_Freq, Peak,
%   PowerLF, PowerHF, PowerTot, SpectralEntropy
%
% Signals available:
%   Acceleration, Angular_Velocity,
%   Module_Acceleration, Module_Angular_Velocity
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

% 4 segments in columns 4-7
Segments = {'Shoulder', 'Arm', 'Forearm', 'Hand'};

% G1 Liszt participants
G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];

%% -----------------------------------------------------------------------
%  MAIN LOOP
%  -----------------------------------------------------------------------
% Storage for results table
VarNames    = {};
SegNames    = {};
Rho_median  = [];
Pval        = [];
N_valid     = [];

fprintf('%-45s | %-8s | %-8s | %-4s | %s\n', 'Variable', 'Median rho', 'p-value', 'N', 'Sig');
fprintf('%s\n', repmat('-', 1, 80));

for iF = 1:length(Features)
    feat = Features{iF};
    for iS = 1:length(Signals)
        sig = Signals{iS};

        % Build filename
        fname = ['datpart1_' feat '_' sig '.mat'];
        fpath = fullfile(PathData, fname);

        if ~exist(fpath, 'file')
            fprintf('  MISSING: %s\n', fname);
            continue;
        end

        load(fpath);  % -> dat_part1

        % Loop over 4 segments (columns 4-7)
        for iSeg = 1:4
            col = iSeg + 3;  % columns 4,5,6,7

            % Compute Spearman rho per participant
            rho_all = NaN(length(G1_Liszt), 1);
            for iG = 1:length(G1_Liszt)
                iP  = G1_Liszt(iG);
                idx = dat_part1(:,1) == iP;
                if sum(idx) < 3, continue; end

                RPE   = dat_part1(idx, 2);
                Score = dat_part1(idx, col);

                % Remove NaN
                valid = ~isnan(RPE) & ~isnan(Score);
                if sum(valid) < 3, continue; end

                rho_all(iG) = corr(RPE(valid), Score(valid), ...
                    'Type', 'Spearman');
            end

            % Group statistics
            rho_clean = rho_all(~isnan(rho_all));
            if length(rho_clean) < 3
                med_rho = NaN; p_val = NaN; n = 0;
            else
                med_rho = median(rho_clean);
                p_val   = signrank(rho_clean);
                n       = length(rho_clean);
            end

            % Significance marker
            if isnan(p_val),      sig_str = 'NaN';
            elseif p_val < 0.001, sig_str = '***';
            elseif p_val < 0.01,  sig_str = '**';
            elseif p_val < 0.05,  sig_str = '*';
            else,                  sig_str = 'ns';
            end

            % Store
            varname = [feat '_' sig '_' Segments{iSeg}];
            VarNames{end+1}   = varname;
            SegNames{end+1}   = Segments{iSeg};
            Rho_median(end+1) = med_rho;
            Pval(end+1)       = p_val;
            N_valid(end+1)    = n;

            % Print
            if ~strcmp(sig_str, 'ns')
                fprintf('%-45s | %+8.3f | %8.4f | %3d | %s  ←\n', ...
                    varname, med_rho, p_val, n, sig_str);
            end
        end
    end
end

%% -----------------------------------------------------------------------
%  SAVE RESULTS
%  -----------------------------------------------------------------------
T = table(VarNames', Rho_median', Pval', N_valid', ...
    'VariableNames', {'Variable', 'Median_rho', 'p_value', 'N'});

% Sort by absolute rho descending
[~, sortIdx] = sort(abs(Rho_median), 'descend');
T_sorted = T(sortIdx, :);

% Save
writetable(T_sorted, fullfile(PathSave, 'Correlation_Goubault_Features_RPE.xlsx'));
save(fullfile(PathSave, 'Correlation_Goubault_Features_RPE.mat'), ...
    'T_sorted', 'VarNames', 'Rho_median', 'Pval', 'N_valid');

fprintf('\n\nTop 20 variables by |rho| :\n');
fprintf('%-45s | %+8s | %8s | %s\n', 'Variable', 'rho', 'p-value', 'Sig');
fprintf('%s\n', repmat('-',1,75));
for i = 1:min(20, height(T_sorted))
    p = T_sorted.p_value(i);
    if isnan(p),      s='NaN';
    elseif p<0.001,   s='***';
    elseif p<0.01,    s='**';
    elseif p<0.05,    s='*';
    else,             s='ns';
    end
    fprintf('%-45s | %+8.3f | %8.4f | %s\n', ...
        T_sorted.Variable{i}, T_sorted.Median_rho(i), p, s);
end

%% -----------------------------------------------------------------------
%  HEATMAP — rho par feature × segment (meilleur signal)
%  -----------------------------------------------------------------------
% Reshape pour heatmap : Features x Segments x Signals
N_feat = length(Features);
N_sig  = length(Signals);
N_seg  = 4;

RhoMatrix = NaN(N_feat * N_sig, N_seg);
PvalMatrix = NaN(N_feat * N_sig, N_seg);
RowNames = {};

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
            end
        end
    end
end

% Remove rows with all NaN
validRows = any(~isnan(RhoMatrix), 2);
RhoMatrix  = RhoMatrix(validRows, :);
PvalMatrix = PvalMatrix(validRows, :);
RowNames   = RowNames(validRows);

figure('Name', 'Correlation Heatmap - Goubault Features', ...
    'NumberTitle', 'off', 'Position', [50 50 900 1000]);

imagesc(RhoMatrix);
colormap(redblue_cmap());
colorbar;
caxis([-1 1]);

% Add values and significance stars
[nr, nc] = size(RhoMatrix);
for r = 1:nr
    for c = 1:nc
        if ~isnan(RhoMatrix(r,c))
            val = RhoMatrix(r,c);
            pv  = PvalMatrix(r,c);
            if ~isnan(pv) && pv < 0.05
                str = sprintf('%.2f*', val);
                fw  = 'bold';
            else
                str = sprintf('%.2f', val);
                fw  = 'normal';
            end
            text(c, r, str, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 8, ...
                'FontWeight', fw, 'Color', 'k');
        end
    end
end

set(gca, 'XTick', 1:4, 'XTickLabel', Segments, 'FontSize', 9, ...
    'YTick', 1:nr, 'YTickLabel', RowNames, 'TickLabelInterpreter', 'none');
xtickangle(0);
title('Spearman rho — Goubault Features vs RPE (* p<0.05)', ...
    'FontSize', 12, 'FontWeight', 'bold');
xlabel('Segment', 'FontSize', 11);
ylabel('Feature — Signal', 'FontSize', 11);

saveas(gcf, fullfile(PathSave, 'Fig_Correlation_Goubault_Heatmap.png'));
fprintf('\nHeatmap saved. Results saved to %s\n', PathSave);

%% -----------------------------------------------------------------------
%  LOCAL FUNCTION: red-white-blue colormap
%  -----------------------------------------------------------------------
function cmap = redblue_cmap()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    cmap = [r', g', b'];
end