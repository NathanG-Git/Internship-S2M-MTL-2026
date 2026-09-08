%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%        Spearman Correlation — Workload Scores vs RPE              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For each variable (Total, ByAxis, BySegMod, BySegAxis),
% for each signal (Acceleration, Jerk),
% for each score (Mean, Median):
%   1. Compute Spearman rho between RPE and score for each participant
%   2. Compute group median rho
%   3. Test if rho significantly != 0 (Wilcoxon signed-rank test)
%
% Output:
%   Results struct + Excel table saved to PathSave
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

% Valid participants
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
%  HELPER: compute Spearman rho for all participants for a given variable
%  Returns rho vector (N_valid x 1), median_rho, p_value (signrank)
%  -----------------------------------------------------------------------
function [rho_all, med_rho, p_val] = computeCorr(Workload, validIdx, ...
                                                   ScoreName, SigName, field, col)
    rho_all = NaN(length(validIdx), 1);
    for i = 1:length(validIdx)
        iG = validIdx(i);
        if ~isfield(Workload(iG), ScoreName), continue; end
        RPE   = Workload(iG).RPE_values;
        Score = Workload(iG).(ScoreName).(SigName).(field)(:, col);
        % Need at least 3 points for meaningful correlation
        if length(RPE) < 3, continue; end
        rho_all(i) = corr(RPE(:), Score(:), 'Type', 'Spearman', 'Rows', 'complete');
    end
    % Remove NaN participants
    rho_clean = rho_all(~isnan(rho_all));
    if length(rho_clean) < 3
        med_rho = NaN;
        p_val   = NaN;
    else
        med_rho = median(rho_clean);
        % Wilcoxon signed-rank test: H0 = rho median == 0
        p_val = signrank(rho_clean);
    end
end

%% -----------------------------------------------------------------------
%  COMPUTE CORRELATIONS
%  -----------------------------------------------------------------------
Results = struct();
VarNames  = {};
MedRho_Mean_Accel  = [];
MedRho_Mean_Jerk   = [];
MedRho_Med_Accel   = [];
MedRho_Med_Jerk    = [];
Pval_Mean_Accel    = [];
Pval_Mean_Jerk     = [];
Pval_Med_Accel     = [];
Pval_Med_Jerk      = [];

for iScore = 1:2
    ScoreName = Scores{iScore};
    for iSig = 1:2
        SigName = Signals{iSig};

        fprintf('\n=== %s %s ===\n', ScoreName, SigName);

        % --- TOTAL ---
        [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'Total', 1);
        Results.(ScoreName).(SigName).Total.rho    = rho;
        Results.(ScoreName).(SigName).Total.median = med;
        Results.(ScoreName).(SigName).Total.pval   = p;
        fprintf('  Total          : rho=%.3f, p=%.4f\n', med, p);
        if iScore==1 && iSig==1, VarNames{end+1}='Total'; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
        if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
        if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
        if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end

        % --- BY AXIS (3) ---
        for iAx = 1:3
            [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'ByAxis', iAx);
            Results.(ScoreName).(SigName).ByAxis(iAx).rho    = rho;
            Results.(ScoreName).(SigName).ByAxis(iAx).median = med;
            Results.(ScoreName).(SigName).ByAxis(iAx).pval   = p;
            fprintf('  Axis %s        : rho=%.3f, p=%.4f\n', Axes_lbl{iAx}, med, p);
            varname = ['Axis_' Axes_lbl{iAx}];
            if iScore==1 && iSig==1, VarNames{end+1}=varname; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
            if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
            if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
            if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end
        end

        % --- BY SEGMENT MODULE (7) ---
        for iSeg = 1:7
            [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'BySegMod', iSeg);
            Results.(ScoreName).(SigName).BySegMod(iSeg).rho    = rho;
            Results.(ScoreName).(SigName).BySegMod(iSeg).median = med;
            Results.(ScoreName).(SigName).BySegMod(iSeg).pval   = p;
            fprintf('  Seg %s module  : rho=%.3f, p=%.4f\n', Segments{iSeg}, med, p);
            varname = ['Mod_' Segments{iSeg}];
            if iScore==1 && iSig==1, VarNames{end+1}=varname; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
            if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
            if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
            if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end
        end

        % --- BY SEGMENT+AXIS (21) ---
        BySegAxis_cols = Workload_Labels.BySegAxis_cols;
        for iCol = 1:21
            [rho, med, p] = computeCorr(Workload, validIdx, ScoreName, SigName, 'BySegAxis', iCol);
            Results.(ScoreName).(SigName).BySegAxis(iCol).rho    = rho;
            Results.(ScoreName).(SigName).BySegAxis(iCol).median = med;
            Results.(ScoreName).(SigName).BySegAxis(iCol).pval   = p;
            fprintf('  %s     : rho=%.3f, p=%.4f\n', BySegAxis_cols{iCol}, med, p);
            varname = BySegAxis_cols{iCol};
            if iScore==1 && iSig==1, VarNames{end+1}=varname; MedRho_Mean_Accel(end+1)=med; Pval_Mean_Accel(end+1)=p; end
            if iScore==1 && iSig==2, MedRho_Mean_Jerk(end+1)=med; Pval_Mean_Jerk(end+1)=p; end
            if iScore==2 && iSig==1, MedRho_Med_Accel(end+1)=med; Pval_Med_Accel(end+1)=p; end
            if iScore==2 && iSig==2, MedRho_Med_Jerk(end+1)=med; Pval_Med_Jerk(end+1)=p; end
        end
    end
end

%% -----------------------------------------------------------------------
%  SUMMARY TABLE
%  -----------------------------------------------------------------------
% Significance markers
function sig = sigMarker(p)
    if isnan(p),     sig = 'NaN';
    elseif p < 0.001, sig = '***';
    elseif p < 0.01,  sig = '**';
    elseif p < 0.05,  sig = '*';
    else,              sig = 'ns';
    end
end

fprintf('\n\n========== SUMMARY TABLE ==========\n');
fprintf('%-20s | Mean_Accel rho (p) | Mean_Jerk rho (p) | Med_Accel rho (p) | Med_Jerk rho (p)\n', 'Variable');
fprintf('%s\n', repmat('-', 1, 100));
for iV = 1:length(VarNames)
    fprintf('%-20s | %+.3f %-3s          | %+.3f %-3s         | %+.3f %-3s         | %+.3f %-3s\n', ...
        VarNames{iV}, ...
        MedRho_Mean_Accel(iV), sigMarker(Pval_Mean_Accel(iV)), ...
        MedRho_Mean_Jerk(iV),  sigMarker(Pval_Mean_Jerk(iV)), ...
        MedRho_Med_Accel(iV),  sigMarker(Pval_Med_Accel(iV)), ...
        MedRho_Med_Jerk(iV),   sigMarker(Pval_Med_Jerk(iV)));
end

%% -----------------------------------------------------------------------
%  SAVE RESULTS + EXCEL TABLE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'Correlation_Workload_RPE.mat'), 'Results', 'VarNames', ...
    'MedRho_Mean_Accel', 'MedRho_Mean_Jerk', 'MedRho_Med_Accel', 'MedRho_Med_Jerk', ...
    'Pval_Mean_Accel',   'Pval_Mean_Jerk',   'Pval_Med_Accel',   'Pval_Med_Jerk');

% Excel table
T = table(VarNames', ...
    MedRho_Mean_Accel', Pval_Mean_Accel', ...
    MedRho_Mean_Jerk',  Pval_Mean_Jerk', ...
    MedRho_Med_Accel',  Pval_Med_Accel', ...
    MedRho_Med_Jerk',   Pval_Med_Jerk', ...
    'VariableNames', {'Variable', ...
        'Mean_Accel_rho', 'Mean_Accel_p', ...
        'Mean_Jerk_rho',  'Mean_Jerk_p', ...
        'Median_Accel_rho','Median_Accel_p', ...
        'Median_Jerk_rho', 'Median_Jerk_p'});

writetable(T, fullfile(PathSave, 'Correlation_Workload_RPE.xlsx'));
fprintf('\nResults saved to %s\n', PathSave);

%% -----------------------------------------------------------------------
%  VISUALIZATION — Heatmap of median rho
%  -----------------------------------------------------------------------
figure('Name', 'Correlation Heatmap', 'NumberTitle', 'off', ...
    'Position', [50 50 1400 600]);

data_heatmap = [MedRho_Mean_Accel; MedRho_Mean_Jerk; ...
                MedRho_Med_Accel;  MedRho_Med_Jerk]';

% Significance mask
sig_mask = [Pval_Mean_Accel; Pval_Mean_Jerk; ...
            Pval_Med_Accel;  Pval_Med_Jerk]' < 0.05;

imagesc(data_heatmap);
colormap(redblue_colormap());
colorbar;
caxis([-1 1]);

% Add significance stars
[rows, cols] = size(data_heatmap);
for r = 1:rows
    for c = 1:cols
        if ~isnan(sig_mask(r,c)) && sig_mask(r,c)
            text(c, r, '*', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', 'white');
        end
        text(c, r-0.3, sprintf('%.2f', data_heatmap(r,c)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontSize', 6, 'Color', 'black');
    end
end

set(gca, 'XTick', 1:4, ...
    'XTickLabel', {'Mean Accel', 'Mean Jerk', 'Median Accel', 'Median Jerk'}, ...
    'YTick', 1:length(VarNames), 'YTickLabel', VarNames, ...
    'FontSize', 9, 'TickLabelInterpreter', 'none');
title('Spearman rho — Workload scores vs RPE (* p<0.05)', ...
    'FontSize', 13, 'FontWeight', 'bold');
xlabel('Score type', 'FontSize', 11);
ylabel('Variable', 'FontSize', 11);

saveas(gcf, fullfile(PathSave, 'Fig8_Correlation_Heatmap.png'));
fprintf('Heatmap saved.\n');

%% -----------------------------------------------------------------------
%  LOCAL FUNCTION: red-blue colormap centered on 0
%  -----------------------------------------------------------------------
function cmap = redblue_colormap()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    cmap = [r', g', b'];
end