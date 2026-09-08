%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%           Visualization - Workload Analysis Li (G1)               %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Style: individual curves transparent in background + group mean±SD
%
% Fig 1 : Total agrégé       — 2x2 (Mean/Median x Accel/Jerk)
% Fig 2 : Par axe agrégé     — 2x3 (Mean/Median x X/Y/Z) — Accel
% Fig 3 : Par axe agrégé     — 2x3 (Mean/Median x X/Y/Z) — Jerk
% Fig 4 : Par segment module — 2x4 (Mean/Median x 7 segs) — Accel
% Fig 5 : Par segment module — 2x4 (Mean/Median x 7 segs) — Jerk
% Fig 6 : Par segment+axe    — 2x7 (Mean/Median x 7 segs, 3 axes superposés) — Accel
% Fig 7 : Par segment+axe    — 2x7 (Mean/Median x 7 segs, 3 axes superposés) — Jerk
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
Segments = Workload_Labels.BySegMod_cols;  % 7 segments
Axes_lbl = Workload_Labels.ByAxis_cols;    % X Y Z

% Valid participants
validPart = false(N_participants, 1);
SubjNames = cell(N_participants, 1);
for iG = 1:N_participants
    if isfield(Workload(iG), 'SubjectID') && ~isempty(Workload(iG).SubjectID)
        validPart(iG) = true;
        SubjNames{iG} = Workload(iG).SubjectID;
    end
end
N_valid = sum(validPart);
fprintf('%d valid participants\n', N_valid);

% Colors
cmap      = lines(N_valid);   % one color per valid participant
SegColors = lines(7);         % one color per segment
AxColors  = {[0.2 0.4 0.8], [0.8 0.2 0.2], [0.1 0.7 0.3]};  % X=blue Y=red Z=green
RPE_grid  = 1:10;

%% -----------------------------------------------------------------------
%  HELPER: compute group mean±SD on RPE_grid from Workload data
%  Returns GroupMean, GroupSD, RPE_valid (where N>=3 participants)
%  -----------------------------------------------------------------------
function [RPE_v, Mu, SD] = groupStats(Workload, validPart, RPE_grid, ScoreName, SigName, field, col)
    N = length(Workload);
    GroupMatrix = NaN(N, length(RPE_grid));
    for iG = 1:N
        if ~validPart(iG), continue; end
        if ~isfield(Workload(iG), ScoreName), continue; end
        RPE   = Workload(iG).RPE_values;
        Score = Workload(iG).(ScoreName).(SigName).(field)(:, col);
        for iR = 1:length(RPE_grid)
            if RPE_grid(iR) >= min(RPE) && RPE_grid(iR) <= max(RPE)
                GroupMatrix(iG,iR) = interp1(RPE, Score, RPE_grid(iR), 'linear');
            end
        end
    end
    N_contrib = sum(~isnan(GroupMatrix), 1);
    Mu_all    = nanmean(GroupMatrix, 1);
    SD_all    = nanstd(GroupMatrix, 0, 1);
    validIdx  = N_contrib >= 3;
    RPE_v = RPE_grid(validIdx);
    Mu    = Mu_all(validIdx);
    SD    = SD_all(validIdx);
end

%% -----------------------------------------------------------------------
%  HELPER: plot individual curves + group mean±SD on current axes
%  -----------------------------------------------------------------------
function plotGroupStyle(Workload, validPart, cmap, RPE_grid, ...
                        ScoreName, SigName, field, col, lineColor, SubjNames)
    N      = length(Workload);
    validG = find(validPart);
    GroupMatrix = NaN(N, length(RPE_grid));
    for iG = 1:N
        if ~validPart(iG), continue; end
        if ~isfield(Workload(iG), ScoreName), continue; end
        RPE   = Workload(iG).RPE_values;
        Score = Workload(iG).(ScoreName).(SigName).(field)(:, col);
        % Individual curve — transparent
        cidx = find(validG == iG, 1);
        if isempty(cidx), cidx = 1; end
        plot(RPE, Score, '-', 'Color', [cmap(cidx,:) 0.2], 'LineWidth', 0.8);
        for iR = 1:length(RPE_grid)
            if RPE_grid(iR) >= min(RPE) && RPE_grid(iR) <= max(RPE)
                GroupMatrix(iG,iR) = interp1(RPE, Score, RPE_grid(iR), 'linear');
            end
        end
    end
    % Group mean ± SD
    N_contrib = sum(~isnan(GroupMatrix), 1);
    Mu_all    = nanmean(GroupMatrix, 1);
    SD_all    = nanstd(GroupMatrix, 0, 1);
    validIdx  = N_contrib >= 3;
    RPE_v = RPE_grid(validIdx);
    Mu    = Mu_all(validIdx);
    SD    = SD_all(validIdx);
    if ~isempty(RPE_v)
        fill([RPE_v fliplr(RPE_v)], [Mu+SD fliplr(Mu-SD)], ...
            lineColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(RPE_v, Mu, 'o-', 'LineWidth', 2.5, 'MarkerSize', 7, ...
            'Color', lineColor, 'MarkerFaceColor', lineColor);
    end
end

%% -----------------------------------------------------------------------
%  FIG 1 : TOTAL — 2x2 (Mean/Median x Accel/Jerk)
%  -----------------------------------------------------------------------
figure('Name','Fig1 - Total','NumberTitle','off','Position',[50 50 1200 700]);
for iScore = 1:2
    for iSig = 1:2
        subplot(2, 2, (iScore-1)*2 + iSig);
        hold on;
        plotGroupStyle(Workload, validPart, cmap, RPE_grid, ...
            Scores{iScore}, Signals{iSig}, 'Total', 1, [0 0 0], SubjNames);
        xlabel('RPE (Borg CR-10)', 'FontSize', 10);
        ylabel([Scores{iScore} ' Total'], 'FontSize', 10);
        title([Scores{iScore} ' — Total ' Signals{iSig}], 'FontSize', 11);
        xlim([0.5 10.5]); grid on; box on;
        hold off;
    end
end
sgtitle('Total Workload Score — Mean & Median', 'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, fullfile(PathSave, 'Fig1_Total.png'));

%% -----------------------------------------------------------------------
%  FIG 2-3 : PAR AXE — 2x3 (Mean/Median x X/Y/Z) — une figure par signal
%  -----------------------------------------------------------------------
for iSig = 1:2
    figure('Name', ['Fig' num2str(iSig+1) ' - ByAxis - ' Signals{iSig}], ...
        'NumberTitle', 'off', 'Position', [50 50 1400 700]);
    for iScore = 1:2
        for iAx = 1:3
            subplot(2, 3, (iScore-1)*3 + iAx);
            hold on;
            plotGroupStyle(Workload, validPart, cmap, RPE_grid, ...
                Scores{iScore}, Signals{iSig}, 'ByAxis', iAx, AxColors{iAx}, SubjNames);
            xlabel('RPE', 'FontSize', 10);
            ylabel([Scores{iScore} ' Axis ' Axes_lbl{iAx}], 'FontSize', 10);
            title([Scores{iScore} ' — Axis ' Axes_lbl{iAx}], 'FontSize', 11, ...
                'Color', AxColors{iAx});
            xlim([0.5 10.5]); grid on; box on;
            hold off;
        end
    end
    sgtitle([Signals{iSig} ' By Axis (sum over 7 segments) — Mean & Median'], ...
        'FontSize', 13, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig' num2str(iSig+1) '_ByAxis_' Signals{iSig} '.png']));
end

%% -----------------------------------------------------------------------
%  FIG 4-5 : PAR SEGMENT MODULE — 2x4 (Mean/Median x 7 segs)
%  une figure par signal
%  -----------------------------------------------------------------------
for iSig = 1:2
    figure('Name', ['Fig' num2str(iSig+3) ' - BySegment - ' Signals{iSig}], ...
        'NumberTitle', 'off', 'Position', [50 50 1600 700]);
    for iScore = 1:2
        for iSeg = 1:7
            % Layout: 2 rows x 4 cols, seg 5-7 in row2 cols 1-3, col4 empty
            if iSeg <= 4
                spIdx = (iScore-1)*4 + iSeg;
            else
                spIdx = (iScore-1)*4 + iSeg;
            end
            % Use 2x7 layout to fit all 7 segments cleanly
            subplot(2, 7, (iScore-1)*7 + iSeg);
            hold on;
            plotGroupStyle(Workload, validPart, cmap, RPE_grid, ...
                Scores{iScore}, Signals{iSig}, 'BySegMod', iSeg, SegColors(iSeg,:), SubjNames);
            xlabel('RPE', 'FontSize', 8);
            ylabel([Scores{iScore} ' module'], 'FontSize', 8);
            title([Scores{iScore} ' — ' Segments{iSeg}], 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', SegColors(iSeg,:));
            xlim([0.5 10.5]); grid on; box on;
            hold off;
        end
    end
    sgtitle([Signals{iSig} ' By Segment (module) — Mean & Median'], ...
        'FontSize', 13, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig' num2str(iSig+3) '_BySegment_' Signals{iSig} '.png']));
end

%% -----------------------------------------------------------------------
%  FIG 6-7 : PAR SEGMENT+AXE — 2x7 subplots, 3 axes superposés
%  une figure par signal
%  -----------------------------------------------------------------------
for iSig = 1:2
    figure('Name', ['Fig' num2str(iSig+5) ' - BySegAxis - ' Signals{iSig}], ...
        'NumberTitle', 'off', 'Position', [30 30 1800 700]);
    for iScore = 1:2
        for iSeg = 1:7
            subplot(2, 7, (iScore-1)*7 + iSeg);
            hold on;
            for iAx = 1:3
                col_idx = (iSeg-1)*3 + iAx;
                plotGroupStyle(Workload, validPart, cmap, RPE_grid, ...
                    Scores{iScore}, Signals{iSig}, 'BySegAxis', col_idx, ...
                    AxColors{iAx}, SubjNames);
            end
            xlabel('RPE', 'FontSize', 8);
            ylabel(Scores{iScore}, 'FontSize', 8);
            title([Scores{iScore} ' ' Segments{iSeg}], 'FontSize', 9, ...
                'FontWeight', 'bold', 'Color', SegColors(iSeg,:));
            xlim([0.5 10.5]); grid on; box on;
            if iSeg == 7 && iScore == 1
                legend({'','X','','Y','','Z'}, 'Location', 'best', 'FontSize', 7);
            end
            hold off;
        end
    end
    sgtitle([Signals{iSig} ' By Segment & Axis — Mean & Median (Blue=X, Red=Y, Green=Z)'], ...
        'FontSize', 12, 'FontWeight', 'bold');
    saveas(gcf, fullfile(PathSave, ['Fig' num2str(iSig+5) '_BySegAxis_' Signals{iSig} '.png']));
end

fprintf('All figures saved as PNG to %s\n', PathSave);