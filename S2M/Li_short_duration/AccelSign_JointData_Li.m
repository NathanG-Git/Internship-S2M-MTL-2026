%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%     Accélération propulsive vs freinatrice — JointData (G1 Li)   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Idée de l'encadrant :
%   Séparer l'accélération articulaire selon son signe par rapport
%   à la vitesse articulaire :
%     - Propulsive  : sign(accel) == sign(velocity)  → accélère le mvt
%     - Freinatrice : sign(accel) != sign(velocity)  → freine le mvt
%
% Pour chaque articulation et chaque cycle musical :
%   - Mean |accel| propulsive
%   - Mean |accel| freinatrice
%   - Ratio propulsif / (propulsif + freinateur)  [0-1]
%
% Données : S10X_Li_JointData.mat (format v7.3 → load via hdf5/matfile)
%   JointData.<joint>.jointAngle    : 3 x N_frames
%   JointData.<joint>.jointVelocity : 3 x N_frames  (rad/s)
%
% Articulations analysées : rGH, rEL, rWR (bras droit — côté joueur)
%
% Output : AccelSign_Li_G1.mat + corrélations Spearman vs RPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  PATHS & PARAMETERS
%  -----------------------------------------------------------------------
PathJoint = 'J:\Piano_Fatigue\Data_Angles\Compute_JointData\CalculatedAngles\';
PathCycle = 'J:\Piano_Fatigue\Data_Exported\';
PathInfo  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
PathSave  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Li_short_duration\';

fs = 60;          % Hz
dt = 1 / fs;

% G1 Liszt participants
G1_Liszt = [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49];

% Articulations d'intérêt (bras droit principal + gauche optionnel)
Joints = {'rGH', 'rEL', 'rWR'};
Axes   = {'x', 'y', 'z'};

%% -----------------------------------------------------------------------
%  LOAD CYCLE TIMESTAMPS & RPE
%  -----------------------------------------------------------------------
load(fullfile(PathCycle, 'Cycle_Li_Felipe.mat'));       % -> cycles
load(fullfile(PathInfo,  'Info_participants_corrected.mat'));  % -> Info_participants

%% -----------------------------------------------------------------------
%  MAIN LOOP
%  -----------------------------------------------------------------------
AccelSign = struct();

for iG = 1:length(G1_Liszt)
    iP = G1_Liszt(iG);
    SubjID = sprintf('S1%02d', iP);

    fprintf('Processing %s ...\n', SubjID);

    %% --- Load JointData (v7.3 → h5read) ---
    fname = fullfile(PathJoint, sprintf('%s_Li_JointData.mat', SubjID));
    if ~exist(fname, 'file')
        fprintf('  MISSING: %s\n', fname);
        continue;
    end

    %% --- RPE vector ---
    RPE_raw = Info_participants(iP).Liszt;

    % Correction S122 (participant index 22 dans G1 = iP=22)
    if iP == 22 && length(RPE_raw) >= 4
        RPE_raw(end) = [];
        fprintf('  S122: RPE corrected (last value removed)\n');
    end
    N_RPE = length(RPE_raw);

    %% --- Cycle timestamps ---
    t_cycles = cycles(iP).seq(:, 1);   % début de chaque cycle (secondes)
    N_cycles = length(t_cycles);

    %% --- Loop over joints ---
    for iJ = 1:length(Joints)
        jnt = Joints{iJ};

        % Load velocity via h5read (fichiers v7.3)
        % h5read en MATLAB inverse les dims HDF5 → retourne directement (N_frames x 3)
        try
            vel_raw = h5read(fname, ['/JointData/' jnt '/jointVelocity']);  % h5read inverse déjà les dims → (N_frames x 3), pas besoin de transposer
        catch
            fprintf('  Joint %s not found in %s\n', jnt, SubjID);
            continue;
        end

        N_frames = size(vel_raw, 1);

        %% --- Compute joint acceleration (numerical derivative) ---
        % Forward difference + padding first row
        accel_raw = [vel_raw(2,:) - vel_raw(1,:); diff(vel_raw, 1, 1)] / dt;
        % Units: rad/s²

        %% --- Low-pass filter (Butterworth 15 Hz, same as pipeline) ---
        [b_filt, a_filt] = butter(2, 2*15/fs);
        min_len = 3 * 2 + 1;  % filtfilt requires at least 3*order+1 = 7 points
        if N_frames < min_len
            fprintf('  Joint %s: signal trop court (%d frames), ignore\n', jnt, N_frames);
            continue;
        end
        vel   = filtfilt(b_filt, a_filt, vel_raw);
        accel = filtfilt(b_filt, a_filt, accel_raw);

        %% --- Propulsive / braking decomposition ---
        % Same sign  → propulsive (segment accelerates in its direction)
        % Opposite   → braking    (segment decelerates)
        same_sign = sign(accel) == sign(vel);   % N_frames x 3, logical

        % Weighted by |accel|
        accel_prop  = abs(accel) .* same_sign;         % propulsive component
        accel_brake = abs(accel) .* ~same_sign;        % braking component

        %% --- Aggregate by cycle ---
        PropMean  = NaN(N_cycles, 3);   % mean |accel| propulsive, per axis
        BrakeMean = NaN(N_cycles, 3);   % mean |accel| braking, per axis
        Ratio     = NaN(N_cycles, 3);   % propulsive ratio [0-1]
        RPE_cycle = NaN(N_cycles, 1);

        for iC = 1:N_cycles
            % Frame indices for this cycle
            t_start = t_cycles(iC);
            if iC < N_cycles
                t_end = t_cycles(iC + 1);
            else
                t_end = (N_frames - 1) / fs;
            end

            idx_start = max(1, round(t_start * fs) + 1);
            idx_end   = min(N_frames, round(t_end * fs));

            if idx_end <= idx_start, continue; end

            % Features
            PropMean(iC, :)  = mean(accel_prop(idx_start:idx_end, :), 1);
            BrakeMean(iC, :) = mean(accel_brake(idx_start:idx_end, :), 1);

            total = PropMean(iC,:) + BrakeMean(iC,:);
            Ratio(iC, :) = PropMean(iC,:) ./ max(total, eps);

            % Associate RPE: window of 30s
            idx_rpe = ceil(t_start / 30);
            idx_rpe = max(1, min(idx_rpe, N_RPE));
            RPE_cycle(iC) = RPE_raw(idx_rpe);
        end

        %% --- Aggregate by unique RPE value ---
        RPE_unique = unique(RPE_cycle(~isnan(RPE_cycle)));
        N_rpe_u = length(RPE_unique);

        PropByRPE  = NaN(N_rpe_u, 3);
        BrakeByRPE = NaN(N_rpe_u, 3);
        RatioByRPE = NaN(N_rpe_u, 3);

        for iR = 1:N_rpe_u
            idx = RPE_cycle == RPE_unique(iR);
            PropByRPE(iR,:)  = nanmean(PropMean(idx,:),  1);
            BrakeByRPE(iR,:) = nanmean(BrakeMean(idx,:), 1);
            RatioByRPE(iR,:) = nanmean(Ratio(idx,:),     1);
        end

        %% --- Store ---
        AccelSign(iG).SubjectID        = SubjID;
        AccelSign(iG).RPE_values       = RPE_unique;
        AccelSign(iG).(jnt).PropMean   = PropByRPE;    % N_RPE x 3 (x,y,z)
        AccelSign(iG).(jnt).BrakeMean  = BrakeByRPE;
        AccelSign(iG).(jnt).Ratio      = RatioByRPE;
    end

    fprintf('  %s done. RPE range: [%d-%d], %d cycles\n', ...
        SubjID, min(RPE_raw), max(RPE_raw), N_cycles);
end

%% -----------------------------------------------------------------------
%  SAVE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'AccelSign_Li_G1.mat'), 'AccelSign', 'Joints', 'Axes');
fprintf('\nAccelSign_Li_G1.mat saved.\n');

%% -----------------------------------------------------------------------
%  SPEARMAN CORRELATION vs RPE
%  -----------------------------------------------------------------------
fprintf('\n--- Spearman rho — AccelSign features vs RPE ---\n');
fprintf('%-35s | %+8s | %8s | %s\n', 'Variable', 'Rho', 'p_raw', 'p_FDR');
fprintf('%s\n', repmat('-', 1, 70));

Features_names = {'PropMean', 'BrakeMean', 'Ratio'};
Axis_names     = {'X', 'Y', 'Z'};

% Collect all rho and p for FDR (3 joints x 3 features x 3 axes = 27 variables)
All_VarNames = {};
All_Rho      = [];
All_Pval     = [];

validG = find(arrayfun(@(s) isfield(s,'SubjectID') && ~isempty(s.SubjectID), AccelSign));

for iJ = 1:length(Joints)
    jnt = Joints{iJ};
    for iF = 1:length(Features_names)
        feat = Features_names{iF};
        for iAx = 1:3
            rho_all = NaN(length(validG), 1);
            for ii = 1:length(validG)
                iG = validG(ii);
                if ~isfield(AccelSign(iG), jnt), continue; end
                if ~isfield(AccelSign(iG).(jnt), feat), continue; end
                RPE   = AccelSign(iG).RPE_values;
                Score = AccelSign(iG).(jnt).(feat)(:, iAx);
                valid = ~isnan(RPE) & ~isnan(Score);
                if sum(valid) < 3, continue; end
                rho_all(ii) = corr(RPE(valid), Score(valid), 'Type', 'Spearman');
            end
            rho_clean = rho_all(~isnan(rho_all));
            if length(rho_clean) < 3
                med_rho = NaN; p_val = NaN;
            else
                med_rho = median(rho_clean);
                p_val   = signrank(rho_clean);
            end
            varname = [jnt '_' feat '_' Axis_names{iAx}];
            All_VarNames{end+1} = varname;
            All_Rho(end+1)      = med_rho;
            All_Pval(end+1)     = p_val;
        end
    end
end

% FDR correction
valid_idx = ~isnan(All_Pval);
[adj_p, h_fdr] = bh_fdr(All_Pval(valid_idx), 0.05);
All_Padj      = NaN(size(All_Pval));
All_H         = false(size(All_Pval));
All_Padj(valid_idx) = adj_p;
All_H(valid_idx)    = h_fdr;

fprintf('\n  FDR: %d / %d variables significatives (q<0.05)\n', ...
    sum(All_H), sum(valid_idx));

for i = 1:length(All_VarNames)
    if isnan(All_Pval(i)), continue; end
    raw_s = sigMarker(All_Pval(i));
    fdr_s = sigMarker(All_Padj(i));
    if ~strcmp(raw_s,'ns') || ~strcmp(fdr_s,'ns')
        fprintf('%-35s | %+8.3f | %8.4f | %s\n', ...
            All_VarNames{i}, All_Rho(i), All_Pval(i), fdr_s);
    end
end

%% -----------------------------------------------------------------------
%  HEATMAP — rho par articulation x feature (axes moyennés)
%  -----------------------------------------------------------------------
% Reshape: rows = joint x feature, cols = axes X Y Z
n_rows = length(Joints) * length(Features_names);
RhoMat  = NaN(n_rows, 3);
PvalMat = NaN(n_rows, 3);
AdjMat  = NaN(n_rows, 3);
RowLbls = {};

iRow = 0;
for iJ = 1:length(Joints)
    for iF = 1:length(Features_names)
        iRow = iRow + 1;
        RowLbls{iRow} = [Joints{iJ} ' — ' Features_names{iF}];
        for iAx = 1:3
            varname = [Joints{iJ} '_' Features_names{iF} '_' Axis_names{iAx}];
            idx = strcmp(All_VarNames, varname);
            if any(idx)
                RhoMat(iRow, iAx)  = All_Rho(idx);
                PvalMat(iRow, iAx) = All_Pval(idx);
                AdjMat(iRow, iAx)  = All_Padj(idx);
            end
        end
    end
end

figure('Name','AccelSign Heatmap','NumberTitle','off','Position',[50 50 700 600]);
imagesc(RhoMat);
colormap(redblue_cmap());
cb = colorbar; cb.Label.String = 'Spearman \rho'; cb.Label.FontSize = 10;
caxis([-1 1]);

[nr, nc] = size(RhoMat);
for r = 1:nr
    for c = 1:nc
        if isnan(RhoMat(r,c)), continue; end
        val  = RhoMat(r,c);
        praw = PvalMat(r,c);
        padj = AdjMat(r,c);
        if ~isnan(padj) && padj < 0.05
            fw = 'bold'; fc = 'k'; marker = '†';
        elseif ~isnan(praw) && praw < 0.05
            fw = 'bold'; fc = [0.3 0.3 0.3]; marker = '*';
        else
            fw = 'normal'; fc = [0.5 0.5 0.5]; marker = '';
        end
        text(c, r-0.15, sprintf('%.2f', val), ...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'FontSize',9,'FontWeight',fw,'Color',fc);
        if ~isempty(marker)
            text(c, r+0.25, marker, ...
                'HorizontalAlignment','center','VerticalAlignment','middle',...
                'FontSize',11,'FontWeight','bold','Color',fc);
        end
    end
end

set(gca,'XTick',1:3,'XTickLabel',Axis_names,'FontSize',10,...
    'YTick',1:nr,'YTickLabel',RowLbls,'TickLabelInterpreter','none');
title({'Spearman \rho — Accel propulsive/freinatrice vs RPE',...
    '\bf†\rm = FDR q<0.05   |   \bf*\rm = raw p<0.05 only'},...
    'FontSize',11,'FontWeight','bold');
xlabel('Axe','FontSize',11);
ylabel('Articulation — Feature','FontSize',11);

saveas(gcf, fullfile(PathSave, 'Fig_AccelSign_Heatmap.png'));
fprintf('\nHeatmap saved.\n');

%% -----------------------------------------------------------------------
%  VISUALIZATION — Ratio propulsif moyen par articulation vs RPE
%  -----------------------------------------------------------------------
RPE_grid  = 1:10;
JntColors = {[0.2 0.5 0.9], [0.9 0.3 0.2], [0.1 0.7 0.3]};  % rGH rEL rWR

figure('Name','Propulsive Ratio vs RPE','NumberTitle','off',...
    'Position',[50 50 1100 400]);

for iJ = 1:length(Joints)
    jnt = Joints{iJ};
    subplot(1, 3, iJ); hold on;

    GroupMatrix = NaN(length(validG), length(RPE_grid));
    for ii = 1:length(validG)
        iG = validG(ii);
        if ~isfield(AccelSign(iG), jnt), continue; end
        RPE   = AccelSign(iG).RPE_values;
        % Mean ratio over 3 axes
        Ratio_mean = nanmean(AccelSign(iG).(jnt).Ratio, 2);
        if length(RPE) < 2, continue; end
        plot(RPE, Ratio_mean, '-', 'Color', [JntColors{iJ} 0.2], 'LineWidth', 0.8);
        for iR = 1:length(RPE_grid)
            if RPE_grid(iR) >= min(RPE) && RPE_grid(iR) <= max(RPE)
                GroupMatrix(ii, iR) = interp1(RPE, Ratio_mean, RPE_grid(iR), 'linear');
            end
        end
    end

    N_contrib = sum(~isnan(GroupMatrix), 1);
    Mu = nanmean(GroupMatrix, 1);
    SD = nanstd(GroupMatrix, 0, 1);
    validIdx_r = N_contrib >= 3;
    RPE_v = RPE_grid(validIdx_r);
    Mu_v  = Mu(validIdx_r);
    SD_v  = SD(validIdx_r);

    if ~isempty(RPE_v)
        fill([RPE_v fliplr(RPE_v)], [Mu_v+SD_v fliplr(Mu_v-SD_v)], ...
            JntColors{iJ}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(RPE_v, Mu_v, 'o-', 'LineWidth', 2.5, 'MarkerSize', 7, ...
            'Color', JntColors{iJ}, 'MarkerFaceColor', JntColors{iJ});
    end

    yline(0.5, '--k', 'LineWidth', 1, 'Alpha', 0.4);  % ligne équilibre prop/brake
    xlabel('RPE (Borg CR-10)', 'FontSize', 10);
    ylabel('Ratio propulsif [0-1]', 'FontSize', 10);
    title(jnt, 'FontSize', 12, 'FontWeight', 'bold', 'Color', JntColors{iJ});
    xlim([0.5 10.5]); ylim([0 1]); grid on; box on;
    hold off;
end

sgtitle('Ratio propulsif (accel dans sens du mvt) vs RPE — G1 Liszt', ...
    'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf, fullfile(PathSave, 'Fig_PropRatio_vsRPE.png'));
fprintf('Figures saved to %s\n', PathSave);

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

function [adj_p, h] = bh_fdr(p_vals, alpha)
% Benjamini-Hochberg FDR — aucun toolbox requis
    n = numel(p_vals);
    [p_sorted, sort_idx] = sort(p_vals(:)');
    adj_sorted = p_sorted;
    for i = n-1:-1:1
        adj_sorted(i) = min(p_sorted(i) * n / i, adj_sorted(i+1));
    end
    adj_sorted = min(adj_sorted, 1);
    adj_p = NaN(1, n);
    adj_p(sort_idx) = adj_sorted;
    h = adj_p < alpha;
end

function cmap = redblue_cmap()
    n = 256;
    r = [linspace(0,1,n/2), ones(1,n/2)];
    g = [linspace(0,1,n/2), linspace(1,0,n/2)];
    b = [ones(1,n/2), linspace(1,0,n/2)];
    cmap = [r', g', b'];
end