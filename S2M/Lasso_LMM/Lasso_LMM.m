%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  LASSO + LMM - Niveau individuel                                  %%%%
%%%%  Variables : toutes celles significatives (q<0.05) en Spearman   %%%%
%%%%  G1 (26 participants) | G2 (23 participants)                      %%%%
%%%%                                                                   %%%%
%%%%  Pipeline :                                                        %%%%
%%%%  1. Charger les .mat pipeline                                      %%%%
%%%%  2. Construire tableau long individuel (participant x bin)         %%%%
%%%%  3. LASSO separement G1 et G2 avec CV 5-fold                     %%%%
%%%%  4. LMM sur variables selectionnees (effet aleatoire participant)  %%%%
%%%%  5. Figure coefficients beta G1 vs G2 avec sc annote             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;
rng(42);  % Graine fixe pour reproductibilite des resultats LASSO

%% -----------------------------------------------------------------------
%  PATHS
%  -----------------------------------------------------------------------
PathSave   = fileparts(mfilename('fullpath'));
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

% Participants valides
Vi1     = find(d1.valid_subj);    % brut ~14 G1
Vi2     = find(d2.valid_subj);    % brut ~15 G2
Vi1_26  = find(d1.valid_26);      % Features3 ~25 G1
Vi2_26  = find(d2.valid_26);      % Features3 ~23 G2
Vi1g    = find(dG1.valid_subj);   % Goubault G1
Vi2g    = find(dG2.valid_subj);   % Goubault G2
Vi1e    = find(eG1.valid_subj);   % EMG G1
Vi2e    = find(eG2.valid_subj);   % EMG G2

N_bins   = 10;
Segments = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Axes_lbl = {'X','Y','Z'};
N_seg    = 7;
Channels = eG1.Channels;
VarNames_emg = eG1.VarNames;
N_ch     = length(Channels);
N_max_G1 = 26;
N_max_G2 = 23;

fprintf('G1: brut n=%d | Features3 n=%d | Goubault n=%d | EMG n=%d\n',...
    length(Vi1), length(Vi1_26), length(Vi1g), length(Vi1e));
fprintf('G2: brut n=%d | Features3 n=%d | Goubault n=%d | EMG n=%d\n',...
    length(Vi2), length(Vi2_26), length(Vi2g), length(Vi2e));

%% -----------------------------------------------------------------------
%  DEFINITION DES BLOCS DE VARIABLES
%  Bloc A : Features3 BySegMod + Total (n=25/23) - Vi1_26 / Vi2_26
%  Bloc B : BySegAxis brut (n=14/15)             - Vi1    / Vi2
%  Bloc C : Goubault (n=25/23)                   - Vi1g   / Vi2g
%  Bloc D : EMG sans Amplitude (n=26/23)          - Vi1e   / Vi2e
%  -----------------------------------------------------------------------

%% -----------------------------------------------------------------------
%  SECTION 1 : CONSTRUCTION DES TABLEAUX LONGS PAR BLOC
%  -----------------------------------------------------------------------
fprintf('\n=== Construction des tableaux longs ===\n');

%% --- Bloc A : Features3 BySegMod + Total ---
var_names_A = {}; X_A_G1 = []; X_A_G2 = [];
RPE_A_G1 = []; RPE_A_G2 = []; subj_A_G1 = []; subj_A_G2 = [];

for iSig = 1:2
    if iSig==1
        sn='Total'; m1=d1.Mat_Total_Accel_26; m2=d2.Mat_Total_Accel_26;
        jn='TotalJerk'; j1=d1.Mat_Total_Jerk_26; j2=d2.Mat_Total_Jerk_26;
    end
end
% Total Accel
[x1,r1,s1] = buildBlock(Vi1_26, d1.Mat_Total_Accel_26, d1.Mat_RPE_26, N_bins);
[x2,r2,s2] = buildBlock(Vi2_26, d2.Mat_Total_Accel_26, d2.Mat_RPE_26, N_bins);
var_names_A{end+1}='Total_Accel'; X_A_G1=[X_A_G1 x1]; X_A_G2=[X_A_G2 x2];
RPE_A_G1=r1; RPE_A_G2=r2; subj_A_G1=s1; subj_A_G2=s2;

% Total Jerk
[x1,~,~] = buildBlock(Vi1_26, d1.Mat_Total_Jerk_26, d1.Mat_RPE_26, N_bins);
[x2,~,~] = buildBlock(Vi2_26, d2.Mat_Total_Jerk_26, d2.Mat_RPE_26, N_bins);
var_names_A{end+1}='Total_Jerk'; X_A_G1=[X_A_G1 x1]; X_A_G2=[X_A_G2 x2];

% BySegMod
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_Seg_Accel_26; m2=d2.Mat_Seg_Accel_26;
    else,        sn='Jerk';  m1=d1.Mat_Seg_Jerk_26;  m2=d2.Mat_Seg_Jerk_26; end
    for iSeg = 1:N_seg
        [x1,~,~] = buildBlock(Vi1_26, squeeze(m1(:,:,iSeg)), d1.Mat_RPE_26, N_bins);
        [x2,~,~] = buildBlock(Vi2_26, squeeze(m2(:,:,iSeg)), d2.Mat_RPE_26, N_bins);
        var_names_A{end+1} = [sn '_' Segments{iSeg}];
        X_A_G1=[X_A_G1 x1]; X_A_G2=[X_A_G2 x2];
    end
end
fprintf('  Bloc A (Features3) : %d variables, G1 n=%d obs, G2 n=%d obs\n',...
    length(var_names_A), size(X_A_G1,1), size(X_A_G2,1));

%% --- Bloc B : BySegAxis brut ---
var_names_B = {}; X_B_G1 = []; X_B_G2 = [];
RPE_B_G1 = []; RPE_B_G2 = []; subj_B_G1 = []; subj_B_G2 = [];
first_B = true;
for iSig = 1:2
    if iSig==1, sn='Accel'; m1=d1.Mat_SegAxis_Accel; m2=d2.Mat_SegAxis_Accel;
    else,        sn='Jerk';  m1=d1.Mat_SegAxis_Jerk;  m2=d2.Mat_SegAxis_Jerk; end
    for iSeg = 1:N_seg
        for iAx = 1:3
            [x1,r1,s1] = buildBlock(Vi1, squeeze(m1(:,:,iSeg,iAx)), d1.Mat_RPE, N_bins);
            [x2,r2,s2] = buildBlock(Vi2, squeeze(m2(:,:,iSeg,iAx)), d2.Mat_RPE, N_bins);
            var_names_B{end+1} = [sn '_Axe' Axes_lbl{iAx} '_' Segments{iSeg}];
            X_B_G1=[X_B_G1 x1]; X_B_G2=[X_B_G2 x2];
            if first_B
                RPE_B_G1=r1; RPE_B_G2=r2; subj_B_G1=s1; subj_B_G2=s2;
                first_B=false;
            end
        end
    end
end
fprintf('  Bloc B (BySegAxis brut) : %d variables, G1 n=%d obs, G2 n=%d obs\n',...
    length(var_names_B), size(X_B_G1,1), size(X_B_G2,1));

%% --- Bloc C : Goubault ---
var_names_C = {}; X_C_G1 = []; X_C_G2 = [];
RPE_C_G1 = []; RPE_C_G2 = []; subj_C_G1 = []; subj_C_G2 = [];
first_C = true;
for iV = 1:size(dG1.Vars, 1)
    [x1,r1,s1] = buildBlock(Vi1g, squeeze(dG1.Mat_Vars(:,:,iV)), dG1.Mat_RPE, N_bins);
    [x2,r2,s2] = buildBlock(Vi2g, squeeze(dG2.Mat_Vars(:,:,iV)), dG2.Mat_RPE, N_bins);
    vn = strrep(dG1.Vars{iV,1},' ','_');
    vn = strrep(vn,'-','_'); vn = strrep(vn,'—','_');
    var_names_C{end+1} = ['Goub_' vn];
    X_C_G1=[X_C_G1 x1]; X_C_G2=[X_C_G2 x2];
    if first_C, RPE_C_G1=r1; RPE_C_G2=r2; subj_C_G1=s1; subj_C_G2=s2; first_C=false; end
end
fprintf('  Bloc C (Goubault) : %d variables, G1 n=%d obs, G2 n=%d obs\n',...
    length(var_names_C), size(X_C_G1,1), size(X_C_G2,1));

%% --- Bloc D : EMG sans Amplitude ---
var_names_D = {}; X_D_G1 = []; X_D_G2 = [];
RPE_D_G1 = []; RPE_D_G2 = []; subj_D_G1 = []; subj_D_G2 = [];
first_D = true;
for iV = 1:length(VarNames_emg)
    if strcmp(VarNames_emg{iV},'Amplitude'), continue; end
    for iC = 1:N_ch
        [x1,r1,s1] = buildBlock(Vi1e, squeeze(eG1.Mat_EMG(:,:,iC,iV)), eG1.Mat_RPE, N_bins);
        [x2,r2,s2] = buildBlock(Vi2e, squeeze(eG2.Mat_EMG(:,:,iC,iV)), eG2.Mat_RPE, N_bins);
        var_names_D{end+1} = [VarNames_emg{iV} '_' Channels{iC}];
        X_D_G1=[X_D_G1 x1]; X_D_G2=[X_D_G2 x2];
        if first_D, RPE_D_G1=r1; RPE_D_G2=r2; subj_D_G1=s1; subj_D_G2=s2; first_D=false; end
    end
end
fprintf('  Bloc D (EMG) : %d variables, G1 n=%d obs, G2 n=%d obs\n',...
    length(var_names_D), size(X_D_G1,1), size(X_D_G2,1));

%% -----------------------------------------------------------------------
%  SECTION 2 : LASSO PAR BLOC
%  Standardisation + LASSO CV 5-fold pour chaque bloc et chaque groupe
%  -----------------------------------------------------------------------
fprintf('\n=== LASSO par bloc ===\n');

blocs = {'A','B','C','D'};
bloc_data = {
    X_A_G1, RPE_A_G1, subj_A_G1, var_names_A;
    X_B_G1, RPE_B_G1, subj_B_G1, var_names_B;
    X_C_G1, RPE_C_G1, subj_C_G1, var_names_C;
    X_D_G1, RPE_D_G1, subj_D_G1, var_names_D;
};
bloc_data_G2 = {
    X_A_G2, RPE_A_G2, subj_A_G2, var_names_A;
    X_B_G2, RPE_B_G2, subj_B_G2, var_names_B;
    X_C_G2, RPE_C_G2, subj_C_G2, var_names_C;
    X_D_G2, RPE_D_G2, subj_D_G2, var_names_D;
};

sel_vars_G1 = {}; sel_coef_G1 = []; sel_bloc_G1 = {};
sel_vars_G2 = {}; sel_coef_G2 = []; sel_bloc_G2 = {};

for iB = 1:4
    for iGrp = 1:2
        if iGrp==1
            X = bloc_data{iB,1}; Y = bloc_data{iB,2};
            vnames = bloc_data{iB,4}; gname = 'G1';
        else
            X = bloc_data_G2{iB,1}; Y = bloc_data_G2{iB,2};
            vnames = bloc_data_G2{iB,4}; gname = 'G2';
        end

        % Enlever NaN et Inf
        ok = ~any(isnan(X)|isinf(X),2) & ~isnan(Y) & ~isinf(Y);
        if sum(ok) < 10
            fprintf('  Bloc %s %s : pas assez de donnees (%d obs), skip\n', blocs{iB}, gname, sum(ok));
            continue;
        end
        Xok = X(ok,:); Yok = Y(ok);

        % Standardiser
        mu_X = mean(Xok); sd_X = std(Xok); sd_X(sd_X==0)=1;
        Xs = (Xok - mu_X) ./ sd_X;

        fprintf('  LASSO Bloc %s %s (%d obs, %d vars)...\n', blocs{iB}, gname, sum(ok), size(Xs,2));
        try
            [B, FitInfo] = lasso(Xs, Yok, 'CV', 5, 'Alpha', 1);
            idx_opt = FitInfo.Index1SE;
            coef    = B(:, idx_opt);
            sel     = coef ~= 0;
            fprintf('    => %d variables selectionnees (lambda=%.4f)\n', sum(sel), FitInfo.Lambda(idx_opt));

            if iGrp==1
                sel_vars_G1 = [sel_vars_G1; vnames(sel)'];
                sel_coef_G1 = [sel_coef_G1; coef(sel)];
                sel_bloc_G1 = [sel_bloc_G1; repmat({blocs{iB}}, sum(sel), 1)];
            else
                sel_vars_G2 = [sel_vars_G2; vnames(sel)'];
                sel_coef_G2 = [sel_coef_G2; coef(sel)];
                sel_bloc_G2 = [sel_bloc_G2; repmat({blocs{iB}}, sum(sel), 1)];
            end

            for j=1:length(vnames)
                if sel(j)
                    fprintf('    [%s] %-45s coef=%.4f\n', blocs{iB}, vnames{j}, coef(j));
                end
            end
        catch ME
            fprintf('    ERREUR: %s\n', ME.message);
        end
    end
end

fprintf('\n=== Recapitulatif LASSO (avant top 8) ===\n');
fprintf('G1 : %d variables selectionnees\n', length(sel_vars_G1));
fprintf('G2 : %d variables selectionnees\n', length(sel_vars_G2));

% Limiter au top 8 variables par bloc par groupe (|coef| le plus eleve)
N_top = 8;
for iGrp = 1:2
    if iGrp==1, sv=sel_vars_G1; sc=sel_coef_G1; sb=sel_bloc_G1;
    else,        sv=sel_vars_G2; sc=sel_coef_G2; sb=sel_bloc_G2; end

    sv_new = {}; sc_new = []; sb_new = {};
    for iB = 1:4
        mask = strcmp(sb, blocs{iB});
        sv_b = sv(mask); sc_b = sc(mask);
        if length(sv_b) > N_top
            [~, idx_sort] = sort(abs(sc_b), 'descend');
            sv_b = sv_b(idx_sort(1:N_top));
            sc_b = sc_b(idx_sort(1:N_top));
        end
        sv_new = [sv_new; sv_b(:)];
        sc_new = [sc_new; sc_b(:)];
        sb_new = [sb_new; repmat({blocs{iB}}, length(sv_b), 1)];
    end
    if iGrp==1, sel_vars_G1=sv_new; sel_coef_G1=sc_new; sel_bloc_G1=sb_new;
    else,        sel_vars_G2=sv_new; sel_coef_G2=sc_new; sel_bloc_G2=sb_new; end
end
fprintf('Apres top %d par bloc : G1=%d vars | G2=%d vars\n', ...
    N_top, length(sel_vars_G1), length(sel_vars_G2));

%% -----------------------------------------------------------------------
%  SECTION 3 : LMM PAR BLOC sur variables selectionnees
%  -----------------------------------------------------------------------
fprintf('\n=== LMM par bloc ===\n');

lmm_G1 = struct(); lmm_G2 = struct();

for iB = 1:4
    bn = blocs{iB};

    for iGrp = 1:2
        if iGrp==1
            X = bloc_data{iB,1}; Y = bloc_data{iB,2};
            S = bloc_data{iB,3}; vnames = bloc_data{iB,4};
            gname = 'G1';
            sel_v = sel_vars_G1(strcmp(sel_bloc_G1, bn));
        else
            X = bloc_data_G2{iB,1}; Y = bloc_data_G2{iB,2};
            S = bloc_data_G2{iB,3}; vnames = bloc_data_G2{iB,4};
            gname = 'G2';
            sel_v = sel_vars_G2(strcmp(sel_bloc_G2, bn));
        end

        if isempty(sel_v)
            fprintf('  LMM Bloc %s %s : aucune variable selectionnee\n', bn, gname);
            continue;
        end

        % Indices des variables selectionnees
        sel_idx = find(ismember(vnames, sel_v));
        Xsel = X(:, sel_idx);

        ok = ~any(isnan(Xsel)|isinf(Xsel),2) & ~isnan(Y) & ~isinf(Y);
        Xsel = Xsel(ok,:); Yok = Y(ok); Sok = S(ok);
        ok_idx = find(ok);  % indices originaux des obs valides

        % Standardiser
        mu_s = mean(Xsel); sd_s = std(Xsel); sd_s(sd_s==0)=1;
        Xstd = (Xsel - mu_s) ./ sd_s;

        % Noms propres pour table
        col_names = cellfun(@(x) regexprep(x,'[^a-zA-Z0-9]','_'), vnames(sel_idx), 'UniformOutput', false);

        % Ajouter RPE_lag1 comme covariable pour controler l'autocorrelation AR(1)
        N_obs_ok = sum(ok);
        RPE_lag  = NaN(N_obs_ok, 1);
        for iSubj = unique(Sok)'
            idx_s = find(Sok == iSubj);
            if length(idx_s) > 1
                RPE_lag(idx_s(2:end)) = Yok(idx_s(1:end-1));
            end
        end
        RPE_lag(isnan(RPE_lag)) = nanmean(Yok);

        T = array2table([Yok Xstd RPE_lag double(Sok)], ...
            'VariableNames', ['RPE' col_names {'RPE_lag1'} {'Participant'}]);
        T.Participant = categorical(T.Participant);

        pred_str = strjoin([col_names {'RPE_lag1'}], ' + ');
        formula  = ['RPE ~ ' pred_str ' + (1|Participant)'];

        fprintf('  LMM Bloc %s %s (%d obs, %d vars)...\n', bn, gname, sum(ok), length(sel_idx));
        try
            lme = fitlme(T, formula);
            [~,~,FEStats] = fixedEffects(lme, 'DFMethod', 'satterthwaite');

            res.coef  = FEStats.Estimate(2:end);
            res.ci_lo = FEStats.Lower(2:end);
            res.ci_hi = FEStats.Upper(2:end);
            res.pval  = FEStats.pValue(2:end);
            res.names = col_names;
            res.R2    = lme.Rsquared.Adjusted;
            res.bloc  = bn;
            res.lme   = lme;
            % Autocorrelation des residus (lag 1)
            resid = residuals(lme);
            r_lag1 = corr(resid(1:end-1), resid(2:end));
            res.r_lag1 = r_lag1;
            fprintf('    R2 ajuste = %.3f  |  r_lag1 = %.3f\n', res.R2, r_lag1);

            % R2 marginal par variable : R2_full - R2_sans_variable
            % Mesure la contribution unique de chaque variable
            R2_marg = NaN(length(col_names), 1);
            for iVar = 1:length(col_names)
                % Modele sans cette variable
                col_reduced = col_names([1:iVar-1, iVar+1:end]);
                col_reduced_lag = [col_reduced {'RPE_lag1'}];
                pred_red = strjoin(col_reduced_lag, ' + ');
                formula_red = ['RPE ~ ' pred_red ' + (1|Participant)'];
                try
                    lme_red = fitlme(T, formula_red);
                    R2_marg(iVar) = res.R2 - lme_red.Rsquared.Adjusted;
                catch
                    R2_marg(iVar) = NaN;
                end
            end
            res.R2_marg = R2_marg;
            fprintf('    R2 marginal par variable :\n');
            for iVar = 1:length(col_names)
                if ~isnan(R2_marg(iVar))
                    fprintf('      %-40s : dR2 = %.4f\n', col_names{iVar}, R2_marg(iVar));
                end
            end

            if iGrp==1, lmm_G1.(bn) = res;
            else,        lmm_G2.(bn) = res; end
        catch ME
            fprintf('    ERREUR LMM: %s\n', ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
%  SECTION 4 : FIGURE - Coefficients beta par bloc
%  -----------------------------------------------------------------------
fprintf('\n=== Figures coefficients LMM ===\n');

col_G1 = [0.85 0.15 0.15];
col_G2 = [0.15 0.45 0.85];

for iB = 1:4
    bn = blocs{iB};
    has_G1 = isfield(lmm_G1, bn);
    has_G2 = isfield(lmm_G2, bn);
    if ~has_G1 && ~has_G2, continue; end

    % Collecter tous les noms de variables du bloc
    all_names = {};
    if has_G1, all_names = [all_names; lmm_G1.(bn).names']; end
    if has_G2, all_names = [all_names; lmm_G2.(bn).names']; end
    all_names = unique(all_names, 'stable');
    N_v = length(all_names);

    fig_name = sprintf('LMM_Bloc%s_G1vsG2', bn);
    figure('Name', fig_name, 'NumberTitle', 'off', ...
        'Position', [50 50 max(800, N_v*80+200) 550], 'Color', 'white');
    hold on;

    offset = 0.2;
    for i = 1:N_v
        vn = all_names{i};

        % G1
        if has_G1
            idx1 = find(strcmp(lmm_G1.(bn).names, vn));
            if ~isempty(idx1)
                b1=lmm_G1.(bn).coef(idx1); lo1=lmm_G1.(bn).ci_lo(idx1);
                hi1=lmm_G1.(bn).ci_hi(idx1); p1=lmm_G1.(bn).pval(idx1);
                bar(i-offset, b1, 0.35, 'FaceColor',col_G1,'EdgeColor','none','FaceAlpha',0.8);
                errorbar(i-offset, b1, b1-lo1, hi1-b1, 'k','LineWidth',1.2,'CapSize',5);
                if p1<0.001, sig='***'; elseif p1<0.01, sig='**'; elseif p1<0.05, sig='*'; else, sig=''; end
                if ~isempty(sig)
                    text(i-offset, hi1+0.05, sig,'HorizontalAlignment','center',...
                        'FontSize',9,'Color',col_G1,'FontWeight','bold');
                end
            end
        end

        % G2
        if has_G2
            idx2 = find(strcmp(lmm_G2.(bn).names, vn));
            if ~isempty(idx2)
                b2=lmm_G2.(bn).coef(idx2); lo2=lmm_G2.(bn).ci_lo(idx2);
                hi2=lmm_G2.(bn).ci_hi(idx2); p2=lmm_G2.(bn).pval(idx2);
                bar(i+offset, b2, 0.35, 'FaceColor',col_G2,'EdgeColor','none','FaceAlpha',0.8);
                errorbar(i+offset, b2, b2-lo2, hi2-b2, 'k','LineWidth',1.2,'CapSize',5);
                if p2<0.001, sig='***'; elseif p2<0.01, sig='**'; elseif p2<0.05, sig='*'; else, sig=''; end
                if ~isempty(sig)
                    text(i+offset, hi2+0.05, sig,'HorizontalAlignment','center',...
                        'FontSize',9,'Color',col_G2,'FontWeight','bold');
                end
            end
        end
    end

    yline(0,'--k','LineWidth',1,'Alpha',0.5);
    xticks(1:N_v);
    xticklabels(cellfun(@normLabel, all_names, 'UniformOutput', false));
    xtickangle(35);
    ylabel('Coefficient β standardise (± IC 95%)','FontSize',10);

    r2_str = '';
    if has_G1, r2_str = [r2_str sprintf('G1 R²=%.2f  ', lmm_G1.(bn).R2)]; end
    if has_G2, r2_str = [r2_str sprintf('G2 R²=%.2f', lmm_G2.(bn).R2)]; end
    title(sprintf('LMM Bloc %s (lag1) - %s', bn, r2_str), 'FontSize',11,'FontWeight','bold');

    % Legende avec couleurs forcees
    h_leg1 = patch(NaN, NaN, col_G1, 'FaceAlpha', 0.8, 'EdgeColor','none');
    h_leg2 = patch(NaN, NaN, col_G2, 'FaceAlpha', 0.8, 'EdgeColor','none');
    legend([h_leg1 h_leg2], {'G1 Short Duration','G2 Long Duration'}, ...
        'Location','best','FontSize',9);
    grid on; box on; hold off;
    xlim([0.5 N_v+0.5]);

    saveas(gcf, fullfile(PathSave, [fig_name '.png']));
    fprintf('  Figure : %s.png\n', fig_name);
end


%% -----------------------------------------------------------------------
%  TABLEAU RECAPITULATIF : beta + R2 marginal par variable
%  Classe par R2 marginal decroissant (max des deux groupes)
%  G1 et G2 sur la meme ligne
%  -----------------------------------------------------------------------
fprintf('\n=== Tableau recapitulatif biomarqueurs ===\n');

blocs_labels = {'A - BySegMod','B - BySegAxis','C - Goubault','D - EMG'};

for iB = 1:4
    bn = blocs{iB};
    has_G1 = isfield(lmm_G1, bn) && isfield(lmm_G1.(bn), 'R2_marg');
    has_G2 = isfield(lmm_G2, bn) && isfield(lmm_G2.(bn), 'R2_marg');
    if ~has_G1 && ~has_G2, continue; end

    % Collecter toutes les variables du bloc
    all_vars = {};
    if has_G1, all_vars = [all_vars; lmm_G1.(bn).names']; end
    if has_G2, all_vars = [all_vars; lmm_G2.(bn).names']; end
    all_vars = unique(all_vars, 'stable');
    N_v = length(all_vars);

    % Construire matrice de donnees : [beta_G1, R2m_G1, beta_G2, R2m_G2]
    tab_data = NaN(N_v, 4);
    tab_pval = NaN(N_v, 2);
    for iv = 1:N_v
        vn = all_vars{iv};
        if has_G1
            idx = find(strcmp(lmm_G1.(bn).names, vn));
            if ~isempty(idx)
                tab_data(iv,1) = lmm_G1.(bn).coef(idx);
                tab_data(iv,2) = lmm_G1.(bn).R2_marg(idx);
                tab_pval(iv,1) = lmm_G1.(bn).pval(idx);
            end
        end
        if has_G2
            idx = find(strcmp(lmm_G2.(bn).names, vn));
            if ~isempty(idx)
                tab_data(iv,3) = lmm_G2.(bn).coef(idx);
                tab_data(iv,4) = lmm_G2.(bn).R2_marg(idx);
                tab_pval(iv,2) = lmm_G2.(bn).pval(idx);
            end
        end
    end

    % Trier par max(R2m_G1, R2m_G2) decroissant
    max_r2 = max(abs(tab_data(:,2)), abs(tab_data(:,4)));
    max_r2(isnan(max_r2)) = 0;
    [~, sort_idx] = sort(max_r2, 'descend');
    all_vars  = all_vars(sort_idx);
    tab_data  = tab_data(sort_idx,:);
    tab_pval  = tab_pval(sort_idx,:);
    N_v = length(all_vars);

    % --- Figure tableau ---
    col_headers = {'Variable', 'beta G1', 'R2m G1', 'sig G1', 'beta G2', 'R2m G2', 'sig G2'};
    col_x = [0.02, 0.40, 0.52, 0.61, 0.67, 0.79, 0.88];
    col_G1_bg = [0.95 0.88 0.88];
    col_G2_bg = [0.88 0.91 0.97];

    row_h  = min(0.055, 0.82/max(N_v,1));
    fig_h  = max(300, 80 + N_v * 28);
    fig_w  = 1000;

    fig = figure('Visible','off','Position',[50 50 fig_w fig_h],'Color','white');
    ax  = axes('Position',[0 0 1 1],'Visible','off');
    ax.XLim = [0 1]; ax.YLim = [0 1];
    hold on;

    % Titre
    rectangle('Position',[0 0.95 1 0.05],'FaceColor',[0.15 0.25 0.45],'EdgeColor','none');
    text(0.5, 0.975, sprintf('LMM - Bloc %s : %s', bn, blocs_labels{iB}), ...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',11,'FontWeight','bold','Color','white');

    % Sous-titres R2 global
    r2_str = '';
    if has_G1, r2_str=[r2_str sprintf('G1 R²=%.2f  ', lmm_G1.(bn).R2)]; end
    if has_G2, r2_str=[r2_str sprintf('G2 R²=%.2f', lmm_G2.(bn).R2)]; end
    rectangle('Position',[0 0.90 1 0.05],'FaceColor',[0.85 0.88 0.93],'EdgeColor','none');
    text(0.5, 0.925, r2_str, 'HorizontalAlignment','center',...
        'FontSize',9,'Color',[0.2 0.2 0.4]);

    % En-tetes colonnes
    head_y = 0.87;
    rectangle('Position',[0 head_y-0.005 1 0.035],'FaceColor',[0.75 0.78 0.85],'EdgeColor','none');
    for iC = 1:length(col_headers)
        text(col_x(iC), head_y+0.012, col_headers{iC},...
            'FontSize',8,'FontWeight','bold','Color',[0.1 0.1 0.3],...
            'VerticalAlignment','middle');
    end
    line([0 1],[head_y-0.005 head_y-0.005],'Color',[0.5 0.5 0.6],'LineWidth',1);

    % Lignes de donnees
    for iv = 1:N_v
        y_row = head_y - 0.005 - iv*row_h;
        y_txt = y_row + row_h*0.5;

        % Fond alterne
        if mod(iv,2)==0
            rectangle('Position',[0 y_row 1 row_h],'FaceColor',[0.97 0.97 0.99],'EdgeColor','none');
        end

        % Fond couleur si R2m > 0.03
        if ~isnan(tab_data(iv,2)) && tab_data(iv,2) > 0.03
            rectangle('Position',[col_x(2) y_row col_x(5)-col_x(2) row_h],...
                'FaceColor',[1.0 0.92 0.88],'EdgeColor','none');
        end
        if ~isnan(tab_data(iv,4)) && tab_data(iv,4) > 0.03
            rectangle('Position',[col_x(5) y_row 1-col_x(5) row_h],...
                'FaceColor',[0.88 0.93 1.0],'EdgeColor','none');
        end

        % Nom variable
        vn_display = normLabel(all_vars{iv});
        text(col_x(1)+0.005, y_txt, vn_display,...
            'FontSize',7,'VerticalAlignment','middle','Interpreter','none',...
            'Color',[0.1 0.1 0.1]);

        % G1 beta
        if ~isnan(tab_data(iv,1))
            fw = 'normal';
            if ~isnan(tab_pval(iv,1)) && tab_pval(iv,1)<0.05, fw='bold'; end
            text(col_x(2)+0.02, y_txt, sprintf('%.3f', tab_data(iv,1)),...
                'FontSize',7.5,'VerticalAlignment','middle',...
                'FontWeight',fw,'HorizontalAlignment','center',...
                'Color',[0.7 0.1 0.1]);
        end
        % G1 R2m
        if ~isnan(tab_data(iv,2))
            col_r2 = [0.4 0.4 0.4];
            if tab_data(iv,2) > 0.05, col_r2=[0.1 0.5 0.1]; end
            text(col_x(3)+0.02, y_txt, sprintf('%.4f', tab_data(iv,2)),...
                'FontSize',7.5,'VerticalAlignment','middle',...
                'HorizontalAlignment','center','Color',col_r2,'FontWeight','bold');
        end
        % G1 sig
        if ~isnan(tab_pval(iv,1))
            p=tab_pval(iv,1);
            if p<0.001, sig='***'; elseif p<0.01, sig='**'; elseif p<0.05, sig='*'; else, sig=''; end
            text(col_x(4)+0.015, y_txt, sig,'FontSize',8,'FontWeight','bold',...
                'VerticalAlignment','middle','HorizontalAlignment','center','Color',[0.7 0.1 0.1]);
        end
        % G2 beta
        if ~isnan(tab_data(iv,3))
            fw = 'normal';
            if ~isnan(tab_pval(iv,2)) && tab_pval(iv,2)<0.05, fw='bold'; end
            text(col_x(5)+0.02, y_txt, sprintf('%.3f', tab_data(iv,3)),...
                'FontSize',7.5,'VerticalAlignment','middle',...
                'FontWeight',fw,'HorizontalAlignment','center',...
                'Color',[0.1 0.2 0.7]);
        end
        % G2 R2m
        if ~isnan(tab_data(iv,4))
            col_r2 = [0.4 0.4 0.4];
            if tab_data(iv,4) > 0.05, col_r2=[0.1 0.5 0.1]; end
            text(col_x(6)+0.02, y_txt, sprintf('%.4f', tab_data(iv,4)),...
                'FontSize',7.5,'VerticalAlignment','middle',...
                'HorizontalAlignment','center','Color',col_r2,'FontWeight','bold');
        end
        % G2 sig
        if ~isnan(tab_pval(iv,2))
            p=tab_pval(iv,2);
            if p<0.001, sig='***'; elseif p<0.01, sig='**'; elseif p<0.05, sig='*'; else, sig=''; end
            text(col_x(7)+0.015, y_txt, sig,'FontSize',8,'FontWeight','bold',...
                'VerticalAlignment','middle','HorizontalAlignment','center','Color',[0.1 0.2 0.7]);
        end

        line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
    end

    % Legende
    y_leg = head_y - 0.005 - (N_v+0.8)*row_h;
    text(0.02, y_leg, '* p<0.05  ** p<0.01  *** p<0.001  |  R2m vert >= 0.05  |  Fond colore = R2m > 0.03  |  G1=rouge  G2=bleu',...
        'FontSize',6.5,'Color',[0.3 0.3 0.3],'FontAngle','italic');
    line([0 1 1 0 0],[0 0 1 1 0],'Color',[0.4 0.4 0.5],'LineWidth',1);
    hold off;

    fig_name = sprintf('LMM_Tab_Bloc%s', bn);
    fpath = fullfile(PathSave, [fig_name '.png']);
    exportgraphics(fig, fpath, 'Resolution', 150);
    close(fig);
    fprintf('  Tableau : %s.png\n', fig_name);
end

%% -----------------------------------------------------------------------
%  FIGURE COMPARATIVE INTER-BLOCS
%  1. Bar chart R2 global par bloc (G1 et G2)
%  2. Bar chart R2 marginal toutes variables triees
%  -----------------------------------------------------------------------

blocs_labels_short = {'BySegMod','BySegAxis','Goubault','EMG'};
bloc_colors = {[0.20 0.65 0.30], [0.90 0.55 0.10], [0.75 0.10 0.10], [0.15 0.35 0.75]};
col_G1 = [0.85 0.15 0.15];
col_G2 = [0.15 0.45 0.85];

%% --- Figure 1 : R2 global par bloc ---
figure('Name','R2_Global_ParBloc','NumberTitle','off',...
    'Position',[50 50 700 450],'Color','white');
hold on;

R2_G1_vec = NaN(4,1);
R2_G2_vec = NaN(4,1);
for iB = 1:4
    bn = blocs{iB};
    if isfield(lmm_G1,bn), R2_G1_vec(iB) = lmm_G1.(bn).R2; end
    if isfield(lmm_G2,bn), R2_G2_vec(iB) = lmm_G2.(bn).R2; end
end

x = 1:4;
b1 = bar(x-0.2, R2_G1_vec, 0.35, 'FaceColor','flat', 'EdgeColor','none');
b2 = bar(x+0.2, R2_G2_vec, 0.35, 'FaceColor','flat', 'EdgeColor','none');

for iB = 1:4
    b1.CData(iB,:) = col_G1;
    b2.CData(iB,:) = col_G2;
    % Valeurs au dessus des barres
    if ~isnan(R2_G1_vec(iB))
        text(iB-0.2, R2_G1_vec(iB)+0.01, sprintf('%.2f', R2_G1_vec(iB)),...
            'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',col_G1);
    end
    if ~isnan(R2_G2_vec(iB))
        text(iB+0.2, R2_G2_vec(iB)+0.01, sprintf('%.2f', R2_G2_vec(iB)),...
            'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',col_G2);
    end
end

% Ligne reference Goubault 2023
yline(0.65, '--k', 'Goubault 2023 max', 'LineWidth', 1.5,...
    'LabelHorizontalAlignment','left','FontSize',8);
yline(0.45, ':k', 'Goubault 2023 min', 'LineWidth', 1,...
    'LabelHorizontalAlignment','left','FontSize',8);

xticks(1:4);
xticklabels(blocs_labels_short);
ylabel('R² ajuste (LMM + lag1)', 'FontSize', 11);
xlabel('Bloc de variables', 'FontSize', 11);
title('R² global par bloc - G1 (rouge) vs G2 (bleu)', 'FontSize', 12, 'FontWeight', 'bold');
legend([b1 b2], {'G1 Short Duration', 'G2 Long Duration'},...
    'Location','southeast','FontSize',10);
ylim([0 1]); grid on; box on; hold off;

saveas(gcf, fullfile(PathSave, 'LMM_R2_ParBloc.png'));
fprintf('Figure : LMM_R2_ParBloc.png\n');

%% --- Figure 2 : R2 marginal toutes variables triees ---
% Collecter toutes les variables de tous les blocs
all_vnames = {};
all_r2m_G1 = [];
all_r2m_G2 = [];
all_beta_G1 = [];
all_beta_G2 = [];
all_pval_G1 = [];
all_pval_G2 = [];
all_bloc_idx = [];

for iB = 1:4
    bn = blocs{iB};
    % Variables G1
    if isfield(lmm_G1,bn) && isfield(lmm_G1.(bn),'R2_marg')
        for iv = 1:length(lmm_G1.(bn).names)
            vn = lmm_G1.(bn).names{iv};
            % Chercher si deja dans la liste
            idx_exist = find(strcmp(all_vnames, strrep(vn,'_',' ')));
            if isempty(idx_exist)
                all_vnames{end+1} = strrep(vn,'_',' ');
                all_r2m_G1(end+1)  = lmm_G1.(bn).R2_marg(iv);
                all_r2m_G2(end+1)  = NaN;
                all_beta_G1(end+1) = lmm_G1.(bn).coef(iv);
                all_beta_G2(end+1) = NaN;
                all_pval_G1(end+1) = lmm_G1.(bn).pval(iv);
                all_pval_G2(end+1) = NaN;
                all_bloc_idx(end+1) = iB;
            end
        end
    end
    % Variables G2
    if isfield(lmm_G2,bn) && isfield(lmm_G2.(bn),'R2_marg')
        for iv = 1:length(lmm_G2.(bn).names)
            vn = lmm_G2.(bn).names{iv};
            label = strrep(vn,'_',' ');
            idx_exist = find(strcmp(all_vnames, label));
            if isempty(idx_exist)
                all_vnames{end+1} = label;
                all_r2m_G1(end+1)  = NaN;
                all_r2m_G2(end+1)  = lmm_G2.(bn).R2_marg(iv);
                all_beta_G1(end+1) = NaN;
                all_beta_G2(end+1) = lmm_G2.(bn).coef(iv);
                all_pval_G1(end+1) = NaN;
                all_pval_G2(end+1) = lmm_G2.(bn).pval(iv);
                all_bloc_idx(end+1) = iB;
            else
                all_r2m_G2(idx_exist)  = lmm_G2.(bn).R2_marg(iv);
                all_beta_G2(idx_exist) = lmm_G2.(bn).coef(iv);
                all_pval_G2(idx_exist) = lmm_G2.(bn).pval(iv);
            end
        end
    end
end

% Trier par max(R2m_G1, R2m_G2) decroissant
max_r2m = max(abs(all_r2m_G1), abs(all_r2m_G2));
max_r2m(isnan(max_r2m)) = 0;
[~, sort_idx] = sort(max_r2m, 'descend');
all_vnames   = all_vnames(sort_idx);
all_r2m_G1   = all_r2m_G1(sort_idx);
all_r2m_G2   = all_r2m_G2(sort_idx);
all_pval_G1  = all_pval_G1(sort_idx);
all_pval_G2  = all_pval_G2(sort_idx);
all_bloc_idx = all_bloc_idx(sort_idx);

N_all = length(all_vnames);
fig_w2 = 1200;
fig_h2 = max(500, 100 + N_all * 30);

fig2 = figure('Name','R2_Marginal_AllVars','NumberTitle','off',...
    'Position',[50 50 fig_w2 fig_h2],'Color','white');
hold on;

x_pos = 1:N_all;
offset = 0.2;

for iv = 1:N_all
    bc = bloc_colors{all_bloc_idx(iv)};

    % G1 barre verticale
    if ~isnan(all_r2m_G1(iv)) && all_r2m_G1(iv) > 0
        bar(x_pos(iv)-offset, all_r2m_G1(iv), 0.35,...
            'FaceColor', col_G1, 'EdgeColor','none', 'FaceAlpha', 0.85);
        if ~isnan(all_pval_G1(iv)) && all_pval_G1(iv) < 0.05
            if all_pval_G1(iv)<0.001, sig='***'; elseif all_pval_G1(iv)<0.01, sig='**'; else, sig='*'; end
            text(x_pos(iv)-offset, all_r2m_G1(iv)+0.001, sig,...
                'FontSize',7,'Color',col_G1,'FontWeight','bold',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end

    % G2 barre verticale
    if ~isnan(all_r2m_G2(iv)) && all_r2m_G2(iv) > 0
        bar(x_pos(iv)+offset, all_r2m_G2(iv), 0.35,...
            'FaceColor', col_G2, 'EdgeColor','none', 'FaceAlpha', 0.85);
        if ~isnan(all_pval_G2(iv)) && all_pval_G2(iv) < 0.05
            if all_pval_G2(iv)<0.001, sig='***'; elseif all_pval_G2(iv)<0.01, sig='**'; else, sig='*'; end
            text(x_pos(iv)+offset, all_r2m_G2(iv)+0.001, sig,...
                'FontSize',7,'Color',col_G2,'FontWeight','bold',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
        end
    end

    % Marqueur bloc en bas (rectangle colore sous l'axe X)
    patch([x_pos(iv)-0.45 x_pos(iv)+0.45 x_pos(iv)+0.45 x_pos(iv)-0.45],...
        [-0.003 -0.003 -0.001 -0.001], bc, 'EdgeColor','none');
end

% Lignes seuil
yline(0.03, '--k', 'LineWidth', 1.2);
yline(0.05, ':k',  'LineWidth', 1.2);
text(N_all+0.5, 0.030, 'R²m=0.03','FontSize',8,'VerticalAlignment','bottom');
text(N_all+0.5, 0.050, 'R²m=0.05','FontSize',8,'VerticalAlignment','bottom');

% Axes et labels
xticks(x_pos);
xticklabels(cellfun(@normLabel, all_vnames, 'UniformOutput', false));
xtickangle(40);
ax2 = gca; ax2.XAxis.FontSize = 7.5;
ylabel('R² marginal (contribution unique au RPE)', 'FontSize', 11);
title('R² marginal par variable - G1 (rouge) vs G2 (bleu)', ...
    'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5 N_all+0.5]);
ylim([-0.005 max(max(all_r2m_G1, [], 'omitnan'), max(all_r2m_G2, [], 'omitnan'))*1.25]);

% Legende
h_g1 = patch(NaN, NaN, col_G1, 'FaceAlpha',0.85);
h_g2 = patch(NaN, NaN, col_G2, 'FaceAlpha',0.85);
h_b  = zeros(4,1);
for iB = 1:4
    h_b(iB) = patch(NaN, NaN, bloc_colors{iB});
end
legend([h_g1 h_g2 h_b'], {'G1','G2', blocs_labels_short{:}},...
    'Location','northeast','FontSize',8);

grid on; box on; hold off;

fpath2 = fullfile(PathSave, 'LMM_R2m_AllVars.png');
exportgraphics(fig2, fpath2, 'Resolution', 150);
close(fig2);
fprintf('Figure : LMM_R2m_AllVars.png\n');


%% -----------------------------------------------------------------------
%  BOUCLE DE STABILITE LASSO
%  -----------------------------------------------------------------------
fprintf('\n=== Boucle de stabilite LASSO (100 graines) ===\n');

seeds = 1:100;
N_seeds = length(seeds);

% Compteurs de selection par variable par groupe
count_G1 = struct('A',{{}},'B',{{}},'C',{{}},'D',{{}});
count_G2 = struct('A',{{}},'B',{{}},'C',{{}},'D',{{}});
for iB = 1:4
    bn = blocs{iB};
    count_G1.(bn) = containers.Map('KeyType','char','ValueType','double');
    count_G2.(bn) = containers.Map('KeyType','char','ValueType','double');
end

for iSeed = 1:N_seeds
    rng(seeds(iSeed));
    if mod(iSeed,10)==0
        fprintf('  Graine %d/%d...\n', iSeed, N_seeds);
    end

    for iB = 1:4
        bn = blocs{iB};
        for iGrp = 1:2
            if iGrp==1
                X = bloc_data{iB,1}; Y = bloc_data{iB,2};
                vnames = bloc_data{iB,4}; cnt = count_G1.(bn);
            else
                X = bloc_data_G2{iB,1}; Y = bloc_data_G2{iB,2};
                vnames = bloc_data_G2{iB,4}; cnt = count_G2.(bn);
            end

            ok = ~any(isnan(X)|isinf(X),2) & ~isnan(Y) & ~isinf(Y);
            if sum(ok) < 10, continue; end
            Xok=X(ok,:); Yok=Y(ok);
            mu_X=mean(Xok); sd_X=std(Xok); sd_X(sd_X==0)=1;
            Xs=(Xok-mu_X)./sd_X;

            try
                [B,FitInfo] = lasso(Xs, Yok, 'CV',5, 'Alpha',1);
                coef = B(:, FitInfo.Index1SE);
                sel  = coef ~= 0;
                for j = find(sel)'
                    vn = vnames{j};
                    if isKey(cnt, vn), cnt(vn) = cnt(vn)+1;
                    else,              cnt(vn) = 1; end
                end
            catch; end

            if iGrp==1, count_G1.(bn) = cnt;
            else,        count_G2.(bn) = cnt; end
        end
    end
end

% Reinitialiser la graine principale
rng(42);

fprintf('\n=== Tableau recapitulatif stabilite LASSO (%d graines) ===\n', N_seeds);

%% Generer un tableau PNG de stabilite par bloc
blocs_labels_stab = {'A - BySegMod','B - BySegAxis','C - Goubault','D - EMG'};

for iB = 1:4
    bn = blocs{iB};

    % Collecter toutes les variables des deux groupes
    all_vn = {};
    cnt1 = count_G1.(bn); cnt2 = count_G2.(bn);
    all_vn = [all_vn keys(cnt1) keys(cnt2)];
    all_vn = unique(all_vn);
    if isempty(all_vn), continue; end

    % Construire vecteurs de frequence
    freq1 = zeros(length(all_vn),1);
    freq2 = zeros(length(all_vn),1);
    for iv = 1:length(all_vn)
        vn = all_vn{iv};
        if isKey(cnt1,vn), freq1(iv) = cnt1(vn)/N_seeds*100; end
        if isKey(cnt2,vn), freq2(iv) = cnt2(vn)/N_seeds*100; end
    end

    % Trier par max(freq1, freq2) decroissant
    [~,si] = sort(max(freq1,freq2),'descend');
    all_vn = all_vn(si); freq1=freq1(si); freq2=freq2(si);
    N_v = length(all_vn);

    % Figure tableau
    row_h = max(0.012, min(0.030, 0.82/max(N_v,1)));
    fig_h = max(350, 80 + N_v*24);
    fig_w = 950;

    fig = figure('Visible','off','Position',[50 50 fig_w fig_h],'Color','white');
    ax  = axes('Position',[0 0 1 1],'Visible','off');
    ax.XLim=[0 1]; ax.YLim=[0 1]; hold on;

    % Titre
    rectangle('Position',[0 0.95 1 0.05],'FaceColor',[0.15 0.25 0.45],'EdgeColor','none');
    text(0.5,0.975,sprintf('Stabilite LASSO - Bloc %s : %s (%d graines)', bn, blocs_labels_stab{iB}, N_seeds),...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'FontSize',10,'FontWeight','bold','Color','white');

    % En-tetes
    head_y = 0.90;
    rectangle('Position',[0 head_y-0.005 1 0.030],'FaceColor',[0.75 0.78 0.85],'EdgeColor','none');
    col_x = [0.02 0.42 0.56 0.70 0.84];
    hdrs = {'Variable','Freq G1 (%)','Freq G2 (%)','Stable G1','Stable G2'};
    for ic=1:5
        text(col_x(ic), head_y+0.010, hdrs{ic},'FontSize',8,'FontWeight','bold',...
            'Color',[0.1 0.1 0.3],'VerticalAlignment','middle');
    end
    line([0 1],[head_y-0.005 head_y-0.005],'Color',[0.5 0.5 0.6],'LineWidth',1);

    thresh_pct = 75;  % stable si selectionne dans >= 50% des runs

    for iv = 1:N_v
        y_row = head_y - 0.005 - iv*row_h;
        y_txt = y_row + row_h*0.5;

        if mod(iv,2)==0
            rectangle('Position',[0 y_row 1 row_h],'FaceColor',[0.97 0.97 0.99],'EdgeColor','none');
        end

        % Fond colore si stable dans au moins un groupe
        is_stab1 = freq1(iv) >= thresh_pct;
        is_stab2 = freq2(iv) >= thresh_pct;
        if is_stab1
            rectangle('Position',[col_x(2) y_row col_x(4)-col_x(2) row_h],...
                'FaceColor',[1.0 0.92 0.88],'EdgeColor','none');
        end
        if is_stab2
            rectangle('Position',[col_x(3) y_row col_x(5)-col_x(3) row_h],...
                'FaceColor',[0.88 0.93 1.0],'EdgeColor','none');
        end

        % Nom variable
        text(col_x(1)+0.005, y_txt, normLabel(all_vn{iv}),...
            'FontSize',7,'VerticalAlignment','middle','Interpreter','none',...
            'Color',[0.1 0.1 0.1]);

        % Freq G1
        col_f1 = [0.4 0.4 0.4]; if freq1(iv)>=thresh_pct, col_f1=[0.7 0.1 0.1]; end
        if freq1(iv)>0
            text(col_x(2)+0.04, y_txt, sprintf('%.0f%%', freq1(iv)),...
                'FontSize',7.5,'VerticalAlignment','middle','HorizontalAlignment','center',...
                'FontWeight','bold','Color',col_f1);
        end

        % Freq G2
        col_f2 = [0.4 0.4 0.4]; if freq2(iv)>=thresh_pct, col_f2=[0.1 0.2 0.7]; end
        if freq2(iv)>0
            text(col_x(3)+0.04, y_txt, sprintf('%.0f%%', freq2(iv)),...
                'FontSize',7.5,'VerticalAlignment','middle','HorizontalAlignment','center',...
                'FontWeight','bold','Color',col_f2);
        end

        % Stable G1
        if is_stab1
            text(col_x(4)+0.04, y_txt, '✓','FontSize',10,'VerticalAlignment','middle',...
                'HorizontalAlignment','center','Color',[0.7 0.1 0.1],'FontWeight','bold');
        end

        % Stable G2
        if is_stab2
            text(col_x(5)+0.04, y_txt, '✓','FontSize',10,'VerticalAlignment','middle',...
                'HorizontalAlignment','center','Color',[0.1 0.2 0.7],'FontWeight','bold');
        end

        line([0 1],[y_row y_row],'Color',[0.88 0.88 0.90],'LineWidth',0.3);
    end

    y_leg = head_y - 0.005 - (N_v+0.8)*row_h;
    text(0.02, y_leg, sprintf('Stable = selectionne dans >= %d%% des runs  |  G1=rouge  G2=bleu', thresh_pct),...
        'FontSize',7,'Color',[0.3 0.3 0.3],'FontAngle','italic');
    line([0 1 1 0 0],[0 0 1 1 0],'Color',[0.4 0.4 0.5],'LineWidth',1);
    hold off;

    fig_name = sprintf('LASSO_Stabilite_Bloc%s.png', bn);
    fpath_stab = fullfile(PathSave, fig_name);
    exportgraphics(fig, fpath_stab, 'Resolution',150);
    close(fig);
    fprintf('  Tableau : %s\n', fig_name);
end

%% -----------------------------------------------------------------------
%  SAVE
%  -----------------------------------------------------------------------
save(fullfile(PathSave, 'LASSO_LMM_Results.mat'), ...
    'sel_vars_G1','sel_vars_G2','sel_coef_G1','sel_coef_G2',...
    'sel_bloc_G1','sel_bloc_G2','lmm_G1','lmm_G2');
fprintf('\nResultats sauvegardes : LASSO_LMM_Results.mat\n');
fprintf('Termine.\n');


% Fonction de normalisation des labels pour l'affichage

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------

function [X_long, RPE_long, subj_long] = buildBlock(Vi, Mat_data, Mat_RPE, N_bins)
% Mat_data : N_all x N_bins
% Retourne format long centre par participant : X(i,k) - mean_k(X(i,:))
% -> capture l'evolution intra-individuelle, coherent avec Spearman groupe
    N_subj = length(Vi);
    N_obs  = N_subj * N_bins;
    X_long    = NaN(N_obs, 1);
    RPE_long  = NaN(N_obs, 1);
    subj_long = NaN(N_obs, 1);
    row = 0;
    for iG = 1:N_subj
        iP = Vi(iG);
        x_subj   = Mat_data(iP, :);          % 1 x N_bins
        rpe_subj = Mat_RPE(iP, :);           % 1 x N_bins
        % Centrage individuel : soustraire la moyenne du participant
        x_cent   = x_subj - nanmean(x_subj);
        rpe_cent = rpe_subj - nanmean(rpe_subj);
        for k = 1:N_bins
            row = row + 1;
            X_long(row)    = x_cent(k);
            RPE_long(row)  = rpe_cent(k);
            subj_long(row) = iG;
        end
    end
end

function lbl = normLabel(lbl)
    lbl = strrep(lbl, '_', ' ');
    lbl = strrep(lbl, ' — ', ' - ');
    lbl = strrep(lbl, '—', ' - ');
    lbl = strrep(lbl, 'TFR ', '');
    lbl = strrep(lbl, 'Goub ', '');
    lbl = strrep(lbl, '  ', ' ');
    lbl = strtrim(lbl);
end