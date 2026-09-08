%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES FEATURES SPECTRALES — JEU DE DONNÉES ERGO (Robin2)    %%%%
%%%%  Pipeline Goubault et al. (2023) — version conforme               %%%%
%%%%                                                                   %%%%
%%%%  Corrections par rapport à la version précédente :                %%%%
%%%%   1. Filtre Butterworth ordre 2 zero-lag 0.5–29 Hz sur acc/gyro   %%%%
%%%%   2. CWT sur le SIGNAL CONTINU ENTIER (et non par fenêtre 1s) :   %%%%
%%%%      features extraites instant par instant du scalogramme puis   %%%%
%%%%      moyennées par tranche d'1 s                                  %%%%
%%%%   3. Morlet complexe, nombre d'onde 8, grille linéaire 0.05 Hz    %%%%
%%%%      (via cwtft — cwtfilterbank impose une grille logarithmique)  %%%%
%%%%                                                                   %%%%
%%%%  Features temporelles : inchangées (fenêtre 1s non chevauchante)  %%%%
%%%%  1260 variables : 7 seg × 3 sig × 4 comp × 15 features     
%%%%                                                                   %%%%
%%%%  SpectralFeatures_Ergo.mat (CWT par fenêtre d'1 s, sans filtre),  %%%%
%%%%  Ergo2 (CWT corrigée mais instants biaisés, RPE linéaire), Ergo3 
%%%% (tout corrigé)                                                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data  = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\';
path_xsens = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin\April28\';
path_save  = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

%% -----------------------------------------------------------------------
%  PARAMÈTRES
%  -----------------------------------------------------------------------
Fs            = 60;
dt            = 1/Fs;
N_bins        = 10;
N_frames_min  = 9000;
baseline_pct  = 0.10;
n_top_vals    = 10;   % facteur d'échelle = moyenne des N valeurs |.| les plus hautes
win_sec       = 1;
win_frames    = win_sec * Fs;      % 60 frames
half_bin_sec  = 2;                 % ±2s autour de chaque instant RPE
half_bin_fr   = half_bin_sec * Fs; % ±120 frames

% ── Niveau d'application de la normalisation Goubault ───────────────────
%  'signal'  : normalisation du signal brut avant calcul des features
%              (comportement des versions précédentes du pipeline)
%  'feature' : normalisation de chaque série temporelle de features
%              (lecture littérale de l'ordre du texte de Goubault 2023)
norm_level = 'feature';

% ── Filtre Butterworth (Goubault) ───────────────────────────────────────
filt_order = 2;
filt_band  = [0.5 29];             % Hz
[z_f, p_f, k_f] = butter(filt_order, filt_band/(Fs/2), 'bandpass');
SOS_filt = zp2sos(z_f, p_f, k_f);  % forme SOS : stable près de Nyquist

% Forme RPE normalisée G1+G2 (Goubault)
RPE_shape = [0, 0.2757, 0.4889, 0.6171, 0.7147, 0.7717, 0.8557, 0.9065, 0.9555, 1.0];

% Segments
Seg_names = {'L5','T8','Head','Shoulder','Arm','Forearm','Hand'};
Seg_idx   = [2, 5, 7, 8, 9, 10, 11];
N_seg     = numel(Seg_names);

% Signaux et composantes
sig_types  = {'Accel','AngVel','Jerk'};
comp_names = {'X','Y','Z','Mod'};
feat_temporal = {'Mean','Std','Median','Max','p10','p25','p75','p90'};
feat_cwt      = {'MedianFreq','SpectralEntropy','Power_below4Hz',...
                 'Power_above4Hz','TotalPower','PeakPower','PeakPowerFreq'};
feat_names    = [feat_temporal, feat_cwt];
N_feat        = numel(feat_names);

% ── CWT — Morlet complexe, nombre d'onde 8, pas linéaire 0.05 Hz ────────
omega0      = 8;                       % nombre d'onde (Goubault)
freq_grid   = (0.05:0.05:24.05)';      % 481 fréquences, pas linéaire exact
% Relation échelle↔fréquence — facteur de Fourier exact de la Morlet
% (convention Torrence & Compo, celle utilisée par cwtft) :
%   f = (omega0 + sqrt(2 + omega0^2)) / (4*pi*s)
% NB : la convention "fréquence de pic" f = omega0/(2*pi*s) introduit un
%      biais multiplicatif de +0.78% pour omega0 = 8.
fourier_factor = (omega0 + sqrt(2 + omega0^2)) / (4*pi);
cwt_scales  = fourier_factor ./ freq_grid;   % en secondes
cwt_thresh  = 4.0;                     % seuil bande basse/haute (Hz)

fprintf('=== Paramètres ===\n');
fprintf('  Fs=%d Hz | N_bins=%d\n', Fs, N_bins);
fprintf('  Filtre : Butterworth ordre %d zero-lag %.1f–%.0f Hz\n', ...
    filt_order, filt_band(1), filt_band(2));
fprintf('  Fenêtre features temporelles : %ds\n', win_sec);
fprintf('  CWT : signal continu | Morlet omega0=%d | %d fréquences (%.2f–%.2f Hz, pas %.2f)\n', ...
    omega0, numel(freq_grid), freq_grid(1), freq_grid(end), freq_grid(2)-freq_grid(1));
fprintf('  Fenêtre RPE : ±%ds\n', half_bin_sec);
fprintf('  Normalisation : niveau "%s" | baseline = %.0f%% initiaux | échelle = moyenne des %d plus hautes |.|\n', ...
    norm_level, baseline_pct*100, n_top_vals);
fprintf('  Segments : %s\n', strjoin(Seg_names,', '));
fprintf('  Total variables : %d seg × 3 sig × 4 comp × %d feat = %d\n\n',...
    N_seg, N_feat, N_seg*3*4*N_feat);

%% -----------------------------------------------------------------------
%  AUTO-TEST CWT — vérifie empiriquement la relation échelle↔fréquence
%  Un sinus pur à 5 Hz doit produire un pic de puissance à 5 Hz.
%  -----------------------------------------------------------------------
f_test = 5.0;
t_test = (0:1/Fs:20-1/Fs)';
x_test = sin(2*pi*f_test*t_test);
S_test = struct('val', x_test, 'period', dt);
try
    cw_test = cwtft(S_test, 'scales', cwt_scales, 'wavelet', {'morl', omega0});
catch ME
    error(['cwtft a échoué avec {''morl'',%d} : %s\n' ...
           'Vérifier la syntaxe de la Morlet paramétrée dans cette version.'], ...
           omega0, ME.message);
end
pow_test = mean(abs(cw_test.cfs).^2, 2);
[~, i_pk] = max(pow_test);
f_pk = freq_grid(i_pk);
err_rel = 100*abs(f_pk - f_test)/f_test;
fprintf('=== Auto-test CWT ===\n');
fprintf('  Sinus %.2f Hz → pic détecté à %.3f Hz (écart %.2f%%)\n', ...
    f_test, f_pk, err_rel);
if err_rel > 2
    warning(['Écart > 2%% entre fréquence imposée et fréquence détectée. ' ...
             'La relation échelle↔fréquence est probablement incorrecte.']);
else
    fprintf('  Relation échelle↔fréquence validée.\n');
end
% Comparaison avec les fréquences renvoyées par cwtft, si disponibles
if isfield(cw_test, 'frequencies') && ~isempty(cw_test.frequencies)
    d_max = max(abs(cw_test.frequencies(:) - freq_grid));
    fprintf('  Écart max grille imposée / grille cwtft : %.4f Hz\n', d_max);
end
fprintf('\n');
clear cw_test pow_test x_test t_test S_test;

%% -----------------------------------------------------------------------
%  FICHIERS XSENS
%  -----------------------------------------------------------------------
xsens_files = {
    fullfile(path_xsens, 'XSens_Fatigue_1.mat');
    fullfile(path_xsens, 'XSens_Fatigue_2.mat');
    fullfile(path_xsens, 'XSens_Fatigue_3.mat');
    fullfile(path_xsens, 'XSens_Fatigue_4.mat');
    fullfile(path_xsens, 'XSens_Fatigue_5.mat');
};

Results_Spectral = struct();

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------

    function sig_norm = normalize_goubault(sig, baseline_pct, n_top_vals)
        sig  = double(sig(:));
        N    = length(sig);
        n_bl = max(1, floor(N * baseline_pct));
        base = mean(sig(1:n_bl), 'omitnan');
        sc   = sig - base;
        % Facteur d'échelle = moyenne des n_top_vals valeurs absolues les plus hautes
        a    = abs(sc);
        a    = a(isfinite(a));
        k    = min(n_top_vals, numel(a));
        if k < 1
            sig_norm = sc;
            return;
        end
        a_sorted = sort(a, 'descend');
        nf   = mean(a_sorted(1:k));
        if isfinite(nf) && nf > 1e-10
            sig_norm = sc / nf;
        else
            sig_norm = sc;
        end
    end

    function feats = compute_temporal(x)
        x = x(isfinite(x));
        if isempty(x)
            feats = struct('Mean',NaN,'Std',NaN,'Median',NaN,'Max',NaN,...
                'p10',NaN,'p25',NaN,'p75',NaN,'p90',NaN);
            return;
        end
        feats.Mean   = mean(x);
        feats.Std    = std(x,0);
        feats.Median = median(x);
        feats.Max    = max(x);
        feats.p10    = prctile(x,10);
        feats.p25    = prctile(x,25);
        feats.p75    = prctile(x,75);
        feats.p90    = prctile(x,90);
    end

    % ── CWT sur signal continu → 7 features par instant → moyenne par 1s ──
    function feats = compute_cwt_timeseries(x, scales, omega0, dt, ...
                                            f_grid, thr_band, win_frames, N_wins)
        % Sortie : struct de 7 champs, chacun (N_wins × 1)
        blank = NaN(N_wins,1);
        feats = struct('MedianFreq',blank,'SpectralEntropy',blank,...
            'Power_below4Hz',blank,'Power_above4Hz',blank,...
            'TotalPower',blank,'PeakPower',blank,'PeakPowerFreq',blank);

        x = double(x(:));
        if ~any(isfinite(x)) || numel(x) < 2*numel(scales)/10
            return;
        end
        x(~isfinite(x)) = 0;

        try
            S  = struct('val', x, 'period', dt);
            cw = cwtft(S, 'scales', scales, 'wavelet', {'morl', omega0});
        catch
            return;
        end

        % Puissance instantanée : (N_freq × N_frames), en single pour la RAM
        P = single(abs(cw.cfs).^2);
        clear cw;

        tot = sum(P, 1);                       % 1 × N_frames
        ok  = tot > 0 & isfinite(tot);

        % --- Fréquence médiane (par instant) ---
        cp   = cumsum(P, 1) ./ tot;
        imed = sum(cp < 0.5, 1) + 1;
        imed = min(max(imed,1), numel(f_grid));
        v_med = f_grid(imed)';                 % 1 × N_frames
        clear cp;

        % --- Entropie spectrale (par instant) ---
        pn = P ./ tot;
        tt = pn .* log2(pn);
        tt(~isfinite(tt)) = 0;
        v_ent = -sum(tt, 1);
        clear pn tt;

        % --- Puissances par bande ---
        m_lo = f_grid <= thr_band;
        m_hi = f_grid >  thr_band;
        v_lo = sum(P(m_lo,:), 1);
        v_hi = sum(P(m_hi,:), 1);

        % --- Pic de puissance ---
        [v_pk, ipk] = max(P, [], 1);
        v_pkf = f_grid(ipk)';

        v_tot = tot;
        clear P;

        % Invalider les instants dégénérés
        v_med(~ok)=NaN; v_ent(~ok)=NaN; v_lo(~ok)=NaN;
        v_hi(~ok)=NaN;  v_pk(~ok)=NaN;  v_pkf(~ok)=NaN; v_tot(~ok)=NaN;

        % --- Moyenne par tranche d'1 s ---
        n_use = N_wins * win_frames;
        pack  = @(v) mean(reshape(double(v(1:n_use)), win_frames, N_wins), ...
                          1, 'omitnan')';

        feats.MedianFreq      = pack(v_med);
        feats.SpectralEntropy = pack(v_ent);
        feats.Power_below4Hz  = pack(v_lo);
        feats.Power_above4Hz  = pack(v_hi);
        feats.TotalPower      = pack(v_tot);
        feats.PeakPower       = pack(v_pk);
        feats.PeakPowerFreq   = pack(v_pkf);
    end

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
fprintf('=== Calcul des features spectrales (méthode Goubault) ===\n');

for iFile = 1:numel(xsens_files)
    fpath = xsens_files{iFile};
    if ~exist(fpath,'file')
        fprintf('  Introuvable : %s\n', fpath); continue;
    end
    fprintf('\nChargement : %s\n', fpath);
    tmp = load(fpath);
    XF  = tmp.XSens_Fatigue;
    subjects = fieldnames(XF);

    for iS = 1:numel(subjects)
        subj = subjects{iS};
        if strcmp(subj,'P07'), continue; end
        conditions = fieldnames(XF.(subj));
        fprintf('  %s — %s\n', subj, strjoin(conditions,' + '));

        for iK = 1:numel(conditions)
            kbd = conditions{iK};
            fprintf('    [%s] ', kbd);
            tic;

            %% --- Extraction signaux bruts ---
            try
                acc_raw = cell(1,N_seg);
                gyr_raw = cell(1,N_seg);
                for iSeg = 1:N_seg
                    acc_raw{iSeg} = XF.(subj).(kbd).segmentData(Seg_idx(iSeg)).acceleration;
                    gyr_raw{iSeg} = XF.(subj).(kbd).segmentData(Seg_idx(iSeg)).angularVelocity;
                end
            catch ME
                fprintf('Erreur : %s — skip\n', ME.message); continue;
            end

            N_frames = size(acc_raw{1},1);
            if N_frames < N_frames_min
                fprintf('Trop court (%d frames) — skip\n', N_frames); continue;
            end

            %% --- Filtrage Butterworth zero-lag (Goubault) ---
            for iSeg = 1:N_seg
                a = double(acc_raw{iSeg});
                g = double(gyr_raw{iSeg});
                a(~isfinite(a)) = 0;
                g(~isfinite(g)) = 0;
                acc_raw{iSeg} = filtfilt(SOS_filt, 1, a);
                gyr_raw{iSeg} = filtfilt(SOS_filt, 1, g);
            end

            %% --- Jerk = dérivée temporelle de l'accélération filtrée ---
            jerk_acc = cell(1,N_seg);
            for iSeg = 1:N_seg
                jerk_acc{iSeg} = [diff(acc_raw{iSeg},1,1)*Fs; zeros(1,3)];
            end

            %% --- Normalisation Goubault (si niveau 'signal') ---
            acc_n = cell(1,N_seg);
            gyr_n = cell(1,N_seg);
            jrk_n = cell(1,N_seg);
            for iSeg = 1:N_seg
                if strcmpi(norm_level,'signal')
                    acc_n{iSeg} = zeros(N_frames,3);
                    gyr_n{iSeg} = zeros(N_frames,3);
                    jrk_n{iSeg} = zeros(N_frames,3);
                    for iAx = 1:3
                        acc_n{iSeg}(:,iAx) = normalize_goubault(acc_raw{iSeg}(:,iAx), baseline_pct, n_top_vals);
                        gyr_n{iSeg}(:,iAx) = normalize_goubault(gyr_raw{iSeg}(:,iAx), baseline_pct, n_top_vals);
                        jrk_n{iSeg}(:,iAx) = normalize_goubault(jerk_acc{iSeg}(:,iAx), baseline_pct, n_top_vals);
                    end
                else
                    acc_n{iSeg} = acc_raw{iSeg};
                    gyr_n{iSeg} = gyr_raw{iSeg};
                    jrk_n{iSeg} = jerk_acc{iSeg};
                end
            end

            %% --- Séries temporelles de features (résolution 1 s) ---
            N_wins  = floor(N_frames / win_frames);
            feat_ts = struct();

            for iSeg = 1:N_seg
                sname = Seg_names{iSeg};
                sig_all = {acc_n{iSeg}, gyr_n{iSeg}, jrk_n{iSeg}};

                for iSt = 1:3
                    sig_name = sig_types{iSt};
                    xyz      = sig_all{iSt};
                    modv     = sqrt(sum(xyz.^2,2));
                    comps    = {xyz(:,1), xyz(:,2), xyz(:,3), modv};

                    for iCo = 1:4
                        comp  = comp_names{iCo};
                        x_sig = comps{iCo};

                        % ── Features temporelles : fenêtre 1s non chevauchante ──
                        tmp_t = struct();
                        for iFt = 1:numel(feat_temporal)
                            tmp_t.(feat_temporal{iFt}) = NaN(N_wins,1);
                        end
                        for iW = 1:N_wins
                            i0 = (iW-1)*win_frames + 1;
                            i1 = i0 + win_frames - 1;
                            tf = compute_temporal(x_sig(i0:i1));
                            for iFt = 1:numel(feat_temporal)
                                fn = feat_temporal{iFt};
                                tmp_t.(fn)(iW) = tf.(fn);
                            end
                        end
                        for iFt = 1:numel(feat_temporal)
                            fn = feat_temporal{iFt};
                            vn = sprintf('%s_%s_%s_%s', sig_name, comp, fn, sname);
                            feat_ts.(vn) = tmp_t.(fn);
                        end

                        % ── Features CWT : signal continu → moyenne par 1s ──
                        cf = compute_cwt_timeseries(x_sig, cwt_scales, omega0, ...
                                dt, freq_grid, cwt_thresh, win_frames, N_wins);
                        for iFt = 1:numel(feat_cwt)
                            fn = feat_cwt{iFt};
                            vn = sprintf('%s_%s_%s_%s', sig_name, comp, fn, sname);
                            feat_ts.(vn) = cf.(fn);
                        end
                    end
                end
            end

            %% --- Normalisation Goubault (si niveau 'feature') ---
            fn_list = fieldnames(feat_ts);
            if strcmpi(norm_level,'feature')
                for iV = 1:numel(fn_list)
                    feat_ts.(fn_list{iV}) = ...
                        normalize_goubault(feat_ts.(fn_list{iV}), baseline_pct, n_top_vals);
                end
            end

            %% --- Extraction aux 10 instants RPE (moyenne ±2s) ---
            rpe_frames = round(linspace(0, 1, N_bins) * N_frames);
            rpe_frames = max(1, min(rpe_frames, N_frames));

            for iV = 1:numel(fn_list)
                vname    = fn_list{iV};
                ts_val   = feat_ts.(vname);
                bin_vals = NaN(1,N_bins);

                for iBin = 1:N_bins
                    t_center = rpe_frames(iBin);
                    t_lo = t_center - half_bin_fr;
                    t_hi = t_center + half_bin_fr;
                    w_lo = max(1,      ceil(t_lo / win_frames));
                    w_hi = min(N_wins, floor(t_hi / win_frames));
                    if w_hi >= w_lo
                        bin_vals(iBin) = mean(ts_val(w_lo:w_hi), 'omitnan');
                    end
                end

                Results_Spectral.(subj).(kbd).(vname) = bin_vals;
            end

            elapsed = toc;
            fprintf('OK (%d features, %.1fs)\n', numel(fn_list), elapsed);
            clear feat_ts;
        end
    end
    clear XF tmp;
end

%% -----------------------------------------------------------------------
%  RÉSUMÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Résumé ===\n');
subjects_done = fieldnames(Results_Spectral);
N_subj_done   = numel(subjects_done);
n_cross       = 0;
for iS = 1:N_subj_done
    subj = subjects_done{iS};
    kbds = fieldnames(Results_Spectral.(subj));
    fprintf('  %-6s  %s\n', subj, strjoin(kbds,' + '));
    if numel(kbds)==2, n_cross=n_cross+1; end
end
fprintf('\n  Total : %d sujets | %d crossover\n', N_subj_done, n_cross);

if N_subj_done > 0
    s1 = subjects_done{1}; k1 = fieldnames(Results_Spectral.(s1));
    fprintf('  Variables : %d\n', numel(fieldnames(Results_Spectral.(s1).(k1{1}))));
end

n_ok=0; n_tot=0;
for iS=1:N_subj_done
    subj=subjects_done{iS};
    for kk=fieldnames(Results_Spectral.(subj))'
        for ff=fieldnames(Results_Spectral.(subj).(kk{1}))'
            v=Results_Spectral.(subj).(kk{1}).(ff{1});
            n_ok=n_ok+sum(~isnan(v)); n_tot=n_tot+numel(v);
        end
    end
end
fprintf('  Valeurs valides : %d/%d (%.1f%%)\n', n_ok, n_tot, 100*n_ok/max(n_tot,1));

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save_path = fullfile(path_save,'SpectralFeatures_Ergo3.mat');
save(save_path,'Results_Spectral','Seg_names','Seg_idx','sig_types',...
    'comp_names','feat_names','feat_temporal','feat_cwt',...
    'freq_grid','omega0','cwt_thresh','N_bins','Fs','baseline_pct','n_top_vals',...
    'win_sec','half_bin_sec','RPE_shape','norm_level',...
    'filt_order','filt_band','-v7.3');
fprintf('\nSauvegardé : SpectralFeatures_Ergo3.mat\n');
fprintf('Structure  : Results_Spectral.(Sujet).(Condition).(Variable) = [1×%d]\n', N_bins);
fprintf('NomFeature : <Signal>_<Comp>_<Feature>_<Segment>\n');
fprintf('\n=== TERMINÉ ===\n');