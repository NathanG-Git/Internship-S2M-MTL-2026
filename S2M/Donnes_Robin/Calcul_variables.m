%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  CALCUL DES FEATURES FEATURES3 DEPUIS SIGNAL BRUT IMU             %%%%
%%%%  Brouillon — S101 — Tâche Li                                      %%%%
%%%%                                                                    %%%%
%%%%  15 features × 6 signaux :                                        %%%%
%%%%  Stats : Peak, Mean, Median, STD, P10, P25, P75, P90             %%%%
%%%%  CWT   : MedianFreq, SpectEnt, PowerLF, PowerHF, PowerTot,       %%%%
%%%%          PeakPower, PeakPowerFreq                                 %%%%
%%%%  Signaux : Accel, AngVel, Jerk (axes + module)                   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS ET FONCTIONS
%  -----------------------------------------------------------------------
path_raw   = 'J:\Piano_Fatigue\Data_Exported\Xsens\';
path_info  = 'J:\Piano_Fatigue\Matlab_matrix\Info_participants\';
path_cyc   = 'J:\Piano_Fatigue\Data_Exported\';
path_save  = fileparts(mfilename('fullpath'));
if isempty(path_save), path_save = pwd; end

addpath('C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Donnees_Robin');

iP_str = 'S101';

%% -----------------------------------------------------------------------
%  PARAMÈTRES CWT
%  -----------------------------------------------------------------------
Fs_acq    = 60;
Fs_cwt    = 60;    % Fréquence CWT — conforme pipeline original (60 Hz)
ond       = 'cmor8-1';
fc_morlet = 1.0;
f         = 0.05:0.05:15;
scale_cwt = Fs_cwt * fc_morlet ./ f;

% Indices fréquences — dynamiques, indépendants de Fs
% LowPower  : 0.1–4 Hz
% HighPower : 6–12 Hz  (pipeline original : idx 120:240 à 40Hz = 6–12 Hz)
% TotalPower: 0.1–12 Hz
idx_lf_lo  = find(f >= 0.1,  1);
idx_lf_hi  = find(f >= 4.0,  1);
idx_hf_lo  = find(f >= 6.0,  1);   % 6 Hz, pas 4 Hz
idx_hf_hi  = find(f >= 12.0, 1);
idx_tot_lo = idx_lf_lo;
idx_tot_hi = idx_hf_hi;

fprintf('Bandes (Fs=%d Hz) : LF=%.1f-%.1f Hz | HF=%.1f-%.1f Hz | Tot=%.1f-%.1f Hz\n', ...
    Fs_cwt, f(idx_lf_lo), f(idx_lf_hi), f(idx_hf_lo), f(idx_hf_hi), f(idx_tot_lo), f(idx_tot_hi));

%% -----------------------------------------------------------------------
%  CHARGEMENT
%  -----------------------------------------------------------------------
fprintf('Chargement S101...\n');

raw   = load(fullfile(path_raw, 'S101_Li_Accel.mat'));
Accel = raw.Accel.data;            % [17067 x 69] — 23 seg x 3 axes

gyr   = load(fullfile(path_raw, 'S101_Li_Free.mat'));
AngVel = gyr.FreeAccel.data;       % [17067 x 51] — 17 seg x 3 axes

info  = load(fullfile(path_info, 'Info_participants_corrected.mat'));
RPE_raw = info.Info_participants(1).Liszt(:);

cyc   = load(fullfile(path_cyc, 'Cycle_Li_Felipe.mat'));
seq   = cyc.cycles(1).seq;         % [97 x 5] timestamps en secondes
N_cyc = size(seq, 1);

fprintf('  Accel : %dx%d | AngVel : %dx%d | Cycles : %d\n', ...
    size(Accel,1), size(Accel,2), size(AngVel,1), size(AngVel,2), N_cyc);

%% -----------------------------------------------------------------------
%  RÉÉCHANTILLONNAGE ET JERK
%  -----------------------------------------------------------------------
fprintf('Calcul Jerk (pas de rééchantillonnage — Fs_cwt = Fs_acq = 60 Hz)...\n');

t_orig = (0:size(Accel,1)-1) / Fs_acq;
t_new  = t_orig;   % même grille temporelle
N_new  = length(t_new);

Accel_rs  = Accel;
AngVel_rs = AngVel;

% Jerk = diff(Accel) * Fs
Jerk_rs = [diff(Accel_rs, 1, 1) * Fs_cwt; zeros(1, size(Accel_rs,2))];

fprintf('  Signal : %d frames à %d Hz\n', N_new, Fs_cwt);

%% -----------------------------------------------------------------------
%  DÉFINITION DES SIGNAUX À TRAITER
%  -----------------------------------------------------------------------
% Pour ce brouillon : on prend les 7 premiers segments disponibles
% (mapping exact à faire sur les vraies données)

% Indices colonnes pour 7 segments (on prend cols non-nulles)
% Accel : colonnes 4-24 pour les premiers segments avec signal
% On prend simplement les colonnes 4,7,10,13,16,19,22 (axe X de 7 seg)
seg_cols_accel  = [4,7,10,13,16,19,22];   % axe X des 7 segments
seg_cols_angvel = [4,7,10,13,16,19,22];   % idem pour AngVel (brouillon)
N_seg = 7;
Seg_names = {'Seg1','Seg2','Seg3','Seg4','Seg5','Seg6','Seg7'};

% Construction des 6 signaux [N_frames x 7] pour axes X
% et [N_frames x 7] pour modules
Signals = struct();

% --- Accélération axe X ---
Signals.Acceleration = Accel_rs(:, seg_cols_accel);

% --- Accélération module ---
Mod_Accel = zeros(N_new, N_seg);
for iSeg = 1:N_seg
    cx = (iSeg-1)*3 + 4;   % col X du segment
    Mod_Accel(:,iSeg) = sqrt(Accel_rs(:,cx).^2 + ...
                              Accel_rs(:,cx+1).^2 + ...
                              Accel_rs(:,cx+2).^2);
end
Signals.Module_Acceleration = Mod_Accel;

% --- Vitesse angulaire axe X ---
Signals.Angular_Velocity = AngVel_rs(:, seg_cols_angvel);

% --- Vitesse angulaire module ---
Mod_AngVel = zeros(N_new, N_seg);
for iSeg = 1:N_seg
    cx = (iSeg-1)*3 + 4;
    cx = min(cx, size(AngVel_rs,2)-2);  % sécurité taille
    Mod_AngVel(:,iSeg) = sqrt(AngVel_rs(:,cx).^2 + ...
                               AngVel_rs(:,cx+1).^2 + ...
                               AngVel_rs(:,cx+2).^2);
end
Signals.Module_Angular_Velocity = Mod_AngVel;

% --- Jerk axe X ---
Signals.Jerk = Jerk_rs(:, seg_cols_accel);

% --- Jerk module ---
Mod_Jerk = zeros(N_new, N_seg);
for iSeg = 1:N_seg
    cx = (iSeg-1)*3 + 4;
    Mod_Jerk(:,iSeg) = sqrt(Jerk_rs(:,cx).^2 + ...
                             Jerk_rs(:,cx+1).^2 + ...
                             Jerk_rs(:,cx+2).^2);
end
Signals.Module_Jerk = Mod_Jerk;

sig_names = fieldnames(Signals);
N_sig = length(sig_names);
fprintf('  %d signaux définis\n', N_sig);

%% -----------------------------------------------------------------------
%  CALCUL DES FEATURES PAR CYCLE ET PAR SIGNAL
%  -----------------------------------------------------------------------
fprintf('\nCalcul des features (%d cycles x %d signaux)...\n', N_cyc, N_sig);

% Noms des features
feat_stats = {'Peak','Mean','Median','STD','P10','P25','P75','P90'};
feat_cwt   = {'MedianFreq','SpectralEntropy','PowerLF','PowerHF', ...
              'PowerTot','PeakPower','PeakPower_Freq'};
feat_all   = [feat_stats, feat_cwt];
N_feat     = length(feat_all);

% Initialisation : une struct Features_Calc.(feat).(signal) = [N_cyc x N_seg]
Features_Calc = struct();
for iF = 1:N_feat
    for iS = 1:N_sig
        Features_Calc.(feat_all{iF}).(sig_names{iS}) = NaN(N_cyc, N_seg);
    end
end

% Boucle principale
for iCyc = 1:N_cyc
    if mod(iCyc, 20) == 0
        fprintf('  Cycle %d/%d...\n', iCyc, N_cyc);
    end

    % Extraire les frames de ce cycle
    t_lo  = seq(iCyc, 1);
    t_hi  = seq(iCyc, 5);
    idx_c = find(t_new >= t_lo & t_new <= t_hi);

    if length(idx_c) < 5, continue; end

    for iS = 1:N_sig
        sname  = sig_names{iS};
        sig_mat = Signals.(sname);   % [N_frames x N_seg]

        for iSeg = 1:N_seg
            x = sig_mat(idx_c, iSeg);

            % Remplacer NaN par interpolation
            if any(isnan(x))
                idx_ok = find(~isnan(x));
                if length(idx_ok) < 2, continue; end
                x = interp1(idx_ok, x(idx_ok), 1:length(x), 'linear', 'extrap');
                x = x(:);
            end

            % --- STATS SIMPLES ---
            Features_Calc.Peak.(sname)(iCyc,iSeg)   = max(abs(x));
            Features_Calc.Mean.(sname)(iCyc,iSeg)   = mean(x);
            Features_Calc.Median.(sname)(iCyc,iSeg) = median(x);
            Features_Calc.STD.(sname)(iCyc,iSeg)    = std(x);
            Features_Calc.P10.(sname)(iCyc,iSeg)    = prctile(x, 10);
            Features_Calc.P25.(sname)(iCyc,iSeg)    = prctile(x, 25);
            Features_Calc.P75.(sname)(iCyc,iSeg)    = prctile(x, 75);
            Features_Calc.P90.(sname)(iCyc,iSeg)    = prctile(x, 90);

            % --- CWT : somme des 3 axes (conforme pipeline Features3) ---
            % TFR = CWT(x) + CWT(y) + CWT(z)
            % Pour les signaux modules : on recalcule depuis les 3 axes
            try
                % Récupérer les 3 axes selon le type de signal
                if contains(sname, 'Module')
                    % Signal module : retrouver les 3 axes depuis le signal brut
                    if contains(sname, 'Acceleration')
                        cx = (iSeg-1)*3 + 4;
                        x_ax = Accel_rs(idx_c, cx);
                        y_ax = Accel_rs(idx_c, cx+1);
                        z_ax = Accel_rs(idx_c, cx+2);
                    elseif contains(sname, 'Angular')
                        cx = (iSeg-1)*3 + 4;
                        cx = min(cx, size(AngVel_rs,2)-2);
                        x_ax = AngVel_rs(idx_c, cx);
                        y_ax = AngVel_rs(idx_c, cx+1);
                        z_ax = AngVel_rs(idx_c, cx+2);
                    else  % Module_Jerk
                        cx = (iSeg-1)*3 + 4;
                        x_ax = Jerk_rs(idx_c, cx);
                        y_ax = Jerk_rs(idx_c, cx+1);
                        z_ax = Jerk_rs(idx_c, cx+2);
                    end
                else
                    % Signal axe unique : dupliquer sur les 3 axes
                    % (approximation pour brouillon — à corriger avec vrai mapping)
                    x_ax = x; y_ax = x; z_ax = x;
                end

                % CWT sur chaque axe puis somme
                coef_x = cwt(x_ax(:), scale_cwt, ond, 'ExtendSignal', 1);
                coef_y = cwt(y_ax(:), scale_cwt, ond, 'ExtendSignal', 1);
                coef_z = cwt(z_ax(:), scale_cwt, ond, 'ExtendSignal', 1);
                TFR_full = abs(coef_x) + abs(coef_y) + abs(coef_z);

                if size(TFR_full,1) < idx_hf_hi, continue; end

                TFR  = TFR_full(idx_lf_lo:idx_hf_hi, :);
                Freq = f(idx_lf_lo:idx_hf_hi);

                if any(isnan(TFR(:))), continue; end

                % Puissances — indices basés sur Freq (vecteur local)
                idx_lf = Freq >= 0.1 & Freq <= 4.0;
                idx_hf = Freq >= 6.0 & Freq <= 12.0;  % 6 Hz, pas 4 Hz
                Features_Calc.PowerLF.(sname)(iCyc,iSeg)  = mean(mean(TFR(idx_lf,:)));
                Features_Calc.PowerHF.(sname)(iCyc,iSeg)  = mean(mean(TFR(idx_hf,:)));
                Features_Calc.PowerTot.(sname)(iCyc,iSeg) = mean(mean(TFR));

                % PeakPower
                mean_t = mean(TFR, 2);
                [pk_val, pk_idx] = max(mean_t);
                Features_Calc.PeakPower.(sname)(iCyc,iSeg)      = pk_val;
                Features_Calc.PeakPower_Freq.(sname)(iCyc,iSeg) = Freq(pk_idx);

                % MedianFreq
                mf = Compute_Median_Frequency(TFR, Freq);
                Features_Calc.MedianFreq.(sname)(iCyc,iSeg) = mean(mf(~isnan(mf)));

                % SpectralEntropy
                se = Compute_Spectral_Entropy(TFR, Freq);
                Features_Calc.SpectralEntropy.(sname)(iCyc,iSeg) = mean(se(~isnan(se)));

            catch
                % CWT échouée — NaN conservé
            end
        end
    end
end

%% -----------------------------------------------------------------------
%  VÉRIFICATION : COMPARAISON AVEC FEATURES3
%  -----------------------------------------------------------------------
fprintf('\n=== Vérification vs Features3 ===\n');
try
    path_f3 = 'J:\Piano_Fatigue\Data_Exported\Features_XSENS_Cycles\Li\';
    f3 = load(fullfile(path_f3, 'Features3_XSENS_Cycles_S101_Li.mat'));
    F3 = f3.Features_XSENS_Cycles;

    % Comparer MedianFreq sur Module_Acceleration (premier segment)
    our_mf = Features_Calc.MedianFreq.Module_Acceleration(:,1);
    ref_mf = F3.MedianFreq.Module_Acceleration(:,1);
    n_comp = min(length(our_mf), length(ref_mf));
    r = corr(our_mf(1:n_comp), ref_mf(1:n_comp), 'Rows', 'complete');
    fprintf('  MedianFreq Module_Accel seg1 : r=%.4f (N=%d cycles)\n', r, n_comp);

    % Comparer Mean
    our_mn = Features_Calc.Mean.Module_Acceleration(:,1);
    ref_mn = F3.Mean.Module_Acceleration(:,1);
    r2 = corr(our_mn(1:n_comp), ref_mn(1:n_comp), 'Rows', 'complete');
    fprintf('  Mean Module_Accel seg1       : r=%.4f\n', r2);

catch ME
    fprintf('  Comparaison impossible : %s\n', ME.message);
end

%% -----------------------------------------------------------------------
%  AFFICHAGE RÉSUMÉ
%  -----------------------------------------------------------------------
fprintf('\n=== Résumé des features calculées ===\n');
fprintf('%-20s', 'Feature');
for iS = 1:N_sig
    fprintf('%-22s', sig_names{iS});
end
fprintf('\n%s\n', repmat('-', 1, 20 + 22*N_sig));

for iF = 1:N_feat
    fprintf('%-20s', feat_all{iF});
    for iS = 1:N_sig
        mat = Features_Calc.(feat_all{iF}).(sig_names{iS});
        n_ok = sum(~isnan(mat(:)));
        n_tot = numel(mat);
        fprintf('%-22s', sprintf('%d/%d valides', n_ok, n_tot));
    end
    fprintf('\n');
end

%% -----------------------------------------------------------------------
%  SAUVEGARDE
%  -----------------------------------------------------------------------
save(fullfile(path_save, sprintf('Features_Calc_%s_Li.mat', iP_str)), ...
    'Features_Calc', 'Seg_names', 'sig_names', 'feat_all');
fprintf('\nSauvegardé : Features_Calc_%s_Li.mat\n', iP_str);