%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  PLOT Accel_Mod_Hand — Expression musicale                         %%%%
%%%%  Conditions : Competition + Extrait (IE/PS, takes A/B/C)          %%%%
%%%%  Chopin ignoré                                                     %%%%
%%%%                                                                    %%%%
%%%%  Couleurs : rouge = Competition | vert = IE | bleu = PS           %%%%
%%%%  Intensité : clair = A | normal = B | foncé = C                   %%%%
%%%%  Chiffre en fin de courbe : ordre de passage (1-6, Comp exclu)    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
path_data    = 'C:\Users\Nathan\OneDrive - Aix-Marseille Université\Documents\MATLAB\Piano_fatigue\Analyse_Descriptive\';
path_p01_p08 = fullfile(path_data, 'XSens_MIDI_cut_P01_P08 1.mat');
path_p09_p14 = fullfile(path_data, 'XSens_MIDI_cut_P09_P14 1.mat');
path_save    = path_data;
%  PARAMÈTRES
%  -----------------------------------------------------------------------
Subjects_p1 = {'P01','P02','P04','P05','P06','P07','P08'};
Subjects_p2 = {'P09','P10','P11','P12','P13','P14'};

Expressions   = {'IE', 'PS'};
Takes_default = {'A', 'B', 'C'};
seg_hand      = 11;   % rHand = segment 11

%% -----------------------------------------------------------------------
%  ORDRE DE PASSAGE — 1er morceau joué (IE ou PS) par participant
%  Competition toujours en 1er, non numérotée
%  Si IE 1er : IE.A=1, IE.B=2, IE.C=3, PS.A=4, PS.B=5, PS.C=6
%  Si PS 1er : PS.A=1, PS.B=2, PS.C=3, IE.A=4, IE.B=5, IE.C=6
%  -----------------------------------------------------------------------
first_expr = struct(...
    'P01','IE', 'P02','PS', 'P04','PS', 'P05','IE', 'P06','PS', 'P07','PS', ...
    'P08','IE', 'P09','IE', 'P10','PS', 'P11','PS', 'P12','IE', 'P13','IE', 'P14','PS');

%% -----------------------------------------------------------------------
%  COULEURS
%  Compétition : rouge
%  IE          : vert  (clair / normal / foncé pour A/B/C)
%  PS          : bleu  (clair / normal / foncé pour A/B/C)
%  -----------------------------------------------------------------------
col_comp = [0.85 0.10 0.10];   % rouge

% Vert : clair / normal / foncé (contraste fort)
col_IE = [0.75 1.00 0.75;   % A — très clair (menthe)
          0.15 0.65 0.15;   % B — normal (vert vif)
          0.00 0.25 0.00];  % C — très foncé (vert forêt)

% Bleu : clair / normal / foncé (contraste fort)
col_PS = [0.65 0.85 1.00;   % A — très clair (bleu ciel)
          0.10 0.35 0.85;   % B — normal (bleu roi)
          0.00 0.05 0.45];  % C — très foncé (bleu marine)

take_labels = {'A','B','C'};

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
%  CHARGEMENT ET FUSION DES DEUX FICHIERS (logique identique à ComputeJointAngles.m)
%  -----------------------------------------------------------------------
fprintf('Chargement P01-P08...\n');
tmp = load(path_p01_p08, 'XSens_final');
XSens_part1 = tmp.XSens_final;
fprintf('  OK\n');

fprintf('Chargement P09-P14...\n');
tmp = load(path_p09_p14, 'XSens_final');
XSens_part2 = tmp.XSens_final;
fprintf('  OK\n');

% Fusion manuelle des deux structs (remplace catstruct non disponible)
XSens_cur = XSens_part1;
f2 = fieldnames(XSens_part2);
for iF = 1:length(f2)
    XSens_cur.(f2{iF}) = XSens_part2.(f2{iF});
end
Subj_cur  = [Subjects_p1, Subjects_p2];
clear XSens_part1 XSens_part2 tmp;
fprintf('Fusion OK — %d sujets\n\n', length(Subj_cur));

for iFile = 1:1   % boucle conservée pour la structure, un seul passage

    for iS = 1:length(Subj_cur)
        subj = Subj_cur{iS};
        fprintf('  %s...\n', subj);

        first = first_expr.(subj);   % 'IE' ou 'PS'

        % Ordre de passage : position 1-6 pour chaque take
        % first = IE : IE.A→1, IE.B→2, IE.C→3, PS.A→4, PS.B→5, PS.C→6
        % first = PS : PS.A→1, PS.B→2, PS.C→3, IE.A→4, IE.B→5, IE.C→6
        passage_order = struct();
        if strcmp(first, 'IE')
            passage_order.IE.A = 1; passage_order.IE.B = 2; passage_order.IE.C = 3;
            passage_order.PS.A = 4; passage_order.PS.B = 5; passage_order.PS.C = 6;
        else
            passage_order.PS.A = 1; passage_order.PS.B = 2; passage_order.PS.C = 3;
            passage_order.IE.A = 4; passage_order.IE.B = 5; passage_order.IE.C = 6;
        end

        %% -----------------------------------------------------------
        %  COLLECTE DES ESSAIS
        %  -----------------------------------------------------------
        trials = struct('label',{},'color',{},'order_num',{},'signal',{});

        % --- Competition ---
        if ~strcmp(subj, 'P13') && isfield(XSens_cur.(subj), 'Competition')
            try
                acc = XSens_cur.(subj).Competition.segmentData(seg_hand).acceleration;
                trials(end+1).label     = 'Competition';
                trials(end).color       = col_comp;
                trials(end).order_num   = 0;   % 0 = Competition, pas numérotée 1-6
                trials(end).signal      = sqrt(sum(acc.^2, 2));
            catch
                fprintf('    Competition manquante pour %s\n', subj);
            end
        end

        % --- Extrait IE et PS ---
        for iE = 1:length(Expressions)
            expr = Expressions{iE};

            if strcmp(subj, 'P14')
                Takes = {'B'};
            elseif ismember(subj, {'P05','P11'})
                Takes = {'A','B'};
            else
                Takes = Takes_default;
            end

            for iT = 1:length(Takes)
                take = Takes{iT};
                take_idx = find(strcmp(take_labels, take));

                try
                    acc = XSens_cur.(subj).Extrait.(expr).(take).segmentData(seg_hand).acceleration;
                    sig = sqrt(sum(acc.^2, 2));

                    if strcmp(expr, 'IE')
                        col = col_IE(take_idx, :);
                    else
                        col = col_PS(take_idx, :);
                    end

                    ord = passage_order.(expr).(take);

                    trials(end+1).label     = sprintf('Extrait %s %s', expr, take);
                    trials(end).color       = col;
                    trials(end).order_num   = ord;
                    trials(end).signal      = sig;
                catch
                    fprintf('    Manquant : %s Extrait %s %s\n', subj, expr, take);
                end
            end
        end

        if isempty(trials)
            fprintf('    Aucun essai valide — skip\n');
            continue;
        end

        %% -----------------------------------------------------------
        %  DÉCOUPAGE EN 100 BINS (moyenne par intervalle de 1%)
        %  -----------------------------------------------------------
        N_bins = 20;
        x_bins = linspace(0.5, 99.5, N_bins);

        for iT = 1:length(trials)
            sig  = trials(iT).signal;
            N_i  = length(sig);
            bins = NaN(1, N_bins);
            for iBin = 1:N_bins
                i0 = round((iBin-1)/N_bins * N_i) + 1;
                i1 = round(iBin/N_bins * N_i);
                i0 = max(1, i0); i1 = min(N_i, i1);
                bins(iBin) = mean(sig(i0:i1), 'omitnan');
            end
            trials(iT).signal_interp = bins;
        end

        %% -----------------------------------------------------------
        %  FIGURE
        %  -----------------------------------------------------------
        fig = figure('Name', subj, 'NumberTitle','off', ...
            'Position',[50 50 1500 520], 'Color','white');
        hold on;

        h_legend = [];
        lbl_legend = {};

        for iT = 1:length(trials)
            sig  = trials(iT).signal_interp;
            col  = trials(iT).color;
            ord  = trials(iT).order_num;
            lbl  = trials(iT).label;

            % Épaisseur de ligne selon le take
            if contains(lbl, ' A') || strcmp(lbl, 'Competition')
                lw = 1.8;
            elseif contains(lbl, ' B')
                lw = 1.3;
            else
                lw = 0.9;
            end

            h = plot(x_bins, sig, '-', 'Color', col, 'LineWidth', lw);

            % Chiffre en fin de courbe
            x_end = x_bins(end);
            y_end = sig(end);
            if ord == 0
                txt = 'C';   % Competition
            else
                txt = num2str(ord);
            end
            text(x_end + 0.3, y_end, txt, 'FontSize', 7, 'Color', col, ...
                'FontWeight', 'bold', 'VerticalAlignment', 'middle');

            % Légende : une entrée par condition/take principale
            if contains(lbl, ' A') || strcmp(lbl, 'Competition')
                h_legend(end+1) = h;
                lbl_legend{end+1} = lbl;
            end
        end

        % Annotation ordre de passage
        ord_txt = sprintf('1er joué : %s  (1=1er passage … 6=dernier)', first);
        annotation(fig,'textbox',[0.01 0.92 0.6 0.05], ...
            'String', ord_txt, 'FontSize', 8, 'EdgeColor','none', 'Color',[0.3 0.3 0.3]);

        legend(h_legend, lbl_legend, 'Location','northeast', 'FontSize',7);
        xlabel('% de la session', 'FontSize',10);
        ylabel('Accel Mod Main (m/s²)', 'FontSize',10);
        title(sprintf('%s — Accel\\_Mod\\_Hand | Comp(rouge) | IE(vert) | PS(bleu) | clair=A normal=B foncé=C | chiffre=ordre de passage', subj), ...
            'FontSize', 9, 'FontWeight','bold');
        xlim([0 100]); grid on; box on;

        fname = fullfile(path_save, sprintf('Fig_AccelHand_%s.png', subj));
        saveas(fig, fname);
        fprintf('    Fig_AccelHand_%s.png sauvegardée\n', subj);
        close(fig);
    end

end

fprintf('\n=== TERMINÉ ===\n');