% ==============================
% Tableau style publication (coefficients)
% ==============================

clear; clc;

% Données
headers = {'Variable', 'Bloc', '\beta_{raw}', '\beta_{std}', 'p'};

data = {
    'Accel Head',               'SegMod',    '+2.303', '+1.093', '***';
    'Accel Hand',               'SegMod',    '-2.359', '-1.120', '***';
    'Hand MedianFreq Accel_Y',  'Goubault',  '+2.013', '+0.956', '***';
    'Hand PeakPower Accel_Module', 'Goubault','-2.072', '-0.984', '***';
    'Head PeakPower AngVel_X',  'Goubault',  '+1.169', '+0.555', '***';
    'Forearm SpectralEntropy AngVel_Mod', 'Goubault', '+1.890', '+0.897', '***';
    'SampleEntropy Biceps',     'EMG',       '-1.945', '-0.879', '**';
    'TFR_MedianFreq Biceps',    'EMG',       '-0.705', '-0.319', '*';
    'TFR_MedianFreq Triceps',   'EMG',       '-1.282', '-0.579', '***';
    'TFR_SpectralEntropy DeltAnt','EMG',     '-1.309', '-0.591', '***';
};

title_txt = 'Coefficients de régression par bloc';

% Taille
nCols = numel(headers);
nRows = size(data, 1);

% Figure
fig = figure('Color', 'w', 'Units', 'pixels', 'Position', [100 100 1200 420]);
ax = axes(fig, 'Position', [0.04 0.08 0.92 0.82]);
axis(ax, [0 1 0 1]);
axis(ax, 'off');
hold(ax, 'on');

% Titre
text(0.5, 0.93, title_txt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'tex');

% Géométrie du tableau
left = 0.04;
right = 0.96;
top = 0.78;
rowH = 0.095;
tableWidth = right - left;
colW = tableWidth / nCols;

% Règles horizontales (style booktabs)
plot([left right], [top top], 'k-', 'LineWidth', 1.2);                    % top rule
plot([left right], [top-rowH top-rowH], 'k-', 'LineWidth', 0.8);          % mid rule (sous l'en-tête)
plot([left right], [top-(nRows+1)*rowH top-(nRows+1)*rowH], ...
    'k-', 'LineWidth', 1.2);                                              % bottom rule

% En-têtes
for c = 1:nCols
    x = left + (c-0.5)*colW;
    y = top - rowH/2;
    text(x, y, headers{c}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'tex');
end

% Corps du tableau
for r = 1:nRows
    for c = 1:nCols
        x = left + (c-0.5)*colW;
        y = top - rowH - (r-0.5)*rowH;
        
        % Alignement : colonne Variable (gauche), autres (centre)
        if c == 1
            halign = 'left';
            x = x + 0.01; % petit décalage vers la droite
        else
            halign = 'center';
        end
        
        text(x, y, data{r,c}, ...
            'HorizontalAlignment', halign, ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Times New Roman', ...
            'FontSize', 12, ...
            'Interpreter', 'tex');
    end
end

% Export
exportgraphics(ax, 'tableau_coefficients_matlab.png', 'Resolution', 300);