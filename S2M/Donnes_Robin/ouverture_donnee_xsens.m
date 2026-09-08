

fprintf('Chargement des fichiers Xsens en cours (cela peut prendre quelques minutes)...\n');

% On charge chaque fichier dans une structure indépendante (F1, F2, etc.)
F1 = load('C:\Users\jerem\Documents\Projet de stage\DATA_Jeremy\April28\XSens_Fatigue_1.mat');
fprintf('Fichier 1 chargé.\n');

F2 = load('C:\Users\jerem\Documents\Projet de stage\DATA_Jeremy\April28\XSens_Fatigue_2.mat');
fprintf('Fichier 2 chargé.\n');

F3 = load('C:\Users\jerem\Documents\Projet de stage\DATA_Jeremy\April28\XSens_Fatigue_3.mat');
fprintf('Fichier 3 chargé.\n');

F4 = load('C:\Users\jerem\Documents\Projet de stage\DATA_Jeremy\April28\XSens_Fatigue_4.mat');
fprintf('Fichier 4 chargé.\n');

F5 = load('C:\Users\jerem\Documents\Projet de stage\DATA_Jeremy\April28\XSens_Fatigue_5.mat');
fprintf('Fichier 5 chargé.\n');

fprintf('Tous les fichiers sont chargés avec succès !\n');