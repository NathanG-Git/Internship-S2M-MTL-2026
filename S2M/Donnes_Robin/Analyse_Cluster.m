% --- Analyse_Cluster_Par_Bloc.m ---
% 1. Filtrage des sujets valides (Common entre Norm & CS60, sans P07)
sujets = fieldnames(Results_IMU);
sujets_communs = {};
for i = 1:length(sujets)
    f = fieldnames(Results_IMU.(sujets{i}));
    if ismember('Norm', f) && ismember('CS60', f) && ~strcmp(sujets{i}, 'P07')
        sujets_communs{end+1} = sujets{i};
    end
end

% 2. Définition des blocs (noms exacts des champs)
blocs.SegMod   = {'Accel_Mod_Head', 'Accel_Mod_Hand'};
blocs.Goubault = {'MedianFreq_Accel_Y_Hand', 'PeakPower_Accel_Mod_Hand', ...
                  'PeakPower_AngVel_X_Head', 'SpectralEntropy_AngVel_Mod_Forearm'};
blocs.EMG      = {'TFR_MedianFreq_Triceps', 'TFR_SpectralEntropy_Deltoid', 'SampleEntropy_Biceps'};
nom_blocs = fieldnames(blocs);

% 3. Analyse par bloc
for b = 1:length(nom_blocs)
    vars = blocs.(nom_blocs{b});
    X = zeros(length(sujets_communs), length(vars));
    
    for i = 1:length(sujets_communs)
        subj = sujets_communs{i};
        % Calcul de la dérive (Norm uniquement pour l'instant)
        for v = 1:length(vars)
            d = Results_IMU.(subj).Norm.(vars{v}); 
            X(i,v) = d(10) - d(1);
        end
    end
    
    % Nettoyage NaN éventuels
    X(any(isnan(X), 2), :) = []; 
    
    % 4. Clustering (K-means)
    X_norm = zscore(X);
    idx = kmeans(X_norm, 2, 'Replicates', 20);
    
    % 5. Corrélation de Spearman (en supposant Y_drpe_norm aligné)
    % Vous devrez extraire Y_drpe pour les sujets_communs ici
    fprintf('\n--- Bloc : %s (N=%d) ---\n', nom_blocs{b}, length(sujets_communs));
    % ... rho = corr(...) ...
end