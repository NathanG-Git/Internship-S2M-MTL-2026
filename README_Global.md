# README — Pipeline d'analyse de la fatigue pianistique

Ce dépôt contient l'ensemble des scripts MATLAB développés lors du stage M1 au S2M Laboratory (Université de Montréal) pour l'identification de biomarqueurs biomécanique et neuromusculaires de la fatigue pianistique.

---

## Structure du dépôt

### `Analyse_descriptive/`
Scripts d'**analyse de passages musicaux ciblés** sur le dataset Expression musicale (13 pianistes, conditions Competition/IE/PS). À partir de passages identifiés manuellement dans les enregistrements, les scripts découpent les fichiers audio (WAV) et MIDI correspondants et génèrent des projets Reaper (`.RPP`) autonomes pour chaque passage, permettant l'écoute synchronisée audio + MIDI dans Reaper.

---

### `Correlation/`
Scripts de **corrélation de Spearman** entre les variables biomécaniques/EMG et le RPE, sur le dataset Goubault. Comprend deux paires de scripts (avec et sans correction FDR Benjamini-Hochberg) pour les variables cinématiques (format cycles) et les features Goubault. Ces scripts ont fourni les premières listes de candidats biomarqueurs avant le passage au PLSR.

---

### `Donnes_Robin/`
Scripts de **calcul de features et d'analyse statistique** sur le dataset Robin2/Ergo (23 sujets, clavier Normal vs CS60 ergonomique). Contient le pipeline LMM de validation externe (`Analyse_Validation.m`), les scripts de calcul de features IMU et EMG (`Analyse_IMU.m`, `Analyse_EMG.m`, `Compute_*`), les analyses RPE et sous-groupes taille de main, ainsi que le pipeline Spearman cross-dataset sur 5 groupes et 93 variables avec filtrage IM et modélisation GAMM.

---

### `LD vs SD/`
Scripts d'**analyse comparative entre petites mains (SD) et grandes mains (LD)**, condition Normal uniquement (dataset Robin2). Comprend les tests Mann-Whitney sur le RPE, les corrélations dérive↔ΔRPE par groupe, les corrélations partielles contrôlant RPE₀, et les forest plots associés. Ces scripts testent l'hypothèse que SD ≈ G1 (fatigue plus rapide) en raison d'une contrainte mécanique relative plus grande sur clavier standard.

---

### `Lasso_LMM/`
Scripts de **sélection de variables et modélisation** sur le dataset Goubault. Contient trois versions successives du pipeline (LASSO classique → Group LASSO FISTA + Mundlak HAC → variante multi-fenêtres temporelles), un script autonome de forest plots, et un script de mise en page pour publication. L'objectif est d'identifier les biomarqueurs prédicteurs du RPE au niveau individuel avec contrôle de l'autocorrélation intra-sujet.

---

### `Li_long_duration/`
Scripts de **calcul des features de référence pour G2** (Long Duration, dataset Goubault) depuis les données brutes XSens et EMG. Produit les trois fichiers `.mat` chargés sous `d2`, `dG2` et `eG2` dans tous les scripts d'analyse : variables cinématiques temporelles, 10 variables Goubault, et 6 variables EMG × 5 canaux.

---

### `Li_short_duration/`
Scripts de **calcul des features de référence pour G1** (Short Duration, dataset Goubault) et scripts de corrélation exploratoires. Produit les trois fichiers `.mat` chargés sous `d1`, `dG1` et `eG1`. Contient également `AccelSign_JointData_Li.m` (accélération propulsive vs freinatrice) et les scripts Spearman sur les variables Workload et Goubault.

---

### `PLSR/`
Scripts du **pipeline PLSR + VIP itératif** sur le dataset Robin2/Ergo. Comprend le calcul des 1260 features spectrales en repère global (`Compute_SpectralFeatures_Ergo.m`, version de référence Ergo3) et en repère local (`Compute_SpectralFeatures_Ergo_Local.m`), les deux scripts PLSR (avec et sans EMG), et le script de test de corrélation pour la réduction de redondance entre variables. Les résultats comparatifs des différents runs sont compilés dans `Dimensions_Variables1.ods`.
