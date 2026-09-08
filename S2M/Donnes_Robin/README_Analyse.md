# README — Scripts d'analyse

Ce groupe contient tous les scripts qui **chargent des features pré-calculées et produisent des résultats statistiques et des figures**. Ils n'écrivent pas de nouvelles variables brutes.

---

## Vue d'ensemble et chaîne d'analyse

```
Étape 1 — Calcul des features (scripts Compute_* / Analyse_EMG / Analyse_IMU)
        ↓
Étape 2 — Validation LMM sur le dataset Robin2 (Normal vs CS60)
        Analyse_EMG.m  → EMG_Features_NewDataset.mat
        Analyse_Validation.m  ← IMU_Features_NewDataset.mat + EMG_Features_NewDataset.mat
              → LMM_Validation_Results.mat (LMM_Results, Index_Results, Derive_Results, RPE)
        ↓
Étape 3 — Analyse RPE et sous-groupes taille de main
        Analyse_RPE.m  ← LMM_Validation_Results.mat + données brutes
              → RPE_Analysis_Results.mat + forest plots
        Analyse_RPE_SDvsLD.m  ← LMM_Validation_Results.mat + IMU + EMG
              → RPE_SD_LD_Results.mat + scatter plots
        ↓
Étape 4 — Analyse Spearman cross-dataset (5 groupes)
        Analyse_Spearman.m  ← Goubault (jeu 1) + Robin2 (jeu 2, tous fichiers features)
              → Derive_DRPE_Results.mat + tableaux PNG
        ↓
Étape 5 — Filtrage IM + modélisation non-linéaire
        Analyse_IM_GAM.m  ← Derive_DRPE_Results.mat
              → IM_GAM_Results.mat + Tableau_IM_FDR.png + Fig_GAM_*.png
        Analyse_IM_GAM_Bins.m  ← Derive_DRPE_Results.mat + IM_GAM_Results.mat
              → GAMM_Bins_Results.mat + Fig_GAMM_*.png / Fig_Drift_*.png

Scripts exploratoires (pas de sortie persistante) :
        Analyse_Cluster.m
```

---

## Abréviations critiques — disambiguïsation SD/LD vs G1/G2

> **Attention — deux usages du même sigle SD :**
>
> Dans `Analyse_RPE_SDvsLD.m` et `Analyse_Spearman.m` (dataset Robin2) :
> - **SD = Small Hand** (petites mains, groupe contraint par le clavier standard)
> - **LD = Large Hand** (grandes mains)
>
> Dans `Spearman_Corr_G1_G2.m` et `Indiv_Spearman_Corr_G1_G2.m` (dataset Goubault) :
> - **G1 = Short Duration** (participants ayant atteint RPE ≥ 7 tôt)
> - **G2 = Long Duration** (participants n'ayant jamais atteint le seuil)
>
> L'analogie entre SD/LD (taille de main) et G1/G2 (régime de fatigue) est mentionnée dans les commentaires d'`Analyse_RPE_SDvsLD.m` comme hypothèse de travail — les petites mains subissent une contrainte mécanique relative plus grande sur clavier standard, d'où une fatigue potentiellement plus rapide comme G1.

---

## Scripts d'analyse

---

### `Analyse_Validation.m` — LMM principal, dataset Robin2

**Ce que fait le script :** pipeline central de la validation externe. Pour chaque variable des 3 blocs (SegMod, Goubault, EMG), ajuste plusieurs modèles LMM et produit des figures de trajectoires.

**Données d'entrée :** `IMU_Features_NewDataset.mat` + `EMG_Features_NewDataset.mat` (produits par `Analyse_IMU.m` et `Analyse_EMG.m`).

**Variables analysées (9) :**
- Bloc SegMod (2) : `Accel_Mod_Head`, `Accel_Mod_Hand`
- Bloc Goubault (4) : `MedianFreq_Accel_Y_Hand`, `PeakPower_Accel_Mod_Hand`, `PeakPower_AngVel_X_Head`, `SpectralEntropy_AngVel_Mod_Forearm`
- Bloc EMG (3) : `TFR_MedianFreq_Triceps`, `TFR_SpectralEntropy_Deltoid`, `SampleEntropy_Biceps`

**Modèles LMM ajustés, dans l'ordre :**

**Modèle principal** (`LMM_Results`) :
`Variable ~ Temps + Condition + Temps:Condition + (1|Participant)`
- Variables z-scorées globalement. Effet Temps : tous les sujets en Normal (N≈23–24). Effet Condition + interaction : sujets crossover uniquement (N=13).

**LMM effet d'ordre** (`LMM_Order_Results`) :
`Variable ~ Temps + Condition + Order + Temps:Condition + Temps:Order + Condition:Order + (1|Participant)`
- N=13 crossover uniquement. Comparaison via LRT (modèle réduit sans Order vs complet) à 3 degrés de liberté. Order = quelle condition a été jouée en premier (Normal ou CS60). Motivation : 6/9 biomarqueurs étaient confondus par l'ordre de passage. Ajusté en ML (pas REML) pour que le LRT soit valide.

**LMM SD vs LD** (`LMM_SDvsLD_Results`) :
`Variable ~ Temps + Group + Temps:Group + (1|Participant)`
- Deux variantes : complète (SD=Subject_Both N=13, LD N=10) et restreinte (SD=sujets ayant joué Normal en premier, N=5, LD N=10). La variante restreinte contrôle l'ordre de passage. LRT à 2 df (Group + Temps:Group).

**Index de bloc** (`Index_Results`, `LMM_Index`) : moyenne pondérée des z-scores signés des variables d'un même bloc. Calculé en 8 variantes : 2 jeux de signes (Goubault = directions attendues d'après Goubault G1 ; observés = direction réelle du β_Temps dans ce dataset) × 4 schémas de pondération (non pondéré, pondéré par β_std Goubault, par fréquence de sélection sur 100 runs PLSR, par β×fréquence).

**Modèles RPE** (sur les crossover, N=13, N_obs=26) :
- **RPE + ordre** : `ΔRPE ~ Condition + Order + Condition:Order + (1|Participant)` — LRT 2df
- **RPE₀ + ordre** : `RPE_début ~ Condition + Order + Condition:Order + (1|Participant)` — LRT 2df
- **Modèle B** : `RPE_fin ~ RPE_début + Condition + (1|Participant)` — RPE₀ influence-t-il le niveau final ? LRT 1df sur tous sujets (N_obs≤36)
- **Modèle C** : `ΔRPE ~ RPE_début + Condition + (1|Participant)` — RPE₀ module-t-il la pente (pas seulement le niveau) ? Complément algébrique du modèle B.
- **Modèle D** : `ΔRPE ~ Dérive + Condition + Order + Dérive:Order + (1|Participant)` — le lien Dérive↔ΔRPE dépend-il de l'ordre ? LRT 2df.

**Figures produites :**
- `Fig_Trajectoires_<bloc>.png` — courbes moyennes par condition (Normal N=13, CS60 N=13, Normal N=23, Normal LD N=10), normalisées par rapport au bin 1.
- `Fig_Index_<bloc>.png` — idem pour les 8 variantes d'index.
- `Fig_RPE_Scatter_<bloc>.png` — scatter dérive ↔ ΔRPE pour 4 séries (Normal 13, CS60, Normal 23, LD).
- `Fig_Tableau_LMM.png` — tableau récapitulatif LMM (9 variables).
- `Fig_Tableau_LMM_Index.png` — tableau récapitulatif LMM index (3 blocs × 8 variantes).

**Sortie principale :** `LMM_Validation_Results.mat` (contient `LMM_Results`, `LMM_Index`, `Index_Results`, `RPE`, `Derive_Results`, `RPE_Scatter_Corr`, et tous les résultats des modèles d'ordre et SD vs LD).

**R² de Nakagawa :** calculé par `computeNakagawaR2` — R²m (variance expliquée par les effets fixes seuls) et R²c (fixes + aléatoires).

---

### `Analyse_RPE.m` — corrélations dérive ↔ ΔRPE, dataset Robin2

**Ce que fait le script :** calcule les corrélations (Pearson et Spearman) entre la dérive de chaque biomarqueur (valeur bin 10 − valeur bin 1, z-scorée) et ΔRPE (RPE_fin − RPE_début) pour les 13 sujets crossover, séparément par condition (Normal et CS60).

**Données d'entrée :** `LMM_Validation_Results.mat` (pour les z-scores de référence `mu_z`, `sd_z` et les `Index_Results`) + `IMU_Features_NewDataset.mat` + `EMG_Features_NewDataset.mat`.

**RPE :** saisies manuelles en dur dans le script (23 sujets, P07 exclu). Format : RPE_début et RPE_fin en CR-100 pour Normal et CS60.

**Ce que produit le script :**
- Corrélations Pearson et Spearman pour chaque variable + index de bloc, par condition.
- Forest plot `Fig_Forest_RPE_Correlations.png` avec IC95% via transformation de Fisher z.
- `RPE_Analysis_Results.mat`.

> **Relation avec `Analyse_Validation.m` :** `Analyse_Validation.m` calcule déjà des corrélations Pearson similaires (section `RPE_Scatter_Corr`) pour 4 séries incluant Normal N=23 et LD. `Analyse_RPE.m` est plus ciblé : uniquement Normal N=13 et CS60 N=13, avec Pearson ET Spearman côte à côte dans un forest plot dédié.

---

### `Analyse_RPE_SDvsLD.m` — RPE et dérive par taille de main

**Ce que fait le script :** analyse les différences de RPE et de dérive biomarqueur entre le groupe **SD (Small Hand, petites mains, N=13)** et le groupe **LD (Large Hand, grandes mains, N=10)**, en condition Normal uniquement.

> **SD ici = taille de main** (Small Hand), pas Short Duration. Voir note d'ambiguïté en tête de ce README.

**Hypothèse testée :** les petites mains subissent une contrainte mécanique relative plus grande sur clavier standard → fatigue plus rapide ou plus marquée, analogue au groupe G1 de Goubault.

**Étapes du script :**

1. **Mann-Whitney SD vs LD** sur RPE_début, RPE_fin, ΔRPE. Produit `Fig_RPE_SD_vs_LD.png` (boxplots).

2. **Dérive par variable** : z-score depuis `LMM_Results`, dérive = z(bin 10) − z(bin 1). Inclut aussi les index de bloc (variante `observed_np`).

3. **Spearman dérive ↔ ΔRPE par groupe** (SD séparément, LD séparément). IC95% via Fisher z. Produit `Fig_Forest_RPE_SD_LD.png`.

4. **Scatter triés par dérive** : `Fig_RPE_SortedScatter_<bloc>.png` — sujets triés par dérive croissante sur l'axe X, ΔRPE sur l'axe Y, avec étiquettes de sujet pour identifier les points influents.

5. **Spearman RPE_début ↔ Dérive par groupe** : vérifie si un RPE_début élevé prédit une dérive biomécanique plus forte. Produit `Fig_Forest_RPE0_Derive_SD_LD.png`.

6. **Corrélation partielle dérive ↔ ΔRPE | RPE_début** : contrôle l'effet du niveau de fatigue initial. Implémentée localement (rang + corrélation partielle de Pearson sur les rangs, test t avec df=n−3).

7. **Comparaison brut vs partiel** : tableau console montrant l'évolution de ρ après contrôle de RPE₀.

**Option `exclude_outliers`** : flag à `true` pour exclure des sujets identifiés visuellement comme atypiques (`Outliers_SD = {'P04','P08','P12','P15'}`, `Outliers_LD = {'P14'}`). Par défaut `false`.

**Sortie :** `RPE_SD_LD_Results.mat`.

---

### `Analyse_Spearman.m` — corrélations dérive ↔ ΔRPE, 5 groupes et 93 variables

**Ce que fait le script :** pipeline Spearman sur **93 variables** réparties en 4 blocs, calculées sur **5 groupes simultanément** : G1 (Short Duration, Goubault), G2 (Long Duration, Goubault), Normal (SD+LD confondus, Robin2), SD (petites mains, Robin2), LD (grandes mains, Robin2).

**Données d'entrée :**
- Jeu 1 (Goubault) : `Workload_Li_TimeNormalised.mat` / `_G2`, `_Goubault` / `_Goubault_G2`, `EMG_Li_TimeNormalised_G1/G2.mat`
- Jeu 2 (Robin2) : `IMU_Features_NewDataset.mat`, `EMG_Features_NewDataset_25vars.mat`, `SegMod_SegAxis_Features_NewDataset.mat`, `Goubault_Features_NewDataset.mat`

**93 variables réparties en 4 blocs :**
- SegMod (16) : Accel/Jerk × 7 segments + 2 totaux
- SegAxis (42) : Accel/Jerk × 7 segments × 3 axes
- Goubault (10) : 10 variables issues de la Table III de Goubault 2023
- EMG (25) : 5 features × 5 canaux (depuis `EMG_Features_NewDataset_25vars.mat`)

**Ce que calcule le script pour chaque variable × groupe :**
- **Dérive X** = X(bin 10) − X(bin 1) par sujet
- **ΔRPE** = RPE(bin 10) − RPE(bin 1) par sujet (jeu 1) ou RPE_fin − RPE_début (jeu 2, valeurs en dur)
- **Corrélation Spearman brute** (dérive X ↔ ΔRPE)
- **Corrélation partielle Spearman** (contrôle RPE₀ = RPE_début), implémentée localement sans toolbox (rang + corrélation partielle de Pearson sur les rangs)

**Champs `raw_*` :** vecteurs bruts par groupe (dX, dRPE, RPE0), stockés dans `Results_corr` pour réutilisation par `Analyse_IM_GAM.m`.

**Champs `bins_*` :** matrices [N_sujets × N_bins] complètes, stockées pour réutilisation par `Analyse_IM_GAM_Bins.m`. Pour le jeu 2, le RPE par bin est **interpolé linéairement** entre RPE_début et RPE_fin (approximation explicitement documentée dans le code — hypothesis d'une montée de fatigue constante, non vérifiée).

**Tableaux PNG produits :** un tableau par bloc × {corrélation simple, corrélation partielle} = 8 PNG. Colonnes : Variable, G1, G2, Normal, SD, LD. Fond coloré et texte en gras si p < 0.05.

**Sortie :** `Derive_DRPE_Results.mat` (array de 93 structs `Results_corr`).

---

### `Analyse_IM_GAM.m` — information mutuelle + GAM

**Ce que fait le script :** filtre les 93 variables de `Derive_DRPE_Results.mat` via l'**information mutuelle de Kraskov (KSG-1)**, puis ajuste un **GAM (spline de lissage)** sur les variables survivantes pour visualiser la forme non-linéaire de la relation dérive↔ΔRPE.

**Entrée :** `Derive_DRPE_Results.mat` (les champs `raw_*` sont obligatoires).

**Étape 1 — Information mutuelle de Kraskov :**
L'estimateur KSG-1 est non paramétrique : pour chaque point, il trouve le k-ième plus proche voisin dans l'espace joint (distance de Chebyshev), puis compte les voisins dans les espaces marginaux. Formule : `I(X;Y) = ψ(k) − ⟨ψ(nx+1) + ψ(ny+1)⟩ + ψ(N)`. k adaptatif : k=3 si N<15, k=5 si N≥15. Tronqué à 0 (l'IM théorique est ≥0, les estimateurs peuvent légèrement dépasser). Test de significativité par **1000 permutations** de Y.

**Correction FDR :** Benjamini-Hochberg sur les 93×5=465 p-values de permutation simultanément.

**Mode de sélection (`SELECTION_MODE`) :** `'fdr'` (conservateur, q<0.05 après BH) ou `'raw'` (exploratoire, p<0.05 brut sans correction — ~23 faux positifs attendus sous H0 sur 465 tests, clairement documenté dans le code).

**Étape 2 — GAM :** pour chaque variable et groupe survivants, `fitrgam` est appelé (nécessite Statistics and Machine Learning Toolbox ≥ R2021a). La forme de la relation ΔRPE ~ s(dérive X) est visualisée sur une grille de 100 points.

**Sorties :** `Tableau_IM_FDR.png` + `Fig_GAM_<source>_<variable>_<groupe>.png` + `IM_GAM_Results.mat` (contient `IM_vals`, `IM_pvals`, `IM_qvals`, `idx_survive`).

**Fonctions locales clés :**
- `kraskov_mi` — estimateur KSG-1, implémenté localement sans dépendance externe
- `psi_local` — fonction digamma approximée (récurrence + approximation d'Euler-Maclaurin)
- `bh_fdr_local` — correction BH locale

---

### `Analyse_IM_GAM_Bins.m` — GAMM par bins avec effet aléatoire sujet

**Ce que fait le script :** remplace le GAM d'`Analyse_IM_GAM.m` par un **GAMM** (Generalized Additive Mixed Model) qui exploite les 10 bins de chaque sujet au lieu d'un seul point dérive. Utilise les champs `bins_*` de `Derive_DRPE_Results.mat`.

**Entrée :** `Derive_DRPE_Results.mat` (champs `bins_*`) + `IM_GAM_Results.mat` (pour `idx_survive`).

**Deux branches selon le groupe :**

**G1 et G2 (RPE mesuré réellement par bin) → GAMM lag-1 :**
Pour chaque bin t, calcule ΔX(t) = X(t)−X(t−1) et ΔRPE(t) = RPE(t)−RPE(t−1). Empile toutes les paires (ΔX, ΔRPE) de tous les sujets, avec un effet aléatoire (1|Sujet) pour respecter la non-indépendance intra-sujet.

Base de spline : B-splines cubiques à 2 nœuds internes (quantiles de ΔX), implémentés localement (algorithme de Cox-de Boor, sans Curve Fitting Toolbox).

Deux modèles comparés :
- **Modèle brut** : `ΔRPE ~ B-splines(ΔX) + (1|Sujet)`, test de Wald model-based
- **Modèle corrigé** : `ΔRPE ~ Time + B-splines(ΔX) + (1|Sujet)` où Time = indice de bin z-scoré (contrôle une dérive temporelle commune à ΔX et ΔRPE). **Test retenu : Wald robuste en cluster (sujet)**, F(q, G−1), qui ne suppose pas le modèle de covariance correctement spécifié (recommandé par Cameron & Miller pour G modéré).

**Normal, SD, LD (RPE interpolé linéairement) → Corrélation drift :**
Le delta lag-1 d'une interpolation linéaire est une constante par sujet → le GAMM lag-1 + (1|Sujet) y est structurellement dégénéré (l'effet aléatoire absorbe exactement cette constante). Retour à un point par sujet : corrélation de Spearman entre dérive totale X et ΔRPE total (même logique qu'`Analyse_IM_GAM.m`).

**FDR final** : BH appliqué séparément par groupe sur les p_global de tous les ajustements de cette exécution.

**Diagnostics inclus :**
- Autocorrélation lag-1 des résidus intra-sujet (paires de bins consécutifs réels, pas à travers NaN)
- AIC des deux modèles (avec et sans contrôle du temps)

**Sorties :** `GAMM_Bins_Results.mat` + `Fig_GAMM_<source>_<variable>_<groupe>.png` (G1/G2) + `Fig_Drift_<source>_<variable>_<groupe>.png` (Normal/SD/LD).

---

### `Analyse_Cluster.m` — exploration de clustering (script exploratoire)

**Ce que fait le script :** tente un clustering K-means (k=2) sur la dérive biomécanique de 9 variables (issues de `Results_IMU`) pour les sujets ayant les deux conditions. Calcule ensuite (section non implémentée) une corrélation de Spearman entre les clusters et ΔRPE.

**Statut :** script exploratoire inachevé. Ne crée aucun fichier de sortie. La section corrélation est commentée (`% ... rho = corr(...) ...`). Les variables de dérive calculées ici (X(bin 10) − X(bin 1)) sont cohérentes avec les autres scripts mais ne sont pas z-scorées.

**Variables analysées :**
- Bloc SegMod : `Accel_Mod_Head`, `Accel_Mod_Hand`
- Bloc Goubault : `MedianFreq_Accel_Y_Hand`, `PeakPower_Accel_Mod_Hand`, `PeakPower_AngVel_X_Head`, `SpectralEntropy_AngVel_Mod_Forearm`
- Bloc EMG : `TFR_MedianFreq_Triceps`, `TFR_SpectralEntropy_Deltoid`, `SampleEntropy_Biceps`

---

## Abréviations communes

| Abréviation | Signification dans ces scripts |
|---|---|
| SD | **Small Hand** (dataset Robin2) — petites mains. Ne pas confondre avec Short Duration (dataset Goubault). |
| LD | **Large Hand** (dataset Robin2) — grandes mains. Ne pas confondre avec Long Duration (dataset Goubault). |
| G1 | **Short Duration** (dataset Goubault) — participants ayant atteint RPE ≥ 7 tôt, arrêtés tôt. |
| G2 | **Long Duration** (dataset Goubault) — participants n'ayant jamais atteint le seuil de fatigue. |
| RPE | Rate of Perceived Exertion — perception de l'effort, échelle CR-100 ici |
| RPE₀ / RPE_début | RPE au début de la session (niveau de fatigue initial) |
| ΔRPE | RPE_fin − RPE_début — progression de la fatigue perçue sur la session |
| Dérive | Valeur au bin 10 − valeur au bin 1 (z-scorée) — proxy de l'évolution d'une variable sur la session |
| LMM | Linear Mixed Model — modèle linéaire mixte (`fitlme` MATLAB) |
| GAMM | Generalized Additive Mixed Model — LMM avec terme non-linéaire (spline) |
| LRT | Likelihood Ratio Test — test du rapport de vraisemblance entre deux modèles imbriqués |
| FDR | False Discovery Rate — taux de faux positifs contrôlé (correction Benjamini-Hochberg) |
| q-value | p-value ajustée après correction FDR |
| R²m | R² marginal de Nakagawa — variance expliquée par les effets fixes seuls |
| R²c | R² conditionnel de Nakagawa — variance expliquée par les effets fixes + aléatoires |
| IM | Information mutuelle — mesure de dépendance non-linéaire entre deux variables |
| KSG-1 | Estimateur de Kraskov-Stögbauer-Grassberger, variante 1 — estimateur non-paramétrique de l'IM |
| GAM | Generalized Additive Model — modèle additif avec termes spline (`fitrgam`) |
| B-splines | Base de splines cubiques — approximation non-linéaire sans itération, implémentée localement |
| Lag-1 | Différence entre deux bins consécutifs : ΔX(t) = X(t) − X(t−1) |
| Order | Ordre de passage — quelle condition (Normal ou CS60) a été jouée en premier par un sujet crossover |
| crossover | Sujets ayant joué les deux conditions Normal et CS60 (N=13 dans le dataset Robin2) |
| β_Temps | Coefficient de l'effet Temps dans le LMM — pente de la variable en fonction du pourcentage de session |
| β_std | Coefficient β standardisé (sur variable z-scorée) |
| ml / REML | Maximum Likelihood / Restricted ML — méthode d'ajustement du LMM (ML obligatoire pour LRT sur effets fixes) |
| Nakagawa R² | Méthode de calcul du R² pour les LMM (variance des effets fixes / variance totale) |
| Normal N=13 | Sous-ensemble crossover en condition normale (sujets ayant aussi joué CS60) |
| Normal N=23 | Tous les sujets en condition normale, y compris ceux sans CS60 |
| Normal LD | Sujets grandes mains en condition normale uniquement (N=10) |

---

## Fichiers d'entrée / sortie récapitulatifs

| Fichier | Produit par | Consommé par |
|---|---|---|
| `IMU_Features_NewDataset.mat` | `Analyse_IMU.m` | `Analyse_Validation.m`, `Analyse_RPE.m`, `Analyse_RPE_SDvsLD.m`, `Analyse_Spearman.m` |
| `EMG_Features_NewDataset.mat` | `Analyse_EMG.m` | `Analyse_Validation.m`, `Analyse_RPE.m`, `Analyse_RPE_SDvsLD.m` |
| `EMG_Features_NewDataset_25vars.mat` | `Compute_EMG_NewDataSet.m` | `Analyse_Spearman.m` |
| `SegMod_SegAxis_Features_NewDataset.mat` | `Compute_SegMod_SegAxis_NewDataSet.m` | `Analyse_Spearman.m` |
| `Goubault_Features_NewDataset.mat` | `Compute_Goubault_NewDataSet.m` | `Analyse_Spearman.m` |
| `LMM_Validation_Results.mat` | `Analyse_Validation.m` | `Analyse_RPE.m`, `Analyse_RPE_SDvsLD.m` |
| `RPE_Analysis_Results.mat` | `Analyse_RPE.m` | — |
| `RPE_SD_LD_Results.mat` | `Analyse_RPE_SDvsLD.m` | — |
| `Derive_DRPE_Results.mat` | `Analyse_Spearman.m` | `Analyse_IM_GAM.m`, `Analyse_IM_GAM_Bins.m` |
| `IM_GAM_Results.mat` | `Analyse_IM_GAM.m` | `Analyse_IM_GAM_Bins.m` (liste `idx_survive`) |
| `GAMM_Bins_Results.mat` | `Analyse_IM_GAM_Bins.m` | — |
