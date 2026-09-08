# README — Scripts LASSO + LMM, dataset Goubault (G1 vs G2)

Ce groupe contient les scripts de **sélection de variables et de modélisation** sur le dataset Goubault (jeu 1, 50 pianistes, tâche Li/Chord), plus deux scripts utilitaires de visualisation associés. L'objectif est d'identifier, parmi les variables cinématiques et EMG, celles qui prédisent le mieux l'évolution du RPE au niveau individuel.

---

## Vue d'ensemble et évolution des scripts

```
Lasso_LMM.m             ← Version 1 (ancienne) — LASSO classique, CV 5-fold
        ↓ remplacé par
Indiv_Lasso_LMM.m       ← Version 2 — Group LASSO FISTA + Mundlak + HAC + 100 graines
        ↓ étendu par
Indiv_Lasso_LMM_V2.m    ← Version 3 — idem + boucle sur 4 variantes de fenêtre temporelle

Sorties :
Stability_Archive.mat   ← Produit par Indiv_Lasso_LMM.m
        ↓ rechargé par
Forest_Plot.m           ← Script autonome de forest plots (sans relancer les 100 graines)

Graphiques_LASSO_LMM.m  ← Mise en page manuelle d'un tableau de coefficients (publication)
```

---

## Dataset commun

| Paramètre | Valeur |
|---|---|
| Dataset | Goubault (50 pianistes experts, S101–S150) |
| G1 (Short Duration) | N≈26, participants ayant atteint RPE ≥ 7 tôt |
| G2 (Long Duration) | N≈23, participants n'ayant jamais atteint le seuil |
| Blocs de variables | A=BySegMod (Features3), B=BySegAxis, C=Goubault, D=EMG |
| Bins | 10 intervalles temporels normalisés |

---

## Scripts

---

### `Lasso_LMM.m` — version initiale (LASSO classique)

**Statut :** ancienne version, remplacée par `Indiv_Lasso_LMM.m`.

**Ce que fait le script :** pipeline LASSO + LMM en 4 étapes, sur les mêmes 4 blocs de variables que les versions suivantes.

**Différences avec `Indiv_Lasso_LMM.m` :**

| | `Lasso_LMM.m` | `Indiv_Lasso_LMM.m` |
|---|---|---|
| Algorithme LASSO | LASSO classique (`lasso`, `fitglmnet`) | **Group LASSO FISTA** (implémenté localement) |
| Validation croisée | CV 5-fold standard | **GroupKFold par sujet** (évite la fuite intra-sujet) |
| Décomposition | Aucune | **Mundlak** (within + between) |
| Erreurs standard | Classiques | **HAC Newey-West** (robustes à l'autocorrélation) |
| Stabilité | 100 graines, comptage binaire | 100 graines, stockage des β_std par run |
| Sortie stabilité | Tableau PNG fréquences | `Stability_Archive.mat` + forest plots |
| Contrôle autocorrélation résiduelle | RPE_lag1 comme covariable | Décomposition Mundlak (pas de lag-1 explicite) |

**Étapes :**
1. Construction des tableaux longs par bloc (`buildBlock` : centrage individuel, format participant × bin)
2. LASSO CV 5-fold par bloc et par groupe, sélection `lambda_1SE`, top 8 variables par bloc
3. LMM avec covariable `RPE_lag1` + calcul R² marginal par variable (contribution unique = R²_full − R²_sans_variable)
4. Figures coefficients β par bloc (barres G1/G2), tableau β + R²m, R² global par bloc

**Sortie :** `LASSO_LMM_Results.mat`.

---

### `Indiv_Lasso_LMM.m` — version principale (Group LASSO FISTA)

**Ce que fait le script :** pipeline complet Group LASSO + LMM avec correction des biais statistiques identifiés dans `Lasso_LMM.m`. C'est la version de référence, utilisée pour les résultats publiés.

**Pipeline en 7 étapes :**

**Étape 1 — Chargement et correction EMG :**
Les valeurs `Inf` de `SampleEntropy` dans `Mat_EMG` sont remplacées par NaN, puis imputées par la moyenne intra-individuelle. Les sujets avec 100 % de valeurs invalides sur un canal sont exclus uniquement du bloc EMG (pas des autres blocs).

**Étape 2 — Construction des tableaux longs :**
`buildBlock` centre chaque variable par participant (within) et calcule sa moyenne (between) — décomposition de Mundlak en amont de l'algorithme. Format long : N_sujets × N_bins observations empilées.

**Étape 3b — Group LASSO FISTA + GroupKFold :**
- **Group LASSO** : les variables within et between d'une même feature cinématique forment un groupe — elles sont sélectionnées ensemble ou pas du tout. Pénalité L1 sur les normes de groupe.
- **FISTA** (Fast Iterative Shrinkage-Thresholding Algorithm) : algorithme de descente proximale accélérée, implémenté localement (`fitGroupLassoFISTA`). Convergence sur `||β(t) − β(t−1)||₂ < TOL_FISTA = 1e-6`.
- **GroupKFold par sujet** (`groupKFoldBySubject`) : K=5 folds définis en entier sur les sujets (pas sur les observations), évitant que les bins d'un même participant se retrouvent en train et en test.
- **Mundlak** : matrice de design `[X_w1, X_b1, X_w2, X_b2, ...]` avec groupes appariés within/between, permettant de séparer l'effet intra-individuel (évolution dans le temps) de l'effet inter-individuel (différences de niveau).
- **HAC Newey-West** (`hac`, bande passante = 2) : erreurs standard robustes à l'autocorrélation résiduelle et à l'hétéroscédasticité.

**Étape 4 — LMM sur variables sélectionnées :**
`Variable ~ variables_sélectionnées + (1|Participant)`, sans terme lag-1 (la décomposition Mundlak gère l'autocorrélation intra-sujet). Correction Cochrane-Orcutt itérative (Prais-Winsten) si l'autocorrélation résiduelle persiste. R² de Nakagawa (R²m marginal, R²c conditionnel). Vérification homoscédasticité via corrélation de Spearman résidus² ~ temps (test de Breusch-Pagan simplifié).

**Étape 5 — Figures LMM :** trajectoires des coefficients β standardisés, chemin de régularisation Group LASSO, comparaison G1/G2.

**Étape 6 — Boucle de stabilité (100 graines) :**
Pour chaque graine 1–100, le Group LASSO est relancé. Pour chaque variable sélectionnée dans un run, le β_std (coefficient standardisé) est stocké. La structure `stab_G1/G2.(bloc).(variable)` accumule les champs : `n_sel` (nombre de sélections), `beta_std` (vecteur des β_std sur les runs où la variable est sélectionnée), `pval`, `dR2`, `R2m_z`, `R2c_z`.

**Étape 7 — Forest plots et sauvegarde :**
Forest plots intégrés (variables avec fréquence ≥ THRESH_STAB = 50 %). Sauvegarde dans `Stability_Archive.mat`.

**Fonctions locales clés :**

| Fonction | Rôle |
|---|---|
| `fitGroupLassoFISTA` | Group LASSO via FISTA, implémenté localement |
| `groupSoftThreshold` | Opérateur proximal du Group LASSO (seuillage par groupe) |
| `computeLambdaMax` | Calcule le λ maximum (toutes variables à zéro) |
| `groupKFoldBySubject` | K-fold par sujet (folds entiers) |
| `buildBlock` | Construction tableau long + décomposition Mundlak within/between |
| `makeLagBySubject` | Lag-1 intra-sujet du RPE (covariable, utilisée dans Lasso_LMM.m mais pas ici) |
| `computeNakagawaR2` | R² marginal et conditionnel (Nakagawa) |
| `fisherTtest` | Corrélation de Spearman individuelle + transformation de Fisher + t-test |
| `buildAR1blockDiag` | Matrice de covariance AR(1) bloc-diagonale par sujet |
| `ar1NegLogLik` | Log-vraisemblance négative AR(1) pour Cochrane-Orcutt |
| `shadedErrorBar_simple` | Courbe avec bande d'erreur (utilitaire de tracé) |
| `normLabel` | Normalise les noms de variables pour l'affichage |
| `p2star` | Convertit une p-value en étoiles (* ** ***) |

**Sorties :** `Stability_Archive.mat` + `ForestPlot_<bloc>_<groupe>.png` + figures LMM.

---

### `Indiv_Lasso_LMM_V2.m` — boucle sur 4 variantes de fenêtre temporelle

**Ce que fait le script :** identique à `Indiv_Lasso_LMM.m` dans sa logique, avec une **boucle externe sur 4 variantes de fenêtre de différence temporelle** pour le calcul de ΔRPE et ΔX.

**Différence principale avec `Indiv_Lasso_LMM.m` :**

| | `Indiv_Lasso_LMM.m` | `Indiv_Lasso_LMM_V2.m` |
|---|---|---|
| Fenêtre de différence | Une seule (lag-1 implicite dans buildBlock) | **4 variantes** : lag1, lag2, lag3, global |
| Graines stabilité | 100 | **10** (réduit pour limiter le temps de calcul × 4 variantes) |
| Sorties | `Stability_Archive.mat` | `lag1_Stability_Archive.mat`, `lag2_...`, etc. |
| Forest plots | `ForestPlot_<bloc>_<groupe>.png` | `lag1_ForestPlot_<bloc>_<groupe>.png`, etc. |

**4 variantes de fenêtre (`lag_modes`) :**
- `lag1` : ΔRPE(t) = RPE(t) − RPE(t−1), ΔX(t) = X(t) − X(t−1) — dynamique fine, 9 observations par sujet, bruit élevé
- `lag2` : fenêtre de 2 bins — lissage modéré, 8 observations par sujet
- `lag3` : fenêtre de 3 bins — lissage plus fort, 7 observations par sujet
- `global` : ΔX = X(bin10) − X(bin1), ΔRPE = RPE_fin − RPE_début — un seul point par sujet (identique à la dérive calculée dans `Analyse_Spearman.m`)

**Objectif :** évaluer la sensibilité des sélections LASSO au choix de la fenêtre de différence. Si les mêmes variables ressortent sur lag1 et global, la sélection est robuste à l'échelle temporelle.

---

### `Forest_Plot.m` — script autonome de forest plots

**Ce que fait le script :** recharge `Stability_Archive.mat` et régénère les forest plots de stabilité LASSO, sans avoir à relancer les 100 graines. Utile pour ajuster les paramètres visuels sans recalcul.

**Entrée :** `Stability_Archive.mat` (produit par `Indiv_Lasso_LMM.m`).

**Ce qu'affiche chaque forest plot (un par bloc × groupe) :**
- Ligne bleue fine : étendue [min, max] de β_std sur les 100 runs où la variable est sélectionnée
- Ligne rouge épaisse : IC 95% de β_std (intervalle de confiance sur la moyenne des β_std)
- Point blanc cerclé rouge : β_std moyen
- Chiffre annoté : fréquence de sélection (%)

Seules les variables avec fréquence ≥ `THRESH_STAB` (défaut = 50 %) sont affichées.

**Sorties :** `ForestPlot_<bloc>_<groupe>.png` pour chaque combinaison bloc × {G1, G2}.

---

### `Graphiques_LASSO_LMM.m` — tableau de coefficients pour publication

**Ce que fait le script :** génère un tableau formaté style publication (style booktabs) avec les coefficients de régression β_raw, β_std et les niveaux de significativité pour les 10 variables retenues, organisées par bloc.

**Les données sont saisies en dur** dans la variable `data` — ce script ne charge aucun fichier `.mat`. Il sert uniquement à produire une figure PNG propre pour un article ou un rapport.

**Sortie :** `tableau_coefficients_matlab.png`.

---

## Abréviations

| Abréviation | Signification dans ces scripts |
|---|---|
| G1 | Short Duration — participants ayant atteint RPE ≥ 7 tôt (N≈26) |
| G2 | Long Duration — participants n'ayant jamais atteint le seuil (N≈23) |
| LASSO | Least Absolute Shrinkage and Selection Operator — régression pénalisée L1 |
| Group LASSO | Variante du LASSO où les variables sont regroupées et sélectionnées ensemble (pénalité L2 intra-groupe) |
| FISTA | Fast Iterative Shrinkage-Thresholding Algorithm — algorithme de descente proximale accélérée pour le Group LASSO |
| GroupKFold | Validation croisée K-fold où les K folds sont définis sur les sujets entiers (pas sur les observations) |
| Mundlak | Décomposition d'une variable en composante within (intra-sujet, centrée) et between (inter-sujet, moyenne) pour séparer les effets longitudinaux et transversaux |
| Within (Xw) | Composante intra-individuelle : X(i,t) − mean_t(X(i,:)) — capture l'évolution dans le temps |
| Between (Xb) | Composante inter-individuelle : mean_t(X(i,:)) — niveau moyen du sujet |
| HAC | Heteroskedasticity and Autocorrelation Consistent — erreurs standard robustes |
| HAC Newey-West | Variante HAC avec noyau de Bartlett-Priestley, bande passante = 2 bins |
| Cochrane-Orcutt | Correction itérative de l'autocorrélation AR(1) des résidus (aussi dite transformation Prais-Winsten) |
| AR(1) | Autocorrélation d'ordre 1 — corrélation entre le résidu au temps t et le résidu au temps t−1 |
| λ_1SE | Lambda sélectionné par la règle « 1 erreur standard » en validation croisée — plus parcimonieux que lambda_min |
| β_std | Coefficient β standardisé (sur variable z-scorée), comparable entre variables d'unités différentes |
| β_raw | Coefficient β non standardisé, dans l'unité originale de la variable |
| R²m | R² marginal de Nakagawa — variance expliquée par les effets fixes seuls |
| R²c | R² conditionnel de Nakagawa — variance expliquée par les effets fixes + aléatoires |
| R²m marginal par variable | ΔR² = R²_full − R²_sans_variable — contribution unique de chaque variable |
| Breusch-Pagan simplifié | Test d'homoscédasticité : corrélation de Spearman entre résidus² et temps |
| THRESH_STAB | Seuil de fréquence de sélection (%) au-delà duquel une variable est considérée stable (défaut : 50 %) |
| n_sel | Nombre de runs (sur 100) où une variable est sélectionnée |
| Bloc A | BySegMod — module Accel/Jerk par segment, Features3 (n≈25/23) |
| Bloc B | BySegAxis — Accel/Jerk par segment et par axe (n≈14/15) |
| Bloc C | Goubault — 10 variables de Goubault 2023 Table III |
| Bloc D | EMG — toutes variables EMG sauf Amplitude |
| lag1/lag2/lag3 | Fenêtre de différence temporelle : ΔX(t) = X(t) − X(t−k), k=1,2,3 |
| global | Différence totale sur la session : X(bin10) − X(bin1) |

---

## Fichiers d'entrée / sortie

| Fichier | Produit par | Consommé par |
|---|---|---|
| `LASSO_LMM_Results.mat` | `Lasso_LMM.m` | — |
| `Stability_Archive.mat` | `Indiv_Lasso_LMM.m` | `Forest_Plot.m` |
| `lag1_Stability_Archive.mat` … `global_Stability_Archive.mat` | `Indiv_Lasso_LMM_V2.m` | — |
| `ForestPlot_<bloc>_<groupe>.png` | `Indiv_Lasso_LMM.m` ou `Forest_Plot.m` | — |
| `lag*_ForestPlot_*.png` | `Indiv_Lasso_LMM_V2.m` | — |
| `tableau_coefficients_matlab.png` | `Graphiques_LASSO_LMM.m` | — |
