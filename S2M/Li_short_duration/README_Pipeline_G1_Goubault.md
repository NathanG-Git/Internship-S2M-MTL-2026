# README — Pipeline G1 : calcul, visualisation et corrélations (dataset Goubault)

Ce groupe contient les scripts de **calcul depuis les données brutes, visualisation et corrélations Spearman** pour le groupe G1 (Short Duration) du dataset Goubault. Ils constituent le pipeline originel G1, dont les scripts G2 (`Visualise_Worload_Li_TimeNorm_LD.m` etc.) sont les équivalents symétriques.

---

## Vue d'ensemble et chaîne de données

```
Données brutes XSens (S10X_Li_Accel.mat, Fs=60 Hz)
        │
        ├── Visualize_Workload_Li_TimeNormalised.m
        │       → Workload_Li_TimeNormalised.mat  (d1 dans les scripts d'analyse)
        │
        ├── Visualize_Workload_Li_TimeNormalised_Goubault.m
        │       → Workload_Li_TimeNormalised_Goubault.mat  (dG1)
        │
        └── Visualize_EMG_Li_TimeNormalised_G1.m
                → EMG_Li_TimeNormalised_G1.mat  (eG1)

Données brutes JointData (S10X_Li_JointData.mat)
        └── AccelSign_JointData_Li.m
                → AccelSign_Li_G1.mat

Workload_Li_G1.mat (format cycles, produit par un script antérieur)
        ├── Visualize_Workload_Li.m             → figures PNG (sans .mat)
        ├── SpearmanCorr_Workload_Li.m          → Correlation_Workload_RPE.mat + .xlsx
        └── SpearmanCorr_Workload_Li_FDR.m      → Correlation_Workload_RPE_FDR.mat + .xlsx

Tables_XSENS_LMM/ (datpart1_*.mat, format par participant)
        ├── SpearmanCorr_GoubaultFeatures_Li.m  → Correlation_Goubault_Features_RPE.mat + .xlsx
        └── SpearmanCorr_GoubaultFeatures_Li_FDR.m → Correlation_Goubault_Features_RPE_FDR.mat + .xlsx
```

---

## Scripts de calcul (producteurs de `.mat`)

---

### `Visualize_Workload_Li_TimeNormalised.m` — variables cinématiques G1

**Équivalent G1 de `Visualise_Worload_Li_TimeNorm_LD.m`.**

**Ce que fait le script :** calcule les variables cinématiques G1 depuis les données brutes XSens frame par frame, avec une différence clé par rapport au script G2 : **coupure RPE**. Chaque participant G1 est tronqué dès que son RPE discret atteint 7 deux fois consécutives (critère d'arrêt). Les intervalles au-delà sont mis à NaN.

**Pipeline :** filtrage Butterworth 0.5–29 Hz → `|accel|` et `jerk = diff(accel)/dt` → agrégation par segment/axe/total (sommes) → 10 bins proportionnels → interpolation RPE → normalisation par baseline individuelle.

**Sections BySegMod_26 :** même logique que dans le script G2 — chargement depuis `Features3_XSENS` pour le sous-ensemble de participants ayant ces fichiers pré-calculés (`valid_26`).

**Sortie `.mat` :** `Workload_Li_TimeNormalised.mat` — matrices `Mat_Total_Accel`, `Mat_Total_Jerk`, `Mat_Axis_Accel`, `Mat_Axis_Jerk`, `Mat_Seg_Accel`, `Mat_Seg_Jerk`, `Mat_SegAxis_Accel`, `Mat_SegAxis_Jerk`, `Mat_RPE`, `valid_subj`, `valid_26`, `Mat_RPE_26`.

> Chargé sous `d1` dans les scripts d'analyse. La présence de `valid_26` (Features3) distingue ce fichier de `Workload_Li_TimeNormalised_G2.mat`.

---

### `Visualize_Workload_Li_TimeNormalised_Goubault.m` — 10 variables Goubault G1

**Équivalent G1 de `Visualise_Workload_Li_TimeNorm_LD_Goubault.m`.**

**Ce que fait le script :** calcule les 10 variables Goubault 2023 depuis `Features3_XSENS_S10X_Li.mat` pour G1. Pipeline et variables identiques au script G2 (même définition des colonnes, même normalisation par la moyenne des bins 1–3, même application de `abs()` sur les variables Mean). Pas de coupure RPE.

**Sortie `.mat` :** `Workload_Li_TimeNormalised_Goubault.mat` — chargé sous `dG1`.

**Figure produite :** `Fig_TimeNorm_Goubault.png` — grille 2 × 5, double axe Y, droite de tendance rouge.

---

### `Visualize_EMG_Li_TimeNormalised_G1.m` — variables EMG G1

**Équivalent G1 de `Visualise_EMG_Li_TimeNormalised_G2.m`.**

**Ce que fait le script :** calcule 6 variables EMG × 5 canaux pour G1 depuis des fichiers pré-calculés par variable. Pipeline identique au script G2 : normalisation Amplitude en % MVC (indices 43–47 du vecteur `EMG_max`), normalisation Activity par médiane bin1 sur les figures, nettoyage SampleEntropy (Inf → NaN, ±3 SD sur figures Median). Même gestion de l'exception S122.

**Sortie `.mat` :** `EMG_Li_TimeNormalised_G1.mat` — chargé sous `eG1`.

---

### `AccelSign_JointData_Li.m` — accélération propulsive vs freinatrice, G1

**Script spécifique G1 sans équivalent G2.** Idée de l'encadrant.

**Ce que fait le script :** pour chaque articulation du bras droit (rGH, rEL, rWR) et chaque cycle musical, sépare l'accélération articulaire selon son signe relatif à la vitesse :
- **Propulsive** : `sign(accel) == sign(velocity)` — accélère le mouvement
- **Freinatrice** : `sign(accel) != sign(velocity)` — freine le mouvement

**Pour chaque articulation × axe × cycle :**
- Mean |accel| propulsive
- Mean |accel| freinatrice
- Ratio propulsif / (propulsif + freinateur) ∈ [0, 1]

**Source des données :** `S10X_Li_JointData.mat` au format v7.3, lu via `h5read` (champ `JointData.<joint>.jointAngle` et `.jointVelocity`, dimensions 3 × N_frames). L'accélération articulaire est calculée comme `diff(velocity) * Fs` (dérivée numérique de la vitesse angulaire).

**Corrélations Spearman** intégrées : pour chaque feature (Mean propulsive, Mean freinatrice, Ratio), corrélation vs RPE par participant, puis test de Wilcoxon sur la médiane de groupe + correction FDR Benjamini-Hochberg.

**Sortie `.mat` :** `AccelSign_Li_G1.mat` (structure `AccelSign`).

---

## Scripts de visualisation (format cycles, sans production de `.mat` new)

---

### `Visualize_Workload_Li.m` — visualisation format cycles, G1

**Ce que fait le script :** visualise les données de `Workload_Li_G1.mat` (format par cycles musicaux, produit par un script antérieur non inclus ici). Courbes individuelles transparentes + moyenne ± SD de groupe en premier plan.

**Format des données d'entrée :** `Workload_Li_G1.mat` contient une structure `Workload` (un élément par participant) avec des champs `Mean` et `Median`, chacun contenant les sous-champs `Acceleration` et `Jerk`, avec les sous-sous-champs `Total`, `ByAxis`, `BySegMod`, `BySegAxis`. L'axe X est le **RPE** (pas le temps), de 1 à 10.

> **Différence fondamentale avec `Visualize_Workload_Li_TimeNormalised.m` :** ici les données sont tracées en fonction du RPE (abscisse = niveau de fatigue perçue), pas du temps normalisé.

**7 figures produites :**
- `Fig1_Total.png` — grille 2×2 (Mean/Median × Accel/Jerk)
- `Fig2_ByAxis_Acceleration.png`, `Fig3_ByAxis_Jerk.png` — grille 2×3 (Mean/Median × X/Y/Z)
- `Fig4_BySegment_Acceleration.png`, `Fig5_BySegment_Jerk.png` — grille 2×7 (segments)
- `Fig6_BySegAxis_Acceleration.png`, `Fig7_BySegAxis_Jerk.png` — grille 2×7, 3 axes superposés par couleur (X=bleu, Y=rouge, Z=vert)

---

## Scripts de corrélations Spearman

Ces 4 scripts forment deux paires (avec et sans correction FDR) sur deux sources de données différentes.

---

### `SpearmanCorr_Workload_Li.m` et `SpearmanCorr_Workload_Li_FDR.m`

**Source :** `Workload_Li_G1.mat` (format cycles).

**Ce que font ces scripts :** pour chaque variable (Total, ByAxis, BySegMod, BySegAxis) × signal (Accel, Jerk) × statistique (Mean, Median), calculent la corrélation de Spearman entre le score et le RPE **par participant**, puis testent si la médiane de groupe est significativement différente de 0 (test de Wilcoxon signed-rank).

**32 variables** au total (1 Total + 3 ByAxis + 7 BySegMod + 21 BySegAxis = 32, × 2 signaux = **64 p-values par statistique**, soit **128 p-values globales** sur les 4 combinaisons Score × Signal).

| | `SpearmanCorr_Workload_Li.m` | `SpearmanCorr_Workload_Li_FDR.m` |
|---|---|---|
| Correction multiple | Non — p-values brutes | **FDR Benjamini-Hochberg** sur les 128 p-values |
| Heatmap | `*` = p < 0.05 brut | `†` = survit au FDR, `*` = p brut seulement |
| Sorties | `Correlation_Workload_RPE.mat` + `.xlsx` | `Correlation_Workload_RPE_FDR.mat` + `.xlsx` |

**Implémentation BH locale** (`bh_fdr`) : sans dépendance à aucun toolbox. Algorithme step-up identique à `fdr_bh(...,'pdep')` et `scipy.stats.false_discovery_control`.

---

### `SpearmanCorr_GoubaultFeatures_Li.m` et `SpearmanCorr_GoubaultFeatures_Li_FDR.m`

**Source :** `Tables_XSENS_LMM/datpart1_<Feature>_<Signal>.mat` — format long par participant (N_obs × 7 colonnes : ID, RPE interpolé, timestamp, 4 segments Shoulder/Arm/Forearm/Hand).

**Ce que font ces scripts :** idem pour les **features Goubault** : 8 features × 4 signaux × 4 segments = **128 variables**. La corrélation est calculée par participant (intra-individuelle), puis la médiane de groupe est testée par Wilcoxon.

**8 features** : MedianFreq, PeakPower, PeakPower_Freq, Peak, PowerLF, PowerHF, PowerTot, SpectralEntropy.
**4 signaux** : Acceleration, Angular_Velocity, Module_Acceleration, Module_Angular_Velocity.
**4 segments** (colonnes 4–7) : Shoulder, Arm, Forearm, Hand.

| | `SpearmanCorr_GoubaultFeatures_Li.m` | `SpearmanCorr_GoubaultFeatures_Li_FDR.m` |
|---|---|---|
| Correction multiple | Non | **FDR BH** sur les 128 p-values |
| Heatmap | Grille Feature×Signal × Segment, `*` brut | `†` FDR, `*` brut seulement, grille avec lignes |
| Sorties | `.mat` + `.xlsx` + heatmap | `.mat` + `.xlsx` + heatmap FDR |

> **Différence de source par rapport à `Visualize_Workload_Li_TimeNormalised_Goubault.m` :** les scripts Goubault visualisent des features issues de `Features3_XSENS` en fonction du temps normalisé. Les scripts `SpearmanCorr_Goubault*` calculent des corrélations sur des tables `datpart1_*` en fonction du RPE (format différent, axe différent, 4 segments vs 7).

---

## Abréviations communes

| Abréviation | Signification dans ces scripts |
|---|---|
| G1 | Short Duration — participants ayant atteint RPE ≥ 7 deux fois de suite (N=26) |
| rGH | Articulation gléno-humérale droite (Right Glenohumeral) = épaule |
| rEL | Coude droit (Right Elbow) |
| rWR | Poignet droit (Right Wrist) |
| Propulsive | Accélération de même signe que la vitesse — accélère le mouvement articulaire |
| Freinatrice | Accélération de signe opposé à la vitesse — freine le mouvement articulaire |
| Ratio propulsif | Mean propulsive / (Mean propulsive + Mean freinatrice) ∈ [0, 1] |
| JointData | Données angulaires articulaires (angle, vitesse angulaire) issues de XSens |
| h5read | Fonction MATLAB pour lire les fichiers .mat au format v7.3 (HDF5) |
| datpart1 | Tables longues par participant (format LMM), stockées dans `Tables_XSENS_LMM/` |
| FDR | False Discovery Rate — taux de faux positifs contrôlé |
| BH | Benjamini-Hochberg — méthode de correction FDR step-up |
| q-value / p_adj | P-value ajustée après correction FDR |
| † (dagger) | Symbole utilisé dans les heatmaps pour indiquer la survie à la correction FDR |
| Wilcoxon signed-rank | Test non-paramétrique sur la médiane d'un seul groupe (H0 : médiane = 0) |
| signrank | Implémentation MATLAB du test de Wilcoxon signed-rank |
| PowerLF | Puissance dans les basses fréquences (0.1–4 Hz) du spectre CWT |
| PowerHF | Puissance dans les hautes fréquences (6–12 Hz) du spectre CWT |
| PowerTot | Puissance totale (0.1–12 Hz) du spectre CWT |
| Cycles musicaux | Répétitions d'un motif musical, utilisées comme unité temporelle (axe X dans `Workload_Li_G1.mat`) |
| Bins temporels | Intervalles de durée égale (1/10e de la session), axe X dans les scripts TimeNormalised |

---

## Fichiers d'entrée / sortie

| Fichier | Produit par | Consommé par (alias) |
|---|---|---|
| `Workload_Li_TimeNormalised.mat` | `Visualize_Workload_Li_TimeNormalised.m` | Scripts d'analyse sous `d1` |
| `Workload_Li_TimeNormalised_Goubault.mat` | `Visualize_Workload_Li_TimeNormalised_Goubault.m` | Scripts d'analyse sous `dG1` |
| `EMG_Li_TimeNormalised_G1.mat` | `Visualize_EMG_Li_TimeNormalised_G1.m` | Scripts d'analyse sous `eG1` |
| `AccelSign_Li_G1.mat` | `AccelSign_JointData_Li.m` | — |
| `Correlation_Workload_RPE.mat` + `.xlsx` | `SpearmanCorr_Workload_Li.m` | — |
| `Correlation_Workload_RPE_FDR.mat` + `.xlsx` | `SpearmanCorr_Workload_Li_FDR.m` | — |
| `Correlation_Goubault_Features_RPE.mat` + `.xlsx` | `SpearmanCorr_GoubaultFeatures_Li.m` | — |
| `Correlation_Goubault_Features_RPE_FDR.mat` + `.xlsx` | `SpearmanCorr_GoubaultFeatures_Li_FDR.m` | — |
