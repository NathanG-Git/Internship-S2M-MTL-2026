# README — Scripts de calcul et visualisation, dataset Goubault G2

Ce groupe contient les trois scripts qui **calculent les variables depuis les données brutes et produisent les fichiers `.mat` de référence** pour le groupe G2 (Long Duration) du dataset Goubault. Ils génèrent également des figures de visualisation temporelle associées.

> Ces scripts sont les **équivalents G2** des scripts G1 (`Workload_Li_Cycles.m`, `EMG_Li_Cycles.m`, etc.). Ils appliquent les mêmes pipelines de traitement sur les participants G2 et produisent des fichiers `.mat` qui sont ensuite chargés par les scripts d'analyse (`Spearman_Corr_G1_G2.m`, `Lasso_LMM.m`, etc.).

---

## Contexte commun

| Paramètre | Valeur |
|---|---|
| Groupe | G2 = Long Duration — participants ayant joué la session complète sans atteindre RPE ≥ 7 |
| N participants | 23–24 selon le script (voir ci-dessous) |
| Tâche | Li (Liszt = tâche Chord du dataset Goubault) |
| Bins | 10 intervalles temporels normalisés |
| RPE | Mesures toutes les 30 s (échelle CR-10 0–10), interpolées linéairement sur la grille temporelle |
| Exception S122 | Dernière mesure RPE supprimée (artefact connu) |

> **G2 = Long Duration** : participants n'ayant jamais atteint RPE ≥ 7 deux fois de suite. Ne pas confondre avec LD (Large Hand) du dataset Robin2.

---

## Scripts

---

### `Visualise_Worload_Li_TimeNorm_LD.m` — variables cinématiques G2

**Ce que fait le script :** calcule depuis les **données brutes XSens frame par frame** (`S10X_Li_Accel.mat`) les variables cinématiques pour G2, les normalise temporellement en 10 bins, puis sauvegarde et visualise.

**Pipeline :**
1. Chargement du signal d'accélération brut (N_frames × 69 colonnes — 7 segments × 3 axes + colonnes additionnelles)
2. Filtre Butterworth passe-bande 0.5–29 Hz (ordre 2, Fs=60 Hz) — paramètres conformes à Goubault
3. Calcul de `|accel|` et `jerk = diff(accel) / dt` frame par frame
4. Agrégation par segment/axe/total (sommes, cohérent avec le script G1 de référence `Workload_Li_Cycles.m`)
5. Découpage en 10 bins proportionnels à la durée de session
6. Interpolation RPE sur les 10 bins
7. Normalisation par la moyenne individuelle (ratio / baseline)

**Convention d'agrégation :**
- `ByAxis` : somme des `|accel|` sur les 7 segments pour chaque axe (X, Y, Z)
- `BySegMod` : norme euclidienne `√(X²+Y²+Z²)` pour chaque segment
- `Total` : somme sur les 21 colonnes (= somme de ByAxis)

**Deux sources de données pour BySegMod :**
Le script contient une section spécifique « 26 participants » qui recharge les `Features3_XSENS_S10X_Li.mat` pour les participants G2 ayant ces fichiers pré-calculés. Cela permet d'avoir un sous-ensemble plus cohérent avec les données G1 (qui utilisaient aussi Features3). Les matrices correspondantes sont `Mat_Seg_Accel_26`, `Mat_Seg_Jerk_26`, `valid_26`.

**Participants G2 :** `[3,8,10,11,12,13,14,16,17,18,19,24,25,26,27,33,34,36,40,41,43,47,48,50]` — N=24 (incluant S117 = iP=17, absent dans `Visualise_Workload_Li_TimeNorm_LD_Goubault.m` et `Visualise_EMG_Li_TimeNormalised_G2.m` qui ont N=23).

**Figures produites :**
- `Fig_TimeNorm_Total_G2.png` — Total Accel et Jerk + RPE
- `Fig_TimeNorm_ByAxis_G2_Acceleration.png` et `_Jerk.png` — 3 axes × 2 signaux
- Figures BySegMod et BySegAxis (même logique)

**Sortie `.mat` :** `Workload_Li_TimeNormalised_G2.mat` — contient `Mat_Total_Accel`, `Mat_Total_Jerk`, `Mat_Axis_Accel`, `Mat_Axis_Jerk`, `Mat_Seg_Accel`, `Mat_Seg_Jerk`, `Mat_SegAxis_Accel`, `Mat_SegAxis_Jerk`, `Mat_RPE`, `valid_subj`, `Mat_Seg_Accel_26`, `Mat_Seg_Jerk_26`, `valid_26`.

> Ce fichier est chargé sous le nom `d2` dans les scripts d'analyse (`Spearman_Corr_G1_G2.m`, `Lasso_LMM.m`, `Analyse_Spearman.m`).

---

### `Visualise_Workload_Li_TimeNorm_LD_Goubault.m` — 10 variables Goubault G2

**Ce que fait le script :** calcule depuis les `Features3_XSENS_S10X_Li.mat` (features pré-calculées, 1 timepoint/seconde) les **10 variables Goubault 2023** (Table III, tâche Chord) pour G2, puis sauvegarde et visualise.

**Pipeline :**
1. Chargement de `Features3_XSENS_S10X_Li.mat` (champ `Features_XSENS`)
2. Délimitation de la portion Liszt via `cycles_Li.mat` (timestamps `t_do`)
3. Découpage en 10 bins proportionnels
4. Interpolation RPE sur les bins
5. Normalisation par la moyenne des bins 1–3 (baseline des 3 premiers intervalles)

**10 variables calculées** (même définition que dans `Visualise_Goubault_Brut_G1vsG2.m` et `Compute_Goubault_NewDataSet.m`) :

| # | Segment | Feature | Signal | Colonne |
|---|---|---|---|---|
| 1 | Hand | MedianFreq | Acceleration Y | 20 |
| 2 | Hand | Percentile90 | Acceleration X | 19 |
| 3 | Hand | PeakPower | Module_Acceleration | 7 |
| 4 | Hand | SpectralEntropy | Acceleration Y | 20 |
| 5 | Head | PeakPower | Angular_Velocity X | 7 |
| 6 | Forearm | Mean (|valeur absolue|) | Acceleration Y | 17 |
| 7 | Forearm | SpectralEntropy | Module_Angular_Velocity | 6 |
| 8 | Forearm | PeakPower_Freq | Module_Angular_Velocity | 6 |
| 9 | Shoulder | Mean (|valeur absolue|) | Angular_Velocity X | 10 |
| 10 | Shoulder | Mean (|valeur absolue|) | Angular_Velocity Y | 11 |

> **Note sur la normalisation :** les variables `Mean` sont passées en valeur absolue avant normalisation, évitant la division par des valeurs proches de 0 qui produirait des ratios aberrants.

**Participants G2 :** `[3,8,10,11,12,13,14,16,18,19,24,25,26,27,33,34,36,40,41,43,47,48,50]` — N=23 (S117 absent, pas de fichier Features3).

**Figure produite :** `Fig_TimeNorm_Goubault_G2.png` — grille 2 lignes × 5 colonnes (10 variables). Double axe Y : variable normalisée (gauche) + RPE pointillés (droite). Droite de tendance linéaire rouge superposée. Code couleur par segment : Hand = rouge, Head = violet, Forearm = bleu, Shoulder = vert.

**Sortie `.mat` :** `Workload_Li_TimeNormalised_Goubault_G2.mat` — contient `Mat_Vars` (N × N_bins × 10), `Mat_RPE`, `valid_subj`, `Vars`, `N_bins`, `G2_Liszt`.

> Ce fichier est chargé sous le nom `dG2` dans les scripts d'analyse.

---

### `Visualise_EMG_Li_TimeNormalised_G2.m` — variables EMG G2

**Ce que fait le script :** calcule les **6 variables EMG** pour G2 depuis des fichiers pré-calculés par variable et par participant, puis sauvegarde et visualise.

**Variables EMG (6) :** Amplitude, Activity, Mobility, SampleEntropy, TFR_MedianFreq, TFR_SpectralEntropy.

**5 canaux musculaires :** Biceps, Triceps, DeltAnt, DeltMed, SupTrap. Les indices MVC sont aux positions 43–47 d'un vecteur de 47 éléments (`EMG_max`, contenant aussi les 42 canaux HD-EMG).

**Source des données :** fichiers pré-calculés par variable dans des dossiers séparés :
- `Amplitude_Li/Amplitude_<canal>_Li_EMG.mat`
- `Activity_Li/Activity<canal>_Li_EMG.mat`
- `Mobility_Li/Mobility<canal>_Li_EMG.mat`
- `SampleEntropy_Li/SampleEntropy_<canal>_Li_EMG.mat`
- `TFR_MedianFreq_Li/TFR_MedianFreq_<canal>_Li_EMG.mat`
- `TFR_SpectralEntropy_Li/TFR_SpectralEntropy_<canal>_Li_EMG.mat`
- Valeurs MVC : `EMG_Max_MVC/EMG_Max_S1XX_MVC.mat` (champ `EMG_max`)

**Traitement spécifique par variable :**
- `Amplitude` : normalisée en **% MVC** (valeur brute / MVC × 100). Les autres variables ne sont pas normalisées par MVC.
- `Activity` : normalisée par la médiane du bin 1 sur les figures (les valeurs brutes restent dans `Mat_EMG`).
- `SampleEntropy` : nettoyage des outliers ±3 SD sur les figures Median uniquement.

**Participants G2 :** `[3,8,10,11,12,13,14,16,18,19,24,25,26,27,33,34,36,40,41,43,47,48,50]` — N=23.

**Figures produites :** une figure par variable EMG (`Fig_EMG_Amplitude.png`, `Fig_EMG_Activity.png`, etc.) — grille 2 lignes (Mean, Median) × 5 colonnes (canaux). Double axe Y variable + RPE.

**Sortie `.mat` :** `EMG_Li_TimeNormalised_G2.mat` — contient `Mat_EMG` (N × N_bins × N_ch × N_vars), `Mat_RPE`, `valid_subj`, `Channels`, `VarNames`, `N_bins`, `G2_Liszt`, `Mat_MVC`.

> Ce fichier est chargé sous le nom `eG2` dans les scripts d'analyse.

---

## Fonctions locales communes

`Visualise_Workload_Li_TimeNorm_LD_Goubault.m` définit deux fonctions locales :

| Fonction | Rôle |
|---|---|
| `gStats(Mat)` | Moyenne et écart-type colonne par colonne, NaN si N < 3 participants valides dans un bin |
| `plotDual(T, Mu_v, SD_v, Mu_r, SD_r, col_v, lbl_v, ttl)` | Subplot avec double axe Y (variable + RPE), bande ±SD, droite de tendance linéaire rouge |

---

## Abréviations

| Abréviation | Signification dans ces scripts |
|---|---|
| G2 | Long Duration — participants ayant joué la session complète sans atteindre RPE ≥ 7 (N=23–24) |
| LD | Dans le nom de fichier (`Workload_Li_TimeNorm_LD`) : Long Duration (= G2). Ne pas confondre avec Large Hand (LD) du dataset Robin2. |
| RPE | Rate of Perceived Exertion, échelle CR-10 (0–10) |
| MVC | Maximum Voluntary Contraction — contraction maximale volontaire, valeur de référence pour normaliser l'Amplitude EMG |
| % MVC | Pourcentage de la contraction maximale — indicateur standardisé de l'effort musculaire |
| Features3 | Fichiers de features pré-calculées (`Features3_XSENS_S10X_Li.mat`), 1 timepoint/seconde |
| valid_26 | Masque booléen des participants G2 ayant un fichier Features3 disponible (sous-ensemble cohérent avec G1) |
| t_do | Timestamps des cycles musicaux (debut-onset) dans `cycles_Li.mat` |
| ByAxis | Somme des modules d'accélération sur les 7 segments pour chaque axe (X, Y, Z) |
| BySegMod | Module (norme euclidienne) par segment |
| Total | Somme sur les 21 colonnes = total cinématique toutes directions et tous segments |
| DeltAnt | Deltoïde antérieur |
| DeltMed | Deltoïde médial |
| SupTrap | Trapèze supérieur |
| TFR_MedianFreq | Fréquence médiane du spectre CWT — diminue avec la fatigue musculaire |
| TFR_SpectralEntropy | Entropie spectrale CWT — mesure l'irrégularité du spectre |
| SampleEntropy | Entropie de l'échantillon — mesure de complexité non-linéaire du signal |

---

## Fichiers d'entrée / sortie

| Fichier | Produit par | Consommé par (alias de chargement) |
|---|---|---|
| `Workload_Li_TimeNormalised_G2.mat` | `Visualise_Worload_Li_TimeNorm_LD.m` | Scripts d'analyse sous `d2` |
| `Workload_Li_TimeNormalised_Goubault_G2.mat` | `Visualise_Workload_Li_TimeNorm_LD_Goubault.m` | Scripts d'analyse sous `dG2` |
| `EMG_Li_TimeNormalised_G2.mat` | `Visualise_EMG_Li_TimeNormalised_G2.m` | Scripts d'analyse sous `eG2` |
