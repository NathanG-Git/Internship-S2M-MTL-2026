# README — Scripts de calcul de variables

Ce groupe contient tous les scripts qui **calculent et sauvegardent des variables** à partir des données brutes XSens et EMG, sans analyse statistique. Ils alimentent les scripts d'analyse.

---

## Vue d'ensemble et liens entre scripts

```
Données brutes XSens (XSens_Fatigue_1..5.mat, EMG_cut.mat)
        │
        ├── Analyse_IMU.m          → IMU_Features_NewDataset.mat         (6 var, validation)
        ├── Compute_SegMod_SegAxis_NewDataSet.m → SegMod_SegAxis_Features_NewDataset.mat (58 var)
        ├── Compute_Goubault_NewDataSet.m       → Goubault_Features_NewDataset.mat      (10 var)
        ├── Analyse_EMG.m          → EMG_Features_NewDataset.mat          (3 var, validation)
        └── Compute_EMG_NewDataSet.m            → EMG_Features_NewDataset_25vars.mat    (25 var)

Données brutes XSens (S101, format Li)
        └── Calcul_variables.m     → Features_Calc_S101_Li.mat            (brouillon)

Fonctions utilitaires (appelées par tous les scripts ci-dessus) :
        ├── Compute_Median_Frequency.m
        ├── Compute_Spectral_Entropy.m
        └── Spectral_Entropy.m
```

---

## Dataset et contexte commun

Tous les scripts de ce groupe opèrent sur le **dataset Robin2/Ergo** : 23–24 sujets pianistes, deux conditions (Normal = clavier standard, CS60 = clavier ergonomique), design crossover pour 13 sujets. P07 est exclu de tous les scripts (session Normal tronquée à ~6784 frames au lieu de ~18400). Les features sont calculées sur **10 bins** (intervalles de 10 % de la durée de session chacun).

---

## Scripts

---

### `Analyse_IMU.m` — 6 variables IMU, validation

**Ce que fait le script :** calcule les 6 variables IMU retenues pour la validation externe, sur les 5 fichiers `XSens_Fatigue_*.mat`. Les variables sont réparties en deux blocs.

**Bloc SegMod** (statistique brute, pas de CWT) :
- `Accel_Mod_Head` — moyenne du module d'accélération de la tête par bin (sensorData, accélération libre sans gravité)
- `Accel_Mod_Hand` — idem pour la main droite

**Bloc Goubault** (CWT Morlet `cmor8-1`, Fs=60 Hz, f=0.05–15 Hz, z-score avant CWT) :
- `MedianFreq_Accel_Y_Hand` — fréquence médiane du spectre CWT, axe Y uniquement
- `PeakPower_Accel_Mod_Hand` — puissance maximale du spectre CWT, somme des 3 axes
- `PeakPower_AngVel_X_Head` — idem, vitesse angulaire tête axe X
- `SpectralEntropy_AngVel_Mod_Forearm` — entropie spectrale CWT, somme des 3 axes de l'avant-bras

**Sortie :** `IMU_Features_NewDataset.mat` — structure `Results_IMU.(Sujet).(Condition).(Variable) = [1 × 10]`

**Lien avec d'autres scripts :** `Analyse_Validation.m` et `Analyse_RPE.m` chargent ce fichier comme entrée principale.

---

### `Compute_SegMod_SegAxis_NewDataSet.m` — 58 variables cinématiques

**Ce que fait le script :** calcule un set élargi de variables cinématiques sur les mêmes 5 fichiers XSens, pour l'analyse Spearman cross-dataset (`Analyse_Spearman.m`). Pas de CWT — uniquement des statistiques de premier ordre.

**Bloc SegMod** (16 variables) — moyenne du module d'accélération ou de jerk par bin, pour 7 segments :
- `Accel_L5`, `Accel_T8`, `Accel_Head`, `Accel_Shoulder`, `Accel_Arm`, `Accel_Forearm`, `Accel_Hand` (+ idem `Jerk_*`)
- `Total_Accel`, `Total_Jerk` — somme (pas moyenne) des modules sur les 7 segments

**Bloc SegAxis** (42 variables) — moyenne du signal brut par axe et par segment : 7 segments × 3 axes (X, Y, Z) × {Accel, Jerk}. Exemple : `Accel_AxeX_Hand`, `Jerk_AxeZ_Shoulder`.

**Jerk :** calculé comme `diff(Accel) * Fs`, complété par une ligne de zéros en fin de signal.

**Module :** norme euclidienne des 3 axes `sqrt(x² + y² + z²)`.

**Mapping segments :** L5=2, T8=5, Head=7, Shoulder=8 (RightShoulder), Arm=9, Forearm=10, Hand=11 dans la hiérarchie XSens 23-segments.

**Sortie :** `SegMod_SegAxis_Features_NewDataset.mat` — deux structures `Results_SegMod` et `Results_SegAxis`, même convention `.(Sujet).(Condition).(Variable) = [1 × 10]`.

> **Différence avec `Analyse_IMU.m` :** `Analyse_IMU.m` calcule 6 variables retenues pour la validation (dont 4 avec CWT). `Compute_SegMod_SegAxis_NewDataSet.m` calcule 58 variables cinématiques sans CWT, pour l'exploration Spearman sur les 5 groupes.

---

### `Compute_Goubault_NewDataSet.m` — 10 variables Goubault

**Ce que fait le script :** recalcule sur le dataset Robin2/Ergo les 10 variables identifiées dans Goubault et al. (2023, Table III) comme biomarqueurs de fatigue pianistique sur la tâche Chord. Ces 10 variables combinent des statistiques simples et des features CWT Morlet.

**Paramètres CWT :** identiques à `Calcul_variables.m` (jeu 1 validé) : Fs=60 Hz, `cmor8-1`, f=0.05–15 Hz. Les 4 segments utilisés : Head=7, Shoulder=8, Forearm=10, Hand=11.

**Variables calculées par bin :**

| # | Segment | Feature | Signal |
|---|---|---|---|
| 1 | Hand | MedianFreq (CWT) | Acceleration Y |
| 2 | Hand | Percentile 90 | Acceleration X |
| 3 | Hand | PeakPower (CWT, somme 3 axes) | Module Acceleration |
| 4 | Hand | SpectralEntropy (CWT) | Acceleration Y |
| 5 | Head | PeakPower (CWT) | Angular Velocity X |
| 6 | Forearm | Mean (valeur absolue) | Acceleration Y |
| 7 | Forearm | SpectralEntropy (CWT, somme 3 axes) | Module Angular Velocity |
| 8 | Forearm | PeakPowerFreq (CWT, somme 3 axes) | Module Angular Velocity |
| 9 | Shoulder | Mean (valeur absolue) | Angular Velocity X |
| 10 | Shoulder | Mean (valeur absolue) | Angular Velocity Y |

**Sortie :** `Goubault_Features_NewDataset.mat` — structure `Results_Goubault.(Sujet).(Condition).(Variable) = [1 × 10]`.

---

### `Analyse_EMG.m` — 3 variables EMG, validation

**Ce que fait le script :** calcule les 3 variables EMG retenues pour la validation externe, depuis `EMG_cut.mat` (données brutes Robin2). Intègre un LMM préliminaire et une section d'index de bloc en fin de script. Ce script produit directement `EMG_Features_NewDataset.mat`, le fichier chargé par `Analyse_Validation.m`.

**Pipeline de traitement :**
1. Filtre Butterworth passe-bande 10–400 Hz (ordre 2, Fs=2000 Hz)
2. Rééchantillonnage 2000 → 1000 Hz (interpolation spline)
3. Découpage en 10 bins proportionnels

**3 variables calculées par bin :**
- `TFR_MedianFreq_Triceps` (colonne 5) — fréquence médiane du spectre CWT Morlet sur le triceps ; la fréquence médiane EMG diminue avec la fatigue musculaire.
- `TFR_SpectralEntropy_Deltoid` (colonne 6) — entropie spectrale CWT sur le deltoïde.
- `SampleEntropy_Biceps` (colonne 4) — entropie de l'échantillon sur le biceps ; sous-échantillonnage à 100 Hz avant calcul, z-score, tolérance r=0.2.

**Sorties :** `EMG_Features_NewDataset.mat` + `LMM_Results_NewDataset.mat` (si le LMM réussit).

> **Différence avec `Compute_EMG_NewDataSet.m` :** même dataset Robin2/EMG_cut, mais `Analyse_EMG.m` calcule uniquement les 3 variables de validation avec un LMM intégré. `Compute_EMG_NewDataSet.m` calcule 25 variables (5 features × 5 canaux) sans LMM, pour l'analyse Spearman cross-dataset.

---

### `Compute_EMG_NewDataSet.m` — 25 variables EMG

**Ce que fait le script :** calcule un set élargi de 25 variables EMG (5 features × 5 canaux musculaires) depuis `EMG_cut.mat`, pour `Analyse_Spearman.m`. Pas de LMM.

**Pipeline de traitement :**
1. Filtre Butterworth passe-bande 10–400 Hz + filtres coupe-bande (notch) à 60, 120 et 180 Hz
2. Enveloppe : redressement (valeur absolue) + filtre passe-bas 9 Hz
3. Sous-échantillonnage 2000 → 1000 Hz avant CWT (factor 2)
4. Découpage en 10 bins

**5 features par canal :**
- `Activity` — variance de l'enveloppe par bin
- `Mobility` — variance de la dérivée de l'enveloppe / variance de l'enveloppe
- `SampleEntropy` — resample 2000→100 Hz, z-score, SampEn(m=2, r=0.2)
- `TFR_MedianFreq` — CWT Morlet 2–450 Hz à Fs=1000 Hz, fréquence médiane
- `TFR_SpectralEntropy` — idem, entropie spectrale

**5 canaux :** Biceps (col 4), Triceps (col 5), DeltAnt (col 6), DeltMed (col 6, dupliqué — un seul capteur deltoïde dans ce dataset), SupTrap (col 7).

> **Attention — DeltAnt et DeltMed dupliqués :** les deux canaux pointent vers la même colonne 6 (Deltoid). C'est une limitation du montage : un seul capteur pour le deltoïde, mappé sur deux noms de canaux par cohérence avec le jeu de données 1 (qui disposait de deux capteurs distincts).

**Sortie :** `EMG_Features_NewDataset_25vars.mat` — structure `Results_EMG.(Sujet).(Condition).(Variable) = [1 × 10]`.

---

### `Calcul_variables.m` — brouillon de vérification, dataset Li (S101)

**Ce que fait le script :** script de **brouillon** servant à recalculer les features du dataset Goubault (jeu 1, Li) depuis les signaux bruts XSens d'un seul participant (S101), afin de vérifier la conformité avec les features de référence `Features3`. Il n'est pas intégré dans le pipeline de validation.

**Ce qu'il calcule (15 features × 6 signaux) :**

- **Statistiques temporelles (8)** : Peak, Mean, Median, STD, P10, P25, P75, P90
- **Features CWT (7)** : MedianFreq, SpectralEntropy, PowerLF (0.1–4 Hz), PowerHF (6–12 Hz), PowerTot, PeakPower, PeakPowerFreq

**Signaux (6)** : Acceleration (axe X), Module_Acceleration, Angular_Velocity (axe X), Module_Angular_Velocity, Jerk (axe X), Module_Jerk.

**Section de vérification :** compare `MedianFreq.Module_Acceleration` et `Mean.Module_Acceleration` calculés ici avec les valeurs de `Features3_XSENS_Cycles_S101_Li.mat` (référence) via une corrélation.

**Sortie :** `Features_Calc_S101_Li.mat` (local, pour inspection).

> **Note :** ce script est un outil de développement. Le mapping des colonnes XSens est simplifié (colonnes 4,7,10... prises directement) et marqué comme approximation à corriger. Ne pas utiliser les sorties dans le pipeline de validation.

---

## Fonctions utilitaires

### `Compute_Median_Frequency.m`

Calcule la fréquence médiane d'un scalogramme CWT frame par frame. Prend en entrée `TFR` (matrice fréquences × temps en valeur absolue) et `Wave_FreqS` (axe des fréquences). Retourne un vecteur de N_frames fréquences médianes : pour chaque instant, la fréquence à laquelle l'intégrale cumulée du spectre atteint 50 % de l'intégrale totale (méthode des trapèzes). Auteur : Robin Mailly.

### `Compute_Spectral_Entropy.m`

Wrapper qui appelle `Spectral_Entropy` pour chaque colonne (instant) d'un scalogramme. Retourne un vecteur de N_frames entropies spectrales.

### `Spectral_Entropy.m`

Calcule l'entropie spectrale d'un spectre de puissance (log₂, normalisée par log₂(N)). Tronque le spectre à 30 Hz avant calcul. Formule : `H = -Σ(p·log₂(p)) / log₂(N)`, avec `p = DSP / Σ(DSP)` (densité de probabilité normalisée). Une valeur proche de 1 indique un spectre plat (irrégularité maximale) ; proche de 0, un spectre concentré sur quelques fréquences. Auteur : Robin Mailly.

### `ouverture_donnee_xsens.m`

Script de chargement des 5 fichiers `XSens_Fatigue_*.mat` dans des variables `F1`–`F5`. Script utilitaire simple, probablement utilisé lors de l'exploration initiale des données. Auteur : Jérémy (chemin `C:\Users\jerem\...`).

---

## Abréviations communes

| Abréviation | Signification |
|---|---|
| IMU | Inertial Measurement Unit — capteur XSens mesurant accélération et vitesse angulaire |
| EMG | Électromyographie — activité musculaire |
| CWT | Continuous Wavelet Transform — transformée en ondelettes continues |
| cmor8-1 | Ondelette de Morlet complexe, nombre d'onde 8, paramètre de bande passante 1 |
| Fs | Fréquence d'échantillonnage |
| Bin | Intervalle temporel normalisé (1/10e de la durée de session) |
| SegMod | Par segment, module (norme euclidienne des 3 axes) |
| SegAxis | Par segment et par axe (X, Y, Z séparément) |
| TFR | Time-Frequency Representation — scalogramme CWT |
| SampEn / SampleEntropy | Sample Entropy — mesure de complexité non-linéaire d'un signal |
| Normal / Norm | Condition clavier standard |
| CS60 | Condition clavier ergonomique |
| Crossover | Sujets ayant joué les deux conditions (N=13) |
| P07 | Sujet exclu — session Normal tronquée (~6784 frames au lieu de ~18400) |
| SD | Small Duration ou Small Hand selon le script — voir README Analyse |
| LD | Long Duration ou Large Hand selon le script — voir README Analyse |
