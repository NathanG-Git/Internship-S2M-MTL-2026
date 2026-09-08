# README — Pipeline PLSR spectral, dataset Robin2/Ergo

Ce groupe contient le pipeline **PLSR + VIP itératif sur features spectrales** du dataset Robin2 (23 sujets, Normal vs CS60). Il comprend deux scripts de calcul de features, deux scripts PLSR (avec et sans EMG), et un script de test de corrélation pour la réduction de variables.

---

## Vue d'ensemble et chaîne de données

```
Données brutes XSens (XSens_Fatigue_1..5.mat, Fs=60 Hz)
        │
        ├── Compute_SpectralFeatures_Ergo.m   ──→ SpectralFeatures_Ergo3.mat  (repère global, référence)
        └── Compute_SpectralFeatures_Ergo_Local.m ──→ SpectralFeatures_Ergo.mat  ⚠️ (repère local, écrase)
                                        ↑
                            (voir note sur la convention de nommage ci-dessous)

SpectralFeatures_Ergo3.mat + EMG_Features_NewDataset_25vars.mat
        └── PLSR_Ergo.m          ──→ PLSR_Results.mat + figures PNG

SpectralFeatures_Ergo.mat (ou Ergo2/3 selon run)
        └── PLSR_Ergo_SansEMG.m  ──→ PLSR_Results_SansEMG.mat + figures _SansEMG

SpectralFeatures_Ergo.mat + EMG_Features_NewDataset_25vars.mat
        └── Test_Correlation.m   ──→ console uniquement (pas de .mat)
```

### Convention de nommage des fichiers SpectralFeatures (historique)

| Fichier | Contenu | Produit par |
|---|---|---|
| `SpectralFeatures_Ergo.mat` | Version initiale (CWT par fenêtre 1s, sans filtre) — aussi écrasée par `Compute_SpectralFeatures_Ergo_Local.m` ⚠️ | `Compute_SpectralFeatures_Ergo_Local.m` (bug : devait sauvegarder sous Ergo4) |
| `SpectralFeatures_Ergo2.mat` | CWT corrigée mais instants RPE biaisés (bug rpe_frames) | version intermédiaire |
| `SpectralFeatures_Ergo3.mat` | **Référence courante** — CWT corrigée + instants équirépartis + repère global | `Compute_SpectralFeatures_Ergo.m` |

> **Attention :** `Compute_SpectralFeatures_Ergo_Local.m` sauvegarde sous `SpectralFeatures_Ergo.mat`, le même nom que la version initiale. C'est un bug identifié — il aurait dû sauvegarder sous `SpectralFeatures_Ergo4.mat`. Vérifier le contenu de `SpectralFeatures_Ergo.mat` avant utilisation. `PLSR_Ergo_SansEMG.m` charge `SpectralFeatures_Ergo.mat` — s'il est chargé après un run de `Ergo_Local`, il utilise les données en repère local.

---

## Scripts de calcul de features

---

### `Compute_SpectralFeatures_Ergo.m` — features spectrales, repère global (référence)

**Ce que fait le script :** calcule les **1260 features spectrales** du pipeline Goubault 2023 sur le dataset Robin2, en repère global (coordonnées XSens telles que fournies). C'est la version de référence pour le PLSR.

**1260 variables :** 7 segments × 3 signaux × 4 composantes × 15 features = 1260.
- Segments : L5, T8, Head, Shoulder, Arm, Forearm, Hand
- Signaux : Accel, AngVel, Jerk
- Composantes : X, Y, Z, Mod (module euclidien)
- 15 features : 8 temporelles + 7 CWT

**8 features temporelles** (calculées sur fenêtres de 1s non chevauchantes) : Mean, Std, Median, Max, p10, p25, p75, p90.

**7 features CWT** (calculées sur le scalogramme du signal continu entier) : MedianFreq, SpectralEntropy, Power_below4Hz, Power_above4Hz, TotalPower, PeakPower, PeakPowerFreq.

**Corrections par rapport aux versions précédentes :**
- Filtre Butterworth ordre 2 zero-lag passe-bande 0.5–29 Hz (forme SOS, stable près de Nyquist)
- CWT sur le **signal continu entier** (pas par fenêtre de 1s) : le scalogramme est calculé une seule fois, puis les features sont extraites instant par instant et moyennées par tranche de 1s
- Ondelette de Morlet complexe, nombre d'onde 8, grille linéaire 0.05 Hz (via `cwtft`, pas `cwtfilterbank` qui impose une grille logarithmique)
- Correction du bug `rpe_frames` : les instants d'extraction sont équirépartis (`linspace(0,1,N_bins) * N_frames`), indépendants de la forme concave RPE_shape

**Normalisation Goubault :** contrôlée par `norm_level` :
- `'feature'` (défaut) : normalisation appliquée à chaque série temporelle de features après calcul — lecture littérale de Goubault 2023
- `'signal'` : normalisation du signal brut avant calcul des features — comportement des versions antérieures

Facteur d'échelle : moyenne des 10 valeurs absolues les plus hautes de chaque série (`n_top_vals = 10`).

**Variable `RPE_shape` :** vecteur de 10 valeurs représentant la forme normalisée du RPE G1+G2 de Goubault `[0, 0.2757, ..., 1.0]`. Utilisée uniquement dans les scripts PLSR pour construire Y — ne sert plus à définir les instants d'extraction depuis la correction du bug.

**Sortie :** `SpectralFeatures_Ergo3.mat` — structure `Results_Spectral.(Sujet).(Condition).(Variable) = [1 × 10]`.

---

### `Compute_SpectralFeatures_Ergo_Local.m` — features spectrales, repère local

**Ce que fait le script :** variante de `Compute_SpectralFeatures_Ergo.m` qui **reprojette les données d'accélération et de vitesse angulaire dans le repère segment** (repère local) avant le calcul des features.

**Différence unique avec la version globale :**
La reprojection applique `a_local = R' * a_global` où R est la matrice de rotation segment→global extraite des quaternions d'orientation XSens (convention scalaire en dernier). Seules les composantes X, Y, Z changent entre repère global et local — le module (norme euclidienne) est invariant par rotation et sert de contrôle.

**Motivation :** Goubault et al. indiquent avoir utilisé les données en repère local, alors que les données XSens sont fournies en repère global. L'écart relatif moyen global/local des composantes est mesuré (~137 %) et affiché à la console en début de run.

**Tous les autres paramètres** (filtre, CWT, normalisation, bugs corrigés) sont identiques à `Compute_SpectralFeatures_Ergo.m`.

**Sortie :** `SpectralFeatures_Ergo.mat` — **même nom que la version initiale** (bug de nommage identifié, devait être `SpectralFeatures_Ergo4.mat`).

**Fonctions locales spécifiques :**
- `rot_to_local(v, R11..R33)` — reprojection d'un vecteur 3D dans le repère local via les 9 composantes de la matrice de rotation

---

## Scripts PLSR

---

### `PLSR_Ergo.m` — PLSR + VIP itératif, 1280 variables (IMU + EMG)

**Ce que fait le script :** ajuste un modèle PLSR avec réduction itérative de variables par VIP, sur 4 modèles de comparaison (taille de main × clavier), avec validation croisée groupée par sujet.

**Entrée :** `SpectralFeatures_Ergo3.mat` (1260 variables IMU) + `EMG_Features_NewDataset_25vars.mat` (25 variables EMG → 20 après suppression de DeltMed, doublon de DeltAnt).

**4 modèles :**
- **M1** : SH (Small Hand) / Norm — petites mains, clavier standard
- **M2** : SH / CS60 — petites mains, clavier ergonomique
- **M3** : LH (Large Hand) / Norm — grandes mains, clavier standard
- **M4** : Pool / Norm — tous sujets confondus, clavier standard

**Variable Y :** RPE reconstruit bin par bin via la forme normalisée Goubault (`RPE_shape`) mise à l'échelle entre `RPE_début` et `RPE_fin` de chaque sujet. L'objectif est de **caractériser les signatures cinématiques**, pas de prédire le RPE absolument.

**Pipeline en 2 passes :**

**Passe 1 — Sélection du meilleur niveau de réduction :**
- Modèle PLSR complet (toutes N variables) pour calculer le VIP initial
- Réduction successive par paliers de VIP (de `VIP_init = 1.0` jusqu'à convergence, pas `VIP_step = 0.1`)
- Pour chaque niveau : 50 répétitions de validation croisée groupée par sujet (LOSO-like), N_latent optimal déterminé en interne
- Métrique : MAE (Mean Absolute Error) sur le jeu de test, moyennée sur les 50 répétitions
- Sélection du niveau avec la MAE minimale

**Passe 2 — Estimation stable des VIP et coefficients :**
- 50 répétitions sur le sous-ensemble optimal retenu en Passe 1
- VIP et coefficients B calculés sur le jeu TRAIN de chaque fold
- B moyenné uniquement sur les folds où la feature est sélectionnée (VIP > seuil)
- Résultat : VIP_mean, VIP_std, VIP_p025, VIP_p975 (IC 95%) par variable

**Validation croisée :** GroupedKFold par sujet — les K folds sont définis sur les sujets entiers (un sujet = un fold dans la variante LOSO), évitant la fuite intra-sujet.

**Figures produites** (nommées avec suffixe standard, pas `_SansEMG`) :
- `Fig_Tableau_Reduction.png` — tableau des niveaux de réduction (MAE par niveau)
- `Fig_VIP_IC95.png` — VIP moyen ± IC 95% des variables retenues
- `Tableau_Top10_PLSR.png` — top 10 variables par VIP
- `Tableau_Reference.png` — variables de référence Goubault/Robin comparées

**Sortie :** `PLSR_Results.mat` — structure `Results_PLSR.(model_name)` avec les champs `reduction`, `feat_used`, `VIP_mean`, `VIP_std`, `VIP_p025`, `VIP_p975`, `feat_full`, `VIP_full`, `src_used`.

---

### `PLSR_Ergo_SansEMG.m` — PLSR + VIP itératif, 1260 variables (IMU uniquement)

**Ce que fait le script :** **identique à `PLSR_Ergo.m`** dans sa logique complète, avec une seule différence : les 20 variables EMG sont exclues. Seules les 1260 variables IMU spectrales sont utilisées.

**Motivation :** tester si l'EMG apporte de l'information prédictive au-delà de la cinématique. Si les performances (MAE, variables sélectionnées) sont similaires entre `PLSR_Ergo.m` et `PLSR_Ergo_SansEMG.m`, l'EMG n'apporte pas de signal supplémentaire.

**Différences avec `PLSR_Ergo.m` :**

| | `PLSR_Ergo.m` | `PLSR_Ergo_SansEMG.m` |
|---|---|---|
| Variables | 1280 (1260 IMU + 20 EMG) | 1260 (IMU uniquement) |
| Entrée spectrale | `SpectralFeatures_Ergo3.mat` | `SpectralFeatures_Ergo.mat` ⚠️ |
| Suffixe des sorties | (standard) | `_SansEMG` |
| Fichier `.mat` | `PLSR_Results.mat` | `PLSR_Results_SansEMG.mat` |

> **Attention :** `PLSR_Ergo_SansEMG.m` charge `SpectralFeatures_Ergo.mat` (pas Ergo3). Si ce fichier a été écrasé par `Compute_SpectralFeatures_Ergo_Local.m`, les résultats sont en repère local et non comparables à `PLSR_Ergo.m` qui utilise Ergo3 en repère global.

**Sorties** (suffixe `_SansEMG`) : `PLSR_Results_SansEMG.mat`, `Fig_Tableau_Reduction_SansEMG.png`, `Fig_VIP_IC95_SansEMG.png`, `Tableau_Top10_PLSR_SansEMG.png`, `Tableau_Reference_SansEMG.png`.

---

### `Test_Correlation.m` — test de redondance entre variables par sous-groupes

**Ce que fait le script :** calcule la matrice de corrélation de Spearman entre les variables d'un même sous-groupe (défini manuellement) pour détecter les redondances avant de relancer le PLSR sur un ensemble élagué.

**Ce que fait le script :** aucune sortie `.mat`, uniquement console. Pour chaque sous-groupe de variables (17 sous-groupes définis manuellement, organisés par segment anatomique), affiche la matrice de corrélation et signale les paires avec `|r| ≥ thresh` (défaut = 0.70).

**17 sous-groupes** couvrant : Tête (2), Épaule (4), Tronc (1), Bras (2), Avant-bras (4), Main (4). Chaque sous-groupe regroupe des variables dont la redondance est suspectée a priori (même feature sur des signaux voisins, mêmes signaux sur des features voisines).

**Sources des données :** `SpectralFeatures_Ergo.mat` (variables IMU) + `EMG_Features_NewDataset_25vars.mat` (variables EMG), via la fonction locale `collect_var` qui agrège toutes les observations (tous sujets × toutes conditions × tous bins).

**Fonctions locales :**
- `collect_var(F_sp, F_emg, vname)` — agrège toutes les valeurs d'une variable sur l'ensemble sujet × condition
- `test_group(F_sp, F_emg, num, gname, vars, thresh)` — affiche la matrice de corrélation et les paires redondantes pour un sous-groupe

**Usage :** ce script est un outil d'exploration interactif. Les sous-groupes et le seuil `thresh` sont modifiés manuellement entre runs pour affiner la liste des variables à éliminer avant de relancer le PLSR.

---

## Abréviations

| Abréviation | Signification dans ces scripts |
|---|---|
| PLSR | Partial Least Squares Regression — régression sur composantes latentes qui maximise la covariance X-Y |
| VIP | Variable Importance in Projection — score de contribution de chaque variable aux composantes PLSR |
| VIP itératif | Réduction successive par paliers de VIP : on retire les variables VIP < seuil, on réajuste, on compare la MAE |
| MAE | Mean Absolute Error — erreur absolue moyenne entre Y prédit et Y observé (en unités de RPE) |
| LOSO | Leave-One-Subject-Out — variante de k-fold où chaque fold contient un sujet entier |
| Passe 1 | Sélection du niveau de réduction optimal (meilleure MAE sur 50 répétitions par niveau) |
| Passe 2 | Estimation stable des VIP et coefficients B sur le sous-ensemble optimal (50 répétitions) |
| N_latent | Nombre de composantes latentes PLSR — optimisé par CV interne dans chaque fold |
| VIP_p025 / VIP_p975 | Percentiles 2.5 % et 97.5 % des VIP sur les 50 répétitions = IC 95 % de stabilité |
| RPE_shape | Vecteur de 10 valeurs `[0, 0.28, ..., 1.0]` — forme normalisée de la courbe RPE G1+G2 Goubault |
| Y reconstruit | RPE bin par bin = `RPE_début + (RPE_fin - RPE_début) × RPE_shape` par sujet |
| SH | Small Hand — petites mains (Subject_Both dans le dataset Robin2, N=13) |
| LH | Large Hand — grandes mains (Subject_LD, N=10) |
| Pool | Tous les sujets Norm confondus (SH + LH) |
| norm_level | `'feature'` (Goubault littéral) ou `'signal'` (versions antérieures) — niveau d'application de la normalisation |
| n_top_vals | Facteur d'échelle Goubault = moyenne des 10 valeurs absolues les plus hautes de la série |
| cwtft | Fonction MATLAB de CWT via FFT — impose une grille linéaire, contrairement à `cwtfilterbank` (grille log) |
| omega0 | Nombre d'onde de l'ondelette de Morlet complexe (= 8 dans Goubault 2023) |
| rpe_frames | Instants d'extraction des features : `round(linspace(0,1,N_bins) * N_frames)` — équirépartis depuis correction du bug |
| SansEMG | Variante du PLSR sans les 20 variables EMG — permet de tester l'apport de l'EMG sur la cinématique seule |
| DeltMed | Deltoïde médial — supprimé du PLSR car doublon confirmé de DeltAnt (même capteur dans le dataset Robin2) |
| Repère global | Coordonnées XSens telles que fournies — référentiel fixe du laboratoire |
| Repère local | Coordonnées reprojetées dans le référentiel du segment (`a_local = R' × a_global`) |
| thresh | Seuil de corrélation (|r| ≥ 0.70) au-delà duquel deux variables sont considérées redondantes |

---

## Fichiers d'entrée / sortie

| Fichier | Produit par | Consommé par |
|---|---|---|
| `SpectralFeatures_Ergo3.mat` | `Compute_SpectralFeatures_Ergo.m` | `PLSR_Ergo.m` |
| `SpectralFeatures_Ergo.mat` | `Compute_SpectralFeatures_Ergo_Local.m` ⚠️ (bug nommage) | `PLSR_Ergo_SansEMG.m`, `Test_Correlation.m` |
| `PLSR_Results.mat` | `PLSR_Ergo.m` | — |
| `PLSR_Results_SansEMG.mat` | `PLSR_Ergo_SansEMG.m` | — |
| `Fig_Tableau_Reduction*.png`, `Fig_VIP_IC95*.png`, `Tableau_Top10*.png`, `Tableau_Reference*.png` | `PLSR_Ergo.m` / `PLSR_Ergo_SansEMG.m` | — |

---

## Résultats PLSR — tableau comparatif (`Dimensions_Variables1.ods`)

Le fichier `Dimensions_Variables1.ods` compile les top variables sélectionnées par le PLSR sur plusieurs runs, organisées par segment anatomique. Il permet de comparer la stabilité des sélections entre jeux de features, modèles et variantes IMU-uniquement.

### Jeux de features comparés

| Jeu | Variables | Source spectrale | EMG |
|---|---|---|---|
| Ergo3 1280 var | 1260 IMU + 20 EMG | `SpectralFeatures_Ergo3.mat` (repère global, référence) | Oui |
| Ergo3 1260 var | 1260 IMU uniquement | `SpectralFeatures_Ergo3.mat` | Non |
| Ergo2 1280 var | 1260 IMU + 20 EMG | `SpectralFeatures_Ergo2.mat` (instants biaisés) | Oui |
| Ergo2 1260 var | 1260 IMU uniquement | `SpectralFeatures_Ergo2.mat` | Non |
| Run sans suffixe (4 modèles) | 1280 var | `SpectralFeatures_Ergo.mat` (repère local probable) | Oui |

### Performances par modèle (Ergo3, référence)

| Modèle | Description | %RPE (1280) | AbsErr (1280) | %RPE (1260) | AbsErr (1260) |
|---|---|---|---|---|---|
| M1 | SH / Norm | 83.91 % | 1.189 (0.068) | 86.9 % | 1.189 (0.056) |
| M2 | SH / CS60 | 86.89 % | 0.977 (0.031) | 96.3 % | 1.228 (0.079) |
| M3 | LH / Norm | 90.0 % | 1.050 (0.055) | — | — |
| M4 | Pool / Norm | 97.9 % | 1.140 (0.071) | — | — |

> **%RPE :** variance de Y expliquée par le modèle PLSR. **AbsErr :** erreur absolue moyenne (± SD sur 50 répétitions) entre RPE prédit et RPE reconstruit. M3 et M4 SansEMG non reportés dans le fichier.

### Variables récurrentes entre jeux et modèles (Ergo3)

Les variables listées ci-dessous apparaissent dans plusieurs modèles avec le même signe — elles constituent les candidats les plus robustes.

**Tête (Head) :**
- `Head_SpectralEntropy_AngVelX` ↑ (M1, M3 Ergo3 1280)
- `Head_PeakPowerFreq_AccelX` ↑ (M1, M3, M4 selon variante)
- `Head_Max_AccelY` ↑ (M4 1280, M3 1260)
- `Head_SpectralEntropy_AccelY` ↑ (M4 1260)

**Épaule (Shoulder / DeltAnt) :**
- `Activity_DeltAnt` ↓ (M1, M2, M4 Ergo3 1280 et 1260)
- `Shoulder_SpectralEntropy_JerkY` ↓ (M3, M4 Ergo3 1260)
- `TFR_SpectralEntropy_DeltAnt` — signe variable selon modèle (↑ M4 1280, ↓ M4 1260)

**Tronc (T8) :**
- `T8_MedianFreq_AccelX` ↑ (M3, M4 Ergo3 1260)
- `T8_PeakPowerFreq_AccelY` ↓ (M3, M4 Ergo3 1260)
- `T8_SpectralEntropy_JerkY` ↓ (M1 Ergo3 1260)

**Pelvis (L5) :**
- `L5_PeakPowerFreq_AccelZ` ↑ (M1 Ergo3 1280 et 1260)
- `L5_Median_AngVelY` ↑ (M1 Ergo3 1280 et 1260)
- `L5_MedianFreq_AccelX` ↑ (M3, M4 Ergo3 1260)
- `L5_MedianFreq_AccelY` ↑ (M3, M4 Ergo3 1260)
- `L5_p25_AngVelY` ↓ (M3 Ergo3 1260)
- `L5_SpectralEntropy_JerkY` ↓ (M1 Ergo3 1280 et 1260)

**Bras (Arm / Biceps / Triceps) :**
- `SampleEntropy_Biceps` ↓ (M1 Ergo3 1280 et 1260)
- `SampleEntropy_Triceps` ↑ (M2, M4 Ergo3 1280 et 1260)
- `TFR_MedianFreq_Biceps` ↓ (M1, M2 Ergo3 1280)
- `TFR_MedianFreq_Triceps` ↓ (M2 1280 ; ↑ M4 1260 — signe instable)
- `Mobility_Triceps` ↑ (M2, M4 Ergo3 1280 et 1260)
- `Arm_SpectralEntropy_JerkZ` ↓ (M4 Ergo3 1280 et 1260)
- `Arm_Median_AngVelZ` ↓ (M3 Ergo3 1260 ; M4 Ergo3 1260)

**Avant-bras (Forearm) :**
- `Forearm_Median_AccelZ` ↓ — la variable la plus récurrente : présente dans M1, M2, M3, M4 de Ergo3 1280 et 1260
- `Forearm_PeakPowerFreq_AngVelMagnitude` ↓ (M2, M3 Ergo3 1280 et 1260)
- `Forearm_SpectralEntropy_AngVelMagnitude` ↑ (M1 Ergo3 1280 et 1260)
- `Forearm_PeakPowerFreq_AccelMagnitude` ↓ (M2, M3 Ergo3 1260)

**Main (Hand) :**
- `Hand_Median_AccelZ` ↓ — présente dans M1, M2, M3, M4 de Ergo3 1280 et 1260
- `Hand_MedianFreq_AccelY` ↑ (M1 LMM Robin et M1 Goubault Chord)
- `Hand_PeakPower_AccelMagnitude` ↓ (M2, M3 Ergo3 1260)
- `Hand_PeakPowerFreq_AngVelMagnitude` ↓ (M3, M4 Ergo3 1260)


