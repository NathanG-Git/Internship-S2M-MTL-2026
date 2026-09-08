# README — Scripts de visualisation, dataset Goubault (G1 vs G2)

Ce groupe contient trois scripts de **visualisation pure** sur le dataset Goubault (jeu 1, 50 pianistes experts, tâche Li/Chord). Aucun ne produit de fichier `.mat` ni de résultat statistique — uniquement des figures PNG.

---

## Dataset commun

| Paramètre | Valeur |
|---|---|
| Dataset | Goubault (50 pianistes, S101–S150) |
| Tâche | Li (Liszt = tâche Chord) |
| G1 (Short Duration) | 26 participants : [1,2,4,5,6,7,9,15,20,21,22,23,28,29,30,31,32,35,37,38,39,42,44,45,46,49] |
| G2 (Long Duration) | 23 participants : [3,8,10,11,12,13,14,16,18,19,24,25,26,27,33,34,36,40,41,43,47,48,50] |
| Bins | 10 intervalles de durée égale, axe x exprimé en % de session (10 %, 20 %... 100 %) |
| Couleurs | G1 = rouge `[0.85 0.15 0.15]` | G2 = bleu `[0.15 0.45 0.85]` | RPE = pointillés (rouge foncé / bleu foncé) |

> **G1 = Short Duration** : participants ayant atteint RPE ≥ 7 deux fois de suite, arrêtés tôt. **G2 = Long Duration** : participants n'ayant jamais atteint ce seuil. Voir note d'ambiguïté dans README_Analyse.md (SD/LD est une abréviation différente pour la taille de main dans le dataset Robin2).

---

## Logique commune aux trois scripts

### Chargement des données

Chaque script charge `Features3_XSENS_S1XX_Li.mat` (features pré-calculées par participant) et aligne les données sur la fenêtre d'intérêt définie par `cycles_Li.mat` (timestamps des cycles musicaux) :

- `tp_s` / `tp_e` = premier et dernier frame du cycle d'intérêt
- Le RPE brut (mesuré toutes les 30 s, échelle CR-10) est **interpolé linéairement** sur la grille temporelle de la session, puis moyenné par bin
- Exception S122 (`iP == 22`) : la dernière mesure RPE est supprimée (artefact connu)

### Deux représentations produites systématiquement

**Figures BRUTES** : valeurs moyennées par bin, sans soustraction. Deux axes Y superposés : axe gauche = variable biomécanique (unité propre), axe droit = RPE (Borg CR-10, 0–10).

**Figures DELTA** : pour chaque participant, `signal(k) − signal(bin1)` avant moyennage de groupe. Tous les participants partent de 0, ce qui élimine les différences de niveau de base et met en évidence l'évolution relative. Le RPE est également centré sur bin 1 (`ΔRPE`). Une ligne pointillée noire horizontale à zéro est ajoutée comme référence.

### Deux statistiques de groupe disponibles

Chaque figure est produite en **deux versions** :
- **Mean** : moyenne des participants valides + écart-type
- **Median** : médiane des participants valides + écart-type (mesure de dispersion robuste, même si la centralité est la médiane)

Un participant est exclu d'un bin si N < 3 participants valides.

---

## Scripts

---

### `RPE.m` — évolution du RPE seul

**Ce que fait le script :** visualise uniquement le profil temporel du RPE pour G1 et G2, sans aucune variable biomécanique.

**Données d'entrée :** `cycles_Li.mat` + `info_participants_corrected.mat` (champ `.Liszt` = série RPE).

**Une seule figure produite :** `Fig_RPE_G1vsG2.png`
- Courbe moyenne G1 (rouge) et G2 (bleu) avec bandes ±1 SD
- Marqueurs circulaires sur chaque bin
- Axe Y : RPE Borg CR-10 (0–10)

Ce script est le plus simple du groupe : il ne contient pas de boucle sur des variables, et pas de double représentation Mean/Median ni de version Delta. Il sert de référence rapide pour la forme des courbes RPE des deux groupes.

---

### `Visualise_Accel_Jerk_Brut_G1vsG2.m` — accélération et jerk par segment

**Ce que fait le script :** visualise l'évolution temporelle du **module d'accélération** et du **jerk** pour les 7 segments corporels, G1 vs G2.

**Données d'entrée :** `Features3_XSENS_S1XX_Li.mat` (champs `feat.Mean.Module_Acceleration` et `feat.Mean.Module_Jerk`).

**7 segments** (dans l'ordre des colonnes de `Features3`) : L5, T8, Head, Shoulder, Arm, Forearm, Hand.

**4 figures produites :**

| Figure | Signal | Représentation |
|---|---|---|
| `Fig_TimeNorm_BySegment_Acceleration_Brut_G1vsG2.png` | Module Accélération | Valeurs brutes |
| `Fig_TimeNorm_BySegment_Jerk_Brut_G1vsG2.png` | Module Jerk | Valeurs brutes |
| `Fig_TimeNorm_BySegment_Acceleration_Delta_G1vsG2.png` | Module Accélération | Évolution depuis bin 1 |
| `Fig_TimeNorm_BySegment_Jerk_Delta_G1vsG2.png` | Module Jerk | Évolution depuis bin 1 |

Chaque figure : grille 2 lignes (Mean, Median) × 7 colonnes (segments). Double axe Y : variable à gauche, RPE à droite (pointillés).

---

### `Visualise_Goubault_Brut_G1vsG2.m` — 10 variables Goubault

**Ce que fait le script :** visualise les 10 variables identifiées dans Goubault et al. (2023, Table III) pour G1 vs G2. Ces mêmes 10 variables sont recalculées sur le dataset Robin2 par `Compute_Goubault_NewDataSet.m`.

**Données d'entrée :** `Features3_XSENS_S1XX_Li.mat` (accès par nom de champ feature + signal + colonne).

**10 variables** (définies dans la table `Vars` du script) :

| # | Nom affiché | Feature | Signal | Colonne |
|---|---|---|---|---|
| 1 | Hand MedianFreq Accel Y | MedianFreq | Acceleration | 20 |
| 2 | Hand Prc90 Accel X | Percentile90 | Acceleration | 19 |
| 3 | Hand PeakPower Accel Mod | PeakPower | Module_Acceleration | 7 |
| 4 | Hand SpectralEntropy Accel Y | SpectralEntropy | Acceleration | 20 |
| 5 | Head PeakPower AngVel X | PeakPower | Angular_Velocity | 7 |
| 6 | Forearm Mean Accel Y | Mean | Acceleration | 17 |
| 7 | Forearm SpectralEntropy AngVel Mod | SpectralEntropy | Module_Angular_Velocity | 6 |
| 8 | Forearm PeakPowerFreq AngVel Mod | PeakPower_Freq | Module_Angular_Velocity | 6 |
| 9 | Shoulder Mean AngVel X | Mean | Angular_Velocity | 10 |
| 10 | Shoulder Mean AngVel Y | Mean | Angular_Velocity | 11 |

> **Note sur `Mean` (variables 6, 9, 10) :** la valeur absolue est appliquée avant moyennage (`abs(blk)`) pour éviter l'annulation de signe entre participants, conformément à la méthode de Goubault.

**4 figures produites** (même logique que `Visualise_Accel_Jerk_Brut_G1vsG2.m`) :

| Figure | Représentation |
|---|---|
| `Fig_Goubault_Brut_G1vsG2_Mean.png` | Valeurs brutes, statistique Mean |
| `Fig_Goubault_Brut_G1vsG2_Median.png` | Valeurs brutes, statistique Median |
| `Fig_Goubault_Delta_G1vsG2_Mean.png` | Évolution depuis bin 1, statistique Mean |
| `Fig_Goubault_Delta_G1vsG2_Median.png` | Évolution depuis bin 1, statistique Median |

Chaque figure : grille 2 lignes × 5 colonnes (10 variables). Double axe Y.

---

## Fonctions locales communes

Les deux scripts de visualisation (`Visualise_*.m`) définissent les mêmes deux fonctions locales :

| Fonction | Rôle |
|---|---|
| `gStats(Mat)` | Retourne la moyenne (`nanmean`) et l'écart-type (`nanstd`) colonne par colonne. Met NaN si N < 3 participants valides dans un bin. |
| `gStatsMedian(Mat)` | Retourne la médiane (`nanmedian`) et l'écart-type colonne par colonne. Même seuil N < 3. |

`RPE.m` n'utilise pas ces fonctions (calcul direct en ligne).

---

## Abréviations

| Abréviation | Signification dans ces scripts |
|---|---|
| G1 | Short Duration — groupe de participants ayant fatigué vite (RPE ≥ 7 atteint deux fois de suite, N=26) |
| G2 | Long Duration — groupe ayant joué jusqu'au bout sans atteindre le seuil (N=23) |
| RPE | Rate of Perceived Exertion — perception de l'effort, échelle Borg CR-10 (0–10) |
| CR-10 | Borg Category-Ratio 10 — échelle de 0 à 10, seuil d'arrêt = 7 dans Goubault |
| ΔRPE | RPE(k) − RPE(bin 1) — évolution centrisée du RPE depuis le début de la session |
| Li / Liszt | Tâche Chord du dataset Goubault — jeu de la pièce de Liszt utilisée comme stimulus |
| t_do | Timestamps des cycles musicaux (debut-onset), utilisés pour délimiter la fenêtre d'analyse |
| S1XX | Identifiant participant Goubault (S101 à S150) |
| Features3 | Fichiers de features pré-calculées par participant (`Features3_XSENS_S1XX_Li.mat`) |
| Module | Norme euclidienne des 3 axes : `sqrt(x² + y² + z²)` |
| Brut | Sans normalisation ni centrage — valeurs dans leur unité physique d'origine |
| Delta | Valeur centrée sur le bin 1 : `signal(k) − signal(bin1)`, tous les participants partent de 0 |
| Mean / Median | Statistique de centralité utilisée pour la courbe de groupe |
| yyaxis | Double axe Y MATLAB — axe gauche pour la variable biomécanique, axe droit pour le RPE |
