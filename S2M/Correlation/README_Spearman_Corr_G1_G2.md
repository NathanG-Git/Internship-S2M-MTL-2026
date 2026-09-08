# README — `Spearman_Corr_G1_G2.m` et `Indiv_Spearman_Corr_G1_G2.m`

## Vue d'ensemble

Ces deux scripts calculent des **corrélations de Spearman entre chaque variable cinématique/EMG et le RPE**, séparément pour G1 et G2, sur le dataset Goubault (50 pianistes experts). Ils produisent des figures PNG de tableaux visuels et une figure de synthèse des variables significatives. Ils couvrent les mêmes variables, les mêmes groupes, et produisent la même structure de sorties — **la seule différence est la méthode de corrélation**.

| | `Spearman_Corr_G1_G2.m` | `Indiv_Spearman_Corr_G1_G2.m` |
|---|---|---|
| **Méthode** | Corrélation sur la courbe moyenne du groupe (1 rho par variable par groupe) | Corrélation individuelle par participant, moyennée via transformation de Fisher (1 rho moyen par variable par groupe) |
| **Test de significativité** | p-value directe de la corrélation de Spearman | t-test à un échantillon sur les z de Fisher (H0 : z_moyen = 0) |
| **Interprétation** | Corrélation entre la forme temporelle moyenne et le RPE moyen | Robustesse du lien RPE-variable au niveau individuel |
| **Préfixe des sorties** | `Spearman_Tab*.png`, `Spearman_Synthese.png` | `Indiv_Spearman_Tab*.png`, `Indiv_Spearman_Synthese.png`, `Indiv_Spearman_Results.mat` |

---

## Dataset

| Paramètre | Valeur |
|---|---|
| Dataset | Goubault (50 pianistes experts) |
| G1 (Short Duration) | Participants ayant atteint RPE ≥ 7 deux fois de suite sur la tâche Chord — arrêtés tôt (n ≈ 26–27) |
| G2 (Long Duration) | Participants n'ayant jamais atteint RPE ≥ 7 — ayant joué jusqu'au bout (n ≈ 23–25) |
| Données IMU | `Workload_Li_TimeNormalised.mat` (G1) / `Workload_Li_TimeNormalised_G2.mat` (G2) |
| Données Goubault | `Workload_Li_TimeNormalised_Goubault.mat` / `_G2.mat` |
| Données EMG | `EMG_Li_TimeNormalised_G1.mat` / `EMG_Li_TimeNormalised_G2.mat` |

> **Attention — ambiguïté de l'abréviation G1/G2 :** dans ce script, G1 = *Short Duration* (fatigue rapide) et G2 = *Long Duration* (fatigue lente). Ces deux groupes sont définis par le **régime de fatigue**, pas par la taille des mains. Dans d'autres scripts du projet (dataset Robin2/Ergo), SD/LH désignent la **taille de main** (Small Hand / Large Hand). Ne pas confondre.

---

## Ce que font les scripts

### 1. Chargement
Six fichiers `.mat` sont chargés : données IMU normalisées temporellement (G1 et G2), variables Goubault (G1 et G2), et données EMG (G1 et G2). Les indices de participants valides (`valid_subj`, `valid_26`) sont extraits pour chaque source.

> `valid_26` correspond à un sous-ensemble de 26 participants (G1) / 25 (G2) pour lequel les variables `Features3` (Total et BySegMod) sont disponibles. Les variables ByAxis et BySegAxis utilisent l'ensemble complet `valid_subj`.

### 2. Quatre tableaux de variables

| Tableau | Variables | Fichier PNG |
|---|---|---|
| Tab1 | Cinématique : Total, BySegMod, ByAxis | `*_Tab1_Cinematique.png` |
| Tab2 | Cinématique BySegAxis — une figure par axe (X, Y, Z) | `*_Tab2_BySegAxis_X/Y/Z.png` |
| Tab3 | EMG (toutes variables sauf Amplitude, tous canaux) | `*_Tab3_EMG.png` |
| Tab4 | Variables Goubault (spectrales, depuis les `.mat` Goubault) | `*_Tab4_Goubault.png` |

### 3. Calcul de la corrélation

**`Spearman_Corr_G1_G2.m` — corrélation sur moyenne de groupe (`spear`)**
- Pour chaque variable et chaque groupe, la moyenne temporelle des participants valides est calculée (vecteur de N_bins valeurs).
- La corrélation de Spearman est calculée entre ce vecteur moyen et le RPE moyen du groupe.
- Résultat : 1 rho et 1 p-value par variable par groupe.
- Minimum : 4 paires de points valides.

**`Indiv_Spearman_Corr_G1_G2.m` — corrélation individuelle + Fisher (`spearIndiv`)**
- Pour chaque participant valide, la corrélation de Spearman est calculée entre son propre profil temporel de variable et son propre profil RPE.
- Chaque rho individuel est transformé en z de Fisher : `z = atanh(rho)`. Les rho sont bornés à ±0.9999 pour éviter les infinis.
- Un t-test à un échantillon est appliqué sur les z (H0 : z_moyen = 0).
- La moyenne des z est retransformée en rho moyen : `rho_moy = tanh(mean(z))`.
- Un participant est exclu si ses données contiennent moins de 8 bins valides.
- Minimum : 3 participants valides par groupe pour produire un résultat.

### 4. Correction FDR (`computeFDR`)
Identique dans les deux scripts. La correction de Benjamini-Hochberg est appliquée **séparément** sur les p-values de G1 et sur celles de G2. Les q-values sont calculées et la monotonie est imposée (cummin depuis le rang le plus élevé). Une variable passe le FDR si `q < 0.05`.

### 5. Score composite (`sc`)
Affiché dans la colonne `sc` de chaque tableau :

```
sc = |rho| × sqrt(n / n_max)
```

où `n_max` est l'effectif maximal du groupe (26 pour G1, 23 pour G2). Ce score pondère la magnitude de la corrélation par la puissance statistique disponible. Un score ≥ 0.8 est mis en vert dans le tableau.

### 6. Visualisation (`plotTable`)
Identique dans les deux scripts (mise en page quasi identique, largeur de figure légèrement différente : 1050 px vs 1300 px). Chaque ligne du tableau affiche : variable, rho, p, q, étoiles, n, score, pour G1 et G2. Le fond des colonnes G1 est coloré en rouge pâle selon le niveau de significativité (q < 0.001 → rouge, q < 0.01 → orange, q < 0.05 → jaune). G2 suit le même code couleur en bleu dans les étoiles.

### 7. Figure de synthèse
Toutes les variables significatives (q < 0.05 pour G1 **ou** G2) sont rassemblées dans une figure unique, triée par ordre anatomique proximal → distal (L5 → T8 → Head → Shoulder → Arm → Forearm → Hand, puis EMG/Goubault en fin).

### 8. Sauvegarde supplémentaire (`Indiv_Spearman_Corr_G1_G2.m` uniquement)
Les variables significatives sont sauvegardées dans `Indiv_Spearman_Results.mat` (structure `Indiv_Sig` contenant labels, data, qvals, npart, fdr_pass). Ce fichier peut être chargé par d'autres scripts pour réutiliser les résultats sans relancer le calcul.

---

## Fonctions locales

| Fonction | Présente dans | Rôle |
|---|---|---|
| `spear(x, y)` | `Spearman_Corr_G1_G2.m` uniquement | Corrélation de Spearman sur deux vecteurs, avec exclusion des NaN |
| `spearIndiv(Mat_var, Mat_RPE, Vi)` | `Indiv_Spearman_Corr_G1_G2.m` uniquement | Corrélation individuelle + transformation de Fisher + t-test |
| `computeFDR(data, alpha)` | Les deux | Correction Benjamini-Hochberg sur les p-values de G1 et G2 |
| `plotTable(...)` | Les deux | Génère et sauvegarde la figure PNG du tableau visuel |
| `getStars(p)` | Les deux | Convertit une p-value en étoiles (* ** ***) |
| `getSigColor(p)` | Les deux | Retourne une couleur de fond selon le niveau de significativité |

---

## Abréviations

| Abréviation | Signification dans ces scripts |
|---|---|
| G1 | Short Duration — participants ayant atteint le seuil de fatigue (RPE ≥ 7) tôt. **Ici : régime de fatigue, pas taille de main.** |
| G2 | Long Duration — participants n'ayant jamais atteint le seuil de fatigue. **Ici : régime de fatigue, pas taille de main.** |
| RPE | Rate of Perceived Exertion — échelle de perception de l'effort (CR-10 dans Goubault, ×10 pour comparaison avec les données Robin en CR-100) |
| rho | Coefficient de corrélation de Spearman (sans unité, entre -1 et 1) |
| z | Transformée de Fisher de rho : `z = atanh(rho)`. Permet de moyenner des corrélations et d'appliquer un t-test. |
| FDR | False Discovery Rate — taux de faux positifs contrôlé par la correction de Benjamini-Hochberg |
| q | Q-value = p-value ajustée après correction FDR |
| sc | Score composite : `|rho| × sqrt(n/n_max)` — combine magnitude et puissance statistique |
| BySegMod | Par segment, norme du vecteur (Modulus) — somme quadratique des composantes X, Y, Z |
| ByAxis | Par axe global (X, Y, Z), tous segments confondus |
| BySegAxis | Par segment ET par axe — la combinaison la plus fine |
| Total | Somme des normes de tous les segments |
| Features3 | Sous-ensemble de variables calculées sur n=26 (G1) / 25 (G2) participants (subset avec données complètes) |
| valid_subj | Masque booléen indiquant les participants avec données valides dans le fichier `.mat` |
| valid_26 | Masque booléen indiquant le sous-ensemble de 26/25 participants pour les Features3 |
| IMU | Inertial Measurement Unit — capteur XSens mesurant accélération et vitesse angulaire |
| EMG | Électromyographie — mesure de l'activité musculaire |
| TFR | Time-Frequency Representation — features spectrales calculées par CWT (voir scripts Compute_SpectralFeatures) |

---

## Fichiers d'entrée / sortie

| Fichier | Rôle |
|---|---|
| `Workload_Li_TimeNormalised.mat` | Données IMU normalisées G1 |
| `Workload_Li_TimeNormalised_G2.mat` | Données IMU normalisées G2 |
| `Workload_Li_TimeNormalised_Goubault.mat` | Variables spectrales Goubault G1 |
| `Workload_Li_TimeNormalised_Goubault_G2.mat` | Variables spectrales Goubault G2 |
| `EMG_Li_TimeNormalised_G1.mat` | Données EMG normalisées G1 |
| `EMG_Li_TimeNormalised_G2.mat` | Données EMG normalisées G2 |
| `Spearman_Tab1–4_*.png` / `Indiv_Spearman_Tab1–4_*.png` | Tableaux visuels par catégorie de variables |
| `Spearman_Synthese.png` / `Indiv_Spearman_Synthese.png` | Figure de synthèse des variables significatives |
| `Indiv_Spearman_Results.mat` | *(Indiv uniquement)* Structure des variables significatives, réutilisable par d'autres scripts |

---

## Lien avec d'autres scripts

- **`Compute_SpectralFeatures_Ergo.m`** : génère les features spectrales (Tableau 4 Goubault) utilisées en entrée.
- **`PLSR_Ergo.m` / `PLSR_Ergo_SansEMG.m`** : utilisent les mêmes groupes G1/G2 et la même définition du RPE comme variable Y.
- **`Analyse_Validation.m`** (dataset Robin2/Ergo) : les groupes SD/LH de ce script correspondent à la **taille de main**, pas au régime de fatigue — ne pas confondre avec G1/G2 ici.
