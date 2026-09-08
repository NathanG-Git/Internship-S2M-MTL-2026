# README — `ComputeJointAngles.m`

## Vue d'ensemble

Script de **calcul des angles articulaires** à partir des données de capture de mouvement XSens, pour le dataset Expression musicale. Il reconstruit les matrices de transformation homogène de chaque segment corporel, puis en dérive les angles articulaires (décomposition d'Euler) pour l'ensemble des articulations du membre supérieur droit et gauche, du tronc et des membres inférieurs.

> **Auteur original : Robin Mailly (ISM / S2M).** Ce script a été fourni en référence dans le cadre du stage. Le chemin `'/Users/robin/...'` correspond à son poste de travail.

---

## Dataset

| Paramètre | Valeur |
|---|---|
| Participants | P01, P02, P04–P14 (13 sujets) |
| Conditions | Competition, Extrait (IE/PS), Chopin |
| Takes par défaut | A, B, C |
| Exceptions de takes | P05, P11, P14 (Extrait) : A et B uniquement ; P10 (Chopin) : B et C uniquement |
| Segments XSens utilisés | 23 segments corporels (Pelvis → orteils) |

---

## Ce que fait le script

### 1. Chargement et fusion des données
Même logique que `Plot_Expression.m` : deux fichiers `.mat` (P01–P08 et P09–P14) sont chargés et fusionnés. Ici la fusion est faite avec `catstruct` (disponible sur le poste de Robin, à remplacer par une fusion manuelle si non disponible).

### 2. Définition de la hiérarchie squelettique
La table `Joints` définit les paires (parent, enfant) pour chaque articulation. Par exemple :

| Nom | Parent | Enfant |
|---|---|---|
| rGH | rShould (8) | rUpperArm (9) |
| rEL | rUpperArm (9) | rForeArm (10) |
| rWR | rForeArm (10) | rHand (11) |

La séquence d'Euler utilisée est `yxz` pour toutes les articulations.

### 3. Correction d'orientation (T-pose)
Pour chaque essai, l'orientation du pelvis en T-pose est extraite (quaternion → matrice de rotation). La composante de rotation autour de l'axe Z (yaw, `theta_z`) est calculée et une matrice de correction `R_correction` est appliquée à tous les segments, de façon à recaler le sujet dans un repère global cohérent quelle que soit son orientation initiale.

### 4. Calcul des matrices de transformation (`compute_joint_angles`)
Pour chaque frame et chaque segment :
- Le quaternion d'orientation XSens est converti en matrice de rotation 3×3.
- `R_correction` est appliquée (produit matriciel).
- La matrice homogène 4×4 `RT` est construite en ajoutant la position 3D du segment.

Les angles articulaires sont ensuite calculés comme la rotation relative entre segment enfant et segment parent (`RJ = R_enfant × R_parent^T`), décomposée selon la séquence `yxz` via `angleRotation`. Les articulations du bras (`rUpperArm`, `rForeArm`, `lUpperArm`, `lForeArm`) sont dépliées (`unwrap`) pour éviter les discontinuités, puis converties en degrés.

### 5. Stockage
Les résultats sont stockés dans trois structures :

| Structure | Contenu |
|---|---|
| `JointAngles` | Angles articulaires calculés (matrice frames × 3 axes × articulations) |
| `XSens_Angles` | Angles articulaires bruts XSens (pour comparaison) |
| `XSens_Segments` | Données brutes de segments (orientation, position) |

L'arborescence suit le schéma `structure.sujet.condition[.expression.take]`.

### 6. Visualisation (optionnelle, inactive)
La fonction `plot_joint_axes` (commentée dans la boucle principale) permet de visualiser les axes locaux de chaque segment et le squelette en 3D pour une frame donnée. Elle n'est pas appelée par défaut.

---

## Fonctions locales

| Fonction | Rôle |
|---|---|
| `compute_joint_angles` | Calcule les matrices RT et les angles articulaires pour un essai |
| `plot_joint_axes` | Visualise les axes locaux et le squelette 3D (inactive par défaut) |

---

## Abréviations

| Abréviation | Signification dans ce script |
|---|---|
| RT | Matrice de transformation homogène 4×4 (Rotation + Translation) |
| RJ | Matrice de rotation articulaire relative (enfant × parent^T) |
| JA | Joint Angles — angles articulaires |
| rGH / lGH | Articulation gléno-humérale droite / gauche (épaule) |
| rEL / lEL | Articulation du coude droit / gauche |
| rWR / lWR | Articulation du poignet droit / gauche |
| rST / lST | Articulation scapulo-thoracique droite / gauche |
| PL | Articulation pelvis-L5 (lombo-sacrée) |
| TH | Articulation thoracique (Pelvis–T8) |
| NK | Articulation cervicale (T8–Head) |
| rCF / lCF | Articulation coxo-fémorale droite / gauche (hanche) |
| rKN / lKN | Articulation du genou droit / gauche |
| rAK / lAK | Articulation de la cheville droite / gauche |
| T-pose | Position de référence debout bras écartés utilisée pour la calibration XSens |
| tpose | Champ de la structure XSens contenant l'orientation et la position de référence |
| yxz | Séquence de décomposition d'Euler (rotation autour de Y, puis X, puis Z) |
| yaw | Rotation autour de l'axe vertical (Z) — composante corrigée par R_correction |
| IE | Interprétation Expressive |
| PS | Performance Standard |
| catstruct | Fonction MATLAB tierce pour concaténer deux structures (à remplacer si non disponible) |

---

## Fichiers d'entrée / sortie

| Fichier | Rôle |
|---|---|
| `XSens_MIDI_cut_P01_P08.mat` | Données XSens P01–P08 |
| `XSens_MIDI_cut_P09_P14.mat` | Données XSens P09–P14 |
| `JointAngles.mat` | Angles articulaires calculés (sauvegarde conditionnelle) |
| `XSens_Angles.mat` | Angles articulaires bruts XSens |
| `XSens_Segments.mat` | Données brutes de segments |

La sauvegarde n'est effectuée que si `save_data = true` (ligne commentée par défaut).

---

## Lien avec d'autres scripts

- **`Plot_Expression.m`** : utilise la même logique de chargement/fusion des deux `.mat` et le même identifiant de segment (`rHand = 11`) pour accéder aux données d'accélération.
