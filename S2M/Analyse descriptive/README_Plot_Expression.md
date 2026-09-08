# README — `Plot_Expression.m`

## Vue d'ensemble

Script de **visualisation** de la variable `Accel_Mod_Hand` (norme de l'accélération de la main droite, segment 11 XSens) pour le dataset Expression musicale. Il produit une figure par participant, exportée en `.png`.

---

## Dataset

| Paramètre | Valeur |
|---|---|
| Participants | P01, P02, P04–P14 (13 sujets, P03 absent) |
| Conditions | Competition, Extrait (IE + PS), Chopin ignoré |
| Takes par défaut | A, B, C |
| Exceptions de takes | P14 : B uniquement ; P05 et P11 : A et B uniquement |
| Segment analysé | rHand = segment 11 (XSens) |

---

## Ce que fait le script

### 1. Chargement et fusion des données
Les données XSens sont réparties en deux fichiers `.mat` (P01–P08 et P09–P14). Le script les charge séparément puis fusionne les structures manuellement champ par champ (`XSens_cur`). Cette fusion manuelle remplace la fonction `catstruct` (non disponible dans toutes les installations MATLAB).

### 2. Calcul du signal
Pour chaque essai, le signal est calculé comme la norme euclidienne de l'accélération 3D :

```
Accel_Mod_Hand = sqrt(ax² + ay² + az²)
```

### 3. Normalisation temporelle (20 bins)
Chaque signal est découpé en **20 intervalles de durée égale** (5 % de la session chacun), et la moyenne est calculée par intervalle. Cela permet de comparer des essais de durées différentes sur un axe commun 0–100 %.

> **Note :** le script utilise 20 bins (axe x : 0.5 à 99.5 %). D'autres scripts du projet utilisent 10 bins — les résultats ne sont pas directement comparables.

### 4. Code couleur et épaisseur de ligne

| Condition | Couleur |
|---|---|
| Competition | Rouge |
| Extrait IE | Vert (clair=A, normal=B, foncé=C) |
| Extrait PS | Bleu (clair=A, normal=B, foncé=C) |

L'épaisseur de ligne varie selon le take : A = 1.8 pt, B = 1.3 pt, C = 0.9 pt.

### 5. Ordre de passage
Chaque participant joue les conditions dans un ordre différent. La structure `first_expr` encode quel extrait (IE ou PS) a été joué en premier. La Competition est toujours en premier et est annotée `C` en fin de courbe. Les extraits sont numérotés de 1 (premier passage) à 6 (dernier passage).

### 6. Sortie
Une figure par participant est sauvegardée sous `Fig_AccelHand_<sujet>.png` dans le dossier `path_save`.

---

## Abréviations

| Abréviation | Signification dans ce script |
|---|---|
| IE | Interprétation Expressive — condition de jeu expressif |
| PS | Performance Standard — condition de jeu neutre/standard |
| A, B, C | Takes successifs d'un même extrait (répétitions) |
| Competition | Condition de jeu en situation de compétition (toujours en premier) |
| rHand | Main droite (right Hand) — segment 11 dans la hiérarchie XSens |
| Accel_Mod | Accélération Modulée = norme de l'accélération 3D |
| Bin | Intervalle temporel normalisé (1/20e de la durée totale) |
| P13 | Participant exclu de la condition Competition (données manquantes) |

---

## Fichiers d'entrée / sortie

| Fichier | Rôle |
|---|---|
| `XSens_MIDI_cut_P01_P08 1.mat` | Données XSens P01–P08 |
| `XSens_MIDI_cut_P09_P14 1.mat` | Données XSens P09–P14 |
| `Fig_AccelHand_<sujet>.png` | Figures de sortie (une par participant) |

---

## Lien avec d'autres scripts

- **`ComputeJointAngles.m`** : utilise la même logique de chargement et fusion des deux `.mat` (P01–P08 + P09–P14). La structure `XSens_cur` est construite de la même façon.
- **Scripts PLSR du dataset Ergo** : utilisent également une normalisation temporelle par bins, mais sur 10 bins et sur un dataset différent (Robin2/Ergo, 23 sujets, clavier ergonomique).
