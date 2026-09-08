# README — `Intervalle_Reaper.m`

## Vue d'ensemble

Script de **génération de projets Reaper (`.RPP`)** pour l'analyse de passages musicaux ciblés. Pour chaque participant, il découpe le fichier audio (WAV microphone) et le fichier MIDI sur des intervalles d'intérêt prédéfinis, puis génère un projet Reaper autonome contenant les deux pistes synchronisées.

---

## Dataset

| Paramètre | Valeur |
|---|---|
| Participants | P01–P14 (13 sujets, P03 absent) |
| Conditions | Competition ou Extrait (IE/PS/take) selon le meilleur essai de chaque sujet |
| Source audio | Fichier microphone `01-mic-*.wav` |
| Synchronisation | Signal TTL (impulsion électrique enregistrée en parallèle) |

---

## Ce que fait le script

### 1. Sélection du meilleur essai (`BEST`)
La table `BEST` définit, pour chaque participant, l'essai de référence (condition + expression + take). Ce choix a été fait manuellement à partir des résultats visuels de `Plot_Expression.m` (essai le plus long ou le plus représentatif).

### 2. Définition des passages d'intérêt (`QUERIES`)
Pour chaque participant, `QUERIES` contient un ou plusieurs intervalles `[t_début, t_fin]` en secondes (t=0 = première note MIDI). Ces passages ont été identifiés manuellement à l'écoute ou par analyse.

> P14 n'a aucun passage défini (`QUERIES.P14 = {}`) — le script le saute automatiquement.

### 3. Détection de l'offset TTL (`detect_ttl_offset`)
Le signal TTL est un déclencheur électrique enregistré simultanément avec l'audio pour marquer le début de l'acquisition XSens. La fonction lit le fichier `*_TTL.wav`, détecte le premier échantillon dépassant 50 % de l'amplitude maximale, et retourne le temps correspondant en secondes. Cet offset permet d'aligner le WAV audio sur la timeline MIDI.

### 4. Découpage audio
Le WAV microphone est découpé entre `t_start + ttl_offset` et `t_end + ttl_offset` (en tenant compte du décalage TTL) et sauvegardé dans un sous-dossier par passage.

### 5. Découpage MIDI (`trim_midi`)
Le fichier MIDI est parsé en binaire (format standard MIDI type 0 ou 1). Les événements compris entre `t_start` et `t_end` sont extraits et réécrits dans un nouveau fichier `.mid`. Les événements de tempo sont toujours conservés (même hors fenêtre) pour garantir la bonne interprétation des ticks. Le premier événement est forcé à delta=0.

### 6. Calcul du silence initial MIDI (`get_midi_first_note_delay`)
Après découpage, le MIDI peut commencer par un silence si la première note du passage n'est pas exactement à t_start. Cette fonction lit le MIDI découpé et retourne le délai en secondes avant la première `note-on` (vélocité > 0).

### 7. Génération du projet Reaper (`write_rpp_local`)
Génère un fichier `.RPP` autonome avec chemins relatifs (le RPP, le WAV et le MIDI sont dans le même dossier). La piste audio est décalée de `midi_offset` secondes via le champ `POSITION` du RPP, de sorte que la première note MIDI et le début de l'audio coïncident à l'ouverture dans Reaper.

---

## Fonctions locales

| Fonction | Rôle |
|---|---|
| `detect_ttl_offset` | Détecte le temps du premier trigger TTL dans le WAV |
| `trim_midi` | Découpe un fichier MIDI entre deux instants (secondes) |
| `get_midi_first_note_delay` | Retourne le délai avant la première note-on dans un MIDI |
| `write_rpp_local` | Génère le RPP avec chemins relatifs et compensation du silence MIDI |
| `read_midi_ticks` | Lit le nombre de ticks par noire dans l'en-tête MIDI |
| `read_varlen` | Lit un entier à longueur variable (format MIDI standard) |
| `write_varlen` | Encode un entier en longueur variable MIDI |
| `event_data_len` | Retourne la longueur des données après un status byte MIDI |
| `basename` | Extrait le nom de fichier depuis un chemin complet |
| `guid_new` | Génère un GUID aléatoire au format Reaper |
| `write_rpp` *(inactive)* | Ancienne version du générateur RPP — chemins absolus, pas de découpage préalable des fichiers. Conservée à titre de référence mais non appelée dans la boucle principale. |

---

## Différence entre `write_rpp_local` et `write_rpp`

| | `write_rpp_local` *(active)* | `write_rpp` *(legacy, inactive)* |
|---|---|---|
| Fichiers référencés | WAV et MIDI pré-découpés, chemins relatifs | WAV et MIDI bruts, chemins absolus |
| Synchronisation audio | `POSITION = midi_offset` (décale la piste audio) | `SOFFS = t_start + ttl_offset` (seek dans le WAV brut) |
| Portabilité du RPP | Oui (dossier autonome déplaçable) | Non (dépend des chemins absolus de la machine) |

---

## Abréviations

| Abréviation | Signification dans ce script |
|---|---|
| RPP | Format de projet Reaper (fichier texte décrivant pistes, items, sources) |
| TTL | Transistor-Transistor Logic — signal électrique impulsionnel utilisé comme déclencheur de synchronisation |
| MIDI | Musical Instrument Digital Interface — fichier d'événements musicaux (notes, tempo, etc.) |
| WAV | Format audio non compressé |
| SOFFS | Source Offset — position de départ dans le fichier source (champ Reaper) |
| POSITION | Position de l'item sur la timeline Reaper (en secondes depuis le début) |
| IE | Interprétation Expressive |
| PS | Performance Standard |
| take | Répétition d'un extrait (A, B ou C) |
| t_start / t_end | Bornes temporelles du passage d'intérêt, en secondes depuis la première note MIDI (t=0) |
| delta | Délai relatif entre deux événements MIDI successifs (en ticks) |
| tpb | Ticks Per Beat — résolution temporelle du fichier MIDI |
| note-on | Événement MIDI indiquant l'appui sur une touche (status 0x90, vélocité > 0) |

---

## Fichiers d'entrée / sortie

| Fichier | Rôle |
|---|---|
| `<sujet>/Audio/<condition>/01-mic-*.wav` | Enregistrement microphone brut |
| `<sujet>/Audio/<condition>/*_TTL.wav` | Signal TTL de synchronisation |
| `<sujet>/MIDI/<condition>.mid` | Fichier MIDI brut |
| `Reaper_Passages/<sujet>/passage_N_Xs_Ys/` | Dossier de sortie par passage (WAV + MIDI découpés + RPP) |

---

## Problème connu

La compensation du silence MIDI via `midi_offset` (champ `POSITION` de l'item audio) ne résout pas toujours la désynchronisation entre audio et MIDI dans les projets générés. Le silence initial varie selon les participants et n'est pas toujours correctement détecté par `get_midi_first_note_delay`.

---

## Lien avec d'autres scripts

- **`Plot_Expression.m`** : les choix dans la table `BEST` sont basés sur les figures produites par ce script.
