%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  GÉNÉRATION DES PROJETS REAPER (.RPP) PAR PARTICIPANT             %%%%
%%%%  Audio (mic WAV) + MIDI découpés sur le passage d'intérêt         %%%%
%%%%  Synchronisation via signal TTL                                    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% -----------------------------------------------------------------------
%  CHEMINS
%  -----------------------------------------------------------------------
ROOT = 'Z:\Piano_Expression_2025\Raw Data';
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
OUT = fullfile(script_dir, 'Reaper_Passages');
if ~exist(OUT, 'dir'), mkdir(OUT); end

%% -----------------------------------------------------------------------
%  ESSAI LE PLUS LONG PAR PARTICIPANT (depuis MATLAB Plot_Expression)
%  Format : {sujet, condition, expression, take}
%  -----------------------------------------------------------------------
BEST = {
    'P01', 'Competition', '',   '';
    'P02', 'Competition', '',   '';
    'P04', 'Extrait',     'IE', 'B';
    'P05', 'Competition', '',   '';
    'P06', 'Extrait',     'IE', 'C';
    'P07', 'Competition', '',   '';
    'P08', 'Competition', '',   '';
    'P09', 'Extrait',     'IE', 'C';
    'P10', 'Competition', '',   '';
    'P11', 'Competition', '',   '';
    'P12', 'Extrait',     'IE', 'C';
    'P13', 'Extrait',     'IE', 'B';
    'P14', 'Competition', '',   '';
};

%% -----------------------------------------------------------------------
%  PASSAGES D'INTÉRÊT PAR PARTICIPANT
%  Format : cell array de [t_debut_s, t_fin_s] en secondes MIDI (t=0 = 1ère note)
%  -----------------------------------------------------------------------
QUERIES = struct();
QUERIES.P01 = {[158, 185]};
QUERIES.P02 = {[80, 84]; [175, 179]};
QUERIES.P04 = {[62, 66]; [215, 219]};
QUERIES.P05 = {[121, 128]};
QUERIES.P06 = {[80, 84]};
QUERIES.P07 = {[118, 149]; [170, 174]};
QUERIES.P08 = {[40, 43]};
QUERIES.P09 = {[34, 37]};
QUERIES.P10 = {[143, 147]};
QUERIES.P11 = {[303, 311]};
QUERIES.P12 = {[17, 22]; [84, 89]; [206, 210]; [245, 250]; [343, 348]};
QUERIES.P13 = {[23, 29]; [217, 223]; [371, 377]};
QUERIES.P14 = {};

%% -----------------------------------------------------------------------
%  BOUCLE PRINCIPALE
%  -----------------------------------------------------------------------
fprintf('\n=== Génération des projets Reaper ===\n\n');

for iS = 1:size(BEST, 1)
    subj = BEST{iS, 1};
    cond = BEST{iS, 2};
    expr = BEST{iS, 3};
    take = BEST{iS, 4};

    queries = QUERIES.(subj);
    if isempty(queries)
        fprintf('%s : aucun passage défini — skip\n', subj);
        continue;
    end

    fprintf('--- %s (%s %s %s) ---\n', subj, cond, expr, take);

    % --- Construire les chemins Audio et MIDI ---
    if strcmp(cond, 'Competition')
        audio_folder = fullfile(ROOT, subj, 'Audio', 'Competition');
        ttl_name     = 'Competition_TTL.wav';
        midi_name    = 'Competition.mid';
    else
        audio_folder = fullfile(ROOT, subj, 'Audio', cond, expr, ...
                                sprintf('%s_%s_%s', cond, expr, take));
        ttl_name     = sprintf('%s_%s_%s_TTL.wav', cond, expr, take);
        midi_name    = sprintf('%s_%s_%s.mid', cond, expr, take);
    end

    midi_path = fullfile(ROOT, subj, 'MIDI', midi_name);
    ttl_path  = fullfile(audio_folder, ttl_name);

    % Trouver le fichier mic (01-mic-*.wav, taille > 1 Mo)
    mic_files = dir(fullfile(audio_folder, '01-mic-*.wav'));
    % Filtrer les fichiers .reapeaks et garder uniquement les vrais WAV (>1Mo)
    mic_files = mic_files(~contains({mic_files.name}, '.reapeaks') & [mic_files.bytes] > 1e6);
    if isempty(mic_files)
        fprintf('  [ERREUR] Fichier mic introuvable dans %s\n', audio_folder);
        continue;
    end
    mic_path = fullfile(audio_folder, mic_files(1).name);

    % Vérifications
    if ~exist(mic_path, 'file')
        fprintf('  [ERREUR] Mic introuvable : %s\n', mic_path); continue;
    end
    if ~exist(ttl_path, 'file')
        fprintf('  [ERREUR] TTL introuvable : %s\n', ttl_path); continue;
    end
    if ~exist(midi_path, 'file')
        fprintf('  [ERREUR] MIDI introuvable : %s\n', midi_path); continue;
    end

    fprintf('  Mic  : %s\n', mic_files(1).name);
    fprintf('  TTL  : %s\n', ttl_name);
    fprintf('  MIDI : %s\n', midi_name);

    % --- Détecter le décalage TTL ---
    ttl_offset = detect_ttl_offset(ttl_path);
    fprintf('  Offset TTL : %.3f s\n', ttl_offset);

    % --- Créer le dossier de sortie ---
    subj_out = fullfile(OUT, subj);
    if ~exist(subj_out, 'dir'), mkdir(subj_out); end

    % --- Générer un RPP par passage ---
    for iP = 1:length(queries)
        t_start = queries{iP}(1);
        t_end   = queries{iP}(2);
        duration = t_end - t_start;

        fprintf('  Passage %d : %.1f s - %.1f s\n', iP, t_start, t_end);

        out_rpp = fullfile(subj_out, ...
            sprintf('%s_passage_%d_%ds_%ds.RPP', subj, iP, round(t_start), round(t_end)));

        try
            % Créer dossier du passage
            passage_folder = fullfile(subj_out, ...
                sprintf('passage_%d_%ds_%ds', iP, round(t_start), round(t_end)));
            if ~exist(passage_folder, 'dir'), mkdir(passage_folder); end

            % Noms des fichiers locaux
            out_wav  = fullfile(passage_folder, sprintf('%s_passage_%d_audio.wav', subj, iP));
            out_mid  = fullfile(passage_folder, sprintf('%s_passage_%d.mid', subj, iP));
            out_rpp2 = fullfile(passage_folder, ...
                sprintf('%s_passage_%d_%ds_%ds.RPP', subj, iP, round(t_start), round(t_end)));

            % 1) Découper le WAV
            info_wav = audioinfo(mic_path);
            sr       = info_wav.SampleRate;
            s1 = max(1,    round((t_start + ttl_offset) * sr) + 1);
            s2 = min(info_wav.TotalSamples, round((t_end + ttl_offset) * sr));
            [y, ~] = audioread(mic_path, [s1, s2]);
            audiowrite(out_wav, y, sr, 'BitsPerSample', info_wav.BitsPerSample);
            fprintf('    WAV découpé : %.1f s\n', (s2-s1)/sr);

            % 2) Découper le MIDI
            trim_midi(midi_path, out_mid, t_start, t_end);
            fprintf('    MIDI découpé\n');

            % 3) Calculer le silence initial du MIDI découpé
            midi_offset = get_midi_first_note_delay(out_mid);
            fprintf('    Silence initial MIDI : %.3f s\n', midi_offset);

            % 4) Générer le RPP avec chemins relatifs
            %    L'audio est avancé de midi_offset pour compenser le silence MIDI
            write_rpp_local(out_rpp2, out_wav, out_mid, subj, iP, t_end - t_start, midi_offset);
            fprintf('    -> %s\n', sprintf('%s_passage_%d_%ds_%ds.RPP', subj, iP, round(t_start), round(t_end)));
        catch ME
            fprintf('    [ERREUR] %s\n    %s\n', ME.message, ME.stack(1).name);
        end
    end
    fprintf('\n');
end

fprintf('=== TERMINÉ ===\n');
fprintf('Projets dans : %s\n', OUT);

%% -----------------------------------------------------------------------
%  FONCTIONS LOCALES
%  -----------------------------------------------------------------------

function offset = detect_ttl_offset(ttl_path)
%DETECT_TTL_OFFSET Retourne le temps (s) du premier trigger TTL dans le WAV.
    try
        info    = audioinfo(ttl_path);
        sr      = info.SampleRate;
        data    = audioread(ttl_path);
        data    = data(:,1);   % mono ou canal 1

        threshold = max(abs(data)) * 0.5;
        idx       = find(abs(data) > threshold, 1, 'first');
        if isempty(idx)
            offset = 0;
        else
            offset = (idx - 1) / sr;
        end
    catch
        offset = 0;
        fprintf('  [WARN] TTL non lisible, offset=0\n');
    end
end


function write_rpp_local(out_path, wav_path, midi_path, subj, passage_idx, duration, midi_offset)
%WRITE_RPP_LOCAL Génère un RPP avec chemins relatifs vers les fichiers locaux découpés.
    guid_twav  = guid_new();
    guid_tmidi = guid_new();
    guid_iwav  = guid_new();
    guid_imidi = guid_new();

    ticks_per_beat = read_midi_ticks(midi_path);

    % Chemins relatifs (le RPP est dans le même dossier que les fichiers)
    wav_rel  = ['.' filesep basename(wav_path)];
    midi_rel = ['.' filesep basename(midi_path)];

    fid = fopen(out_path, 'w', 'n', 'UTF-8');
    if fid == -1, error('Impossible de créer : %s', out_path); end

    fprintf(fid, '<REAPER_PROJECT 0.1 "7.40/win64" 0\r\n');
    fprintf(fid, '  RIPPLE 0 0\r\n');
    fprintf(fid, '  LOOP 0\r\n');
    fprintf(fid, '  TEMPO 120 4 4 0\r\n');
    fprintf(fid, '  PLAYRATE 1 0 0.25 4\r\n');
    fprintf(fid, '  SELECTION 0 0\r\n');
    fprintf(fid, '  SELECTION2 0 0\r\n');

    % Track Audio
    fprintf(fid, '  <TRACK %s\r\n', guid_twav);
    fprintf(fid, '    NAME "Audio - %s passage %d"\r\n', subj, passage_idx);
    fprintf(fid, '    PEAKCOL 16576\r\n');
    fprintf(fid, '    BEAT -1\r\n');
    fprintf(fid, '    VOLPAN 1 0 -1 -1 1\r\n');
    fprintf(fid, '    MUTESOLO 0 0 0\r\n');
    fprintf(fid, '    ISBUS 0 0\r\n');
    fprintf(fid, '    NCHAN 2\r\n');
    fprintf(fid, '    FX 1\r\n');
    fprintf(fid, '    TRACKID %s\r\n', guid_twav);
    fprintf(fid, '    MIDIOUT -1\r\n');
    fprintf(fid, '    MAINSEND 1 0\r\n');
    fprintf(fid, '    <ITEM\r\n');
    fprintf(fid, '      POSITION %.6f\r\n', midi_offset);
    fprintf(fid, '      SNAPOFFS 0\r\n');
    fprintf(fid, '      LENGTH %.6f\r\n', duration);
    fprintf(fid, '      LOOP 0\r\n');
    fprintf(fid, '      FADEIN 1 0.01 0 1 0 0 0\r\n');
    fprintf(fid, '      FADEOUT 1 0.01 0 1 0 0 0\r\n');
    fprintf(fid, '      MUTE 0 0\r\n');
    fprintf(fid, '      IGUID %s\r\n', guid_iwav);
    fprintf(fid, '      IID 1\r\n');
    fprintf(fid, '      VOLPAN 1 0 1 -1\r\n');
    fprintf(fid, '      SOFFS 0\r\n');
    fprintf(fid, '      PLAYRATE 1 1 0 -1 0 0.0025\r\n');
    fprintf(fid, '      CHANMODE 0\r\n');
    fprintf(fid, '      GUID %s\r\n', guid_iwav);
    fprintf(fid, '      <SOURCE WAVE\r\n');
    fprintf(fid, '        FILE "%s"\r\n', wav_rel);
    fprintf(fid, '      >\r\n');
    fprintf(fid, '    >\r\n');
    fprintf(fid, '  >\r\n');

    % Track MIDI
    fprintf(fid, '  <TRACK %s\r\n', guid_tmidi);
    fprintf(fid, '    NAME "MIDI - %s passage %d"\r\n', subj, passage_idx);
    fprintf(fid, '    PEAKCOL 33895697\r\n');
    fprintf(fid, '    BEAT -1\r\n');
    fprintf(fid, '    VOLPAN 1 0 -1 -1 1\r\n');
    fprintf(fid, '    MUTESOLO 0 0 0\r\n');
    fprintf(fid, '    ISBUS 0 0\r\n');
    fprintf(fid, '    NCHAN 2\r\n');
    fprintf(fid, '    FX 1\r\n');
    fprintf(fid, '    TRACKID %s\r\n', guid_tmidi);
    fprintf(fid, '    MIDIOUT -1\r\n');
    fprintf(fid, '    MAINSEND 1 0\r\n');
    fprintf(fid, '    <ITEM\r\n');
    fprintf(fid, '      POSITION 0\r\n');
    fprintf(fid, '      SNAPOFFS 0\r\n');
    fprintf(fid, '      LENGTH %.6f\r\n', duration);
    fprintf(fid, '      LOOP 0\r\n');
    fprintf(fid, '      FADEIN 1 0 0 1 0 0 0\r\n');
    fprintf(fid, '      FADEOUT 1 0 0 1 0 0 0\r\n');
    fprintf(fid, '      MUTE 0 0\r\n');
    fprintf(fid, '      IGUID %s\r\n', guid_imidi);
    fprintf(fid, '      IID 2\r\n');
    fprintf(fid, '      VOLPAN 1 0 1 -1\r\n');
    fprintf(fid, '      SOFFS 0 0\r\n');
    fprintf(fid, '      PLAYRATE 1 1 0 -1 0 0.0025\r\n');
    fprintf(fid, '      CHANMODE 0\r\n');
    fprintf(fid, '      GUID %s\r\n', guid_imidi);
    fprintf(fid, '      <SOURCE MIDI\r\n');
    fprintf(fid, '        FILE "%s"\r\n', midi_rel);
    fprintf(fid, '        IGNTEMPO 0 120 4 4\r\n');
    fprintf(fid, '      >\r\n');
    fprintf(fid, '    >\r\n');
    fprintf(fid, '  >\r\n');
    fprintf(fid, '>\r\n');
    fclose(fid);
end


function delay = get_midi_first_note_delay(midi_path)
%GET_MIDI_FIRST_NOTE_DELAY Retourne le délai (s) avant la première note-on dans le MIDI.
    delay = 0;
    try
        fid = fopen(midi_path, 'rb', 'b');
        if fid == -1, return; end
        % Lire header
        fread(fid, 4, 'uint8');   % MThd
        fread(fid, 1, 'uint32');  % len
        fread(fid, 1, 'uint16');  % format
        n_tracks = fread(fid, 1, 'uint16');
        tpb = fread(fid, 1, 'uint16');

        tempo = 500000;
        first_note_ticks = Inf;

        for iT = 1:n_tracks
            tag = fread(fid, 4, 'uint8=>char')';
            if ~strcmp(tag, 'MTrk'), break; end
            len  = fread(fid, 1, 'uint32');
            data = fread(fid, len, 'uint8')';

            pos = 1; cum = 0; last_s = 0;
            while pos <= length(data)
                [dt, pos] = read_varlen(data, pos);
                cum = cum + dt;
                if pos > length(data), break; end
                s = data(pos);
                if s == 0xFF
                    mt = data(pos+1);
                    [ml, p2] = read_varlen(data, pos+2);
                    if mt == 0x51 && ml == 3
                        tempo = data(p2)*65536 + data(p2+1)*256 + data(p2+2);
                    end
                    pos = p2 + ml;
                elseif s == 0xF0 || s == 0xF7
                    [sl, p2] = read_varlen(data, pos+1);
                    pos = p2 + sl;
                else
                    if s < 0x80
                        type = bitand(last_s, 0xF0);
                        n = event_data_len(last_s, data, pos);
                        if type == 0x90 && n >= 2 && data(pos+1) > 0
                            first_note_ticks = min(first_note_ticks, cum);
                        end
                        pos = pos + n;
                    else
                        last_s = s;
                        type = bitand(s, 0xF0);
                        n = event_data_len(s, data, pos+1);
                        if type == 0x90 && n >= 2 && pos+2 <= length(data) && data(pos+2) > 0
                            first_note_ticks = min(first_note_ticks, cum);
                        end
                        pos = pos + n + 1;
                    end
                end
            end
        end
        fclose(fid);

        if ~isinf(first_note_ticks)
            delay = first_note_ticks / tpb * (tempo / 1e6);
        end
    catch
        delay = 0;
    end
end


function trim_midi(in_path, out_path, t_start, t_end)
%TRIM_MIDI Découpe un fichier MIDI entre t_start et t_end (secondes).
%  Lecture binaire brute du format MIDI standard (type 0 ou 1).
%  Les événements entre t_start et t_end sont extraits et réécrits.

    fid = fopen(in_path, 'rb', 'b');   % big-endian
    if fid == -1, error('MIDI introuvable : %s', in_path); end

    % Lire header MThd
    hdr_tag  = fread(fid, 4, 'uint8=>char')';
    if ~strcmp(hdr_tag, 'MThd'), fclose(fid); error('Pas un fichier MIDI valide'); end
    hdr_len  = fread(fid, 1, 'uint32');   % toujours 6
    fmt      = fread(fid, 1, 'uint16');
    n_tracks = fread(fid, 1, 'uint16');
    tpb      = fread(fid, 1, 'uint16');   % ticks per beat

    % Lire toutes les tracks
    tracks_raw = cell(1, n_tracks);
    for iT = 1:n_tracks
        tag = fread(fid, 4, 'uint8=>char')';
        if ~strcmp(tag, 'MTrk'), break; end
        len  = fread(fid, 1, 'uint32');
        data = fread(fid, len, 'uint8')';
        tracks_raw{iT} = data;
    end
    fclose(fid);

    % Trouver le tempo (defaut 500000 µs/beat = 120 BPM)
    tempo = 500000;
    for iT = 1:n_tracks
        data = tracks_raw{iT};
        pos  = 1;
        cum_ticks = 0;
        while pos <= length(data)
            [dt, pos] = read_varlen(data, pos);
            cum_ticks = cum_ticks + dt;
            if pos > length(data), break; end
            status = data(pos);
            if status == 0xFF  % meta
                meta_type = data(pos+1);
                [meta_len, pos2] = read_varlen(data, pos+2);
                if meta_type == 0x51 && meta_len == 3
                    tempo = data(pos2)*65536 + data(pos2+1)*256 + data(pos2+2);
                end
                pos = pos2 + meta_len;
            else
                pos = pos + 1 + event_data_len(status, data, pos+1);
            end
        end
    end

    ticks_per_sec = 1e6 / tempo * tpb;
    tick_start    = round(t_start * ticks_per_sec);
    tick_end      = round(t_end   * ticks_per_sec);

    % Extraire les événements dans la fenêtre pour chaque track
    new_tracks = cell(1, n_tracks);
    for iT = 1:n_tracks
        data      = tracks_raw{iT};
        pos       = 1;
        cum_ticks = 0;
        out_bytes = uint8([]);
        prev_out_tick = 0;
        last_status = 0;
        first_event = true;   % forcer delta=0 pour le premier événement

        while pos <= length(data)
            [dt, pos] = read_varlen(data, pos);
            cum_ticks = cum_ticks + dt;
            if pos > length(data), break; end

            status = data(pos);
            is_meta = (status == 0xFF);
            is_sysex = (status == 0xF0 || status == 0xF7);

            if is_meta
                meta_type = data(pos+1);
                [meta_len, pos2] = read_varlen(data, pos+2);
                evt_bytes = data(pos:pos2+meta_len-1);
                pos = pos2 + meta_len;
            elseif is_sysex
                [sx_len, pos2] = read_varlen(data, pos+1);
                evt_bytes = data(pos:pos2+sx_len-1);
                pos = pos2 + sx_len;
            else
                if status < 0x80
                    % running status
                    n = event_data_len(last_status, data, pos);
                    evt_bytes = [last_status, data(pos:pos+n-1)];
                    pos = pos + n;
                else
                    last_status = status;
                    n = event_data_len(status, data, pos+1);
                    evt_bytes = data(pos:pos+n);
                    pos = pos + n + 1;
                end
            end

            % Inclure l'événement si dans la fenêtre
            in_window = (cum_ticks >= tick_start && cum_ticks <= tick_end);
            is_tempo_meta = (is_meta && length(evt_bytes) >= 2 && evt_bytes(2) == 0x51);
            is_eot = (is_meta && length(evt_bytes) >= 2 && evt_bytes(2) == 0x2F);

            if in_window || is_tempo_meta
                out_tick = max(0, cum_ticks - tick_start);
                if first_event
                    delta = 0;   % forcer le premier événement à delta=0
                    first_event = false;
                else
                    delta = out_tick - prev_out_tick;
                end
                out_bytes = [out_bytes, write_varlen(uint32(delta)), evt_bytes];
                prev_out_tick = out_tick;
            end

            if is_eot && cum_ticks >= tick_end
                % Ajouter EOT final
                out_bytes = [out_bytes, write_varlen(uint32(0)), uint8([0xFF, 0x2F, 0x00])];
                break;
            end
        end

        if isempty(out_bytes) || ~(out_bytes(end-2) == 0xFF && out_bytes(end-1) == 0x2F)
            out_bytes = [out_bytes, write_varlen(uint32(0)), uint8([0xFF, 0x2F, 0x00])];
        end

        new_tracks{iT} = out_bytes;
    end

    % Écrire le fichier MIDI de sortie
    fid = fopen(out_path, 'wb', 'b');
    if fid == -1, error('Impossible de créer MIDI : %s', out_path); end

    % Header
    fwrite(fid, 'MThd', 'uint8');
    fwrite(fid, uint32(6), 'uint32');
    fwrite(fid, uint16(fmt), 'uint16');
    fwrite(fid, uint16(n_tracks), 'uint16');
    fwrite(fid, uint16(tpb), 'uint16');

    % Tracks
    for iT = 1:n_tracks
        fwrite(fid, 'MTrk', 'uint8');
        fwrite(fid, uint32(length(new_tracks{iT})), 'uint32');
        fwrite(fid, new_tracks{iT}, 'uint8');
    end
    fclose(fid);
end


function [val, pos] = read_varlen(data, pos)
%READ_VARLEN Lit un entier à longueur variable (format MIDI).
    val = 0;
    for k = 1:4
        b   = data(pos); pos = pos + 1;
        val = val * 128 + bitand(b, 127);
        if bitand(b, 128) == 0, break; end
    end
end


function bytes = write_varlen(val)
%WRITE_VARLEN Encode un entier en longueur variable MIDI.
    bytes = uint8(bitand(val, 127));
    val   = bitshift(val, -7);
    while val > 0
        bytes = [uint8(bitor(bitand(val, 127), 128)), bytes];
        val   = bitshift(val, -7);
    end
end


function n = event_data_len(status, data, pos)
%EVENT_DATA_LEN Retourne le nombre d'octets de données après le status byte.
    type = bitand(status, 0xF0);
    if type == 0x80 || type == 0x90 || type == 0xA0 || type == 0xB0 || type == 0xE0
        n = 2;
    elseif type == 0xC0 || type == 0xD0
        n = 1;
    elseif status == 0xF2
        n = 2;
    elseif status == 0xF1 || status == 0xF3
        n = 1;
    else
        n = 0;
    end
end


function b = basename(p)
%BASENAME Retourne le nom de fichier depuis un chemin complet.
    [~, name, ext] = fileparts(p);
    b = [name ext];
end


function write_rpp(out_path, wav_path, midi_path, ttl_offset, t_start, t_end, subj, passage_idx)
%WRITE_RPP Génère un fichier .RPP Reaper avec 2 tracks : audio + MIDI.
%
%  Le WAV est référencé tel quel (chemin complet), le MIDI aussi.
%  L'audio commence à audio_soffs = t_start + ttl_offset dans le WAV brut.
%  Le MIDI commence à t_start dans le fichier MIDI original.

    duration    = t_end - t_start;
    audio_soffs = t_start + ttl_offset;   % position dans le WAV brut

    % GUIDs uniques (format Reaper)
    guid_twav  = guid_new();
    guid_tmidi = guid_new();
    guid_iwav  = guid_new();
    guid_imidi = guid_new();

    % Lire le tempo MIDI pour HASDATA
    ticks_per_beat = read_midi_ticks(midi_path);

    fid = fopen(out_path, 'w', 'n', 'UTF-8');
    if fid == -1
        error('Impossible de créer le fichier : %s', out_path);
    end

    fprintf(fid, '<REAPER_PROJECT 0.1 "7.40/win64" 0\r\n');
    fprintf(fid, '  RIPPLE 0 0\r\n');
    fprintf(fid, '  GROUPOVERRIDE 0 0 0\r\n');
    fprintf(fid, '  AUTOXFADE 129\r\n');
    fprintf(fid, '  LOOP 0\r\n');
    fprintf(fid, '  TEMPO 120 4 4 0\r\n');
    fprintf(fid, '  PLAYRATE 1 0 0.25 4\r\n');
    fprintf(fid, '  SELECTION 0 0\r\n');
    fprintf(fid, '  SELECTION2 0 0\r\n');

    % --- Track Audio ---
    fprintf(fid, '  <TRACK %s\r\n', guid_twav);
    fprintf(fid, '    NAME "Audio - %s passage %d"\r\n', subj, passage_idx);
    fprintf(fid, '    PEAKCOL 16576\r\n');
    fprintf(fid, '    BEAT -1\r\n');
    fprintf(fid, '    AUTOMODE 0\r\n');
    fprintf(fid, '    VOLPAN 1 0 -1 -1 1\r\n');
    fprintf(fid, '    MUTESOLO 0 0 0\r\n');
    fprintf(fid, '    IPHASE 0\r\n');
    fprintf(fid, '    PLAYOFFS 0 1\r\n');
    fprintf(fid, '    ISBUS 0 0\r\n');
    fprintf(fid, '    SHOWINMIX 1 0.6667 0.5 1 0.5 0 0 0\r\n');
    fprintf(fid, '    SEL 0\r\n');
    fprintf(fid, '    REC 0 0 0 0 0 0 0 0\r\n');
    fprintf(fid, '    VU 2\r\n');
    fprintf(fid, '    TRACKHEIGHT 0 0 0 0 0 0 0\r\n');
    fprintf(fid, '    INQ 0 0 0 0.5 100 0 0 100\r\n');
    fprintf(fid, '    NCHAN 2\r\n');
    fprintf(fid, '    FX 1\r\n');
    fprintf(fid, '    TRACKID %s\r\n', guid_twav);
    fprintf(fid, '    PERF 0\r\n');
    fprintf(fid, '    MIDIOUT -1\r\n');
    fprintf(fid, '    MAINSEND 1 0\r\n');
    fprintf(fid, '    <ITEM\r\n');
    fprintf(fid, '      POSITION 0\r\n');
    fprintf(fid, '      SNAPOFFS 0\r\n');
    fprintf(fid, '      LENGTH %.6f\r\n', duration);
    fprintf(fid, '      LOOP 0\r\n');
    fprintf(fid, '      ALLTAKES 0\r\n');
    fprintf(fid, '      FADEIN 1 0.01 0 1 0 0 0\r\n');
    fprintf(fid, '      FADEOUT 1 0.01 0 1 0 0 0\r\n');
    fprintf(fid, '      MUTE 0 0\r\n');
    fprintf(fid, '      SEL 0\r\n');
    fprintf(fid, '      IGUID %s\r\n', guid_iwav);
    fprintf(fid, '      IID 1\r\n');
    fprintf(fid, '      NAME "%s"\r\n', wav_path);
    fprintf(fid, '      VOLPAN 1 0 1 -1\r\n');
    fprintf(fid, '      SOFFS %.6f\r\n', audio_soffs);
    fprintf(fid, '      PLAYRATE 1 1 0 -1 0 0.0025\r\n');
    fprintf(fid, '      CHANMODE 0\r\n');
    fprintf(fid, '      GUID %s\r\n', guid_iwav);
    fprintf(fid, '      <SOURCE WAVE\r\n');
    fprintf(fid, '        FILE "%s"\r\n', wav_path);
    fprintf(fid, '      >\r\n');
    fprintf(fid, '    >\r\n');
    fprintf(fid, '  >\r\n');

    % --- Track MIDI ---
    fprintf(fid, '  <TRACK %s\r\n', guid_tmidi);
    fprintf(fid, '    NAME "MIDI - %s passage %d"\r\n', subj, passage_idx);
    fprintf(fid, '    PEAKCOL 33895697\r\n');
    fprintf(fid, '    BEAT -1\r\n');
    fprintf(fid, '    AUTOMODE 0\r\n');
    fprintf(fid, '    VOLPAN 1 0 -1 -1 1\r\n');
    fprintf(fid, '    MUTESOLO 0 0 0\r\n');
    fprintf(fid, '    IPHASE 0\r\n');
    fprintf(fid, '    PLAYOFFS 0 1\r\n');
    fprintf(fid, '    ISBUS 0 0\r\n');
    fprintf(fid, '    SHOWINMIX 1 0.6667 0.5 1 0.5 0 0 0\r\n');
    fprintf(fid, '    SEL 0\r\n');
    fprintf(fid, '    REC 0 0 0 0 0 0 0 0\r\n');
    fprintf(fid, '    VU 2\r\n');
    fprintf(fid, '    TRACKHEIGHT 0 0 0 0 0 0 0\r\n');
    fprintf(fid, '    INQ 0 0 0 0.5 100 0 0 100\r\n');
    fprintf(fid, '    NCHAN 2\r\n');
    fprintf(fid, '    FX 1\r\n');
    fprintf(fid, '    TRACKID %s\r\n', guid_tmidi);
    fprintf(fid, '    PERF 0\r\n');
    fprintf(fid, '    MIDIOUT -1\r\n');
    fprintf(fid, '    MAINSEND 1 0\r\n');
    fprintf(fid, '    <ITEM\r\n');
    fprintf(fid, '      POSITION 0\r\n');
    fprintf(fid, '      SNAPOFFS 0\r\n');
    fprintf(fid, '      LENGTH %.6f\r\n', duration);
    fprintf(fid, '      LOOP 0\r\n');
    fprintf(fid, '      ALLTAKES 0\r\n');
    fprintf(fid, '      FADEIN 1 0 0 1 0 0 0\r\n');
    fprintf(fid, '      FADEOUT 1 0 0 1 0 0 0\r\n');
    fprintf(fid, '      MUTE 0 0\r\n');
    fprintf(fid, '      SEL 0\r\n');
    fprintf(fid, '      IGUID %s\r\n', guid_imidi);
    fprintf(fid, '      IID 2\r\n');
    fprintf(fid, '      NAME "%s"\r\n', midi_path);
    fprintf(fid, '      VOLPAN 1 0 1 -1\r\n');
    fprintf(fid, '      SOFFS %.6f 0\r\n', t_start);
    fprintf(fid, '      PLAYRATE 1 1 0 -1 0 0.0025\r\n');
    fprintf(fid, '      CHANMODE 0\r\n');
    fprintf(fid, '      GUID %s\r\n', guid_imidi);
    fprintf(fid, '      <SOURCE MIDI\r\n');
    fprintf(fid, '        FILE "%s"\r\n', midi_path);
    fprintf(fid, '        IGNTEMPO 0 120 4 4\r\n');
    fprintf(fid, '      >\r\n');
    fprintf(fid, '    >\r\n');
    fprintf(fid, '  >\r\n');
    fprintf(fid, '>\r\n');

    fclose(fid);
end


function g = guid_new()
%GUID_NEW Génère un GUID aléatoire au format Reaper {XXXXXXXX-XXXX-XXXX-XXXX-XXXXXXXXXXXX}
    hex_chars = '0123456789ABCDEF';
    rnd = randi(16, 1, 32);
    h   = hex_chars(rnd);
    g   = sprintf('{%s-%s-%s-%s-%s}', h(1:8), h(9:12), h(13:16), h(17:20), h(21:32));
end


function ticks = read_midi_ticks(midi_path)
%READ_MIDI_TICKS Lit le nombre de ticks par noire dans l'en-tête MIDI.
    ticks = 960;   % valeur par défaut
    try
        fid = fopen(midi_path, 'rb');
        if fid == -1, return; end
        header = fread(fid, 14, 'uint8')';
        fclose(fid);
        if length(header) >= 14 && isequal(char(header(1:4)), 'MThd')
            ticks = header(13) * 256 + header(14);
        end
    catch; end
end