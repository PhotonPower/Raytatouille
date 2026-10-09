# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/); Versionen nach SemVer.
Neue Einträge stehen bis zum nächsten Release als Fragmente in [`changelog.d/`](changelog.d/README.md)
(`python tools/changelog.py preview` zeigt den vollständigen Abschnitt [Unreleased]).

## [Unreleased]

## [0.6.0] – M4 Multi-Path

### Hinzugefügt
- `rtt-analysis`: Pfadauswertung (`rtt/analysis/paths.hpp`): Transmission je Pfad
  (`path_transmission`) und OPL-Differenz zweier Pfade (`opl_difference`) für vorgegebene
  Startstrahlen, auch in gefalteten Systemen, oder über `make_rays` in rotationssymmetrischen;
  mit Verlustzählung und Abbruch/Fortschritt. Referenzdatei
  `tests/reference/m4/michelson_offset.rtt.json` (#122).
- `rtt-compile`: Ghost-Generator (`rtt/compile/ghosts.hpp`, ADR 0027): `ghost_paths` erzeugt
  aus einem Pfad alle Zweifach-Reflexions-Ghosts als explizite Pfade (N Refract-Ereignisse
  ergeben N(N − 1)/2), `compile_with_ghosts` kompiliert eine Kopie des Systems mit ihnen und
  nennt je Ghost Basispfad und Reflexionsflächen. Keine Formatänderung (#123).
- `rtt-analysis`: Ghost-Ranking (`rtt/analysis/ghosts.hpp`, Nachtrag zu ADR 0027):
  `ghost_ranking` verfolgt alle Ghosts eines Basispfads mit dessen Startstrahlen und sortiert
  sie nach der Bestrahlungsstärke am Bild relativ zum Nutzbild, mit Auflösungsgrenze des
  Detektors (Standard 5 µm); dazu paraxiale Fokuslage und Unschärfe je Ghost (#124).
- Dateiformat: Flächenfeld `diffraction_efficiency` (Leistungsanteil je Beugungsordnung, ADR 0025)
  mit Code `surface.efficiency_invalid`; Python `Surface.diffraction_efficiency` und
  `model.DiffractionEfficiency` (#125).
- `rtt-geom`: Phasenfunktionen der Phasenschichten (`rtt/geom/phase.hpp`, ADR 0025):
  `LinearGratingPhase` und `RadialPhasePolynomial` mit Phase in rad und Gradient in rad/mm,
  Summe der Schichten einer Fläche und `tangential_gradient` für die Projektion in die
  Tangentialebene, geprüft gegen Mansuripur (2007), Gl. (9)–(11) (#126).
- `rtt-trace`: neuer Strahlstatus `RayStatus::Evanescent` (ADR 0025, Punkt 7) für Beugungsordnungen
  ohne reelle Richtung, hinten an das Enum angehängt (`kRayStatusCount` = 8); in Python
  `RayStatus.EVANESCENT`. `TraceStats.rays` und `RayLosses.by_status` haben einen Eintrag mehr,
  Ergebnisformat 0.1.3 (verträglicher Zusatz, ADR 0023) (#127).
- `rtt-trace`: Beugungsordnungen an Flächen mit Phasenschicht (ADR 0025): lokale
  Gittergleichung für `refract`, `reflect` und `transmit` mit `order`, Status `Evanescent` für
  Ordnungen ohne reelle Richtung (Tir der Ordnung 0 zuerst), OPL-Beitrag m φ λ₀/(2π),
  Polarisation P = R(k₀ → k_m)·P₀, Effizienz je Ordnung (`diffraction_efficiency`) im Gewicht.
  `first_order` nimmt Phasenflächen mit Ordnung 0 an, und Ray Aiming geht durch solche Flächen
  vor der Blende (auch mit Effizienzen); für Ordnungen ≠ 0 braucht ein Pfad eigene Startstrahlen
  (ADR 0025, Punkt 4). Neues Referenzsystem `tests/reference/m4/grating_transmission.rtt.json` (#127).
- Dateiformat: einachsige Kristalle als Material `{"ordinary": …, "extraordinary": …}` und
  Elementfeld `optic_axis` (ADR 0026), mit den validate-Codes `crystal.*`; Python
  `Element.crystal`, `Element.optic_axis` und `model.CrystalMaterial` (#128).
- `rtt-material`: `UniaxialMaterial` für einachsige Kristalle aus zwei Materialien für n_O und
  n_E (ADR 0026) und `MaterialLibrary::resolve_uniaxial(ordinary, extraordinary)`, z. B. für das
  Katalogpaar `X` / `X-E`; Wellenlängenbereich als Schnittmenge beider Teile (#129).
- `rtt-polar`: Brechung in ein einachsiges Medium und aus ihm heraus (`rtt/polar/birefringence.hpp`,
  ADR 0026): Wellennormale, Energierichtung mit Walk-off, Brechzahl der o- und e-Mode nach Lam
  (2.11), (2.39), (2.41), Eigenpolarisationen und PRT-Matrizen beim Ein- und Austritt im
  Projektionsmodell von M4 (#130).
- `rtt-compile`: Kristallelemente kompilieren (ADR 0026): eigenes Medium je Kristall mit n_O, n_E
  und globaler optischer Achse, `CompiledEvent::crystal_mode`, Pfadregeln für `ordinary` und
  `extraordinary` mit den Codes `paths.crystal_mode_required`, `paths.mode_without_crystal`,
  `crystal.unsupported`, `crystal.interaction_unsupported` und `crystal.absorbing` (#131).
- `rtt-trace`: Strahlen durch einachsige Kristalle (ADR 0026): Eintritt mit `ordinary` oder
  `extraordinary`, Walk-off des e-Strahls entlang des Poynting-Vektors, optischer Weg nach Lam
  Gl. (2.17), Austritt mit `refract`, auch mit Beugungsordnung. `RayBatch` hat die Spalten
  `wave_x`, `wave_y`, `wave_z` (Wellennormale) und `mode_index`; ein Batch, der nur `dir` setzt,
  bleibt gültig. Abnahme: `tests/reference/m4/calcite_walkoff.rtt.json`. Grenze in M4: `first_order`
  und alles darauf (`make_rays` mit Ray Aiming, OPD, Komfortformen der Analysen) lehnen Pfade in
  einen Kristall ab; solche Pfade brauchen eigene Startstrahlen (#132).
- `rtt-py`: Pfadauswertung und Ghosts aus Python: `rt.analysis.path_transmission` und
  `rt.analysis.opl_difference` mit Startstrahlen (`start=`, auch für gefaltete Pfade) oder über
  `make_rays` (`field`, `wavelength`, `rays`, `aiming`), `rt.compile_with_ghosts`,
  `rt.ghost_paths` und `rt.analysis.ghost_ranking`; Ergebnisse `PathTransmission`,
  `PathOplDifference`, `GhostSystem`, `GhostRanking` mit Spalten als NumPy-Kopien, NaN für
  fehlende Diagnosewerte; Beispiel `examples/python/ghosts.py` (#133).
- Ergebnisformat `raytatouille-result` 0.1.4: `PathTransmission`, `PathOplDifference` und
  `GhostRanking` haben `to_dict()` und `to_json()` (ADR 0023, Nachtrag) (#133).
- `rtt-py`: Beugungsordnungen aus Python bitgleich zu C++ (Gitterbank `m4/grating_transmission`
  mit den Ordnungen −1, 0, +1 und der evaneszenten Ordnung +6, Reflexionsgitter,
  Beugungseffizienzen) und gegen die Gittergleichung geprüft; Beispiel
  `examples/python/grating_orders.py` (#134).
- `rtt-py`: einachsige Kristalle aus Python: `RayBatch.wave_x`/`wave_y`/`wave_z` und
  `mode_index` als beschreibbare NumPy-Ansichten, `CompiledMedium.is_crystal`,
  `reference_extraordinary`, `index_extraordinary` und `optic_axis`; Calcit-Platte bitgleich zu
  C++ und gegen die geschlossenen Formen geprüft; Beispiel
  `examples/python/calcite_double_image.py` (#134).
- Ergebnisformat `raytatouille-result` 0.1.5: `RayBatch` mit `wave_x`, `wave_y`, `wave_z` und
  `mode_index`; ältere Dateien bleiben gültig und laden mit wave = dir, mode_index = 0
  (ADR 0023, Nachtrag) (#134).
- Abnahme M4 Multi-Path (`libs/rtt-analysis/tests/test_m4_acceptance.cpp`, Tag `[m4]`):
  Michelson-Bilanz, Calcit-Walk-off (t·tan ρ), Gittergleichung mit Evaneszenz, Ghost-Ranking an
  zwei Platten und am Singlet sowie die erste
  Beugungsordnung der Feature-Tour gegen Mansuripur (7b), jeweils durch den Tracer mit
  Referenzsystem und analytischem Sollwert. Neue Referenzsysteme
  `tests/reference/m4/ghost_plates.rtt.json` und `ghost_singlet.rtt.json`; Python bitgleich zu C++
  für sie und für die Feature-Tour (#135).

### Geändert
- JSON-Schema: `weight` von Wellenlängen und Feldpunkten muss ≥ 0 sein, wie `validate` es schon
  verlangt; keine neue Schema-Version, weil das Modell unverändert bleibt (#35).
- `rtt-paraxial`: Ein Pfad gilt als afokal, wenn seine Brechkraft bis auf die Rundung null
  ist (|Φ| ≤ 16·N·u·S, S aus demselben y-nu-Trace in Beträgen), statt bei der festen Schwelle
  1e−14 /mm. Kleine, fast afokale Systeme werden damit nicht mehr fälschlich fokal, sehr
  schwache Linsen (EFL > 1e14 mm) nicht mehr fälschlich afokal (#35).
- Dateiformat: Schema 0.3.0 (ADR 0025, ADR 0026). Beugungsordnungen `order` an jedem Ereignis
  einer Fläche mit Phasenschicht; an Flächen ohne Phasenschicht meldet `validate` weiter
  `paths.order_not_allowed`. Dateien mit Schema 0.1 und 0.2 werden gelesen und migriert
  (`diffract` → `transmit` mit derselben Ordnung, die Medien bleiben gleich); `rtt format`
  schreibt sie als 0.3.0 (#125).
- `tests/reference/m0/feature_tour.rtt.json`: Die Beugungsordnungen an der Gitterfläche `G.S1` sind
  jetzt `refract` mit `order` statt `transmit` (ADR 0025, Punkt 4; Freigabe des Maintainers): Der
  Strahl tritt gebeugt in die Gitterplatte ein, statt außen zu bleiben (#135).

### Behoben
- AGF: Eine einzelne TD-Zeile mit mehr als 7 Werten ist jetzt ein `AgfError` mit Datei und Zeile,
  wie schon bei einer TD-Folgezeile; bisher scheiterte sie erst in `CatalogMaterial` ohne Ort
  (#35).

### Entfernt
- Ereignisart `diffract` (`EventKind::Diffract`, Python `EventKind.DIFFRACT`) entfernt; eine
  Beugungsordnung steht jetzt als `order` am Ereignis selbst (ADR 0025). Ältere Dateien lesen
  weiter über die Migration (#125).

## [0.5.0] – GUI-Grundlagen

### Hinzugefügt
- `rtt-material`: Materialbibliothek für eine GUI (#85, G9): `MaterialLibrary::catalogs()` und
  `catalog(name)` listen die geladenen AGF-Kataloge so, wie sie in der Datei stehen (NM-Extras
  status/exclude sub/melt freq, GC, ED, TD, MD, OD mit Resistenzklassen auch als Bereich „a-b“,
  IT); Katalog-Alias `add_catalog(path, name)` für gleichnamige Dateien; `add_catalog_text` lädt
  aus dem Speicher; `index_many` ist bitgleich zum skalaren `index`. Einheiten von ED, MD und IT
  belegt am SCHOTT-Datenblatt N-BK7 (`docs/quellen.md`).
- `rtt-py`: Modul `raytatouille.materials` mit `MaterialLibrary` (Unterklasse der C++-Bibliothek,
  `rt.MaterialLibrary`): `glasses(catalog)`, `glass(reference)` (`GlassInfo`, `ClassRange`,
  `ThermalData`, `MechanicalData`), `glass_map()` mit n_d/ν_d aus den NM-Datensätzen;
  `add_catalog(path, name=None)`, `add_catalog_text(data, name)`, `catalogs()`; `index()` nimmt
  auch NumPy-Arrays von Wellenlängen (bitgleich zum skalaren Aufruf) (#85).
- Stabile Diagnosecodes (ADR 0022, angenommen; #86):
  - Jede `Diagnostic` trägt einen Code in Punktnotation, z. B. `material.unknown`. Die Codes bleiben
    über Versionen gleich; Liste in `docs/diagnostics.md`.
  - Registry in der neuen header-only Bibliothek `rtt-diagnostics` (Schicht Basis). Ein nicht
    registrierter Code kompiliert nicht.
  - Python: `Diagnostic.code`, `CompileError.codes`, `rt.diagnostics.CODES`.
- Diagnosekanal (#86): `CompiledSystem::diagnostics()` hält die Warnungen von `validate()`. Python
  gibt sie zusätzlich als `rt.errors.RaytatouilleWarning` mit `code` und `location` aus.
- Fehlerorte (#86):
  - `CompiledSurface::location` (Python `CompiledSystem.surface_locations`) ist der JSON-Pointer
    jeder Fläche.
  - `ParaxialError` mit `surface` und `location`.
  - `AnalysisError` mit `surface`, `location`, `ray_status`, `field` und `wavelength` des verlorenen
    Strahls.
- `rtt::compile::NoStopError` mit `path_name` und `location` und die einzige Prüfung
  `rtt::compile::require_stop` (`rtt/compile/errors.hpp`). Python: `rt.NoStopError(ParaxialError,
  AnalysisError, ValueError)` (#86).
- `rtt-paraxial`: Prescription-Daten `prescription()` (`prescription.hpp`, G5).
  - Je Event Rand- und Hauptstrahl mit Höhe, Steigung und paraxialem Einfallswinkel sowie die
    Lagrange-Invariante.
  - Systemdaten: Gesamtlänge (bei Spiegeln abgewickelt), Objektabstand, paraxiales Arbeits-F/#,
    paraxiale bildseitige NA, Lagrange-Invariante und die First-Order-Daten.
  - Dieselben Strahlen wie `seidel()`. Ohne Blende oder Pupille sind die betroffenen Werte leer,
    ohne Fehler.
  - Quelle Greivenkamp, OPTI-502 Sec. 9 „Stops and Pupils“, in `docs/quellen.md`.
  - Python `rt.paraxial.prescription` mit NumPy-Arrays je Fläche (NaN für leere Werte), bitgleich
    gegen C++; Tabelle im Beispiel `examples/python/singlet.py` (#84).
- `rtt-trace`: Strahlpfad-Aufzeichnung für die GUI. `SequentialTracer::trace(..., RayPaths&,
  record_rays, max_recorded_rays)` zeichnet Position, Richtung, OPL, `weight` und Status vor dem
  ersten und nach jedem Ereignis auf; verlorene Strahlen sind über `count`, `lost_at` und NaN
  markiert. Ohne Aufzeichnung bleibt das Ergebnis bitgleich (Recorder-Policy), und der letzte
  gültige Slot ist bitgleich mit dem Endzustand. Strahlen werden über `record_rays` ausgewählt,
  ohne Auswahl höchstens `max_recorded_rays` (Standard 10 000), sonst ein Fehler mit dem
  Speicherbedarf (#80).
- `rtt-py`: `rt.trace.trace(..., record_path=True, record_rays=None, max_recorded_rays=10000)` gibt
  `(TraceStats, RayPaths)` zurück, mit NumPy-Views (N, S, 3) bzw. (N, S) ohne Kopie; bitgleich
  gegen C++; Beispiel `examples/python/ray_paths.py` (#80).
- `rtt-compile`: Geometrie-Export für das Layout (`rtt/compile/layout.hpp`).
  - `CompiledSystem::elements()` (Name, Art, Flächen, Medien, `segmented`).
  - `CompiledSurface::element`, `medium_front` und `medium_back`: in der Flächenreihenfolge des Elements, Regeln von ADR 0017.
  - `surface_sag` und `surface_normal`, vektorisiert, lokal oder global.
  - `surface_profile` in Schnittebenen parallel zur lokalen z-Achse, begrenzt durch die Apertur (der Ring ergibt zwei Stücke) und den Formbereich.
  - `element_outlines`: geschlossene Polygone der Glassegmente von Linsen und segmentierten Platten, Randkanten als Stufe (#81).
- `rtt-py`: Modul `rt.layout` mit `surfaces`, `elements`, `sag` und `normal` mit Broadcasting, `profile` und `outlines` (Schnittebene `"yz"`, `"xz"` oder Punkt und Normale); Beispiel `examples/python/layout.py` (Achromat) (#81).
- Tests aus dem Aufräumen: Der Trace lässt Feld, Pupillenkoordinaten und Wellenlänge jedes
  Strahls unverändert, auch bei verlorenen Strahlen. Die Aperturtypen `image_fnumber` und
  `object_na` sind von der Systemangabe bis zum gezielten und verfolgten Strahl geprüft (`object_na`
  paraxial gelesen, tan U = NA/n). Der Schwerpunkt des polychromatischen Spots ist mit
  dispersivem Glas als Σ w·S·c über die Wellenlängen geprüft. Der Brennpunkt der Parabel ist im
  Test hergeleitet (#35).
- `rtt-py`-Tests: Jeder Test läuft unter einem Watchdog (`faulthandler`, Standard 300 s,
  einstellbar mit `RTT_PY_TEST_TIMEOUT`, 0 schaltet ihn ab). Hängt ein Test, gibt der Watchdog die
  Tracebacks aller Threads aus und beendet den Lauf mit Status 1, auch wenn das GIL in C-Code
  gehalten wird, statt die CI bis zum Job-Limit zu blockieren. `rtt_py_pytest` setzt `MPLBACKEND=Agg` (#35).
- Lade-Warnungen der Glaskataloge (ADR 0022, Nachtrag ADR 0008):
  `rtt::material::LoadWarning{code, file, line, message}` in `AgfCatalog::warnings` und
  `MaterialLibrary::load_warnings()`; in Python `MaterialLibrary.load_warnings` und
  `RaytatouilleWarning` beim Laden; Codes `agf.duplicate_glass`, `agf.duplicate_glass_conflict`,
  `agf.preamble_skipped`, `agf.stray_line` (Erzeuger `agf`). `LoadWarning` hat `to_dict()` und
  `to_json()`; das Ergebnisformat `raytatouille-result` wird damit 0.1.1 (ADR 0023) (#71).
- Python liest den ganzen Modellbaum (ADR 0024): `System.root`, `paths`, `fields`, `aperture`
  und `object_space` liefern unveränderliche typisierte Kopien aus dem neuen Modul
  `raytatouille.model` (`Assembly`, `Element`, `Surface`, `Param`, `Pose`, Formen, Aperturen,
  Phasen, Interaktionen, `Path`, `Event` und die Aufzählungen). `Wavelength` hat `__eq__` und
  `__hash__`. Beispiel `examples/python/model_tree.py` (#82).
- `rtt-io`: Systeme mit JSON Patch (RFC 6902) ändern, mit inversen Patches für Undo (ADR 0024):
  `rtt::io::apply_patch(system, patch_json, EditCheck)` auf der Bearbeitungsform
  (`rtt::io::to_edit_json`: alle Werte geschrieben, jedes `Param` als Objekt). Die Operationen
  wendet `rtt-io` selbst an (`nlohmann::json::patch` weicht von RFC 6902 ab). Ein Patch gilt ganz
  oder gar nicht; abgelehnt wird er, wenn er Fehler von `validate()` hinzufügt (`NoNewErrors`).
  Fehler als `rtt::io::EditError` mit Code, Pointer und Operation; neue Codes `edit.*` mit
  Erzeuger `edit` (#82).
- Python ändert Systeme mit JSON Patch und Undo (ADR 0024):
  - `rt.Editor` mit `apply`, `set`, `insert`, `remove`, `move`, `undo` und `redo`; die History
    lässt sich als JSON exportieren und mit `Editor.from_history` wieder abspielen;
  - `rt.apply_patch` und `rt.apply_patch_with_inverse` als Skriptbefehle;
  - `rt.EditError` mit Code, Pointer, Operation und Diagnosen;
  - `System.to_dict()` (Bearbeitungsform), `json_at`, `locate_surface`, `locate_node`,
    `locate_path`;
  - Beispiel `examples/python/edit_undo.py` (#82).
- `rtt-trace`: Abbruch und Fortschritt für lange Rechnungen (`rtt/trace/run_control.hpp`):
  `CancelToken`, `RunControl` (Callback `progress`, Mindestabstand, Blockgröße), Exception
  `Cancelled`. Überladungen von `SequentialTracer::trace` und `make_rays` mit `RunControl`;
  Abbruch nach höchstens einem Block je Worker, Exceptions nur nach dem parallelen Teil.
  Ergebnisse bitgleich und unabhängig von Blockgröße und Threads (ADR 0004, Nachtrag) (#83).
- `rtt-analysis`: Überladungen mit `RunControl` für `spot`, `ray_fan`, `opd_map`, `opd_fan`,
  `distortion`, `field_curvature`, `longitudinal_colour` und `lateral_colour`; die alten
  Signaturen bleiben (#83).
- Benchmark `apps/rtt-bench` (Spot des Cooke-Tripletts mit 1027 Strahlen); in der CI mit
  Zeiten in der Job-Summary und einer Warnung über 30 ms, nie ein Fehler (#83).
- `rtt-py`: Abbruch und Fortschritt aus Python: `rt.CancelToken`, Schlüsselwortargumente
  `cancel=` und `progress=` (`progress(done, total, stage)`, Typ `rt.trace.ProgressCallback`)
  an `rt.trace.trace`, `rt.trace.make_rays` und den Analysen `spot`, `ray_fan`, `opd_map`,
  `opd_fan`, `distortion`, `field_curvature`, `longitudinal_colour` und `lateral_colour`;
  neue Exception `rt.errors.Cancelled`. Die Rechnung läuft ohne GIL, der Callback holt es
  sich; seine Exception wird weitergereicht (#83).
- Ergebnisformat `raytatouille-result` 0.1.0 (ADR 0023):
  - Alle Ergebnisobjekte der Python-API, auch die Geometrie (`SurfaceLayout`, `CompiledElement`),
    haben `to_dict()` und `to_json()`; `rt.results.load_json` liest die Daten bitgleich zurück
    (Arrays mit dtype und shape, NaN/±∞ als Zeichenketten, standardkonformes JSON). Die
    C++-Objekte werden nicht rekonstruiert.
  - JSON-Schema `schema/raytatouille-result.schema.json`, Beispiele je Typ unter
    `tests/reference/results/`, Beispielskript `examples/python/results.py` (#86).
- Strahlverluste als Daten: `SpotDiagram`, `RayFan`, `OpdMap` und `OpdFan` haben `losses`
  (`RayLosses`): gestartete Strahlen, Zählung je `RayStatus` und die Fläche, an der die meisten
  verlorenen Strahlen enden (#86).
- Warnungen der Analysen (ADR 0023): Feld `warnings` in denselben vier Ergebnissen, in Python
  zusätzlich als `RaytatouilleWarning`. `rays.lost`, wenn mehr Strahlen verloren gehen als
  `lost_warning_fraction` (Option, Standard 0,5); `stop.clips_beam`, wenn Strahlen an der Blende
  vignettiert werden. Neue compile-Warnung `stop.not_on_path`, wenn ein Pfad das Blendenelement
  nicht besucht (#86).
- Diagnose-Registry: `CodeInfo::producer` (validate, compile, analysis), in Python
  `rt.diagnostics.CODES[...].producer` (#86).
- Ergebnisformat `raytatouille-result` 0.1.2: `Prescription` (paraxiale Prescription-Daten, #84)
  hat `to_dict()` und `to_json()` wie die übrigen Ergebnistypen (ADR 0023) (#86).
- `rtt-trace`: Objektseitig telezentrische Systeme mit endlichem Objekt (Eintrittspupille im
  Unendlichen) lassen sich tracen.
  - Pupillenkoordinaten sind dann Objektraum-Steigungen, d ∝ (px·u_m, py·u_m, 1); der
    paraxiale Hauptstrahl läuft parallel zur Achse.
  - Beispiel `tests/reference/m1/telecentric_singlet.rtt.json`.
  - u_m = NA/n bei `object_na`, r_Blende/|s| bei `stop_size`; `entrance_pupil_diameter`,
    `image_fnumber` und Feldwinkel ergeben dort einen klaren Fehler.
  - Paraxiale Bildhöhen werden über den achsparallelen Hauptstrahl umgerechnet.
  - `rtt-analysis`: Verzeichnung und Bildfeldwölbung rechnen mit diesem Hauptstrahl; Spot,
    Fans und OPD laufen, solange die Austrittspupille endlich ist.
  - Reale Zielung: Die Schrittweite der Jacobi-Matrix hat eine Untergrenze von 1e−10 des
    Pupillenradius in Pupilleneinheiten (Nachtrag zu ADR 0007). So konvergiert sie auch bei
    sehr weit entfernter EP; bestehende Ergebnisse bleiben bitgleich (#96).
- `rtt-analysis`: OPD bei Austrittspupille im Unendlichen (bild- bzw. beidseitig
  telezentrisch): Referenz ist der Grenzfall R → ∞ der Referenzsphäre um den
  Hauptstrahl-Bildpunkt, `ReferenceSphere::radius` = +∞; OPD-Fächer und -Karte laufen,
  stetig zum endlichen Fall. Referenzsystem `tests/reference/m2/telecentric_4f.rtt.json`
  (4f-Relais aus zwei Descartes-Linsen) (#102).
- `rtt-py`: `rt.trace.make_rays(..., threads=)` begrenzt die Worker-Threads der jetzt parallelen
  Strahlzielung wie bei `rt.trace.trace`. Die Strahlen sind für jede Thread-Anzahl bitgleich
  (#119).

### Geändert
- `rtt-material`: Die AGF-Datensätze NM-Extras, MD, OD und IT werden jetzt streng gelesen statt
  übersprungen: eine falsche Feldanzahl oder ein unbekanntes Token in GC, OD, MD, IT oder den
  NM-Extras ist jetzt ein `AgfError` mit Zeile; nur die Platzhalter der Herstellerdateien (`_`, `-`)
  gelten als „nicht vorhanden“, die melt freq wird als ganze Zahl wie geschrieben gespeichert
  (Nachtrag ADR 0008). Geprüft und weiter ladbar: HOYA, OHARA, Sumita, CDGM, SCHOTT (aktuell und
  IRG), LightPath sowie weitere ältere Kataloge der Zemax-Verteilung; NIKON-HIKARI scheitert wie
  vorher an #71 (#85).
- `rtt validate` und `to_string(Diagnostic)` schreiben den Code mit: `error [code] /pointer:
  Meldung`. Ebenso die Meldung von `CompileError`. Die Exit-Codes bleiben (#86).
- **C++:** Ein Pfad ohne Blende wirft jetzt überall `rtt::compile::NoStopError` (ein
  `std::invalid_argument`) statt `ParaxialError` (`seidel()`) bzw. `AnalysisError` (OPD). Wer nur
  diese Typen fängt, fängt den Fall nicht mehr. Strahlquellen, Spot und Fächer warfen schon
  `std::invalid_argument`. In Python fangen bestehende `except ParaxialError`, `AnalysisError` und
  `ValueError` den Fall weiter, weil `NoStopError` von allen dreien erbt (#86).
- `compile()` verwirft die Warnungen von `validate()` nicht mehr (#86).
- `rtt-polar`: `prt_matrix()` und `geometric_transform()` bilden die Außenprodukte elementweise
  statt als Eigen-Lazy-Ausdruck, mit bitgleichen Ergebnissen. Die falsche GCC-Warnung
  `-Wnull-dereference` bei `-O2` tritt damit nicht mehr auf, und die Unterdrückung in `rtt-py`
  (Hilfsdatei `polar_matrix`) ist entfernt (#35).
- `rtt-material`: NIKON-HIKARI_201911 und die älteren Kataloge des Pakets ZemaxGlass laden jetzt,
  mit Lade-Warnung statt `AgfError`: identische doppelte Gläser werden zusammengelegt, verschiedene
  machen das Glas mehrdeutig (`UnknownMaterial` mit Alias-Hinweis), Textzeilen vor dem ersten
  Datensatz und ein einzelnes Wort vor dem nächsten NM werden übersprungen. Echte Datenfehler
  alter Kopien bleiben `AgfError` (#71).
- **C++:** `rtt::diagnostics::CodeInfo` hat das neue Feld `producer` zwischen `severity` und
  `summary`; eine Brace-Initialisierung von `CodeInfo` außerhalb der Registry muss angepasst
  werden (ADR 0023, #86).
- `compile()` gibt die neue Warnung `stop.not_on_path` aus; Spot, Fächer und OPD haben die
  Option `lost_warning_fraction` (in Python der Parameter gleichen Namens) (#86).
- `rtt-trace`: Pupillenkoordinaten sind an der Referenzwellenlänge orientiert (präzisiert
  #93). Liegt bei endlichem Objekt die EP einer Wellenlänge auf der anderen Seite des Objekts als
  die der Referenzwellenlänge, wird (px, py) gespiegelt verwendet, sodass ein Pupillenpunkt für
  jede Wellenlänge dieselbe Blendenseite meint; Fächer und OPD-Karten bleiben über λ
  vergleichbar. Betrifft nur Systeme mit Seitenwechsel der EP über λ (seit #93 möglich); sonst
  bitgleich (#96).
- `rtt-analysis`: Die OPL bis zur Referenzsphäre wird auslöschungsfrei gerechnet (der
  gemeinsame Term |n′|R entfällt). Bei sehr weit entfernter Austrittspupille war die OPD
  vorher durch Rundung unbrauchbar; bei bestehenden Systemen ändern sich die Werte nur in der
  Rundung (#102).
- `rtt-trace`: `make_rays` zielt die Strahlen parallel (oneTBB), mit denselben
  Zerlegungsregeln wie der Tracer, auch mit Abbruch und Fortschritt. Die Ergebnisse bleiben
  bitgleich und hängen nicht von der Thread-Anzahl ab. Ein Spot des Cooke-Tripletts mit 1027
  Strahlen dauert mit 12 Threads etwa 5 ms statt 20 ms (ADR 0004, Nachtrag) (#119).

### Behoben
- `rtt-trace`: Bei endlichem Objekt und virtueller Eintrittspupille hinter dem Objekt
  (z_EP < z_Objekt, z. B. Blende hinter dem hinteren Brennpunkt der Gruppe davor wie bei
  Mikroskopobjektiven mit Tubuslinse) starteten alle Strahlen in −z und gingen verloren.
  - Jetzt läuft der Strahl auf der Geraden durch Objektpunkt und EP-Punkt ins System, also
    in +z.
  - Die Zuordnung Pupillenkoordinate → Blendenpunkt, die reale Zielung und die
    OPD-Referenzsphäre bleiben unverändert.
  - Damit gilt auch die Feldwinkel-Konvention aus #8 (Hauptstrahl steigt bei θ_y > 0 in +y)
    für diesen Fall.
  - Für z_EP > z_Objekt sind die Startrichtungen bitgleich zu vorher (#93).
- `rtt-paraxial`: `seidel()` und `prescription()` rechnen paraxiale Bildhöhen und Feldwinkel bei
  endlichem Objekt jetzt wie `rtt-trace` bei der Referenzwellenlänge in den Hauptstrahl um, statt
  bei der übergebenen Wellenlänge. Bei einer anderen Wellenlänge als der Referenz stimmt der
  Hauptstrahl damit mit dem paraxial gezielten Strahl von `make_rays` überein; bei der
  Referenzwellenlänge ändert sich nichts (#35).
- `rtt-io`: Eine ganze Zahl außerhalb von `int` (z. B. `order` = 4294967296) ist ein
  `ParseError` mit Pointer statt still abgeschnitten (ADR 0008); das JSON-Schema nennt die
  Grenzen bei `order` (#35).
- Doku: `object_na` wird paraxial gelesen, die Randstrahlsteigung vom axialen Objektpunkt ist
  u = NA / n (Greivenkamp, OPTI-502, S. 9-34), in `system.hpp`, `docs/dateiformat.md` und im
  Python-Docstring von `SystemApertureType`; `docs/dateiformat.md` nennt jetzt Bedeutung und
  Einheit aller Aperturtypen (#35).
- `rtt-io`: Der JSON-Pointer eines `ParseError` maskiert `~` und `/` in Schlüsseln nach RFC 6901
  (`~0`, `~1`), z. B. bei einem unbekannten Schlüssel `a/b` (#82).

## [0.4.0] – M3 Polarisation

### Hinzugefügt
- Abnahme M3 (`libs/rtt-analysis/tests/test_m3_acceptance.cpp`): die M3-Zeilen der Tabelle
  „Validierung und Tests“ Ende zu Ende durch den Tracer mit konventionsfreien Größen: Fresnel an
  N-BK7 und Brewster, Totalreflexion (Grenzwinkel auf 1e−12 in ξ, s/p-Phasensprung),
  λ/4-AR aus MgF₂, verlustfreier Schichtstapel R + T = 1, Malus, λ/4 unter 45° (mit Händigkeit),
  λ/2 unter θ, Metallspiegel unter 45° als ellipsometrisches Verhältnis gegen den idealen Leiter
  (Testmetall ñ = 1,2 + 7,26i), Metall-Periskop und Fresnel-Rhombus. Referenzsysteme unter
  `tests/reference/m3/`, Testkatalog `tests/catalogs/coatings/m3.json` (#63).
- `rtt-polar` (neue Bibliothek, Schicht Physik, nur `rtt-math`): Fresnel-Amplituden r_s, r_p, t_s, t_p
  mit Phase für komplexe Indizes ñ = n + iκ, inklusive Totalreflexion (`fresnel`, kanonische Eingabe
  ξ = Re(ñ_i) sin θ_i; `fresnel_at_angle` für nicht absorbierende Einfallsmedien), Reflexions- und Transmissionsgrade (`fresnel_power`), Normalkomponente mit
  festem Wurzelzweig, Diattenuation und Phasendifferenz Δ; geprüfte Varianten `fresnel_checked`
  und `fresnel_power_checked`. Konvention (Convention A, p = k × s) mit Begründung über den
  idealen Leiter in `docs/architecture.md`; Quelle Byrnes, arXiv:1603.02720 (#56).
- `rtt-coating` (neue Bibliothek, Schicht Physik): Transfermatrix-Methode für Schichtsysteme mit
  komplexen Indizes nach Byrnes (arXiv:1603.02720), Amplituden r_s, r_p, t_s, t_p in der Basis
  [s, p, k] von `rtt-polar` (Convention A, r_p = −r_s bei senkrechtem Einfall), Leistungen R, T, A;
  Winkel als Tangentialinvariante ξ = n sinθ; Dicke in µm oder als QWOT; ideale Coatings ohne
  Retardance (AR, Strahlteiler, Spiegel r_s = −1, r_p = +1) und tabellierte Coatings (bilinear über
  Winkel und Wellenlänge). Konvention „Polarisation und Fresnel“ aus #56, Quelle Byrnes ergänzt
  in `docs/quellen.md` (#58).
- `rtt-polar`: 3D-Polarisations-Raytracing-Matrix P = O_out·diag(a_s, a_p, 1)·O_inᵀ (`prt_basis`,
  `prt_matrix`) mit deterministischer s-Richtung bei senkrechtem Einfall, geometrische Transformation
  Q (`geometric_transform`), Diattenuation aus den Singulärwerten des Transversalanteils und
  Retardance über die Polarzerlegung, gesamt und physikalisch (Q⁻¹P); Golden-Test nach Lam
  (drei Prismen); Quelle Lam, Dissertation Arizona (#57).
- `rtt-coating`: Coating-Kataloge als eigene JSON-Dateien (ADR 0019, Format
  `raytatouille-coatings` 0.1.0, `schema/raytatouille-coatings.schema.json` mit pytest):
  `CoatingDesign` mit Schichten (Materialverweis, `thickness_um` oder `qwot`) und
  Designwellenlänge, strenger Parser mit Datei und JSON-Pointer, `CoatingLibrary` mit
  `add_catalog` (Datei oder Verzeichnis, alles oder nichts) und `resolve("CATALOG:NAME")`;
  Demo-Katalog `tests/catalogs/coatings/demo.json` (#59).
- `rtt-compile`: `compile(system, materials, coatings)` löst `CoatingRef` auf, wertet die
  Schichten je Systemwellenlänge aus (`CompiledSystem::coatings()`) und legt das Substrat je
  Fläche fest (`CompiledSurface::coating`); Kompilierfehler mit Pointer für unbekannte
  Coatings und Schichtmaterialien, Wellenlängenbereiche und mehrdeutige Substratseiten.
  Beispiel `tests/reference/m3/ar_singlet.rtt.json` (#59).
- `rtt-polar`: ideale Elemente (`ideal.hpp`): allgemeine Jones-Matrix in 3D eingebettet,
  linearer Polarisator mit Extinktionsverhältnis als Leistung, linearer und zirkularer Retarder mit
  symmetrischer Phase ±δ/2, Stokes-Parameter; Händigkeit und S3-Vorzeichen nach Lam (Tab. 2.1)
  in `docs/architecture.md`; geprüfte Varianten mit `std::invalid_argument` (#60).
- `rtt-json` (neue header-only Bibliothek, Schicht Basis, nur nlohmann-json): `parse_strict` mit
  `StrictParseError { pointer, message }` als gemeinsamer strikter JSON-Leser von `rtt-io` und
  `rtt-coating` (ADR 0020, #68).
- `rtt-material`: AGF-Formeln 6 (Sellmeier 3, K₁ L₁ … K₄ L₄), 12 (Extended 2, a₀ … a₆; a₇ ≠ 0 bleibt
  ein Fehler) und 13 (Extended 3, a₀ … a₈); Reihenfolge belegt an N(d) und V(d) des freien
  NIKON-HIKARI-Katalogs. Neue Formeln `Extended2Coefficients` und `Extended3Coefficients`. CD- und
  TD-Datensätze dürfen auf Folgezeilen umbrechen, die mit einer Zahl beginnen; CD-Werte jenseits
  der von der Formel genutzten Positionen müssen 0 sein (Formeln 1, 2, 6, 12, 13). Testkatalog
  `tests/catalogs/nikon/nikon-hikari.agf`. Formeln 3, 4, 5, 7–11 bleiben unbelegt und ein klarer
  Fehler; die vollständige NIKON-HIKARI_201911.AGF lädt wegen doppelter Glasnamen noch nicht (#42).
- `rtt-trace`: Interaktionen wirken auf `prt` und `weight` (ADR 0021): Fresnel mit komplexen
  Indizes, ideale Elemente (Spiegel, AR, Strahlteiler, Polarisator, Retarder), Coatings aus dem
  Katalog (umgekehrter Stapel aus dem Substrat), Absorber und Volumenabsorption
  exp(−4πκt/λ). P ist leistungsnormiert (|P E|² = Leistungsanteil), `weight` ist die Leistung
  für eine unpolarisierte Quelle; Beispiel `tests/reference/m3/absorbing_ar_plate.rtt.json`. `apply_event` und `sequential_step` nehmen `EventMedia`
  (komplexe Indizes, Wellenlänge, Schichten); die Formen mit reellen Indizes bleiben (#61).
- `rtt-polar`: leistungsnormierte Grenzflächen-PRT-Matrizen (`interface.hpp`:
  `fresnel_prt`, `interface_prt`, `transmission_power_factors`, `tangential_invariant`) (#61).
- `rtt-compile`: `CompiledEvent::medium_beyond` (Medium hinter der Fläche, auch bei Reflect),
  `CompiledEvent::from_inside` (Strahl kommt aus dem Element, bei Coatings die Substratseite) und
  `CompiledSurface::ideal_axis` (Achse von Polarisator und Retarder global) (#61).
- `rtt-py`: Coatings und Polarisation in Python.
  - `rt.CoatingLibrary` (`add_catalog`, `in`), `rt.compile(system, materials=None,
    coatings=None)` und `rt.CoatingCatalogError`.
  - Neues Modul `rt.polar` mit `initial_directions`, `transverse_polarization`, `transmission`,
    `diattenuation`, `retardance` und `stokes` für eine verfolgte `RayBatch`, dazu
    Einzelmatrix-Funktionen (`prt_matrix`, `geometric_transform`, `diattenuation`, `retardance`,
    `physical_retardance`, `stokes`).
  - Docstrings von `prt` und `weight` nach ADR 0021.
  - Bitgleich gegen C++ mit fünf neuen Fällen.
  - Referenz `tests/reference/m3/polarizer_qwp.rtt.json` und Beispiel
    `examples/python/polarization.py` (#62).

### Geändert
- `rtt-analysis`: Spot-Gewichte enthalten jetzt die Fresnel-Verluste der Strahlen; gewichtete
  Statistiken unvergüteter Systeme verschieben sich leicht (Cooke-Triplett relativ bis 2,8·10⁻³)
  (#61).
- `rtt-compile`: `compile(system, materials)` ohne Coating-Bibliothek meldet für `CoatingRef` jetzt
  einen `CompileError` (bisher ignoriert); `rt.compile` nimmt seit #62 eine `CoatingLibrary`
  (#59, ADR 0019).
- `rtt-io`, `rtt-coating`: Doppelte JSON-Schlüssel auf jeder Ebene sind ein Fehler mit Pointer auf
  das zweite Vorkommen (bisher galt still der letzte Wert); betrifft `rtt validate`,
  `rt.System.from_json` und Coating-Kataloge (Nachtrag ADR 0008, #68).
- `rtt-polar`: Doku von `diattenuation()`: Für das leistungsnormierte P des Tracers ist sie die
  Leistungs-Diattenuation (ADR 0021, #62).

### Behoben
- `rtt-trace`: Beim Cooke-Triplet warfen Feld 1 (und `spot()`) „no start plane before surface 'L2.S1'“.
  Die Schranke für die Startebene konvergiert dort nur linear und war nach 50 Runden nicht fertig;
  jetzt gilt dann das Minimum über die Kreisscheibe mit dem Außenradius der Apertur bzw. dem
  Formgebiet der Fläche. Ein Fehler bleibt nur bei unbegrenzten Flächen ohne Apertur oder mit einer
  Apertur ab 1 km (#72).
- `rtt-io`: Eine Zahl jenseits des `double`-Bereichs (z. B. `1e400`) ergibt einen `ParseError`
  statt einer unübersetzten `nlohmann::json::out_of_range` (#68).

## [0.3.0] – M2 Materialien, Analyse

### Hinzugefügt
- `rtt-material`: Dispersionsformeln Schott, Sellmeier 1–5, Herzberger, Conrady und Cauchy als
  Funktions-Templates über `Real` mit `DispersionMaterial`; `TabulatedMaterial` (n und κ linear in
  λ); `MaterialLibrary::add` registriert Materialien unter einem Namen. Quellen in
  `docs/quellen.md` (#23).
- `rtt-material`: Import von Glaskatalogen im AGF-Format (`MaterialLibrary::add_catalog` für Datei
  oder Verzeichnis, Katalogname = Dateiname in Großbuchstaben, Verweise `KATALOG:NAME`);
  UTF-16LE mit BOM und ANSI/UTF-8, Fehler mit Datei und Zeile; Datensätze NM, CD, LD, TD, ED;
  Formelnummern 1 (Schott) und 2 (Sellmeier 1), andere mit klarem Fehler (#42). Testkatalog
  `tests/catalogs/schott.agf` mit N-BK7 und F2 aus dem freien SCHOTT-Katalog (#24).
- `rtt-analysis` (neue Bibliothek, Schicht Auswertung): Spot-Diagramm je Feld und Wellenlänge
  oder polychromatisch mit den Wellenlängengewichten des Modells (Schwerpunkt, RMS und GEO um
  Schwerpunkt und Hauptstrahl, Anteil vignettierter Strahlen) und Ray Fans tangential/sagittal
  relativ zum Hauptstrahl der Referenzwellenlänge; Datenobjekte, keine Plots. Referenzsystem
  `tests/reference/m2/paraboloid_stop.rtt.json` (#28).
- `rtt-compile`: `CompiledSystem::wavelength_weights()` mit den Wellenlängengewichten des
  Modells (#28).
- `rtt-paraxial`: Seidel-Summen S_I–S_V pro Fläche und für den Pfad aus paraxialem Rand- und
  Hauptstrahl (`seidel()`, `surface_seidel()`), mit Beiträgen von Konik und A4 und den chromatischen
  Termen C_L, C_T für ein wählbares Wellenlängenpaar; Quellen (Sasian, OPTI 517/518; arXiv:2203.02302,
  arXiv:2202.11378) in `docs/quellen.md` (#30).
- `rtt-analysis`: OPD gegen die Referenzsphäre (Zentrum am Hauptstrahl der Referenzwellenlänge,
  Radius bis zur paraxialen Austrittspupille) als Karte über die Pupille und als Fans, in Wellen
  bei der Referenzwellenlänge, mit RMS (ohne Piston) und PV; Vorzeichen nach Wyant & Creath
  (W > 0 voreilend) (#29).
- `rtt-material`: Ciddor-Luft `ciddor_air_index` und `AirMaterial` (`rtt/material/air.hpp`),
  Schott-Temperaturmodell `schott_delta_n_abs` (`rtt/material/thermal.hpp`), beide als
  Templates über `Real`; Quellen NIST (Gl. A21–A41) und SCHOTT TIE-19 in `docs/quellen.md` (#25).
- `rtt-trace`: `aim_ray` für Feldwerte außerhalb der Feldliste des Modells (für Feldverläufe,
  #31).
- `rtt-analysis`: Verzeichnung (real gegen paraxial in der Ebene des Bildflächenscheitels),
  sagittale und tangentiale Bildfeldwölbung (Nachbarstrahlen, O(δ²)), Farblängsfehler paraxial
  und real für ein Wellenlängenpaar, Farbquerfehler als Hauptstrahl je Wellenlänge; Feldverläufe
  über das relative Feld (#31).
- `rtt-py` (neue Bibliothek, Schicht Schnittstellen): Python-Paket `raytatouille` mit nanobind
  (ADR 0002) und scikit-build-core: `load`, `save`, `validate`, `System` (Name und Umgebung
  änderbar), `MaterialLibrary` (Kataloge, `index`), `compile`, `paraxial.first_order` und
  `trace` (`make_rays` mit allen Pupillenverteilungen, `trace` mit `threads`); `RayBatch`-Spalten
  als NumPy-Views ohne Kopie, die den Batch am Leben halten; C++-Fehler als Klassen aus
  `raytatouille.errors`; generierte Typ-Stubs. pytest mit Roundtrip Datei → Python → Datei und
  bitgleichem Vergleich mit C++ (`rtt_py_reference`, 1 und 4 Threads), mypy --strict; Presets
  `ci-linux-python` und `ci-windows-python`, CI-Job „Python“ (#32).
- ADR 0018: Python-Build- und Testwerkzeuge über `pyproject.toml` statt vcpkg; AGENTS.md
  Regel 4 verweist darauf (#32).
- `rtt-py`: Analysen aus Python (`raytatouille.analysis`): Spot-Diagramm, Ray Fans, OPD-Karte
  und -Fans, Farblängs- und Farbquerfehler, Verzeichnung und Bildfeldwölbung (Verlauf und
  einzelner Feldwert) sowie Seidel-Summen (`rt.paraxial.seidel`, Alias `rt.analysis.seidel`);
  `System` oder `CompiledSystem` als erstes Argument, Kurzform für Pupillenverteilungen wie
  `rays="hexapolar:12"`, `threads`; Ergebnislisten als read-only NumPy-Kopien; `rt.Field`;
  `AnalysisError`. Plots in `raytatouille.plot` mit optionalem matplotlib (`raytatouille[plot]`,
  Nachtrag zu ADR 0018). Bitgleiche Tests gegen C++ für jede Analyse, auch mit Vignettierung;
  Beispiel `examples/python/analysis.py` (#33).
- Abnahmetests M2 (`rtt-analysis`, `test_m2_acceptance.cpp`, Tag `[m2]`, #34):
  - Plankonvexlinse: S_I gleich Fit an reale Strahlen (Querabweichung und OPD), relativ 1e−3.
  - Achromat: realer Farblängsfehler gegen einen unabhängigen meridionalen Realstrahl; C_L gegen
    die paraxiale Farbe eines engen Wellenlängenpaars.
  - Cooke-Triplet als Golden-Test nach Sasian (OPTI 517 L20), verglichen mit ray-optics
    (erste Ordnung, Seidel je Fläche, 36 reale Strahlen) und Sasian.
  - Neues Referenzsystem `tests/reference/m2/cooke_triplet.rtt.json`, Glaskatalog
    `tests/catalogs/m2/schott.agf` (N-LAK9, N-SF5). Quellen in `docs/quellen.md`.

### Geändert
- `cmake`: `rtt_add_test` setzt unter MinGW den `PATH` für Testerkennung und Testläufe auf das
  Compiler-Verzeichnis (`DL_PATHS`); `ctest` startet GCC/MSYS2-Tests damit unabhängig von der
  Shell (Git Bash brachte eine falsche `libstdc++-6.dll` mit) (#47).
- `rtt-material`: `Material::index(wavelength_um, temperature_c, pressure_atm)` mit Luftdruck
  (ADR 0014, Punkt 3) und `wavelength_range_um()`; `rtt-compile` übergibt
  `environment.pressure_atm` und meldet Systemwellenlängen außerhalb des gültigen Bereichs eines
  benutzten Mediums als `CompileError` mit JSON-Pointer (#23).
- Dateiformat Schema 0.2.0: `Lens` und `Plate` mit N Flächen tragen ein Material je Segment
  (Segment i zwischen Fläche i und i + 1); `"material"` ist ein String (Kurzform für alle
  Segmente) oder ein Array mit N − 1 Einträgen. `rtt-model`: `Element::segment_materials`,
  `Element::segment_material(i)` und Validierung der Listenlänge mit JSON-Pointer; `rtt-io` liest
  und schreibt beide Formen bitgleich und migriert Dateien mit Schema 0.1; JSON-Schema angepasst;
  Referenzdateien mit `rtt format` neu geschrieben; ADR 0017; neues Beispiel
  `tests/reference/m2/achromat.rtt.json` (#26).
- `rtt-compile`: Kittglieder auf automatischen und expliziten Pfaden. Kurzform und Materialliste
  werden ausgewertet, Fehlerorte zeigen auf den Listeneintrag (`.../material/<i>`). Die Regel für
  brechende Events hängt nur vom Element ab: `Lens` (einheitlich oder gemischt) und `Plate` mit
  verschiedenen Segmentmaterialien folgen der Segmentregel (`Refract` an Fläche i: aus Segment
  i − 1 ins Segment i bzw. zurück, an der ersten und letzten Fläche ins Umgebungsmedium); eine
  `Plate` mit einheitlichem Material behält die Umschaltregel innen ↔ Umgebung an jeder Fläche.
  Der `CompileError` für eine `Lens` mit mehr als 2 Flächen auf dem automatischen Pfad entfällt;
  neu ist ein `CompileError` für eine einheitliche `Plate` mit mehr als 2 Flächen auf dem
  automatischen Pfad und, unter der Segmentregel, für eine innere Fläche zwischen verschiedenen
  Materialien, die ein Pfad von außen erreicht. Festgelegt: Elemente schachteln nicht (Element B
  betreten, während der Strahl in A ist, verlässt A) (#27).
- `rtt-model`: Test, dass ein `ThinElement` keine Materialliste annimmt (#26).
- `rtt-material`: `AIR` ist trockene Luft nach Ciddor (0 % Feuchte, 450 ppm CO₂) bei Temperatur
  und Druck des Systems statt n = 1; Katalogglas liefert die absolute Brechzahl (Daten relativ zu
  Luft bei T_ref und 1 atm, Wellenlänge in Luft) mit dn/dT nach SCHOTT TIE-19. Für Systeme in
  `AIR` sind EFL = 1/Φ und die EPD aus `image_fnumber` damit um den Faktor 1/n_Luft kleiner als
  in M1 (`docs/architecture.md`, Engine 1). Tests, deren Sollwerte für n_außen = 1 gelten, laufen
  in `VACUUM`; der Achromat-Test aus #41 bleibt in `AIR` und vergleicht f′ = n_Luft/Φ mit den
  Designwerten (#25).

### Behoben
- `rtt-trace`: Feldwerte werden bei der Referenzwellenlänge in Richtung bzw. Objektpunkt
  umgerechnet. Bisher geschah das je Strahl-Wellenlänge, sodass bei paraxialer Bildhöhe und bei
  Winkelfeldern mit endlichem Objekt jede Farbe einen anderen Feldpunkt bekam (Farbquerfehler,
  polychromatischer Spot). Fehler aus #8, gefunden in #31.
- `rtt-trace`: Aperturränder sind inklusiv mit `kApertureTolerance` = 1e−9 mm (= Aiming-Toleranz),
  außen wie innen. Randstrahlen mit |p| = 1, die exakt auf den Blendenrand gezielt werden, wurden
  bisher durch Rundung teils als vignettiert markiert (am Referenz-Singlet 7 von 36 je Feld).

## [0.2.0] – M1 Sequenzieller Kern

### Dokumentation
- `docs/quellen.md`: geprüfte Quellen mit Gleichungsnummern; Regel 6 in `AGENTS.md`, Reviewer und PR-Vorlage verweisen darauf.

### Geändert
- `rtt::math::Real` umfasst nur noch `double` (ADR 0015); der `float`-Test in `rtt-geom` entfällt.

### Hinzugefügt
- Anleitung `docs/agenten.md` und Review-Subagent `.claude/agents/physik-reviewer.md`.
- `rtt-geom`: `Shape<T>`-Schnittstelle, `Plane<T>` und `Conic<T>` (Sag, Gradient, Definitionsbereich)
  sowie analytischer Strahlschnitt mit Status `Missed` statt NaN (#3).
- `rtt-material`: Schnittstelle `Material` (komplexer Index n + iκ über Vakuumwellenlänge in µm
  und Temperatur in °C) und thread-sichere `MaterialLibrary` mit den Verweisen `VACUUM`, `AIR`
  (vorläufig n = 1), `CONST:<n>` und `CONST:<n>,<kappa>`; Kataloge folgen in M2 (#2).
- ADR 0014: Laufzeit-Interface für Materialien in `double`; Schnitt-ε = 1e−9 mm in `docs/architecture.md`.
- `rtt-geom`: gerade Asphäre `EvenAsphere<T>` und allgemeiner Newton-Schnitt für beliebige
  `Shape<T>` (Start an der Basis-Konik, Toleranz 1e−12 mm oder Rundungsgrenze des Schritts,
  höchstens 30 Schritte, Status `NoConvergence`, Iterationszahl im Ergebnis) (#4).
- `rtt-geom`: Quelle der Konik-Sag-Formel ist Forbes, Opt. Express 19(10) (2011), Gl. (2.1);
  Begründung der stabilen Form im Kommentar (#13).
- `rtt-compile`: `compile()` erzeugt ein unveränderliches `CompiledSystem` (globale Flächenlagen,
  Formen `Plane`/`Conic`, Medien je Wellenlänge, Pfade mit Medium vor und nach jedem Event,
  automatische Pfade) mit `CompileError` samt JSON-Pointer. Referenzsystem
  `tests/reference/m1/singlet_const.rtt.json` (#5).
- `rtt-trace`: `RayBatch` als Structure-of-Arrays (#5).
- ADR 0016: `CompiledSystem` in eigener Bibliothek `rtt-compile`; feste Reihenfolge
  `rtt-compile` < `rtt-paraxial` < `rtt-trace` in der Schicht Tracing.
- `rtt-trace`: sequenzieller Tracer `SequentialTracer` aus pfadunabhängigen Bausteinen
  `intersect_surface()`, `inside_aperture()` und `apply_event(ray, hit, kind)` (Schnitt im lokalen
  KS, Aperturen → `Vignetted` im sequenziellen Schritt, vektorielle Brechung und Reflexion, TIR → `Tir`,
  `Absorber` → `Absorbed`, M4-Events → `EventImpossible`, OPL mit Re(n)), `TraceStats` je Status,
  Parallelisierung mit oneTBB (`static_partitioner`, bitgleich zu einem Thread). Referenzsystem
  `tests/reference/m1/paraboloid_mirror.rtt.json` (#6).
- `rtt-compile`: `EvenAsphere` in `CompiledShape`; Mangin-Spiegel auf dem automatischen Pfad
  ergibt einen `CompileError` (#6).
- Quelle für Brechung und Reflexion in Vektorform: de Greve, „Reflections and Refractions in Ray
  Tracing“ (2006), Gl. (13), (22)–(25), Abschnitt 6; Born & Wolf §3.2.2 als ergänzende Angabe (#19).
- `rtt-paraxial`: y-nu-Trace (`trace_ray`) und Daten erster Ordnung (`first_order`): EFL, vordere
  und hintere Brennweite, BFL, FFL, Brenn- und Hauptpunkte, paraxiales Bild, Abbildungsmaßstab,
  Winkelvergrößerung, Eintritts- und Austrittspupille; Spiegel über Vorzeichenwechsel des Index;
  `ParaxialError` für nicht rotationssymmetrische Pfade (#7).
- `rtt-trace`: Strahlen aus der Systemdefinition (`sources.hpp`): Feldtypen Winkel
  (Tangens-Konvention), Objekthöhe und paraxiale Bildhöhe; Aperturtypen über die paraxiale
  Eintrittspupille; Ray Aiming real (Newton mit Dämpfung auf den Zielpunkt der Blende,
  1e−9 mm, sonst `NoConvergence`) oder paraxial; Startebene als ebene Welle vor dem System bei
  Objekt im Unendlichen; Verteilungen hexapolar, Gitter, Fächer x/y, Zufall mit Seed
  (plattformunabhängig), Einzelstrahl; `make_rays` und `aim_ray`. Referenzsystem
  `tests/reference/m1/two_lenses_stop_between.rtt.json` (Aperturtyp `stop_size`) (#8).
- Feldwinkel, Objekthöhe und normierte Pupillenkoordinaten als Konvention in
  `docs/architecture.md` (#8).

### Abhängigkeiten
- oneTBB (`tbb` in `vcpkg.json`, `libtbb-dev` für den apt-Weg), Apache-2.0, kompatibel mit MIT
  (ADR 0004, ADR 0010).

## [0.1.0] – M0 Fundament

### Hinzugefügt
- CMake-Projekt mit Presets, vcpkg-Manifest, Compiler-Warnungen, clang-tidy und Sanitizer-Optionen.
- `rtt-math`: Vektor- und Matrixtypen, Einheiten, `Isometry3` mit Pose-Konvention, `Real`-Concept.
- `rtt-model`: Datenmodell (System, Baugruppen, Elemente, Flächen-Stacks, Pfade, Parameter) und semantische Validierung.
- `rtt-io`: strikter JSON-Parser mit JSON-Pointer-Fehlern, kanonischer Writer, Dateiformat Schema 0.1.0.
- `rtt` CLI mit `validate` und `format`.
- JSON-Schema `schema/raytatouille.schema.json` mit pytest-Prüfung.
- Referenzsysteme: Singlet, Michelson-Interferometer, Feature-Tour.
- CI: Format, Schema, Linux GCC/Clang+tidy/Sanitizer, Windows MSVC.
- Dokumentation: Architektur, 13 ADRs, Agentenregeln.
