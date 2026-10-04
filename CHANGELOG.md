# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/); Versionen nach SemVer.

## [Unreleased]

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
