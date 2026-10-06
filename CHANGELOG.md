# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/); Versionen nach SemVer.

## [Unreleased]

### Hinzugefügt
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
