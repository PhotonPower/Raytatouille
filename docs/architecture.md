# Raytatouille – Architektur

> Maßgebliche Fassung. Änderungen nur per PR; Architekturentscheidungen zusätzlich als ADR unter `docs/adr/`.

## Ziel, Scope und Leitprinzipien

Raytatouille ist eine headless Raytracing-Engine für Optikdesign in C++20 mit Python-API, lizenziert unter MIT. Die GUI wird später in Python auf dieser API gebaut. Dieses Dokument ist die verbindliche Vorgabe für Programmieragenten: Es legt Module, Schnittstellen, Konventionen und Abnahmekriterien fest.

**Scope v1:** bildgebende und polarisationsoptische Systeme aus Linsen (sphärisch, asphärisch, Freiform), Spiegeln, Beschichtungen, Strahlteilern, Polarisatoren und Retardern. Tracing paraxial, sequenziell und multi-sequenziell (mehrere benannte Strahlpfade in einem Modell). Analysen, lokale Optimierung und Toleranzierung.

**Später (v2+):** als erste Priorität nicht-sequenzielles Tracing mit Streuung, danach Python-GUI, globale Optimierung, Beamlet-Wellenoptik, GPU-Backend.

**Nicht geplant:** rigorose EM-Solver (RCWA, FDTD), Beleuchtungsdesign, Struktur-FEA. Dafür sind Import-Schnittstellen (Jones-Tabellen, Coating-Dateien) vorgesehen.

**Leitprinzipien**

1. **Polarisation ist nativ.** Jeder reale Strahl trägt ab dem ersten Release eine 3×3-Polarisations-Raytracing-Matrix. Es gibt keinen separaten „Polarisationsmodus“.
2. **Flächen sind Stacks.** Eine Fläche = Form + Apertur + optionale Phase + Interaktion (Fresnel, Coating, Polarisator, Retarder, Strahlteiler). Neue Eigenschaften sind neue Layer, keine neuen Flächentypen.
3. **Pfade sind Daten.** Ein Strahlpfad ist eine explizite Liste von Flächen-Ereignissen (brechen, reflektieren, o-/e-Strahl, Beugungsordnung). Strahlteiler und Doppelbrechung erzeugen mehrere Pfade, keine Sonderfälle im Tracer.
4. **Hierarchie statt Flächenliste.** System → Baugruppe → Element → Fläche, jede Ebene mit eigenem Koordinatensystem und Pivot. Toleranzen wirken auf diese Hierarchie.
5. **Headless und deterministisch.** Keine UI- und keine Zufallsabhängigkeit ohne expliziten Seed. Gleiche Eingabe → bitgleiche Ausgabe auf derselben Plattform.
6. **Physik vor Features.** Kein Modul gilt als fertig ohne Referenztest gegen eine analytische Lösung oder ein publiziertes Ergebnis.
7. **Agentenfreundlich.** Kleine Bibliotheken mit schmalen Schnittstellen, damit Agenten parallel arbeiten können.

## Technologie-Entscheidungen

Details und Begründungen stehen in `docs/adr/`. Kurzfassung:

| ADR | Entscheidung |
| --- | --- |
| 0001 | Kern in C++20, CMake ≥ 3.25, Abhängigkeiten über vcpkg (Manifest `vcpkg.json`) |
| 0002 | Python-Bindings mit nanobind und scikit-build-core, NumPy-Views ohne Kopie |
| 0003 | Eigen (≥ 3.4) für 3-Vektoren und 3×3-Matrizen, `std::complex<double>`, durchgängig `double` |
| 0004 | Parallelität mit oneTBB, `static_partitioner`, deterministisch |
| 0005 | Strahlen als Structure-of-Arrays |
| 0006 | Physik-Bibliotheken als Templates über den Skalartyp (`rtt::math::Real`) |
| 0007 | Ableitungen in v1 per zentraler finiter Differenz |
| 0008 | Dateiformat JSON (nlohmann/json); strenger C++-Parser ist Referenz, JSON-Schema für externe Werkzeuge |
| 0009 | Exceptions nur an API-Grenzen, Strahlprobleme als Status-Flags |
| 0010 | Lizenz MIT, nur kompatible Abhängigkeiten |
| 0011 | Linux und Windows; getestet mit GCC 13, Clang 18, MSVC 2022; kein `-ffast-math` |
| 0012 | GUI später in Python (PySide6, pyvista/VTK, matplotlib) nur über die öffentliche API |
| 0013 | Keine GPU in v1 |

**Konventionen (verbindlich für alle Bibliotheken)**

- **Einheiten:** Längen in mm, Wellenlängen in µm, Winkel intern in rad (API-Parameter und Dateifelder in Grad tragen das Suffix `_deg`), Temperatur in °C, Druck in atm. Hilfsfunktionen in `rtt/math/units.hpp`.
- **Koordinaten:** rechtshändig, optische Achse +z, y ist die meridionale Richtung.
- **Radius:** positiv, wenn der Krümmungsmittelpunkt auf der +z-Seite des Flächenscheitels liegt (übliche Konvention in der Optik). Im Modell steht der Radius; der Tracer rechnet mit c = 1/R. Ebene Flächen sind `Plane`, nie Radius 0 oder unendlich.
- **Transformationen:** Translation, dann Rotation um den Pivot; Rotationen intrinsisch X → Y → Z, also R = Rx·Ry·Rz. p_parent = position + pivot + R·(p_child − pivot). Implementiert in `rtt::math::Isometry3::from_pose`.
- **Felder und Phase:** Ebene Wellen ~ exp(i(k·r − ωt)). Komplexer Index ñ = n + iκ mit κ ≥ 0 für Absorption. Positive Retardance heißt: die langsame Achse ist phasenverzögert.
- **Größen:** OPL in mm, OPD und Wellenfront in Wellen bei der Referenzwellenlänge, Strahlgewicht als Leistung (Quelle normiert auf 1).

## Systemübersicht

14 CMake-Bibliotheken bzw. -Programme in sieben Schichten. Jede darf nur Bibliotheken aus tieferen Schichten verwenden.

| Schicht | Bibliotheken | Status |
| --- | --- | --- |
| Schnittstellen | `rtt-py` (Python-API), `apps/rtt-cli` (Kommandozeile), `rtt-io` (Dateien) | `rtt-cli`, `rtt-io`: M0 |
| Workflows | `rtt-optim`, `rtt-tolerance` | M5, M7 |
| Auswertung | `rtt-analysis` | M2, M6 |
| Tracing | `rtt-paraxial`, `rtt-trace` | M1 |
| Modell | `rtt-model` | M0 |
| Physik | `rtt-geom`, `rtt-material`, `rtt-coating`, `rtt-polar` | M1–M4 |
| Basis | `rtt-math` | M0 |

Die vier Physik-Bibliotheken kennen weder Modell noch Tracer und können daher parallel gebaut werden. `rtt-model` beschreibt nur, `rtt-trace` rechnet, `rtt-analysis` wertet aus. CMake-Targets heißen `rtt_<name>` mit Alias `rtt::<name>`.

## Datenmodell (`rtt-model`, umgesetzt in M0)

Das Modell ist reine Datenstruktur ohne Tracing-Logik. Header: `libs/rtt-model/include/rtt/model/`.

**Hierarchie**

- **`System`:** Name, Umgebung (Temperatur, Druck, Umgebungsmedium), Objektraum (unendlich oder Abstand), Wellenlängen (genau eine Referenz), Systemapertur, Feldpunkte, Wurzel-Baugruppe, Pfade.
- **`Assembly`:** Name, `Pose`, Kinder (`Node` = Assembly oder Element).
- **`Element`:** physischer Körper mit `ElementKind`: `Lens` (≥ 2 Flächen, Material), `Mirror` (≥ 1 Fläche), `Plate` (≥ 2 plane Flächen, Material; Platten, Würfel, Prismen), `ThinElement` (1 Fläche ohne Dicke), `Stop` (1 Fläche mit Apertur, höchstens einer pro System), `Detector` (1 Fläche).
- **`Surface`:** eindeutige `SurfaceId` (lesbarer String wie `"L1.S1"`), `Pose` im Element, Layer-Stack.

**Layer-Stack einer Fläche**

| Layer | Anzahl | Umgesetzt (Schema 0.1) | Geplant |
| --- | --- | --- | --- |
| `shape.base` | genau 1 | `Plane`, `Conic`, `EvenAsphere` | Q-con, Q-bfs, Biconic, Toroid, Axicon, XY-Polynom, ungerade Asphäre |
| `shape.terms` | 0..n, additiv | `ZernikeSag` (Noll) | Grid-Sag, Zernike Fringe |
| `aperture` | 0..1 | kreisförmig (auch Ring), rechteckig, elliptisch | Polygon, Obskurationen |
| `phases` | 0..n, additiv | `LinearGrating`, `RadialPhase` | Zernike-Phase, Metaflächen-Tabellen |
| `interaction` | genau 1 | Fresnel (Default), ideal: Spiegel, AR, Strahlteiler, Polarisator, Retarder; Absorber; Coating-Referenz | Jones-/Mueller-Tabellen |
| `scatter` | 0..1 | – | Lambert, Gauß, ABg (M9) |

**Medien:** innerhalb eines Elements das Element-Material, außerhalb das Umgebungsmedium. Material ist ein Katalogverweis `"KATALOG:NAME"`.

**Pfade:** `Path` = Name + Liste von `Event { SurfaceId surface; EventKind kind; int order; }` oder `automatic` (alle Flächen in Baumreihenfolge). `EventKind`: `Refract`, `Reflect`, `Transmit`, `Ordinary`, `Extraordinary`, `Diffract` (mit `order`). Dieselbe Fläche darf mehrfach vorkommen (Doppeldurchgang, Interferometer). Ein `PathGenerator` (M4) erzeugt z. B. alle Zweifach-Reflexions-Ghosts.

**Parameter:** Jeder optimierbare Wert ist ein `Param { double value; bool variable; std::optional<std::string> pickup; }`. Pickup-Ausdrücke werden ab M5 ausgewertet. Multi-Konfigurationen kommen mit M5.

**Validierung:** `rtt::model::validate(system)` liefert `Diagnostic`s mit JSON-Pointer auf die betroffene Stelle (Wellenlängen, Apertur, Element-Regeln, eindeutige IDs und Namen, gültige Radien und Aperturen, Pfadverweise).

## Dateiformat (`rtt-io`, umgesetzt in M0)

Ein System ist eine Datei `*.rtt.json`. Die vollständige Struktur steht in `schema/raytatouille.schema.json`; Referenz ist aber der C++-Parser.

- **Strikt:** unbekannte Schlüssel, falsche Typen, andere Einheiten und inkompatible `schema_version` sind Fehler mit JSON-Pointer.
- **Versionierung:** `schema_version` SemVer; vor 1.0 müssen Major und Minor exakt passen. Jede Formatänderung erhöht die Version und bringt eine getestete Migration mit.
- **Kanonisch:** `rtt::io::to_json` schreibt immer dieselben Bytes: 2 Leerzeichen Einzug, LF, abschließender Zeilenumbruch, Standardwerte weggelassen, kleine Objekte aus Skalaren auf einer Zeile. Für jede kanonische Datei gilt `to_json(parse(text)) == text`. `rtt format` bringt Dateien in diese Form.
- **Parameter:** als Zahl (`51.68`) oder als Objekt (`{"value": 30.0, "variable": true}`).
- **Import v1 (M8):** AGF-Glaskataloge der Hersteller; ZMX-Dateien als reines Austauschformat (sequenzieller Teil), nicht unterstützte Flächentypen ergeben eine klare Fehlermeldung. **Export v1:** CSV für Analysedaten. **Später:** STEP-Export der Flächen, ISO-10110-Zeichnung.

Beispiel (gekürzt aus `tests/reference/m0/singlet.rtt.json`):

```json
{
  "schema_version": "0.1.0",
  "name": "Plano-convex singlet f = 100 mm",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [
    {"um": 0.4861},
    {"um": 0.5876, "reference": true},
    {"um": 0.6563}
  ],
  "aperture": {"type": "epd", "value": 20.0},
  "fields": {
    "points": [
      {},
      {"y": 5.0}
    ]
  },
  "root": {
    "type": "assembly",
    "name": "system",
    "children": [
      {
        "type": "lens",
        "name": "L1",
        "pose": {"position": [0.0, 0.0, 5.0]},
        "material": "SCHOTT:N-BK7",
        "surfaces": [
          {
            "id": "L1.S1",
            "shape": {
              "base": {"type": "conic", "radius": 51.68}
            },
            "aperture": {"type": "circular", "radius": 12.7}
          },
          {
            "id": "L1.S2",
            "pose": {"position": [0.0, 0.0, 4.0]},
            "aperture": {"type": "circular", "radius": 12.7}
          }
        ]
      }
    ]
  },
  "paths": [
    {"name": "main", "events": "auto"}
  ]
}
```

Alle Modellelemente zeigt `tests/reference/m0/feature_tour.rtt.json`, Mehrfachpfade `michelson.rtt.json`.

## Strahlen und Trace-Engines (ab M1)

Tracing läuft immer über ein unveränderliches `CompiledSystem`: Das Modell wird einmal kompiliert, dann tracen beliebig viele Threads lesend darauf.

**Kompilierung (`rtt-trace`)**

1. Baum flach machen: globale Transformation jeder Fläche, Pickups und Konfigurationen auflösen.
2. Materialien bei allen Wellenlängen, Temperatur und Druck auswerten (komplexer Index; bei Doppelbrechung n_o, n_e und optische Achse).
3. Beschichtungen vorberechnen: wahlweise Tabellen über Einfallswinkel je Wellenlänge oder exakt pro Strahl.
4. Flächentypen in `std::variant` auflösen, damit der Hot Path ohne virtuelle Aufrufe auskommt.
5. Ergebnis ist unveränderlich, thread-sicher lesbar, hält keine Zeiger ins Modell und hat einen Hash für Caching.

**Strahl-Batch (Structure-of-Arrays)**

| Feld | Typ je Strahl | Bedeutung |
| --- | --- | --- |
| `pos`, `dir` | 2 × 3 × `double` | Ort und Einheitsrichtung, global |
| `wl` | `uint16_t` | Index der Systemwellenlänge |
| `opl` | `double` | akkumulierter optischer Weg in mm |
| `prt` | 9 × `std::complex<double>` | akkumulierte 3×3-Polarisations-Raytracing-Matrix P |
| `weight` | `double` | Leistungsgewicht |
| `status` | `uint8_t` (enum) | `Alive`, `Missed`, `Vignetted`, `Tir`, `NoConvergence`, `Absorbed`, `EventImpossible` |
| `field`, `pupil` | `uint16_t`, 2 × `double` | Herkunft für Analysen |
| `last_surface` | `uint32_t` | letzte getroffene Fläche |

Der Eingangs-Jones-Vektor wird nicht gespeichert: E_out = P · E_in gilt für jeden Eingangszustand, damit sind polarisiert und unpolarisiert mit demselben Trace auswertbar.

**Engine 1, paraxial (`rtt-paraxial`):** y-nu-Trace für rotationssymmetrische Systeme: EFL, BFL, Hauptebenen, Pupillen, Vergrößerung, Seidel-Koeffizienten pro Fläche. v2: parabasal für dezentrierte Systeme.

**Engine 2, sequenziell/multi-sequenziell (`rtt-trace`):** pro Pfad und Event: lokal transformieren → Schnitt → Apertur → Normale → Phasengradient → Event (Brechung, Reflexion, Beugung, o-/e-Strahl) → Interaktion auf P und `weight` → global transformieren. Ist ein Event physikalisch unmöglich (z. B. TIR statt Brechung), endet der Strahl mit Status. Ray Aiming per Newton auf die Blendenmitte; Aperturtypen EPD, Bildraum-F/#, objektseitige NA, Blendengröße; Feldtypen Winkel, Objekthöhe, paraxiale Bildhöhe; Verteilungen hexapolar, Gitter, Zufall mit Seed, Fans, Einzelstrahl.

**Engine 3, nicht-sequenziell (erste Priorität nach v1.0):** BVH, Ray Splitting, Strahlbaum mit Abbruch nach Tiefe und Energie, Quellen, Detektoren. Gleiche Physik wie Engine 2; deshalb ist die Event-Ausführung als `apply_event(ray, hit, kind)` gekapselt und jedes Review ab M1 prüft, dass Schnitt, Interaktion und Event-Ausführung pfadunabhängig bleiben.

```cpp
class Tracer {
 public:
  virtual ~Tracer() = default;
  virtual TraceStats trace(const CompiledSystem& sys, PathId path, RayBatch& rays) const = 0;
};
```

## Physik-Module (ab M1)

Jede Physik-Bibliothek ist eine reine Funktionssammlung ohne Kenntnis von Modell oder Tracer und wird isoliert gegen Referenzwerte getestet.

**Geometrie (`rtt-geom`)**

```cpp
template <rtt::math::Real T>
class Shape {
 public:
  virtual ~Shape() = default;
  virtual T sag(T x, T y) const = 0;
  virtual std::pair<T, T> grad(T x, T y) const = 0;  // dz/dx, dz/dy
  virtual std::pair<T, T> base_conic() const = 0;    // (c, k) als Startwert
  virtual std::optional<T> max_radius() const = 0;
};
```

- Schnitt: analytisch mit der Basis-Konik, dann Newton auf F(t) = z(t) − sag(x(t), y(t)); Toleranz 1e−12 mm, höchstens 30 Iterationen, sonst `NoConvergence`. Die nächste Lösung mit t > ε in Ausbreitungsrichtung gilt; ε = 1e−9 mm (`kDefaultTMin` in `rtt-geom`), überall im Schnitt einheitlich zu verwenden.
- Normale aus dem Gradienten: n = (−∂z/∂x, −∂z/∂y, 1) normiert.
- Flächentypen v1: Ebene, Sphäre/Konik, gerade und ungerade Asphäre, Forbes Q-con und Q-bfs, Biconic, Toroid/Zylinder, Axicon, Zernike-Sag (Standard, Fringe), XY-Polynom, Grid-Sag (bikubisch).

**Materialien (`rtt-material`):** Sellmeier (1–5), Schott, Herzberger, Conrady, Cauchy, tabellierte n/k. AGF-Import; Herstellerkataloge werden nicht mitgeliefert. Luft nach Ciddor, Glasdaten relativ zu Luft, dn/dT nach Schott-Modell, CTE. Metalle als n + iκ. Einachsig doppelbrechend in v1, zweiachsig in v2.

**Beschichtungen (`rtt-coating`):** Transfermatrix-Methode (Abelès) für s und p mit komplexen Indizes, Ausgabe r_s, r_p, t_s, t_p. Dicke physikalisch oder als QWOT. Ideale und tabellierte Coatings. Invariante R + T + A = 1.

**Polarisation (`rtt-polar`):** 3D-Polarisations-Raytracing nach Chipman: P_q = O_out · J · O_inᵀ, O_in = [s, p_in, k_in], O_out = [s, p_out, k_out], s = k_in × n normiert, J = diag(a_s, a_p, 1). Bei senkrechtem Einfall s aus stabiler Referenzrichtung. P_gesamt = P_q · P_gesamt. Polarisator (Achse projiziert senkrecht zu k, Extinktionsverhältnis), Retarder (linear, zirkular, allgemeine Jones-Matrix). Depolarisierende Mueller-Elemente erst v2. Auswertung: Diattenuation, Retardance, Jones-Pupille, Stokes.

**Doppelbrechung:** o- und e-Wellenvektor an der Grenzfläche; der e-Strahl läuft entlang des Poynting-Vektors, der OPL über den Wellenvektor. Amplituden aus den Eigenpolarisationen (McClain/Chipman). Der Pfad wählt `Ordinary` oder `Extraordinary`.

**Phasen und Beugung:** lokale Gittergleichung k_t,aus = k_t,ein + m·∇φ; evaneszente Ordnungen enden mit Status. Effizienz pro Ordnung als Wert oder Tabelle.

**Strahlteiler:** kein eigenes Modul, sondern Kombination aus Element und Interaktion (Platte, Würfel, Pellicle, PBS). Der Pfad wählt `Reflect` oder `Refract`/`Transmit`.

**Literatur** (vor Nutzung prüfen): Chipman, Lam, Young: *Polarized Light and Optical Systems* (2018); Yun, McClain, Chipman: *Three-dimensional polarization ray-tracing calculus I/II*, Applied Optics 2011; Forbes: *Shape specification for axially symmetric optical surfaces*, Optics Express 2007; Macleod: *Thin-Film Optical Filters*.

## Analyse, Optimierung, Toleranzierung

Analysen liefern Datenobjekte, niemals Plots. Optimierung und Toleranzierung arbeiten nur über Parameter-Pfade und Analyse-Operanden.

| Gruppe | v1 | v2 |
| --- | --- | --- |
| Erste Ordnung | Systemdaten, Pupillen, Vergrößerung, Seidel pro Fläche | parabasal |
| Geometrisch | Spot (RMS, GEO, Schwerpunkt), Ray Fans, Footprint, Verzeichnung, Bildfeldwölbung, Farbfehler, Vignettierung | Bildsimulation |
| Wellenfront | OPD gegen Referenzsphäre, Zernike-Fit, Strehl-Näherung | – |
| Beugung | FFT-PSF, FFT-MTF, Through-Focus, Encircled Energy | Huygens-PSF, Beamlets, Faserkopplung |
| Polarisation | Jones-Pupille, Diattenuations- und Retardance-Karten, Transmission je Pfad | Stokes-Detektoren |
| Pfade | Transmission je Pfad, Ghost-Ranking | Interferogramm zweier Pfade |

**Optimierung (`rtt-optim`, M5):** Variablen sind alle `Param` mit `variable = true`. Merit-Funktion als gewichtete Summe von Operanden über beliebige Pfade und Konfigurationen; Generatoren für RMS-Spot und RMS-Wellenfront. Levenberg-Marquardt mit paralleler zentraler Differenz; Grenzen per Variablentransformation. v2: globale Suche, Glassubstitution, exakte Gradienten.

**Toleranzierung (`rtt-tolerance`, M7):** Toleranzen auf Fläche, Element und Baugruppe (Kippung um den jeweiligen Pivot), Kompensatoren, Sensitivität, inverse Sensitivität, Monte Carlo mit Seed. Eine Baugruppen-Toleranz bewegt alle Kinder starr; Doppeldurchgänge sehen dieselbe Störung, weil Pfade auf `SurfaceId` verweisen.

## Python-API und CLI

**Python (`rtt-py`, Paket `raytatouille`, ab M2):** dünne Spiegelung der C++-Typen, Ergebnisse als NumPy-Arrays, vollständige Typ-Stubs (`.pyi`). Plots in `raytatouille.plot`. Die spätere GUI (`raytatouille.gui`, PySide6) nutzt ausschließlich die öffentliche API; jede GUI-Aktion ist als Skriptbefehl reproduzierbar.

```python
import raytatouille as rt

sys = rt.load("doublet.rtt.json")
spot = rt.analysis.spot(sys, path="main", field=1, rays="hexapolar:12")
pupil = rt.analysis.jones_pupil(sys, path="main", field=0, grid=64)
mf = rt.optim.MeritFunction.rms_wavefront(sys, rings=6, arms=8)
mf.add("EFL", target=100.0, weight=1.0)
rt.optim.local(sys, mf, max_iter=200)
sys.save("doublet_opt.rtt.json")
```

**CLI (`rtt`):** vorhanden: `rtt validate`, `rtt format [--check]`, `rtt --version`. Geplant: `rtt trace`, `rtt analyze`, `rtt optimize`, `rtt import <zmx> <rtt.json>`.

## Validierung und Tests

Ein Feature ist erst fertig, wenn es mindestens einen Referenztest gegen eine unabhängige Lösung besteht. Pflichtfälle (Testsysteme unter `tests/reference/`):

| Fall | Erwartung | Toleranz |
| --- | --- | --- |
| Paraboloid, Objekt im Unendlichen, auf Achse | perfekter Fokus bei R/2 | RMS-Spot < 1e−9 mm, OPD < 1e−6 λ |
| Dicke Einzellinse | EFL, BFL, Hauptebenen nach Linsenmacherformel | relativ 1e−12 (paraxial) |
| Kleine Apertur, real vs. paraxial | Konvergenz gegen paraxial | Fehler ∝ NA² |
| Sphärische Aberration, Plankonvexlinse | Seidel = Fit an reale Strahlen bei kleiner NA | relativ 1e−3 |
| Fresnel an N-BK7 (n = 1,5168) | senkrecht R ≈ 4,22 %; Brewster arctan(n) mit R_p = 0 | 1e−12 |
| Totalreflexion | Grenzwinkel arcsin(1/n), Phasensprung s/p | 1e−12 |
| λ/4-AR MgF₂ (1,38) auf 1,52 | R ≈ 1,26 % bei Designwellenlänge | 1e−10 |
| Verlustfreier Schichtstapel | R + T = 1 für s und p | 1e−12 |
| Zwei Polarisatoren | Malus: T = cos²θ | 1e−12 |
| λ/4-Platte unter 45° | linear → zirkular | 1e−12 |
| λ/2-Platte unter θ | Drehung um 2θ | 1e−12 |
| Metallspiegel unter 45° | Retardance, Diattenuation aus Fresnel mit ñ | 1e−10 |
| Calcit-Platte, Achse 45° | Walk-off nach analytischer Formel | 1e−10 rad |
| Lineares Gitter | Gittergleichung, evaneszent → Status | 1e−12 |
| Beugungsbegrenztes System | PSF-Nullstelle 1,22 λ F/#, analytische MTF | relativ 1e−3 |
| Michelson | R + T = 1, OPD-Differenz = 2 × Armdifferenz | 1e−10 |

**Teststufen:** Unit-Tests mit Catch2 v3 in jeder Bibliothek; Property-Tests (RapidCheck, ab M1) für Energieerhaltung, |dir| = 1, P·k_ein = k_aus, Umkehrbarkeit, Invarianz bei starrer Bewegung; Golden-Tests mit Standardsystemen aus der Literatur; NSC-Tests (M9) gegen analytische Fälle (Lambert-Strahler, Etendue, Integrationskugel, Leistungsbilanz); Benchmarks mit Google Benchmark (Ziel ≥ 5 Mio. Strahl-Flächen-Schnitte pro Sekunde und Kern bei Sphären inkl. Polarisation, Regression > 10 % blockiert); Schema-Tests mit pytest.

**CI (`.github/workflows/ci.yml`):** clang-format, JSON-Schema, Linux GCC (Release), Linux Clang mit clang-tidy (Release), Linux Clang mit ASan/UBSan (Debug), Windows MSVC (Release); Warnungen sind Fehler.

## Roadmap und Abnahmekriterien

Nach M1 können bis zu drei Agenten parallel arbeiten. M4 und M6 brauchen M2 und M3; M5 braucht nur M2. Ein Meilenstein ist erreicht, wenn alle Kriterien in der CI grün sind.

| Meilenstein | Umfang | Abnahme |
| --- | --- | --- |
| **M0 Fundament** ✅ | CMake mit vcpkg, `rtt-math`, `rtt-model`, `rtt-io` (JSON mit Schema), `rtt` CLI, CI Linux und Windows | Roundtrip Datei → Modell → Datei bitgleich; CI grün |
| M1 Sequenzieller Kern | Ebene, Sphäre, Konik, gerade Asphäre; Schnitt, Brechung, Reflexion; `CompiledSystem`; Ray Aiming; paraxialer Trace | Paraboloid-Fokus, Linsenmacherformel, Konvergenz real → paraxial |
| M2 Materialien, Analyse | Dispersionsformeln, AGF-Import, Luft, dn/dT; Spot, Ray Fans, OPD, Verzeichnung, Seidel; Python-Bindings dafür | Seidel-Fall, Farblängsfehler eines Achromaten; Cooke-Triplet als Golden-Test |
| M3 Polarisation | P-Matrix, Fresnel mit ñ, Transfermatrix-Coatings, ideale Coatings, Polarisator, Retarder | Fresnel, Brewster, TIR, AR, Malus, λ/4, λ/2, Metallspiegel |
| M4 Multi-Path | explizite Pfade, Strahlteiler, einachsige Kristalle, Gitter, Ghost-Generator | Michelson-Bilanz, Calcit-Walk-off, Gittergleichung, Ghost-Ranking |
| M5 Optimierung | Parameter-Pfade, Operanden, Merit-Generatoren, LM, Grenzen, Pickups, Konfigurationen | Singlet auf EFL + minimalen RMS-Spot; reproduzierbar |
| M6 Beugung | Zernike-Fit, FFT-PSF/MTF, Encircled Energy, Jones-Pupille, Retardance-Karten | Airy-Nullstelle, analytische MTF, Jones-Pupille eines Faltspiegels |
| M7 Toleranzierung | Toleranzen auf allen Ebenen, Kompensatoren, Sensitivität, Monte Carlo | Baugruppen-Kippung bewegt alle Kinder; Monte Carlo reproduzierbar |
| M8 Release v1.0 | Forbes-, Zernike-, XY-, Grid-Flächen, ZMX-Import, CLI komplett, Doku, Beispiele | alle Referenzfälle grün; Benchmark-Ziel; ZMX-Testdateien korrekt importiert |
| M9 Nicht-sequenziell | BVH, Quellen, Detektoren, Ray Splitting, Streumodelle, Strahlpfad-Reports | analytische NSC-Testfälle |
| M10 Python-GUI | PySide6-Anwendung auf der öffentlichen API | – |
