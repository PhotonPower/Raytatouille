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
| 0006 | Physik-Bibliotheken als Templates über den Skalartyp (`rtt::math::Real`; nur `double` und später Dual-Zahlen, ADR 0015) |
| 0007 | Ableitungen in v1 per zentraler finiter Differenz |
| 0008 | Dateiformat JSON (nlohmann/json); strenger C++-Parser ist Referenz, JSON-Schema für externe Werkzeuge |
| 0009 | Exceptions nur an API-Grenzen, Strahlprobleme als Status-Flags |
| 0010 | Lizenz MIT, nur kompatible Abhängigkeiten |
| 0011 | Linux und Windows; getestet mit GCC 13, Clang 18, MSVC 2022; kein `-ffast-math` |
| 0012 | GUI später in Python (PySide6, pyvista/VTK, matplotlib) nur über die öffentliche API |
| 0013 | Keine GPU in v1 |
| 0014 | Laufzeit-Interface `Material` in `double`; Formelkerne als Templates |
| 0015 | `rtt::math::Real` umfasst nur `double` (später Dual-Zahlen) |
| 0016 | `CompiledSystem` in eigener Bibliothek `rtt-compile`; feste Reihenfolge in der Schicht Tracing |
| 0017 | Kittglieder: Material je Segment (`Lens`, `Plate`), Schema 0.2 |
| 0018 | Python-Build- und Testwerkzeuge (nanobind, scikit-build-core, numpy, pytest, mypy) über `pyproject.toml` statt vcpkg |

**Konventionen (verbindlich für alle Bibliotheken)**

- **Einheiten:** Längen in mm, Wellenlängen in µm, Winkel intern in rad (API-Parameter und Dateifelder in Grad tragen das Suffix `_deg`), Temperatur in °C, Druck in atm. Hilfsfunktionen in `rtt/math/units.hpp`.
- **Koordinaten:** rechtshändig, optische Achse +z, y ist die meridionale Richtung.
- **Radius:** positiv, wenn der Krümmungsmittelpunkt auf der +z-Seite des Flächenscheitels liegt (übliche Konvention in der Optik). Im Modell steht der Radius; der Tracer rechnet mit c = 1/R. Ebene Flächen sind `Plane`, nie Radius 0 oder unendlich.
- **Transformationen:** Translation, dann Rotation um den Pivot; Rotationen intrinsisch X → Y → Z, also R = Rx·Ry·Rz. p_parent = position + pivot + R·(p_child − pivot). Implementiert in `rtt::math::Isometry3::from_pose`.
- **Felder und Phase:** Ebene Wellen ~ exp(i(k·r − ωt)). Komplexer Index ñ = n + iκ mit κ ≥ 0 für Absorption. Positive Retardance heißt: die langsame Achse ist phasenverzögert.
- **Polarisation und Fresnel** (festgelegt für #56 und #58; Quelle S. J. Byrnes, arXiv:1603.02720v5, siehe `docs/quellen.md`):
  - **Tangentialinvariante:** ξ = n sin θ ist reell und an jeder Grenzfläche erhalten; bei absorbierendem Einfallsmedium ist ξ = Re(ñ_i) sin θ_i, wie bei der Brechung im Tracer. ξ ist die kanonische Eingabe von `rtt-polar` und `rtt-coating`; Varianten mit cos θ_i (`fresnel_at_angle`) gelten nur für nicht absorbierende Einfallsmedien.
  - **Normalkomponente:** q_j = ñ_j cos θ_j = √(ñ_j² − ξ²) mit Im q ≥ 0 und, falls Im q = 0, Re q ≥ 0. Die Welle klingt ab bzw. trägt Leistung vorwärts (Byrnes, Anhang D). Das gilt auch im Einfallsmedium. Totalreflexion: q_t = i√(ξ² − n_t²), |r_s| = |r_p| = 1.
  - **Vorzeichen (Convention A, Byrnes Gl. (6) und Anhang A; Anhang A, Fußnote 6, nennt Lehrbücher mit dieser Konvention):** r_s = (ñ_i cos θ_i − ñ_t cos θ_t)/(ñ_i cos θ_i + ñ_t cos θ_t), r_p = (ñ_t cos θ_i − ñ_i cos θ_t)/(ñ_t cos θ_i + ñ_i cos θ_t), t_s = 2ñ_i cos θ_i/(ñ_i cos θ_i + ñ_t cos θ_t), t_p = 2ñ_i cos θ_i/(ñ_t cos θ_i + ñ_i cos θ_t). Bei senkrechtem Einfall ist r_p = −r_s.
  - **Basis:** Mit der Flächennormale N ins Medium t sind die p-Einheitsvektoren in Byrnes Gl. (5) für die einfallende, die reflektierte und die transmittierte Welle genau p = k × s mit demselben s = k_in × N/|k_in × N|. Kehrt man N um, wechseln s und alle p gemeinsam das Vorzeichen, r und t bleiben gleich. r und t sind damit direkt a_s, a_p in J = diag(a_s, a_p, 1) der PRT-Matrix (Abschnitt Polarisation).
  - **Begründung über den idealen Leiter:** Er reflektiert E_r = −E_tan + E_norm, also E_r = (−I + 2NNᵀ)E_i. Daraus folgt (−I + 2NNᵀ)s = −s und (−I + 2NNᵀ)p_in = p_out für jeden Einfallswinkel, also a_s = −1 und a_p = +1. Auf transversalen Feldern ist P = −R mit der geometrischen Spiegelung R = I − 2NNᵀ: Beide Eigenpolarisationen bekommen dieselbe Phase π, es gibt weder Diattenuation noch Retardance. Fresnel mit ñ = 1 + iκ, κ → ∞, ergibt in Convention A genau r_s → −1 und r_p → +1 (Test in #56). Convention B gäbe r_p → −1, also scheinbar die Retardance π.
  - **Energie:** R = |r|², T_s = |t_s|² Re(q_t)/Re(q_i), T_p = |t_p|² Re(ñ_t cos θ_t*)/Re(ñ_i cos θ_i*) (Byrnes Gl. (21)–(23), mit Konjugation). Bei verlustfreiem Einfallsmedium gilt R + T = 1, auch wenn das zweite Medium absorbiert. Bei absorbierendem Einfallsmedium kann R + T ≠ 1 sein; das ist kein Fehler (Byrnes Abschnitt 4.4.1, Anhang B).
  - **Phasendifferenz (#56):** bei Reflexion Δ = arg(−r_p/r_s) (Ellipsometrie, Byrnes Gl. (16); 0 bei senkrechtem Einfall), bei Transmission Δ = arg(t_p/t_s), jeweils in (−π, π]. Δ > 0 heißt p verzögert; bei exp(−iωt) bedeutet ein größeres arg eine Verzögerung. Die Retardance mit langsamer Achse aus P folgt in #57.
  - **Fehlerbehandlung:** Die Kerne sind `noexcept`, mit dokumentierten Vorbedingungen und `assert` im Debug-Build. An der API-Grenze gibt es geprüfte Varianten (`fresnel_checked`), die `std::invalid_argument` werfen.
- **Größen:** OPL in mm, OPD und Wellenfront in Wellen bei der Referenzwellenlänge, Strahlgewicht als Leistung (Quelle normiert auf 1).
- **OPD-Vorzeichen** (festgelegt in #29): W = OPL_ref − OPL_Strahl, beide bis zur Referenzsphäre gemessen; **W > 0 heißt voreilend**, die Welle krümmt sich stärker als die Referenzsphäre (Wyant & Creath, Abschnitt I, siehe `docs/quellen.md`). Ein zur Austrittspupille hin verschobener Detektor sieht W < 0; unterkorrigierte sphärische Aberration mit Referenz im paraxialen Fokus ergibt W > 0 am Rand.
- **Feldwinkel und Pupille** (festgelegt in #8): Ein Feldwinkel (θx, θy) in Grad gibt die Richtung des Hauptstrahls im Objektraum d ∝ (tan θx, tan θy, 1) an; θy > 0 heißt, der Strahl steigt in +y (der Objektpunkt liegt bei −y). Eine Objekthöhe (x, y) in mm ist der Objektpunkt (x, y, −Objektabstand), positiv = +y. Eine paraxiale Bildhöhe wird über den paraxialen Hauptstrahl in Winkel bzw. Objekthöhe umgerechnet. Feldwerte werden bei der Referenzwellenlänge in Richtung bzw. Objektpunkt umgerechnet (paraxiale Daten der Referenzwellenlänge); ein Feldpunkt ist damit für alle Wellenlängen derselbe physikalische Punkt, während Pupille, Blendenziel und Startebene zur Wellenlänge des Strahls gehören (#31). Normierte Pupillenkoordinaten (px, py): der Einheitskreis ist der Rand der paraxialen Eintrittspupille, +y meridional; Ray Aiming zielt auf (px, py)·R_s auf der Blendenfläche (R_s paraxialer Blendenradius zur Eintrittspupille, vorzeichenbehaftet bei invertierter Abbildung). Ein Feldwinkel bei endlichem Objekt legt den Objektpunkt auf den Hauptstrahl durch die EP-Mitte. Bei Objekt im Unendlichen starten die Strahlen eines Feldes auf einer Ebene senkrecht zur Feldrichtung (ebene Welle, OPL 0), so platziert, dass das Bündel bis 3 EP-Radien um den Hauptstrahl mindestens 1 mm vor der EP und vor jeder erreichbaren Fläche des Pfades startet; bei endlichem Objekt im Objektpunkt.

## Systemübersicht

15 CMake-Bibliotheken bzw. -Programme in sieben Schichten. Jede darf nur Bibliotheken aus tieferen Schichten verwenden. Ausnahmen: die Schicht Tracing mit der festen Reihenfolge `rtt-compile` < `rtt-paraxial` < `rtt-trace`, dort darf eine Bibliothek zusätzlich die in dieser Reihenfolge vor ihr stehenden verwenden (ADR 0016); in der Schicht Schnittstellen dürfen `rtt-py` und `apps/rtt-cli` `rtt-io` verwenden (Dateien lesen und schreiben).

| Schicht | Bibliotheken | Status |
| --- | --- | --- |
| Schnittstellen | `rtt-py` (Python-API), `apps/rtt-cli` (Kommandozeile), `rtt-io` (Dateien) | `rtt-cli`, `rtt-io`: M0; `rtt-py`: M2 |
| Workflows | `rtt-optim`, `rtt-tolerance` | M5, M7 |
| Auswertung | `rtt-analysis` | M2, M6 |
| Tracing | `rtt-compile` < `rtt-paraxial` < `rtt-trace` | M1 |
| Modell | `rtt-model` | M0 |
| Physik | `rtt-geom`, `rtt-material`, `rtt-coating`, `rtt-polar` | M1–M4 |
| Basis | `rtt-math` | M0 |

Die vier Physik-Bibliotheken kennen weder Modell noch Tracer und können daher parallel gebaut werden. `rtt-model` beschreibt nur, `rtt-compile` macht daraus ein unveränderliches `CompiledSystem`, `rtt-paraxial` und `rtt-trace` rechnen, `rtt-analysis` wertet aus. CMake-Targets heißen `rtt_<name>` mit Alias `rtt::<name>`.

## Datenmodell (`rtt-model`, umgesetzt in M0)

Das Modell ist reine Datenstruktur ohne Tracing-Logik. Header: `libs/rtt-model/include/rtt/model/`.

**Hierarchie**

- **`System`:** Name, Umgebung (Temperatur, Druck, Umgebungsmedium), Objektraum (unendlich oder Abstand), Wellenlängen (genau eine Referenz), Systemapertur, Feldpunkte, Wurzel-Baugruppe, Pfade.
- **`Assembly`:** Name, `Pose`, Kinder (`Node` = Assembly oder Element).
- **`Element`:** physischer Körper mit `ElementKind`: `Lens` (≥ 2 Flächen, Material je Segment), `Mirror` (≥ 1 Fläche, optional ein Substratmaterial), `Plate` (≥ 2 plane Flächen, Material je Segment; Platten, Würfel, Prismen), `ThinElement` (1 Fläche ohne Dicke), `Stop` (1 Fläche mit Apertur, höchstens einer pro System), `Detector` (1 Fläche).
- **Segmente (ADR 0017):** Ein Element mit N Flächen hat N − 1 Segmente; Segment i liegt zwischen Fläche i und i + 1 (nullbasiert). `Lens` und `Plate` geben entweder ein Material für alle Segmente an (`material`, Kurzform) oder eine Liste mit genau N − 1 Einträgen (`segment_materials`, Kittglieder); `Element::segment_material(i)` liefert das Material von Segment i in beiden Formen.
- **`Surface`:** eindeutige `SurfaceId` (lesbarer String wie `"L1.S1"`), `Pose` im Element, Layer-Stack.

**Layer-Stack einer Fläche**

| Layer | Anzahl | Umgesetzt (Schema 0.2) | Geplant |
| --- | --- | --- | --- |
| `shape.base` | genau 1 | `Plane`, `Conic`, `EvenAsphere` | Q-con, Q-bfs, Biconic, Toroid, Axicon, XY-Polynom, ungerade Asphäre |
| `shape.terms` | 0..n, additiv | `ZernikeSag` (Noll) | Grid-Sag, Zernike Fringe |
| `aperture` | 0..1 | kreisförmig (auch Ring), rechteckig, elliptisch | Polygon, Obskurationen |
| `phases` | 0..n, additiv | `LinearGrating`, `RadialPhase` | Zernike-Phase, Metaflächen-Tabellen |
| `interaction` | genau 1 | Fresnel (Default), ideal: Spiegel, AR, Strahlteiler, Polarisator, Retarder; Absorber; Coating-Referenz | Jones-/Mueller-Tabellen |
| `scatter` | 0..1 | – | Lambert, Gauß, ABg (M9) |

**Medien:** innerhalb eines Elements das Material des jeweiligen Segments (bei der Kurzform überall dasselbe), außerhalb das Umgebungsmedium. Material ist ein Katalogverweis `"KATALOG:NAME"`. Die Kittfläche i (0 < i < N − 1) ist genau ein Übergang von Segment i − 1 zu Segment i.

**Pfade:** `Path` = Name + Liste von `Event { SurfaceId surface; EventKind kind; int order; }` oder `automatic` (alle Flächen in Baumreihenfolge). `EventKind`: `Refract`, `Reflect`, `Transmit`, `Ordinary`, `Extraordinary`, `Diffract` (mit `order`). Dieselbe Fläche darf mehrfach vorkommen (Doppeldurchgang, Interferometer). Ein `PathGenerator` (M4) erzeugt z. B. alle Zweifach-Reflexions-Ghosts.

**Parameter:** Jeder optimierbare Wert ist ein `Param { double value; bool variable; std::optional<std::string> pickup; }`. Pickup-Ausdrücke werden ab M5 ausgewertet. Multi-Konfigurationen kommen mit M5.

**Validierung:** `rtt::model::validate(system)` liefert `Diagnostic`s mit JSON-Pointer auf die betroffene Stelle (Wellenlängen, Apertur, Element-Regeln inkl. Länge der Materialliste, eindeutige IDs und Namen, gültige Radien und Aperturen, Pfadverweise).

## Dateiformat (`rtt-io`, umgesetzt in M0)

Ein System ist eine Datei `*.rtt.json`. Die vollständige Struktur steht in `schema/raytatouille.schema.json`; Referenz ist aber der C++-Parser.

- **Strikt:** unbekannte Schlüssel, falsche Typen, andere Einheiten und inkompatible `schema_version` sind Fehler mit JSON-Pointer.
- **Versionierung:** `schema_version` SemVer; aktuell 0.2.0. Vor 1.0 müssen Major und Minor exakt passen oder zu einer älteren Version gehören, die `rtt-io` migriert (derzeit 0.1 → 0.2; geschrieben wird immer die aktuelle Version). Jede Formatänderung erhöht die Version und bringt eine getestete Migration mit. Das JSON-Schema beschreibt nur die aktuelle Version.
- **Material:** `"material"` ist ein String (ein Material für alle Segmente) oder ein Array mit einem Eintrag je Segment, z. B. `"material": ["SCHOTT:N-BK7", "SCHOTT:F2"]` (ADR 0017). Die gelesene Form wird unverändert geschrieben.
- **Kanonisch:** `rtt::io::to_json` schreibt immer dieselben Bytes: 2 Leerzeichen Einzug, LF, abschließender Zeilenumbruch, Standardwerte weggelassen, kleine Objekte aus Skalaren auf einer Zeile. Für jede kanonische Datei gilt `to_json(parse(text)) == text`. `rtt format` bringt Dateien in diese Form.
- **Parameter:** als Zahl (`51.68`) oder als Objekt (`{"value": 30.0, "variable": true}`).
- **Import v1:** AGF-Glaskataloge der Hersteller (M2, `rtt-material`); ZMX-Dateien (M8) als reines Austauschformat (sequenzieller Teil), nicht unterstützte Flächentypen ergeben eine klare Fehlermeldung. **Export v1:** CSV für Analysedaten. **Später:** STEP-Export der Flächen, ISO-10110-Zeichnung.

Beispiel (gekürzt aus `tests/reference/m0/singlet.rtt.json`):

```json
{
  "schema_version": "0.2.0",
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

Alle Modellelemente zeigt `tests/reference/m0/feature_tour.rtt.json`, Mehrfachpfade `michelson.rtt.json`, ein Kittglied mit Material je Segment `tests/reference/m2/achromat.rtt.json`.

## Strahlen und Trace-Engines (ab M1)

Tracing läuft immer über ein unveränderliches `CompiledSystem`: Das Modell wird einmal kompiliert, dann tracen beliebig viele Threads lesend darauf.

**Kompilierung (`rtt-compile`, ADR 0016)**

1. Baum flach machen: globale Transformation jeder Fläche, Pickups und Konfigurationen auflösen.
2. Materialien bei allen Wellenlängen, Temperatur und Druck auswerten (komplexer Index; bei Doppelbrechung n_o, n_e und optische Achse).
3. Beschichtungen vorberechnen: wahlweise Tabellen über Einfallswinkel je Wellenlänge oder exakt pro Strahl.
4. Flächentypen in `std::variant` auflösen, damit der Hot Path ohne virtuelle Aufrufe auskommt.
5. Ergebnis ist unveränderlich, thread-sicher lesbar, hält keine Zeiger ins Modell und hat einen Hash für Caching.

**Medien entlang eines Pfads** (festgelegt in #5, für Segmente erweitert in #27): Der Strahl startet im Umgebungsmedium. `Refract`, `Ordinary` und `Extraordinary` an einer Fläche eines Elements mit Material wechseln zwischen dem Inneren dieses Elements und der Umgebung (bei mehreren Segmenten siehe unten); `Reflect`, `Transmit`, `Diffract` und alle Events an Elementen ohne Material behalten das Medium. Der automatische Pfad besucht alle Flächen in Baumreihenfolge: `Lens`/`Plate` → `Refract`, `Mirror` → `Reflect`, `Stop`/`Detector`/`ThinElement` → `Transmit`. Kittglieder (ADR 0017, festgelegt in #27): Eine `Lens`/`Plate` mit N Flächen hat N − 1 Segmente; Segment i liegt zwischen Fläche i und i + 1, die Kurzform `material` gilt für alle Segmente. Ein `Mirror` mit Substratmaterial und mehr als einer Fläche (Mangin-Spiegel) braucht einen expliziten Pfad (`Refract`, `Reflect`, `Refract`); auf dem automatischen Pfad ist er ein Kompilierfehler (#6). Stand M1: Formen `Plane`, `Conic` und `EvenAsphere`; Pickups werden noch nicht ausgewertet; der Hash folgt, wenn Caching gebraucht wird.

Für `Refract`, `Ordinary` und `Extraordinary` hängt die Regel nur von der Art des Elements und seinen Materialien ab, nie von der Art des Pfads; dieselben Events ergeben also auf automatischen und expliziten Pfaden dieselben Medien. „Einheitlich“ heißt Kurzform oder Liste mit lauter gleichen Einträgen.

| Element | automatischer Pfad | expliziter Pfad |
| --- | --- | --- |
| `Lens`, einheitlich oder gemischt | Segmentregel | Segmentregel |
| `Plate`, einheitlich (Prisma, Würfel) | 2 Flächen: Umschaltregel; mehr als 2: Kompilierfehler | Umschaltregel |
| `Plate`, gemischt (z. B. gekitteter Würfel) | Segmentregel | Segmentregel |
| `Mirror` mit Substrat | 1 Fläche: Umschaltregel; mehr: Kompilierfehler (#6) | Umschaltregel |

*Umschaltregel* (Regel aus #5): An jeder Fläche wechselt der Strahl zwischen dem Inneren und der Umgebung; damit lässt sich ein Prisma durch jede Fläche betreten und verlassen. Bei 2 Flächen ist sie gleich der Segmentregel. Eine einheitliche `Plate` mit mehr als 2 Flächen braucht einen expliziten Pfad, weil die Baumreihenfolge ihrer Flächen keinen sinnvollen Strahlweg ergibt.

*Segmentregel*: An Fläche i führt das Event aus Segment i − 1 ins Segment i (an der letzten Fläche ins Umgebungsmedium) und aus Segment i ins Segment i − 1 (an der ersten Fläche ins Umgebungsmedium); aus einem anderen Segment desselben Elements, das die Fläche nicht begrenzt, ins Umgebungsmedium. Von außen führt die erste Fläche ins Segment 0, die letzte ins Segment N − 2 und eine innere Fläche i in das Material der Segmente i − 1 und i, wenn beide gleich sind; sonst ist das Event ein Kompilierfehler (Seite mehrdeutig, z. B. die Kittfläche eines Achromaten von außen). Eine einheitliche `Lens` mit 3 Flächen ergibt damit Luft → Glas → Glas → Luft, ein Achromat Luft → Glas A → Glas B → Luft.

*Mehrere Elemente:* Elemente schachteln nicht. Ein brechendes Event an Element B, während der Strahl in Element A ist, verlässt A und betritt B wie aus der Umgebung (Medium vor dem Event: das von A); das Verlassen von B führt ins Umgebungsmedium, nicht zurück nach A. Folgen Austrittsfläche von A und Eintrittsfläche von B direkt aufeinander, läuft der Strahl über das Umgebungsmedium, auch bei Nullabstand; Kittglieder werden ausschließlich als ein Element mit Segmenten beschrieben.

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

**Engine 1, paraxial (`rtt-paraxial`):** y-nu-Trace für rotationssymmetrische Systeme: EFL, BFL, Hauptebenen, Pupillen, Vergrößerung, Seidel-Koeffizienten pro Fläche (M2, umgesetzt in #30: S_I–S_V je Event und als Summe aus Randstrahl (axialer Objektpunkt, Rand der paraxialen Eintrittspupille bei +y) und Hauptstrahl des maximalen Feldes (radialer Feldwert, Richtung nach der Feldkonvention aus #8 (Winkel: Hauptstrahl steigt in +y, bei endlichem Objekt liegt der Objektpunkt dann bei −y; Höhen: Punkt bei +y), durch die EP-Mitte; S_II und S_V wechseln mit dieser Richtung das Vorzeichen), Beiträge von Konik und A4, chromatische Terme C_L und C_T für ein wählbares Wellenlängenpaar; Normierung W040 = S_I/8 usw. in mm, Quellen in docs/quellen.md). v2: parabasal für dezentrierte Systeme. Konventionen (festgelegt in #7): u = dy/dz global; der Index trägt ein Vorzeichen (positiv bei Ausbreitung in +z) und wechselt es an jedem Spiegel, sodass n·u der optische Richtungskosinus ist und Brechung n′u′ = nu − yφ mit φ = c(n′ − n) sowie Transfer mit globalem Δz auch nach Spiegeln gelten. Brennweiten sind positiv für sammelnde Systeme (Hohlspiegel: EFL = |R|/2); alle Lagen als globale z-Koordinaten, FFL/BFL ab dem ersten bzw. letzten strahlverändernden Scheitel. Nicht rotationssymmetrische Pfade (dezentriert, gekippt, Phasenschichten, Beugung, Doppelbrechung) ergeben einen `ParaxialError`. Da Brechzahlen absolut sind (`AIR` = Ciddor-Luft, n ≈ 1,00027), ist EFL = 1/Φ im absoluten Sinn und `rear_focal_length` = |n′|/Φ der physikalische Abstand; Programme, die relativ zu Luft rechnen (z. B. OpticStudio), geben für ein System in Luft eine um den Faktor n_Luft größere EFL an. Dasselbe gilt für die Systemapertur `image_fnumber`: EPD = EFL/F# mit der absoluten EFL, also bei gleichem F# eine um den Faktor n_Luft kleinere EPD als in OpticStudio. Literaturwerte mit Luft = 1 vergleicht man daher mit `VACUUM` als Umgebung oder rechnet um (#25).

**Engine 2, sequenziell/multi-sequenziell (`rtt-trace`):** pro Pfad und Event: lokal transformieren → Schnitt → Apertur → Normale → Phasengradient → Event (Brechung, Reflexion, Beugung, o-/e-Strahl) → Interaktion auf P und `weight` → global transformieren. Ist ein Event physikalisch unmöglich (z. B. TIR statt Brechung), endet der Strahl mit Status. Aperturränder sind inklusiv mit der Toleranz `kApertureTolerance` = 1e−9 mm (gleich der Aiming-Toleranz), am Außen- wie am Innenrand; ein auf den Blendenrand gezielter Strahl wird damit nicht durch Rundung vignettiert (#50). Ray Aiming per Newton auf die Blendenmitte; Aperturtypen EPD, Bildraum-F/#, objektseitige NA, Blendengröße; Feldtypen Winkel, Objekthöhe, paraxiale Bildhöhe; Verteilungen hexapolar, Gitter, Zufall mit Seed, Fans, Einzelstrahl.

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

- Schnitt: analytisch mit der Basis-Konik, dann Newton auf F(t) = z(t) − sag(x(t), y(t)); konvergiert bei |F| ≤ 1e−12 mm oder |Δt| ≤ 8·ε_double·(1 + |t| + ‖o‖) (Rundungsgrenze bei fernem Strahlursprung), höchstens 30 Iterationen, sonst `NoConvergence`. Die nächste Lösung mit t > ε in Ausbreitungsrichtung gilt; ε = 1e−9 mm (`kDefaultTMin` in `rtt-geom`), überall im Schnitt einheitlich zu verwenden.
- Normale aus dem Gradienten: n = (−∂z/∂x, −∂z/∂y, 1) normiert.
- Flächentypen v1: Ebene, Sphäre/Konik, gerade und ungerade Asphäre, Forbes Q-con und Q-bfs, Biconic, Toroid/Zylinder, Axicon, Zernike-Sag (Standard, Fringe), XY-Polynom, Grid-Sag (bikubisch).

**Materialien (`rtt-material`):** Sellmeier (1–5), Schott, Herzberger, Conrady, Cauchy, tabellierte n/k. AGF-Import; Herstellerkataloge werden nicht mitgeliefert. Luft nach Ciddor, Glasdaten relativ zu Luft, dn/dT nach Schott-Modell, CTE. Umgesetzt in #25: `AIR` ist trockene Luft (0 % Feuchte, 450 ppm CO₂) nach Ciddor bei Systemtemperatur und -druck; Katalogglas wird mit λ_rel = λ_vac/n_Luft(λ_vac, T_ref, 1 atm) ausgewertet und als n_abs = n_rel·n_Luft + Δn_abs(T − T_ref) (SCHOTT TIE-19) absolut geliefert; `VACUUM` bleibt 1, `CONST:` und eigene Formeln bleiben absolut und temperaturunabhängig. Metalle als n + iκ. Einachsig doppelbrechend in v1, zweiachsig in v2.

**Beschichtungen (`rtt-coating`):** Transfermatrix-Methode (Abelès) für s und p mit komplexen Indizes, Ausgabe r_s, r_p, t_s, t_p. Dicke physikalisch oder als QWOT. Ideale und tabellierte Coatings. Invariante R + T + A = 1. Umgesetzt in #58: Transfermatrix in der Grenzflächen-Formulierung nach Byrnes (Gl. 8–15, nicht die Admittanzen nach Macleod, weil diese für p meist die andere Vorzeichenkonvention haben); Schichten vom Einfallsmedium zum Substrat, Dicke in µm (`_um`) oder als QWOT d = q·λ₀/(4 Re n(λ₀)); Kerne als Templates über `Real`, noexcept, Prüfung an der API-Grenze mit `check_stack`. Ideale Coatings mit Reflexionsgrad R ∈ [0, 1] ohne Retardance (r_s = −√R, r_p = +√R, Transmission verlustfrei). Bewusste Modellannahme (festgelegt in #58): Unter Totalreflexion verhält sich ein ideales Coating mit R < 1 wie die unbeschichtete Grenzfläche, einschließlich der TIR-Retardance; der ideale Spiegel (R = 1) bleibt auch dann bei r_s = −1, r_p = +1. Tabellierte Coatings über Einfallswinkel (Gitter muss [0, π/2] abdecken) und Wellenlänge, bilinear in Real- und Imaginärteil, Auswertung noexcept; der Wellenlängenbereich wird bei `compile()` geprüft (#61).

**Polarisation (`rtt-polar`):** 3D-Polarisations-Raytracing nach Chipman: P_q = O_out · J · O_inᵀ, O_in = [s, p_in, k_in], O_out = [s, p_out, k_out], s = k_in × n normiert, J = diag(a_s, a_p, 1). Bei senkrechtem Einfall s aus stabiler Referenzrichtung. P_gesamt = P_q · P_gesamt. Polarisator (Achse projiziert senkrecht zu k, Extinktionsverhältnis), Retarder (linear, zirkular, allgemeine Jones-Matrix). Depolarisierende Mueller-Elemente erst v2. Auswertung: Diattenuation, Retardance, Jones-Pupille, Stokes.

**Doppelbrechung:** o- und e-Wellenvektor an der Grenzfläche; der e-Strahl läuft entlang des Poynting-Vektors, der OPL über den Wellenvektor. Amplituden aus den Eigenpolarisationen (McClain/Chipman). Der Pfad wählt `Ordinary` oder `Extraordinary`.

**Phasen und Beugung:** lokale Gittergleichung k_t,aus = k_t,ein + m·∇φ; evaneszente Ordnungen enden mit Status. Effizienz pro Ordnung als Wert oder Tabelle.

**Strahlteiler:** kein eigenes Modul, sondern Kombination aus Element und Interaktion (Platte, Würfel, Pellicle, PBS). Der Pfad wählt `Reflect` oder `Refract`/`Transmit`.

**Literatur** (vor Nutzung prüfen): Chipman, Lam, Young: *Polarized Light and Optical Systems* (2018); Yun, McClain, Chipman: *Three-dimensional polarization ray-tracing calculus I/II*, Applied Optics 2011; Forbes: *Shape specification for axially symmetric optical surfaces*, Optics Express 2007; Macleod: *Thin-Film Optical Filters*.

## Analyse, Optimierung, Toleranzierung

Analysen liefern Datenobjekte, niemals Plots. Optimierung und Toleranzierung arbeiten nur über Parameter-Pfade und Analyse-Operanden.

**Spot und Ray Fans (`rtt-analysis`, festgelegt in #28):** Bildfläche ist die Fläche des letzten Events im Pfad; Koordinaten x, y in mm im lokalen KS dieser Fläche (Achsen nach ihrer Pose). Bezug ist der Hauptstrahl (Pupillenmitte) bei der Referenzwellenlänge, für alle Spots und Fans, damit der Farbquerfehler sichtbar bleibt. Polychromatisch zählen die Wellenlängengewichte des Modells, normiert auf Summe 1, mal dem Leistungsgewicht des Strahls; Feldgewichte spielen erst bei feldübergreifenden Größen (Merit-Funktion) eine Rolle. Kennzahlen: Schwerpunkt, RMS und GEO jeweils um Schwerpunkt und Hauptstrahl, GEO nur über Punkte mit Gewicht > 0; Anteil vignettierter Strahlen ungewichtet (alles, was nicht mit Status `Alive` auf der Bildfläche ankommt). Ray Fans: tangential ε(p_y) bei p_x = 0, sagittal ε(p_x) bei p_y = 0, je Punkt ε_x, ε_y und Status; nicht angekommene Punkte haben ε = 0 und ihren Status. Kein angekommener Hauptstrahl → `AnalysisError`.

**OPD (`rtt-analysis`, festgelegt in #29):** Referenzsphäre um den realen Durchstoßpunkt des Hauptstrahls der Referenzwellenlänge auf der Bildfläche, Radius bis zur Mitte der paraxialen Austrittspupille; die OPL eines Strahls bis zur Sphäre ist seine OPL auf der Bildfläche minus |n′|·s (s vorzeichenbehafteter Weg von der Sphäre zur Bildfläche entlang des Strahls). Bezugsstrahl ist der Hauptstrahl derselben Wellenlänge (W(0, 0) = 0 je Wellenlänge); Piston und Tilt werden nicht entfernt, der Farbquerfehler erscheint als Tilt. Einheit: Wellen bei der Referenzwellenlänge. OPD-Karte auf einem gleichmäßigen Gitter im Einheitskreis und OPD-Fans; RMS = Standardabweichung über die angekommenen Punkte (Piston entfernt, **Tilt nicht**), PV = max − min; vignettierte Punkte haben W = 0, ihren Status und zählen nicht mit. Bildseitig telezentrische Systeme (Austrittspupille im Unendlichen) → `AnalysisError`.

**Feldabhängige Analysen (`rtt-analysis`, festgelegt in #31):** Feldverläufe laufen über das relative Feld 0 … 1 entlang +y bis zum größten Feldpunkt (linear im Feldwert, Grad bzw. mm). Höhen sind vorzeichenbehaftet: Abstand vom Bildflächenscheitel, projiziert auf die Feldrichtung in der globalen x-y-Ebene. *Verzeichnung* D = (h_real − h_par)/h_par · 100 % mit dem realen Hauptstrahl auf der Bildfläche und dem paraxialen Hauptstrahl **in der Ebene des Bildflächenscheitels** (nicht in der paraxialen Bildebene); D = 0 auf der Achse. Realer und paraxialer Hauptstrahl starten am selben Feldpunkt: Der Feldwert wird wie in `rtt-trace` bei der Referenzwellenlänge umgerechnet (#50), beide Strahlen laufen bei der ausgewerteten Wellenlänge. *Bildfeldwölbung* numerisch: je Feld zwei Nachbarstrahlen bei ±δ (δ = 1e−3 normierte Pupille) tangential (zum Feldpunkt hin) und sagittal; Fokus = Mitte der kürzesten Verbindung der beiden Geraden im Bildraum, angegeben als z-Komponente des Abstands vom Bildflächenscheitel (parallel zur Achse) mit dem Vorzeichen der Ausbreitung im Bildraum; Fehler O(δ²), ungerade Terme (Koma) heben sich auf. *Farblängsfehler* für ein Wellenlängenpaar (Standard erste minus letzte Systemwellenlänge): paraxial aus den Brennpunkten von `rtt-paraxial`, real über den Achsenschnitt des Zonenstrahls (0, zone) des Achsfeldes, beide längs der Ausbreitung. *Farbquerfehler*: realer Hauptstrahl jeder Wellenlänge auf der Bildfläche und sein Versatz gegen den der Referenzwellenlänge.

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

**Stand M2 (#32):** Paket `raytatouille` unter `libs/rtt-py` (C++ in `src/`, Paket in `python/raytatouille/`, Tests in `tests/`), Erweiterung `raytatouille._core` mit nanobind, Build mit scikit-build-core (`pyproject.toml` im Wurzelverzeichnis, `pip install .`) oder in der CMake-Baumstruktur mit `RTT_BUILD_PYTHON=ON` (Presets `ci-linux-python`, `ci-windows-python`; Werkzeuge vorher mit `pip install -r libs/rtt-py/tests/requirements.txt`, ADR 0018). `pip install .` braucht die C++-Abhängigkeiten aus vcpkg über `CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows"` (Linux: Triplet `x64-linux` oder Systempakete). Unter Windows nicht aus einer Visual-Studio-Entwicklerumgebung für x86 bauen: scikit-build-core übernimmt `Platform=x86` bzw. `VSCMD_ARG_TGT_ARCH=x86` und konfiguriert dann Win32, vcpkg baut alle Abhängigkeiten für x86 neu, und das Linken gegen die 64-Bit-Python scheitert; vorher `Platform` und `VSCMD_ARG_TGT_ARCH` auf `x64` setzen oder eine normale Shell nehmen. Vorhanden: `rt.load`, `rt.save`, `rt.validate`, `rt.System` (Name und Umgebung änderbar, volles Editieren folgt), `rt.MaterialLibrary` (`add_catalog`, `index`), `rt.compile`, `rt.paraxial.first_order`, `rt.trace.make_rays` und `rt.trace.trace` (Pfad als Name oder Index, Wellenlänge `None` = Referenz, `threads` ohne Einfluss auf das Ergebnis). `RayBatch`-Spalten sind beschreibbare NumPy-Views ohne Kopie; jede View hält den Batch am Leben, die Größe ist aus Python fest (kein `resize`), `prt_matrices()` ist eine dokumentierte Kopie; während `trace` ist das GIL freigegeben. C++-Fehler kommen als Klassen aus `raytatouille.errors` (`RaytatouilleError` als Basis; `ParseError` mit `pointer`, `CompileError` mit `diagnostics`, `ParaxialError`, `UnknownMaterial`, `AgfError` mit `file` und `line`) mit der C++-Meldung; Dateifehler als `OSError`, sonst ValueError/IndexError wie in nanobind. Die Typ-Stubs `_core.pyi` erzeugt `nanobind.stubgen` beim Build. Tests: Roundtrip Datei → Python → Datei bitgleich für alle Referenzdateien; Ergebnisse aus Python bitgleich zu C++ über das Testprogramm `rtt_py_reference` aus demselben CMake-Baum (1 und 4 Threads); mypy --strict auf die Tests.

**Stand M2 (#33):** `rt.analysis` mit `spot`, `ray_fan`, `opd_map`, `opd_fan`, `longitudinal_colour`, `lateral_colour`, `distortion`, `distortion_at`, `field_curvature`, `field_curvature_at` und `seidel` (dasselbe Objekt wie `rt.paraxial.seidel`). Erstes Argument ist ein `System` (bei jedem Aufruf mit `materials=` kompiliert) oder ein `CompiledSystem`; für mehrere Analysen einmal `rt.compile()` aufrufen. Wellenlänge `None` heißt Referenzwellenlänge, nur bei `spot` polychromatisch mit den Wellenlängengewichten des Modells (wie in C++). Pupillenverteilungen als Objekt oder Kurzform `"hexapolar:N"`, `"grid:N"`, `"fan_x:N"`, `"fan_y:N"`, `"random:N[:SEED]"`, `"single:PX,PY"` (ValueError mit der Kurzform bei Fehlern). Ergebnisse sind die gebundenen C++-Structs; Listen darin (Spotpunkte, Fans, OPD-Punkte, Fokuslagen, Feldverläufe, Seidel-Beiträge je Fläche) kommen als read-only NumPy-Kopien, je Größe ein Array. `AnalysisError` gehört zu `raytatouille.errors`. `raytatouille.plot` zeichnet die Ergebnisse mit matplotlib ohne eigene Auswertung; matplotlib ist optional (`raytatouille[plot]`, ADR 0018 Nachtrag). Jede Analyse ist gegen C++ bitgleich getestet (`rtt_py_reference`, 1 und 4 Threads, auch mit Vignettierung).

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

**CI (`.github/workflows/ci.yml`):** clang-format, JSON-Schema, Linux GCC (Release), Linux Clang mit clang-tidy (Release), Linux Clang mit ASan/UBSan (Debug), Windows MSVC (Release); Warnungen sind Fehler. Job Python (ADR 0018): Paket `raytatouille` in der CMake-Baumstruktur (Presets `ci-linux-python` mit Python 3.10 und 3.13, `ci-windows-python` mit 3.13), pytest einschließlich Bitvergleich mit C++ und mypy --strict über ctest; je OS ein Smoke-Test mit `pip install .`.

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
