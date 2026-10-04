# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/); Versionen nach SemVer.

## [Unreleased]

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
