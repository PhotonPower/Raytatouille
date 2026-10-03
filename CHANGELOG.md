# Changelog

Format nach [Keep a Changelog](https://keepachangelog.com/de/1.1.0/); Versionen nach SemVer.

## [Unreleased]

### Hinzugefügt
- Anleitung `docs/agenten.md` und Review-Subagent `.claude/agents/physik-reviewer.md`.
- `rtt-geom`: `Shape<T>`-Schnittstelle, `Plane<T>` und `Conic<T>` (Sag, Gradient, Definitionsbereich)
  sowie analytischer Strahlschnitt mit Status `Missed` statt NaN (#3).

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
