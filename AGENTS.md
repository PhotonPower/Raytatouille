# Regeln für Programmieragenten

Diese Datei liest jeder Agent zu Beginn jeder Sitzung. Sie gilt zusammen mit
`docs/architecture.md` (Architektur, Konventionen, Roadmap) und `docs/adr/`.

## Projekt in einem Satz

Raytatouille ist eine headless Raytracing-Engine für Optikdesign in C++20 (MIT-Lizenz) mit
Python-API; Linux und Windows.

## Bauen und testen

```bash
# Linux, Systempakete (Ubuntu 24.04):
sudo apt-get install cmake ninja-build g++ clang clang-format clang-tidy \
     libeigen3-dev nlohmann-json3-dev catch2 libtbb-dev
cmake --preset dev && cmake --build --preset dev && ctest --preset dev

# Mit vcpkg (wie die CI; Windows: Preset ci-windows-msvc):
export VCPKG_ROOT=/pfad/zu/vcpkg
cmake --preset ci-linux-gcc && cmake --build --preset ci-linux-gcc && ctest --preset ci-linux-gcc

# Formatierung, Schema, Systemdateien
clang-format -i $(git ls-files '*.cpp' '*.hpp')
pip install -r tests/schema/requirements.txt && pytest -q tests/schema
./build/dev/apps/rtt-cli/rtt validate tests/reference/m0/*.rtt.json
./build/dev/apps/rtt-cli/rtt format tests/reference/m0/*.rtt.json
```

Windows/MSYS2 (Preset `dev` mit GCC aus `C:\msys64\ucrt64`): Tests über `ctest` starten; CMake setzt
dort den `PATH` der Tests selbst (Compiler-Verzeichnis vorn). Direkt gestartete `.exe` aus Git Bash
brauchen `export PATH=/c/msys64/ucrt64/bin:$PATH`, weil Git for Windows eine inkompatible
msvcrt-`libstdc++-6.dll` in `/mingw64/bin` mitbringt (Einsprungpunkt-Fehler beim Start).

Vor jedem PR müssen lokal grün sein: Build mit `-DRTT_WARNINGS_AS_ERRORS=ON`, alle Tests,
`clang-format --dry-run -Werror`, die Schema-Tests. Clang-tidy läuft mit
`-DRTT_ENABLE_CLANG_TIDY=ON`, Sanitizer mit `-DRTT_ENABLE_SANITIZERS=ON`.

## Aufbau

```text
libs/rtt-<name>/include/rtt/<name>/*.hpp   öffentliche Header (Namespace rtt::<name>)
libs/rtt-<name>/src/                         Implementierung
libs/rtt-<name>/tests/                       Catch2-Tests
apps/rtt-cli/                                Programm `rtt`
schema/raytatouille.schema.json              JSON-Schema des Dateiformats
tests/reference/                             Referenzsysteme (*.rtt.json) und Sollwerte
tests/schema/                                pytest für das Schema
docs/                                        Architektur und ADRs
```

Neue Bibliothek: Ordner wie oben, `CMakeLists.txt` nach dem Muster von `libs/rtt-model`
(`rtt_target_defaults(<target> LIBRARY)`, Tests mit `rtt_add_test`), in der Wurzel-
`CMakeLists.txt` eintragen. Abhängigkeiten nur zu tieferen Schichten (siehe Architektur).

## Arbeitsregeln

1. **Ein Auftrag = eine Bibliothek + ein Meilenstein-Ziel.** Abnahmekriterien stehen im Issue
   und in der Roadmap. Nichts darüber hinaus bauen.
2. **Erst Test, dann Code.** Referenztest mit analytischem Sollwert zuerst; er muss vorher
   fehlschlagen. Referenzwerte und bestehende Tests nur mit Freigabe des Maintainers ändern.
   Einen Test anzupassen, damit er grün wird, ist verboten.
3. **Keine Exceptions in Tracing-Schleifen.** Strahlprobleme sind Status-Flags, Eingabefehler
   Exceptions an der API-Grenze. Kein undefiniertes Verhalten; Sanitizer-Läufe sauber.
4. **Modernes C++20:** RAII, kein besitzendes `new`/`delete`, `const` wo möglich, kein
   veränderlicher globaler Zustand, keine C-Casts. Neue Abhängigkeiten nur über `vcpkg.json`
   mit Lizenzprüfung (ADR 0010) und Begründung im PR; Ausnahme für reine Python-Build- und
   Testwerkzeuge siehe ADR 0018.
5. **Konventionen sind Gesetz** (docs/architecture.md, Abschnitt Konventionen): mm, µm, rad
   intern, Suffix `_deg` für Grad, +z optische Achse, Radiusvorzeichen, Rotationsreihenfolge.
   Jede öffentliche Funktion dokumentiert Einheiten, Koordinatensystem und Vorzeichen im
   Doxygen-Kommentar.
6. **Physik mit Quelle:** jede Formel mit Literaturstelle (Quelle, Abschnitt, Gleichungsnummer)
   im Kommentar. Zuerst [docs/quellen.md](docs/quellen.md) verwenden. Fehlt die Formel dort:
   frei zugängliche Quelle mit Gleichungsnummer suchen, im Volltext lesen, Symbol für Symbol mit
   dem Code vergleichen (Vorzeichen, Normalenrichtung, Einheiten) und im selben PR in
   docs/quellen.md eintragen. Gleichungsnummern nie aus dem Gedächtnis oder aus Zitaten
   übernehmen. Ein Frage-Issue erst, wenn keine prüfbare Quelle auffindbar ist.
7. **Determinismus:** Zufall nur mit explizitem Seed; Parallelisierung darf Ergebnisse nicht
   ändern.
8. **Öffentliche Header** anderer Bibliotheken nicht ändern. Brauchst du eine Änderung,
   beschreibe sie im Issue; sie kommt als eigener PR mit ADR.
9. **Dateiformat:** Jede Änderung an Modell oder Format erhöht `kSchemaVersion`, passt
   `rtt-io`, das JSON-Schema und die Referenzdateien an (`rtt format`) und bringt eine Migration.
10. **Unklarheit = Frage, nicht Annahme.** Bei offenen physikalischen Konventionen oder
    widersprüchlichen Vorgaben: Frage im Issue stellen und auf Antwort warten.

## Sprache und Stil

- Code, Kommentare, Commit-Messages und Bezeichner auf Englisch; Projektdokumentation auf Deutsch.
- Stil: `.clang-format` (Google-basiert, 100 Spalten). Dateinamen `snake_case`, Typen
  `PascalCase`, Funktionen und Variablen `snake_case`, Konstanten `kName`, Member mit `_`-Suffix.
- Commits im Imperativ, ein Thema pro Commit, z. B. `geom: add even asphere sag and gradient`.

## Definition of Done (PR-Checkliste)

- [ ] CI grün (Linux GCC/Clang/Sanitizer, Windows MSVC, Format, Schema)
- [ ] Referenztests für jedes neue physikalische Verhalten, mit Quelle
- [ ] Neue Quellen in `docs/quellen.md` eingetragen
- [ ] Öffentliche API mit Doxygen-Kommentaren inkl. Einheiten
- [ ] `CHANGELOG.md` ergänzt
- [ ] Bei Formatänderung: Schema-Version, JSON-Schema, Referenzdateien, Migration
- [ ] Beispiel unter `tests/reference/` oder `examples/` für neue Features
