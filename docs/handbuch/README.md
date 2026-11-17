# Benutzerhandbuch

`raytatouille-handbuch.tex` ist das Benutzerhandbuch für die Python-API (Deutsch, LaTeX). Es
beschreibt genau eine Version, derzeit **0.10.0**; jedes Beispiel wurde mit dieser Version
ausgeführt, die abgedruckten Ausgaben sind unverändert übernommen.

**PDF zum Herunterladen:** an jedem GitHub-Release, für 0.10.0
[raytatouille-handbuch-0.10.0.pdf](https://github.com/PhotonPower/Raytatouille/releases/download/v0.10.0/raytatouille-handbuch-0.10.0.pdf).
Das PDF liegt nicht im Repository.

## Bauen

LuaLaTeX aus TeX Live (mit babel-german, fontspec, listings, tcolorbox, booktabs, longtable,
siunitx, hyperref; Schriften Latin Modern und DejaVu Sans Mono, beide in TeX Live enthalten):

```bash
cd docs/handbuch
latexmk -lualatex raytatouille-handbuch.tex      # oder zweimal lualatex
```

Die Abbildungen liegen in `abbildungen/`; sie entstehen aus den Beispielen `ex15_layout.py` und
`ex16_berichte.py` (siehe unten, Option `--keep`).

## Beispiele prüfen

Die Beispiele stehen im Handbuch als `pycode`-Blöcke mit Dateinamen, ihre Ausgaben darunter.
`tools/handbuch_beispiele.py` schreibt sie in ein Arbeitsverzeichnis, führt jedes `ex*.py` mit
dem angegebenen Python aus und vergleicht die Ausgaben:

```bash
python tools/handbuch_beispiele.py --python /pfad/zur/venv/bin/python
python tools/handbuch_beispiele.py --python … --update         # neue Ausgaben übernehmen
python tools/handbuch_beispiele.py --python … --keep arbeit/    # Dateien und Abbildungen behalten
```

Die venv braucht die beschriebene Version von `raytatouille` (am besten das Wheel von PyPI) und
matplotlib. Zahlen in den letzten Stellen können zwischen Plattformen abweichen (libm); die
abgedruckten Ausgaben stammen von macOS 26 arm64 mit Python 3.14.

## Bei einem neuen Release

1. Version auf Titelseite, im Kopfkommentar, in „Stand und Gültigkeit“ und in der
   Installation anpassen; Neues beschreiben, mit lauffähigem Beispiel.
2. `tools/handbuch_beispiele.py` mit der neuen Version laufen lassen, Abweichungen prüfen und
   mit `--update` übernehmen; Abbildungen mit `--keep` neu erzeugen.
3. Funktionsreferenz, Klassen, Aufzählungen und Diagnosecodes neu erzeugen (sie sind aus der
   installierten Version abgeleitet, siehe die Einleitungen der Kapitel).
4. PDF bauen und an das GitHub-Release hängen.
