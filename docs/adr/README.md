# Architecture Decision Records

Neue Entscheidungen bekommen die nächste Nummer. Bestehende ADRs werden nicht umgeschrieben, sondern durch ein neues ADR mit Status „ersetzt durch“ abgelöst.

Präzisierungen und Ergänzungen einer angenommenen Entscheidung kommen als datierter Nachtrag ans Ende des ADR (mit Issue- oder PR-Nummer). Für Fehler gilt:

- **Falsche Formel, falsches Vorzeichen oder falscher Sollwert:** Korrektur direkt im Text, mit Verweis auf einen Nachtrag „Erratum“, der die Herleitung und den Beleg (Rechnung, Messung) nennt. So baut niemand nach der falschen Fassung (Beispiel: ADR 0033, #225).
- **Ungenaue oder unvollständige Begründung bei richtiger Entscheidung:** Der Text bleibt stehen und bekommt den Hinweis „(Korrigiert im Nachtrag unten.)“ (Beispiel: ADR 0031, #220).

| Nr. | Titel | Status |
| --- | --- | --- |
| [0001](0001.md) | Kern in C++20 mit CMake und vcpkg | angenommen |
| [0002](0002.md) | Python-Bindings mit nanobind | angenommen |
| [0003](0003.md) | Eigen und double | angenommen |
| [0004](0004.md) | Parallelität mit oneTBB | angenommen |
| [0005](0005.md) | Strahlen als Structure-of-Arrays | angenommen |
| [0006](0006.md) | Physik-Templates über den Skalartyp | angenommen, präzisiert durch [0014](0014.md), [0015](0015.md) |
| [0007](0007.md) | Ableitungen per finiter Differenz in v1 (mit Nachträgen zur Zielung, zur Pupillenkonvention und zum Bezug eines anderen Systems, #96, #35 B, #252, #255) | angenommen |
| [0008](0008.md) | JSON-Dateiformat mit strengem Parser | angenommen |
| [0009](0009.md) | Fehlerbehandlung | angenommen |
| [0010](0010.md) | Lizenz MIT | angenommen |
| [0011](0011.md) | Plattformen und Compiler | angenommen |
| [0012](0012.md) | GUI in Python | angenommen |
| [0013](0013.md) | Keine GPU in v1 | angenommen |
| [0014](0014.md) | Laufzeit-Interface für Materialien in double | angenommen, ergänzt 0006 |
| [0015](0015.md) | Real umfasst nur double | angenommen, präzisiert 0006 |
| [0016](0016.md) | CompiledSystem in eigener Bibliothek rtt-compile | angenommen, präzisiert die Schichten in `docs/architecture.md` (kein früheres ADR) |
| [0017](0017.md) | Kittglieder: Material je Segment | angenommen |
| [0018](0018.md) | Python-Build-Abhängigkeiten über pyproject.toml statt vcpkg | angenommen |
| [0019](0019.md) | Coating-Kataloge als eigene JSON-Dateien | angenommen |
| [0020](0020.md) | Gemeinsamer strikter JSON-Parser in der Basisschicht | angenommen |
| [0021](0021.md) | Interaktionen im Tracer, Semantik von prt und weight | angenommen |
| [0022](0022.md) | Fehlercodes, Fehlerorte und Diagnosekanal | angenommen |
| [0023](0023.md) | Ergebnisformat raytatouille-result und Strahlverluste als Daten | angenommen |
| [0024](0024.md) | Modell lesen und ändern, Undo über JSON Patch | angenommen |
| [0025](0025.md) | Gitter als Phasenschicht (Ordnung, Effizienz, Evaneszenz) | angenommen |
| [0026](0026.md) | Kristallmodell (einachsig, optische Achse, Strahlzustand) | angenommen |
| [0027](0027.md) | Ghost-Generator (Zweifachreflexions-Pfade als erzeugte explizite Pfade) | angenommen |
| [0028](0028.md) | Relative Platzierung (Pose.reference, Pose.order; Nachtrag #252: Starrkörper-Überlagerung für Toleranzen) | angenommen |
| [0029](0029.md) | Parametertabelle und Konfigurationen | angenommen |
| [0030](0030.md) | Optimierung (Merit-Funktion, Operanden, Levenberg-Marquardt; Nachtrag #256: Merit-Funktion eines gestörten Prüflings; Nachtrag #281: Hinweise zur Anwendung) | angenommen |
| [0031](0031.md) | Ideale Linse und ideale Zylinderlinse | angenommen |
| [0032](0032.md) | Pupillenfunktion und Wellenfront-Auswertung (Zernike-Fit, Strehl, Through-Focus; Nachtrag #265: nominaler Bezug für Spot und OPD) | angenommen |
| [0033](0033.md) | FFT-Beugung: PSF, MTF, Encircled Energy | angenommen |
| [0034](0034.md) | Polarisation in der Pupille (P je Ereignis, akkumuliertes Q, Jones-Pupille) | angenommen |
| [0035](0035.md) | Toleranzsatz im Systemfile (Abschnitt `tolerances`, Klassen, Kompensatoren, Kriterien; Nachtrag #254: einheitliche Segmentliste ist ein Material) | angenommen |
| [0036](0036.md) | Toleranzanalyse (Störungsmodell, nominaler Bezug, Kompensatoren, Zufall, Statistik; Nachtrag #254: Linien in Normalluft, Prüfungen in compile; Nachtrag #257: ein Prüfling; Nachtrag #258: Sensitivität bei ungültigen Seiten, Abbruch, mehreren Kriterien; Nachtrag #259: Monte Carlo, ungültige Prüflinge, Seed und Warnungen; Nachtrag #293: Radien über Tabellenzeilen und Kompensatoren) | angenommen |
| [0037](0037.md) | Medium im Objektraum (`object.medium`, Schema 0.7; Startmedium aller Pfade, nur endliches Objekt, NA-Prüfung in validate und compile) | angenommen |
| [0038](0038.md) | Bezug des Objektabstands (`object.reference`, Schema 0.8; Abstand vor dem nominalen Scheitel einer Fläche, `CompiledSystem::object_z()`, nominaler Bezug bei der Toleranzierung) | angenommen |
| [0039](0039.md) | Blende auf einer Fläche (`aperture.stop`, Schema 0.8; Verweis auf eine Fläche einer Linse, eines Spiegels, einer Platte oder eines dünnen Elements, höchstens eine Blende; eine Stelle für die Blende in compile; Austrittspupille durch die Blendenfläche selbst; die Blende folgt bei der Toleranzierung der gestörten Fläche) | angenommen |
| [0040](0040.md) | Objektseitige NA real gelesen (`object_na` mit `reading`, Schema 0.8; realer axialer Randstrahl bestimmt das Ziel der Zielung, paraxiale Größen bleiben bei NA/n) | angenommen |
| [0041](0041.md) | Grundlagen für neue Flächentypen (M8 S0: Scheitelentwicklung statt `base_conic` für die Paraxialrechnung, Schnitt mit statischem Dispatch, keine stillen Rückfälle, Schemaplan Welle 3) | angenommen |

## Vorlage

```markdown
# ADR NNNN: Titel

**Status:** vorgeschlagen | angenommen | ersetzt durch NNNN

## Kontext

## Entscheidung

## Verworfen

## Folgen
```
