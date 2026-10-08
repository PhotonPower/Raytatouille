# Architecture Decision Records

Neue Entscheidungen bekommen die nächste Nummer. Bestehende ADRs werden nicht umgeschrieben, sondern durch ein neues ADR mit Status „ersetzt durch“ abgelöst.

| Nr. | Titel | Status |
| --- | --- | --- |
| [0001](0001.md) | Kern in C++20 mit CMake und vcpkg | angenommen |
| [0002](0002.md) | Python-Bindings mit nanobind | angenommen |
| [0003](0003.md) | Eigen und double | angenommen |
| [0004](0004.md) | Parallelität mit oneTBB | angenommen |
| [0005](0005.md) | Strahlen als Structure-of-Arrays | angenommen |
| [0006](0006.md) | Physik-Templates über den Skalartyp | angenommen |
| [0007](0007.md) | Ableitungen per finiter Differenz in v1 | angenommen |
| [0008](0008.md) | JSON-Dateiformat mit strengem Parser | angenommen |
| [0009](0009.md) | Fehlerbehandlung | angenommen |
| [0010](0010.md) | Lizenz MIT | angenommen |
| [0011](0011.md) | Plattformen und Compiler | angenommen |
| [0012](0012.md) | GUI in Python | angenommen |
| [0013](0013.md) | Keine GPU in v1 | angenommen |
| [0014](0014.md) | Laufzeit-Interface für Materialien in double | angenommen |
| [0015](0015.md) | Real umfasst nur double | angenommen |
| [0016](0016.md) | CompiledSystem in eigener Bibliothek rtt-compile | angenommen |
| [0017](0017.md) | Kittglieder: Material je Segment | angenommen |
| [0018](0018.md) | Python-Build-Abhängigkeiten über pyproject.toml statt vcpkg | angenommen |
| [0019](0019.md) | Coating-Kataloge als eigene JSON-Dateien | angenommen |
| [0020](0020.md) | Gemeinsamer strikter JSON-Parser in der Basisschicht | angenommen |
| [0021](0021.md) | Interaktionen im Tracer, Semantik von prt und weight | angenommen |
| [0022](0022.md) | Fehlercodes, Fehlerorte und Diagnosekanal | angenommen |
| [0023](0023.md) | Ergebnisformat raytatouille-result und Strahlverluste als Daten | angenommen |

## Vorlage

```markdown
# ADR NNNN: Titel

**Status:** vorgeschlagen | angenommen | ersetzt durch NNNN

## Kontext

## Entscheidung

## Verworfen

## Folgen
```
