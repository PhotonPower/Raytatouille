# User-Feedback: Anforderungen an die Bibliothek für eine GUI

Stand der Beobachtungen: 2026-10-06, Commit `4221f2d` (Version 0.3.0 mit M3-Vorarbeiten im
Unreleased-Block). Vorgeschlagener Ablageort: `docs/feedback/gui-anforderungen.md`.

## Kontext und Methode

Für das Ausprobieren der Engine wurde eine Streamlit-Oberfläche (`rtt_explorer.py`) auf der
öffentlichen Python-API gebaut. Sie lädt und baut Systeme, zeigt paraxiale Daten und alle
Analysen aus `rtt.analysis` und `rtt.paraxial` an und startet Strahlbündel mit `rtt.trace`.
Dabei traten Lücken auf, die eine echte GUI (M10, PySide6, pyvista) ebenfalls treffen würden.
Dieses Dokument sammelt sie als Aufträge, die zusätzlich zur Roadmap M4 bis M10 sinnvoll sind.

Getestet wurde unter Linux (Ubuntu 24.04, Systempakete statt vcpkg) mit dem aus dem Repo
gebauten Python-Paket. Alle Beispielsysteme unter `tests/reference` wurden geladen und
durch alle Analysen geschickt.

Alle Punkte sind so geschnitten, dass sie nach `AGENTS.md` ein eigener Auftrag sein können
(eine Bibliothek, ein Ziel, ein Referenztest mit analytischem Sollwert).

## Beobachteter Ist-Stand (Python-API 0.3.0)

| Beobachtung | Folge für die GUI |
| --- | --- |
| `trace()` liefert nur den Endzustand und `last_surface` | Strahlen im Layout lassen sich nicht zeichnen |
| `CompiledSystem` kennt `surface_ids`, `path_names`, `media`, `field_count` usw., aber keine Posen, Formen oder Aperturen | Kein Linsenschnitt, keine 3D-Ansicht |
| `System` ist bis auf `name` und `environment` schreibgeschützt; Apertur, Feldtyp und Feldwerte sind nicht abrufbar | Die App musste die JSON-Datei selbst parsen; kein Editor möglich |
| `MaterialLibrary` hat nur `add_catalog` und `index` (skalar) | Kein Glaskatalog-Browser, langsame Dispersionskurven |
| Zwei Kataloge mit gleichem Namen (`tests/catalogs/schott.agf` und `tests/catalogs/m2/schott.agf`) lassen sich nicht gemeinsam laden | Das Cooke-Triplet braucht den zweiten, der Wurzelkatalog den ersten |
| Analysen können nicht abgebrochen werden und melden keinen Fortschritt | Bedienung friert bei großen Rechnungen ein |
| `AnalysisError` und `ParaxialError` tragen weder JSON-Pointer noch Flächen-ID (`CompileError.diagnostics` ist strukturiert, aber ohne Fehlercode) | Fehler lassen sich im GUI nicht an der betroffenen Zeile markieren |

## Priorität 1: ohne diese Punkte wird das GUI nur ein Formular

### G1 Strahlpfad aufzeichnen (`rtt-trace`, `rtt-py`)

- **Problem:** Siehe oben, Strahlen anzeigen ist nicht möglich.
- **Vorschlag:** Option `record_path=True` in `trace()`. Sie liefert pro Strahl und Ereignis
  Position, Richtung, OPL und Status als NumPy-Arrays der Form `(N, n_events, …)`. Bei
  verlorenen Strahlen wird die ID der Fläche mitgeliefert, an der der Strahl endete. Der
  Speicherbedarf muss begrenzbar sein (Chunking oder Maximum).
- **Abnahme:**
  - Der letzte aufgezeichnete Punkt ist bitgleich mit dem normalen Trace-Ergebnis.
  - Das Ergebnis ist unabhängig von der Thread-Anzahl.
  - Ein Strahl auf der Achse durch ein Singlet trifft alle Flächenscheitel.

### G2 Geometrie-Export für Layout und 3D (`rtt-compile`, `rtt-geom`, `rtt-py`)

- **Vorschlag:**
  - Flächenliste mit ID, Element, Typ (Linse, Blende, Detektor), globaler Pose und
    Transformation, Formparametern, Apertur sowie Material davor und danach.
  - `sag(x, y)` und Normale je Fläche.
  - Profil-Polylinien in beliebiger Schnittebene.
  - Element-Umrisse (Rand und Kitt-Flächen).
  - Tessellierbares Mesh für 3D.
- **Abnahme:**
  - Sag-Werte stimmen mit der analytischen Formel überein (Sphäre, Konik, gerade Asphäre).
  - Der Linsenumriss schließt für Singlet, Dublett und Meniskus.
  - Bei Spiegel und gekipptem Element stimmen die globalen Koordinaten.

### G3 Modell lesen und ändern (`rtt-model`, `rtt-py`)

- **Vorschlag:**
  - Lesezugriff auf den ganzen Baum (Wellenlängen, Apertur, Felder, Pfade, Events, Elemente,
    Flächen).
  - Änderungen über JSON-Pointer oder Objektreferenz; Knoten einfügen, löschen, verschieben.
  - Undo und Redo.
  - Eine Flächentabelle (wie im Lens Data Editor) als Sicht auf das Modell, in beide
    Richtungen.
  - Jede Änderung ist als Skriptbefehl wiedergebbar (Vorgabe in `docs/architecture.md`).
- **Abnahme:**
  - Roundtrip Modell → Tabelle → Modell ist für sequenzielle Systeme verlustfrei.
  - Undo stellt den bitgleichen Zustand her.
  - Ungültige Änderungen werfen eine Exception mit Pointer.
- **Hinweis:** Format- und Modelländerung. Schema-Version, JSON-Schema, Referenzdateien,
  Migration und gegebenenfalls ADR gemäß Regel 8 und 9 in `AGENTS.md`.

### G4 Interaktive Performance und Abbruch (`rtt-trace`, `rtt-analysis`, `rtt-py`)

- **Vorschlag:**
  - Cancel-Token und Fortschritts-Callback für alle Analysen.
  - Günstiges Neukompilieren nach kleinen Änderungen (nur betroffene Teile).
  - Benchmarks in der CI.
- **Abnahme:**
  - Ein abgebrochener Aufruf kehrt innerhalb einer definierten Zeit mit einer eigenen
    Exception zurück.
  - Zielwert zum Festlegen: Spot mit etwa 1000 Strahlen am Cooke-Triplet unter 10 ms.

### G5 Systemdaten wie im Prescription Report (`rtt-paraxial`)

- **Vorschlag:**
  - Paraxialer Marginal- und Hauptstrahl pro Fläche (Höhe, Winkel, Einfallswinkel).
  - Gesamtlänge, Arbeits-f/#, bildseitige NA, Vergrößerungen, Lagrange-Invariante.
- **Abnahme:** Analytische Werte für dünne Linse, Singlet und Teleskop; Quellen nach
  Regel 6 in `docs/quellen.md`.

## Priorität 2: macht das GUI komfortabel

### G6 Footprint und freie Öffnung (`rtt-analysis`, braucht G1)

- Strahlfußabdruck pro Fläche und Feld, automatische Berechnung des nötigen
  Halbdurchmessers, Vignettierungsfaktoren pro Feld.
- **Abnahme:** Ein kollimiertes Bündel trifft die erste Fläche als Kreis mit Radius EPD/2.
  Die Halbdurchmesser des Cooke-Triplets stimmen mit einem Referenzwert überein.

### G7 Weitere Standardanalysen (`rtt-analysis`)

- Through-Focus-Spot, longitudinale sphärische Aberration über die Pupillenzonen,
  Spot-Feldraster, Verzeichnungsgitter, relative Beleuchtung, Hauptstrahlwinkel am Bild (CRA).
- PSF, MTF und Encircled Energy bleiben bei M6.
- **Abnahme je Analyse:** Analytischer Fall oder bekannte Referenz, z. B. Parabolspiegel
  ohne sphärische Aberration.

### G8 Fokus- und Skalierungswerkzeuge ohne Optimierer (`rtt-model`, `rtt-analysis`, braucht G3)

- Bildebene auf paraxialen Fokus setzen, Best Focus (minimaler RMS-Spot), System auf
  Ziel-EFL oder f/# skalieren, eine Dicke so lösen, dass die Randstrahlhöhe stimmt.
- M5 kann diese Funktionen später als Bausteine nutzen.
- **Abnahme:** Nach dem Setzen auf den paraxialen Fokus beträgt die Abweichung weniger als
  1e-9 mm. Nach dem Skalieren stimmt die EFL exakt.

### G9 Materialbibliothek durchsuchbar machen (`rtt-material`, `rtt-py`)

- **Vorschlag:**
  - Kataloge und Gläser auflisten, mit n_d, ν_d, Formel, Wellenlängenbereich, dn/dT, Dichte,
    Preis und Transmission.
  - Vektorisiertes `index()` für NumPy-Arrays.
  - Laden eines Katalogs aus dem Speicher statt aus einer Datei.
  - Alias oder Namensraum für Katalognamen, damit gleichnamige Dateien nicht kollidieren.
  - Suche und Daten für die Glaskarte (n_d gegen ν_d).
- **Abnahme:** Das Listing stimmt mit den AGF-Zeilen überein. Vektorisierte und skalare
  Werte sind bitgleich.

### G10 Strukturierte Ergebnisse und Diagnosen (`rtt-py`, `rtt-model`)

- **Vorschlag:**
  - Stabiler Fehlercode pro Diagnostic.
  - Pointer oder Flächen-ID in allen Fehlerklassen, damit das GUI die betroffene Zeile
    markieren kann.
  - `to_dict()` und `to_json()` für alle Ergebnisobjekte mit versioniertem Format (auch
    nützlich, falls die GUI in einem getrennten Prozess läuft).
  - Zusätzliche Warnungen, z. B. Blende außerhalb des Strahlengangs oder viele verlorene
    Strahlen.

## Priorität 3: Roadmap anpassen

### G11 Abhängigkeiten zur bestehenden Roadmap

- **ZMX-Import (M8) vorziehen:** Ohne Import bestehender Designs ist die GUI für echte Arbeit
  kaum nutzbar. Dazu ein Export als Flächentabelle (CSV).
- **Python-Zugriff auf Polarisation und Coatings (#61):** Die GUI braucht ihn für
  Retardance- und Diattenuation-Anzeigen. Der Punkt ist bereits geplant und steht hier nur
  als Abhängigkeit.
- Python-API-Dokumentation und Typ-Stubs für alle neuen Funktionen, wie in
  `docs/architecture.md` gefordert.

## Reihenfolge und Hinweise für die Agenten

- **Parallel startbar:** G1, G2, G5 und G9 sind voneinander unabhängig.
- **Danach:** G6 braucht G1 (und sinnvollerweise G2). G8 braucht G3. G3 und G4 sollten vor
  dem Start von M10 abgeschlossen sein.
- **Pro Auftrag** gilt `AGENTS.md`: Referenztest mit analytischem Sollwert zuerst, Physik mit
  Quelle in `docs/quellen.md`, Einheiten im Doxygen-Kommentar, Determinismus unabhängig von
  der Thread-Anzahl.
- **Bei Formatänderungen** (G3, eventuell G9): Schema-Version erhöhen, JSON-Schema,
  Referenzdateien und Migration anpassen.

## Weitere Beobachtungen (kein Auftrag zur GUI, aber dokumentationsrelevant)

- Die `README.md` nennt als Stand "Meilenstein M0 … noch kein Raytracing". Laut
  `CHANGELOG.md` sind M1 und M2 abgeschlossen und M3 ist weit fortgeschritten.
- In der Roadmap-Tabelle in `docs/architecture.md` ist nur M0 mit ✅ markiert.
- `tests/reference/m0/feature_tour.rtt.json` und `tests/reference/m3/ar_singlet.rtt.json`
  lassen sich in Python 0.3.0 nicht kompilieren (Zernike-Terme erst ab M8, Coatings nicht
  übergebbar). Ein Hinweis dazu in der Dokumentation oder in den Dateien würde Nutzern die
  Fehlersuche ersparen.
- `tests/reference/m1/paraboloid_mirror.rtt.json` und `tests/reference/m0/michelson.rtt.json`
  haben keine Blende bzw. sind gekippt. Spot, OPD und Ray Fans scheitern dort mit
  `ValueError: sources: the path has no stop to aim at`. Eine eigene Fehlerklasse dafür
  wäre für die GUI hilfreich (siehe G10).
