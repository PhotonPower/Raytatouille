# Roadmap-Erweiterungen

> **Status: entschieden** (2026-10-09, Koordinator im Auftrag des Maintainers; siehe [Entscheidungen](#entscheidungen)). Stand nach 0.6.0; M5 läuft (ADR 0028 bis 0030 und Schema 0.4.0 Teil A sind auf `main`, die Einträge sind dagegen abgeglichen). Verbindlich ist ein Eintrag erst, wenn er unter „Entscheidungen“ als freigegeben steht; Einträge mit „Ideenliste“ werden beim Planen des jeweiligen Meilensteins neu bewertet. Die Roadmap-Tabelle in [architecture.md](architecture.md) bleibt die Quelle für die Meilensteine.

## Herkunft und Umgang mit der Quelle

Die Liste entstand aus einem Vergleich mit QUADOA Optical CAD (Handbuch V 26.6). Wie bei ADR 0029 dient das Handbuch nur als Ideenquelle: Die Kapitelnummern in der Spalte „Vorbild“ sind Verweise zum Nachlesen, aus dem Handbuch steht nichts im Repo. Texte, Tabellen, Operandennamen und Dialoge werden nicht übernommen; Formeln stammen aus Literatur und kommen vor der Umsetzung nach [quellen.md](quellen.md) (AGENTS.md, Regel 6).

## Leseregeln

- Jeder Eintrag ist als ein Auftrag gedacht (AGENTS.md, Regel 1: eine Bibliothek, ein Ziel). Ist ein Eintrag zu groß, teilt ihn der Koordinator beim Anlegen des Issues.
- Die Abnahme ist ein Vorschlag mit analytischem Sollwert. Sie ist vor Beginn des Auftrags gegen die Literatur zu prüfen.
- **F** = ändert Modell oder Dateiformat (AGENTS.md, Regel 9: `kSchemaVersion`, Schema, Referenzdateien, Migration). **ADR** = braucht vorher eine Architekturentscheidung.
- Bereits geplante Punkte (M5 bis M10, G6 bis G8 aus [feedback/gui-anforderungen.md](feedback/gui-anforderungen.md)) stehen hier nicht noch einmal.

## Paket P0: kleine Aufträge parallel zu M5

Keine neuen Konzepte, keine Abhängigkeit von M5.

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| R1 | **Reports als Daten:** Raytrace-Report (Strahl × Ereignis mit Ort, Richtung, OPL, Gewicht, Status), Systemdaten-Report (erste Ordnung, Prescription), Maßreport (Mittendicke, Randdicke, Durchmesser je Element); Ausgabe als CSV und im Ergebnisformat (ADR 0023) | `rtt-analysis`, `rtt-py` | Raytrace-Report bitgleich zu `RayPaths` (#80); ein Achsenstrahl durch ein Singlet trifft alle Flächenscheitel; Mittendicke gleich dem Scheitelabstand im Modell | 17.13 |
| R2 | **Ideale Linse und ideale Zylinderlinse** als `ThinElement`-Interaktion (F) | `rtt-model`, `rtt-compile`, `rtt-trace`, `rtt-paraxial` | Punktobjekt im Abstand s wird bei s′ mit 1/s′ = 1/f + 1/s scharf abgebildet (RMS-Spot < 1e−9 mm; OPD < 1e−6 λ nur, wenn die ideale Linse als perfekte Abbildung zwischen konjugierten Punkten definiert ist, nicht als paraxiale Phasenfläche; die Definition von Richtung und optischem Weg an der idealen Linse legt ein kurzer ADR vor dem Code fest); paraxiale EFL = f; die Zylinderlinse fokussiert nur in einer Meridianrichtung und lässt die andere unverändert; Format und Zeitpunkt: siehe Entscheidungen | 7.6 |
| R3 | **Vignettierungsfaktoren automatisch** je Feld, auf dem Footprint (G6) aufbauend | `rtt-analysis` | Ein kreisrundes Fenster mit Radius a in der Eintrittspupillenebene und Bündelradius R lässt den Leistungsanteil (a/R)² durch; ohne Beschnitt ist der Faktor 1 | 14.1.11.1, 17.1.8 |

## Paket P1: Ergänzungen zu M5 (Optimierung)

Setzt ADR 0029 (Parametertabelle) und ADR 0030 (Merit-Funktion im Abschnitt `optimization`) voraus. Die Einträge ergänzen beide und ersetzen nichts.

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| O1 | **Solves:** Abstand, Krümmung und Kippung so lösen, dass paraxiale Randstrahlhöhe oder -steigung ein Ziel erreicht; „Auf Fokus setzen“ paraxial und als bester Fokus (minimaler RMS-Spot). Zusammenlegen mit G8. (ADR) | `rtt-model`, `rtt-analysis` | Nach dem Lösen weicht die Zielgröße um weniger als 1e−9 mm ab; Krümmungs-Solve auf ein Ziel-F/# ergibt EFL/EPD = F# | 18, 3.6.5.2 |
| O2 | **Benutzerdefinierte Operanden:** zusätzlicher Operandentyp `expression` in der Registry von ADR 0030, mit Ausdrücken nach den Regeln von ADR 0029 (kein Code in der Datei); Python-Callback nur in der Python-API, nicht im Dateiformat (ADR: Namen für Analysegrößen in der Ausdrucksgrammatik) | `rtt-optim`, `rtt-py` | Der Ausdruck „EFL − 100“ liefert dieselbe Merit-Funktion, bitgleich, wie der vorgefertigte EFL-Operand mit Ziel 100 | 19.4.5 |
| O3 | **Randbedingungen und Optionen:** Mindest- und Höchstwerte auf abgeleitete Größen (Mittendicke, Randdicke, Baulänge; Grenzen auf Parameter selbst gibt es schon über `min`/`max`, ADR 0029); zu prüfen: zählen die Feld- und Wellenlängengewichte des Modells in den Generatoren mit (Gewichte je Operand und Generator gibt es schon, ADR 0030); Option, den Farbquerfehler im Spot auszuklammern (Schwerpunkt je Wellenlänge); Operanden auf Vignettierung (aus R3) | `rtt-optim` | Eine Singlet-Optimierung mit Mindest-Randdicke hält die Grenze ein (aktive Randbedingung, Abweichung unter der Optimierertoleranz); mit der Option ist der polychromatische RMS gleich dem Mittel der monochromatischen RMS, wenn die Wellenlängen sich nur durch einen Versatz unterscheiden | 19.4.2.2, 19.5 |
| O4 | **Modellglas** aus n_d und ν_d, optimierbar (F) | `rtt-material`, `rtt-model` | n(λ_d) = n_d exakt; (n_F − n_C) = (n_d − 1)/ν_d; ein dünnes, verkittetes Dublett aus zwei Modellgläsern mit φ₁/ν₁ + φ₂/ν₂ = 0 hat paraxial keinen Farblängsfehler zwischen F und C (für dicke oder getrennte Linsen gilt die Bedingung nur näherungsweise) | 8.2, 19.9.2 |
| O5 | **Thermische Längen:** Abstände und Dicken folgen der linearen Ausdehnung als Konfigurationsparameter (Athermalisierung), über die Ausdrücke aus ADR 0029 | `rtt-compile` | L(T) = L₀·(1 + α·ΔT) exakt; die Luftdicke folgt dem Gehäuse-CTE | 16.1.0.4 |
| O6 | **Setup speichern (Rest):** Die Merit-Funktion liegt mit ADR 0030 schon im Abschnitt `optimization`. Dieser Eintrag ergänzt einen gleichartigen optionalen Abschnitt in `*.rtt.json` für den Toleranzsatz (mit M7). Analyse-Voreinstellungen gehören nicht in die Systemdatei (sie beschreibt das System, nicht den Arbeitszustand der GUI); Lösereinstellungen bleiben Laufparameter (ADR 0030, Verworfen) (F, ADR); ein System ohne den Abschnitt bleibt gültig | `rtt-io`, `rtt-optim` | Roundtrip bitgleich; eine Optimierung aus der gespeicherten Datei reproduziert das Ergebnis (ADR 0004) | 3.4.2 |

## Paket P2: Ergänzungen zu M6 (Beugung, Polarisation)

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| W1 | **Interferogramm als Daten:** Intensität I = 1 + cos(2π·W) und entfaltete Phase aus der OPL-Differenz zweier Pfade (Fortsetzung von #122) | `rtt-analysis` | Michelson mit einem um den kleinen Winkel α gekippten Arm: Streifenperiode λ/(2α) in der Pupille | 17.4 |
| W2 | **Feldraster-Analysen:** RMS-Spot und RMS-Wellenfront über ein 2D-Feldgitter, MTF über das Feld, Verzeichnung 2D, longitudinale Pupillenaberration, Wellenfront-Gradient | `rtt-analysis` | Werte am Gitterpunkt bitgleich zur Einzelanalyse; bei rotationssymmetrischem System ist das Raster symmetrisch (bis 1e−9) | 17.1.5, 17.2, 17.3 |
| W3 | **Polarisationsdaten:** Stokes-Punkte für die Poincaré-Kugel, Polarisationskarte über die Pupille (Azimut, Elliptizität), Transmissions-Fan | `rtt-polar`, `rtt-analysis` | Horizontal, vertikal, ±45°, rechts- und linkszirkular liegen auf den Achsen der Kugel (Konvention aus architecture.md); λ/4-Platte unter 45° bringt linear auf einen Pol | 17.11 |
| W4 | **Apodisation und Gaußstrahl-Quelle:** Pupillengewicht (Gauß u. a.); Gaußstrahl paraxial über ABCD und q-Parameter, im Strahlbündel als Äquivalent aus wenigen Strahlen (F, falls die Quelle ins Format kommt) | `rtt-trace`, `rtt-paraxial` | Leistung innerhalb ρ bei Gaußprofil: 1 − exp(−2ρ²/w²); Waist w₀ in der vorderen Brennebene einer Linse f ergibt hinten w₀′ = λ·f/(π·w₀) | 14.1.8, 27 |
| W5 | **Ghost-Erweiterungen:** Ghosts mit Beugungsordnungen, mehr als zwei Reflexionen mit wählbarer Tiefe, Ghost-Footprint, Ghost-Bestrahlungsstärke-Karte (ADR 0027 ergänzen) | `rtt-compile`, `rtt-analysis` | Planplatte bei senkrechtem Einfall: Ghost zu Nutzbild = R₁·R₂ (Fresnel, unbeschichtet R ≈ 0,04, also 0,0016); Pfadzahl bei Tiefe 4 gleich der Kombinatorik | 3.6.4, 17.1.7, 17.7.3 |

## Paket P3: Ergänzung zu M7 (Toleranzierung)

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| T1 | **Toleranzsätze und weitere Toleranzarten:** Default-Satz nach Fertigungsklasse als Daten (eigene Tabelle, keine fremden Werte), Flächenfehler als additive Zernike-Terme, Toleranz auf Brechzahl und Abbe-Zahl, Toleranz auf Konfigurationsparameter | `rtt-tolerance` | Δn auf eine dünne Linse: ΔEFL/EFL = −Δn/(n − 1); ein Zernike-Defokusterm ändert die OPD wie analytisch erwartet | 20.2.8 bis 20.2.12, 3.6.6.1 |

## Paket P4: Ergänzungen zu M8 (Release v1.0)

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| F1 | **Aperturformen und Arrays:** hexagonal, Polygon, Obskurationen (Kreis, Rechteck, Ellipse), Array- und Hex-Array-Operator für Aperturen und Flächen (Linsenarrays), Transformationsoperator (F) | `rtt-geom`, `rtt-model` | Gleichmäßige Kreisausleuchtung mit dem Umkreisradius R einer hexagonalen Apertur: durchgelassener Anteil 3√3/(2π) ≈ 0,827; Linsenarray mit Raster p: Spot-Schwerpunkte auf den Rasterpunkten | 11, 9.2.4 |
| F2 | **Fresnel-Operator** für Flächen (F) | `rtt-geom` | Paraxiale EFL gleich der Mutterfläche; Strahlablenkung innerhalb einer Zone gleich der Mutterfläche | 9.2.4.4 |
| F3 | **Benutzerdefinierte Form und Phase** per Python-Callback (nur zur Laufzeit, nicht im Dateiformat); weitere Phasen: Spiral, 2D-Gitter, Zernike, XYZ (F) | `rtt-geom`, `rtt-py` | Ein Callback mit der Sphärenformel reproduziert den Spot der `Conic`-Fläche (< 1e−9 mm); Spiralphase l·θ hat den Gradienten l/r in azimutaler Richtung und folgt der lokalen Gittergleichung | 9.1.23, 10 |
| F4 | **Import und Export:** eigene CSV-Flächentabelle in beide Richtungen, ZMX-Export des sequenziellen Teils (Austauschformat, kein Abgleich der Ergebnisse); CodeV-`.seq` nur, wenn eine frei zugängliche Formatbeschreibung vorliegt (AGENTS.md, Regel 6; Lizenz) | `rtt-io` | Cooke-Triplett über Export und Import: EFL, BFL und Spot-RMS gleich (< 1e−9) | 26.3.2 |
| F5 | **Konstruktionshelfer** `singlet`, `doublet`, `triplet` (Python) und ein JSON-Format für Linsenkataloge samt Loader, ohne Herstellerdaten (F) | `rtt-py`, `rtt-io` | `doublet(...)` ist bitgleich zu `tests/reference/m2/achromat.rtt.json` | 3.6.3, 22 |

## Paket P5: v1.x nach M8 (braucht ADR)

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| X1 | **GRIN-Medien** (F, ADR): axial, radial, sphärisch (Luneburg), Gitter. Der ADR legt Strahlintegration, OPL = ∫ n ds, Schnitt mit Flächen, Schrittweitenregelung und Determinismus fest. | `rtt-material`, `rtt-trace`, `rtt-model` | Radial n(r) = n₀·(1 − g²r²/2): paraxial y(z) = y₀·cos(gz) + (u₀/g)·sin(gz); Luneburg n = √(2 − (r/R)²): kollimiertes Bündel trifft den gegenüberliegenden Randpunkt (RMS-Grenze im ADR) | 8.4 |
| X2 | **Dämpfer und pupillenabhängige Elemente** (F): konstant, Gauß, periodisch, weiche Apertur; Coatings und Retarder abhängig von Pupillenposition; tabellierte Jones- und Mueller-Elemente | `rtt-coating`, `rtt-polar` | Gaußdämpfer T(r) = exp(−2r²/w²) bei gleichmäßiger Kreisausleuchtung R: Anteil (w²/(2R²))·(1 − exp(−2R²/w²)) | 12.4, 13.2, 13.3 |
| X3 | **Mechanik als Obstruktion** (F, ADR): Quader, Zylinder und Rohr als absorbierende Körper (Prüfung zwischen den Ereignissen), Export als STL und Punktwolke | `rtt-model`, `rtt-trace`, `rtt-io` | Kollimiertes Bündel mit Radius R durch ein Rohr mit Innenradius a < R: durchgelassener Anteil (a/R)² | 21 |

## Paket P6: mit M9 (nicht-sequenziell)

| ID | Auftrag | Bibliothek | Abnahme (Vorschlag) | Vorbild |
| --- | --- | --- | --- | --- |
| N1 | **Rayfile-Export und -Import** für Quellen, zuerst im eigenen dokumentierten Austauschformat (Spalten wie `RayBatch`, Kopf mit Einheiten, Wellenlänge und Gesamtleistung, blockweises Lesen); Fremdformate nur nach Bedarf und nur mit frei zugänglicher Spezifikation (F, ADR) | `rtt-io`, `rtt-trace` | Roundtrip bitgleich in Ort, Richtung und Fluss; Summe der Flüsse bleibt erhalten | 26.3.1 |
| N2 | **STL- und STEP-Import** als Absorber oder Streukörper (braucht die BVH von M9) | `rtt-io`, `rtt-trace` | Würfel als STL: Strahlen treffen die sechs Seiten an den analytischen Orten | 21.0.2 |

## Paket P7: v2 und später

Huygens-PSF und -MTF (17.5.2, 17.6.2), Beamlet-Propagation (17.9), Faserkopplung single- und multimode (17.10), kohärente Bestrahlungsstärke (17.7.2), Bildsimulation (17.8), Glassubstitution und globale Optimierung (19.7), ISO-10110-Zeichnung (3.6.1.4) und STEP-Export der Flächen. Alle stehen schon in architecture.md unter „Später“ oder v2; hier nur zur Vollständigkeit.

## Bewusst nicht aufgenommen

- **MATLAB-Schnittstelle:** Die Python-API ist die einzige Skriptschnittstelle (ADR 0012).
- **Mitgelieferte Herstellerkataloge:** Lizenzfrage (siehe AGF-Regel in architecture.md).
- **Schieberegler und Lookup-Editor:** Das ist GUI-Funktion (M10) auf der Parametertabelle aus ADR 0029, kein Kernthema.

## Reihenfolge und Abhängigkeiten

1. P0 sofort, parallel zu M5; R3 braucht G6.
2. P1 mit M5; O1 und O3 setzen ADR 0029 und 0030 voraus, O2 erweitert die Operandenregistry von ADR 0030, O6 braucht je einen ADR für den Toleranz- und den Analyseabschnitt (Name, Inhalt, Schema-Version, Migration; der Toleranzabschnitt mit M7).
3. P2 und P3 mit M6 und M7; W1 setzt #122 voraus (erledigt), W3 die Jones-Pupille aus M6.
4. P4 mit M8; R2 läuft als eigenes Issue schon früher, das Format kommt mit dem nächsten Schemasprung (siehe Entscheidungen).
5. P5 erst nach v1.0; X1 ist der größte Eingriff in den Tracer und beginnt mit dem ADR, der schon vor M8 vorliegen soll.

## Entscheidungen

Getroffen am 2026-10-09 vom Koordinator; der Maintainer hat die Entscheidung über diesen Entwurf delegiert. „Issue“ heißt: wird jetzt angelegt. „Mit Mx“ heißt: wird beim Planen dieses Meilensteins als Issue geschnitten und dort gegen dessen Umfang geprüft. „Ideenliste“ heißt: bleibt hier stehen, ohne Zusage.

**Grundsatz:** M5 wird nicht erweitert. 0.7.0 bleibt bei dem in #159 bis #170 festgelegten Umfang; Ergänzungen zur Optimierung folgen danach.

| ID | Entscheidung |
| --- | --- |
| R1 | **Freigegeben, #177.** Parallel zu M5, keine Formatänderung. |
| R2 | **Freigegeben, #178.** Zuerst ein kurzer ADR (Richtung und optischer Weg an der idealen Linse, ideale Zylinderlinse), dann die Umsetzung in der vorgeschlagenen Bauabfolge; das Format kommt mit dem nächsten Schemasprung nach 0.4.0. |
| R3 | **Mit G6** (Footprint), wenn G6 eingeplant wird. |
| O1, O2 | **Nach 0.7.0**, je mit eigenem ADR (Solves hängen von Analysen ab und passen nicht in die reinen Tabellenausdrücke von ADR 0029: Reihenfolge, Zyklen, Ort der Auswertung; O2 braucht Namen für Analysegrößen in der Ausdrucksgrammatik). O1 zusammen mit G8. |
| O3 | **Nach 0.7.0**, als erste Ergänzung der Optimierung (Randbedingungen auf Mittendicke, Randdicke, Baulänge). |
| O4, O5 | **Nach 0.7.0**, Ideenliste bis zur Planung der Optimierungs-Ergänzungen. |
| O6 | **Mit M7**, nur der Abschnitt für den Toleranzsatz; Analyse-Voreinstellungen werden nicht in die Systemdatei aufgenommen. |
| W1 bis W5 | **Mit M6.** |
| T1 | **Mit M7.** |
| F1 bis F5 | **Mit M8.** F4 ohne CodeV-`.seq`, solange keine frei zugängliche Formatbeschreibung vorliegt. |
| X1 | **ADR früh** (nach 0.7.0, vor M8), damit Ray Aiming, Paraxial und Analysen nicht die Annahme „Strahlen sind zwischen Flächen Geraden“ fest einbauen; die Umsetzung kommt in v1.x nach v1.0 und vor M10. |
| X2, X3 | **v1.x**, Ideenliste. |
| N1, N2 | **Mit M9.** N1 zuerst im eigenen Austauschformat, Export von Anfang an, Fremdformate nur als Import nach Bedarf. |
| P7 | Wie in architecture.md (v2 und später). |

**Zu den vier Grundsatzfragen des Entwurfs:**

1. Setup in der Systemdatei: optionale Abschnitte ja, aber nur für Daten, die das System oder seine Prüfung beschreiben (Merit-Funktion, Toleranzsatz); keine Analyse- oder GUI-Voreinstellungen.
2. Rayfile: eigenes, dokumentiertes Austauschformat zuerst (siehe N1).
3. Python-Callbacks (O2, F3): nur über die Python-API, nie im Dateiformat, nur vektorisiert (ein Aufruf je NumPy-Array), in der Doku als Prototyping ausgewiesen. In der Datei stehen nur Ausdrücke (ADR 0029) und Daten; ein fehlender Callback beim Laden ergibt einen Diagnosecode.
4. Ideale Linse: siehe R2.