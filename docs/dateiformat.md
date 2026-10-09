# Systemdateien `.rtt.json` erzeugen

Kurzanleitung für das Dateiformat, das die Engine, die CLI (`rtt`) und der Raytatouille Explorer
laden. Alle Beispiele hier wurden gegen Schema `0.4.x` geprüft. Dateien mit Schema 0.1, 0.2 und
0.3 liest `rtt` weiter und schreibt sie mit `rtt format` als 0.4; ein `pickup` aus einer Datei vor
0.4 fällt dabei weg, mit der Warnung `io.pickup_dropped` (siehe „Parametertabelle“).
Die maßgebliche Beschreibung sind `schema/raytatouille.schema.json` und `docs/architecture.md`.

> Hinweis: Der Raytatouille Explorer (Streamlit-Oberfläche `rtt_explorer.py`) ist nicht Teil dieses
> Repos. Verweise darauf gelten für Nutzer dieses Werkzeugs; alles andere funktioniert mit der
> Bibliothek und der CLI allein.

## Grundregeln

- Einheiten sind fest: Längen in **mm**, Wellenlängen in **µm**, Winkel in **Grad**.
- Koordinaten sind rechtshändig, die optische Achse ist **+z**.
- Das System ist ein Baum: **Baugruppe → Element → Fläche**. Jede Ebene hat eine optionale Pose.
  Ohne weitere Angabe gilt die Position einer Fläche relativ zu ihrem Element, die eines Elements
  relativ zur Baugruppe (siehe „Relative Platzierung“ für die anderen Bezüge).
- Die Pose besteht aus `position` [x, y, z] in mm und `rotation_deg` [rx, ry, rz]. Die Rotationen
  laufen intrinsisch X → Y → Z, mit optionalem `pivot`. Standard ist: erst verschieben, dann um
  den Pivot drehen (`"order": "translate_first"`).
- Radien haben das Vorzeichen der Standardkonvention: positiv, wenn der Krümmungsmittelpunkt
  bei +z liegt. Ein fehlendes `shape` bedeutet eine Planfläche.
- Jede Fläche braucht eine eindeutige `id`. Auf sie beziehen sich Pfade, Ergebnisse und
  Fehlermeldungen.
- Ganzzahlige Felder (z. B. `events[].order`) ohne Dezimalpunkt schreiben; `1.0` ist ein Fehler,
  obwohl das JSON-Schema es zulässt.

## Aufbau auf oberster Ebene

| Feld | Pflicht | Inhalt |
| --- | --- | --- |
| `schema_version` | ja | `"0.4.0"` (Muster `0.4.x`) |
| `name` | nein | Anzeigename |
| `units` | ja | genau `{"length": "mm", "wavelength": "um"}` |
| `environment` | nein | `temperature_c`, `pressure_atm`, `medium` (Standard: Luft nach Ciddor, 20 °C, 1 atm) |
| `object` | nein | `{"at_infinity": true}` ist der Standard; sonst `distance` |
| `wavelengths` | ja | Liste von `{"um": …, "weight": …, "reference": true}`. **Genau eine** Wellenlänge ist `reference`. `weight` ≥ 0. |
| `aperture` | ja | `{"type": "epd" \| "image_fnumber" \| "object_na" \| "stop_size", "value": …}`. `epd`: Durchmesser der Eintrittspupille in mm. `image_fnumber`: F-Zahl bei unendlichen Konjugierten, EPD = |EFL| / F#, auch bei endlichem Objekt. `object_na`: objektseitige numerische Apertur, nur bei endlichem Objekt, **paraxial gelesen**: die paraxiale Randstrahlsteigung vom axialen Objektpunkt ist u = NA / n, n der Brechungsindex des Objektraums (Greivenkamp, OPTI-502, Abschn. 9, S. 9-34: NA = n sin U ≈ n u; ein realer Randstrahl mit n sin U = NA weicht bei großer NA ab). Bei paraxialer Zielung hat der Randstrahl tan U = NA / n; bei realer Zielung trifft er den paraxialen Blendenrand R_s, und sein tan U weicht um die Pupillenaberration davon ab. Bei objektseitiger Telezentrie legen nur `object_na` und `stop_size` das Bündel fest. `stop_size` braucht kein `value` (die Blendengröße kommt aus der Apertur der Blende) |
| `fields` | ja | `{"type": "angle_deg" \| "object_height" \| "paraxial_image_height", "points": [{"x":…, "y":…, "weight":…}]}`. `{}` ist der Achspunkt. `weight` ≥ 0. |
| `configurations` | nein | Konfigurationen (ADR 0029), z. B. `[{"name": "wide"}, {"name": "tele"}]`; fehlt der Abschnitt, gibt es nur die Nominalkonfiguration. Nie leer. |
| `parameters` | nein | Parametertabelle (ADR 0029), siehe unten |
| `root` | ja | die oberste Baugruppe |
| `paths` | ja | mindestens ein Pfad, meist `{"name": "main", "events": "auto"}` |

## Elemente (`type`)

| Typ | Verwendung | Material |
| --- | --- | --- |
| `lens` | Linse; N Flächen ergeben N − 1 Glassegmente (Kittglieder) | `material` |
| `plate` | Platte oder Prisma (Prisma: expliziter Pfad nötig, siehe Pfade) | `material` |
| `mirror` | Spiegel; Standard-Interaktion `fresnel` (siehe unten) | optional `material` (Substrat) |
| `thin_element` | dünnes Element ohne Dicke (Polarisator, Retarder, Strahlteiler) | – |
| `stop` | Aperturblende | – |
| `detector` | Bildebene | – |

**Material** ist ein Text `KATALOG:NAME`. Ein Text gilt für alle Segmente, eine Liste enthält
einen Eintrag pro Segment. Gültig sind `AIR`, `VACUUM`, `CONST:1.5168` und Katalogreferenzen
wie `SCHOTT:N-BK7`. Der Katalogname ist der Dateiname der AGF-Datei in Großbuchstaben (oder
ein Alias), die Datei lädt man beim Kompilieren oder im Explorer.

**Kristalle** (einachsig, ADR 0026, ab Schema 0.3): Das Material ist dann ein Objekt mit zwei
Verweisen für die Hauptbrechzahlen n_O und n_E, dazu die optische Achse in Elementkoordinaten
(nur `lens` und `plate`, ein Kristall für alle Segmente):

```json
"material": {"ordinary": "BIREFRINGENT:CALCITE", "extraordinary": "BIREFRINGENT:CALCITE-E"},
"optic_axis": [0.0, 1.0, 1.0]
```

Der Zemax-Katalog `birefringent.agf` legt jeden Kristall so als zwei Gläser `X` und `X-E` ab.
`rtt validate` prüft Kristall und Achse. compile übersetzt Kristalle seit #131 nach den Pfadregeln
in `docs/architecture.md` (Abschnitt „Kristalle im Pfad“): Eintritt nur mit `ordinary` oder
`extraordinary`, Austritt mit `refract`. Verfolgen lassen sich die Moden erst mit #132, bis dahin
endet der Strahl dort mit `EventImpossible`.

**Spiegel.** Ohne Angabe hat jede Fläche die Interaktion `fresnel`. Ein Spiegel **ohne** `material`
wirkt dann als idealer Leiter (r_s = −1, r_p = +1, verlustfrei), dasselbe wie `ideal_mirror`. Ein
Spiegel **mit** `material` (Substrat, z. B. `CONST:1.2,7.26` für ein Metall mit komplexem Index)
reflektiert nach Fresnel am Substrat. Ein Coating auf einem Spiegel (`{"type": "coating", …}`)
braucht das Substrat, ohne `material` ist es ein Kompilierfehler.

## Flächen

| Feld | Inhalt |
| --- | --- |
| `shape.base` | `{"type": "plane"}`, `{"type": "conic", "radius": R, "conic": k}` oder `{"type": "even_asphere", "radius": R, "conic": k, "coefficients": [A4, A6, …]}` |
| `aperture` | `circular` (`radius`, optional `inner_radius`), `rectangular` (`half_width_x`, `half_width_y`) oder `elliptical` (`semi_axis_x`, `semi_axis_y`) |
| `phases` | Phasenschichten: `{"type": "linear_grating", "lines_per_mm": G, "orientation_deg": ψ}` oder `{"type": "radial_phase", "normalization_radius": R, "coefficients": [c1, c2, …]}` (ADR 0025). Einheiten: G in 1/mm, ψ in Grad von der lokalen x- zur lokalen y-Achse (ψ = 0: Rillen parallel zur lokalen y-Achse), φ = 2πG(x cos ψ + y sin ψ); R in mm, Koeffizienten c_k in **rad** (nicht in Wellen, kein Faktor 2π), c1 gehört zu ρ², φ = Σ c_k ρ^{2k} mit ρ = r/R. Die Werte gelten unverändert für jede Wellenlänge; das Vorzeichen von φ und der Ordnung nach ADR 0025, Punkte 1 bis 3 |
| `diffraction_efficiency` | nur an Flächen mit `phases`: Leistungsanteil je Beugungsordnung, z. B. `[{"order": 1, "efficiency": 0.8}]`; fehlt das Feld, hat jede Ordnung 1, sonst haben nicht aufgeführte Ordnungen 0 |
| `interaction` | `fresnel` (Standard), `ideal_mirror`, `ideal_anti_reflection`, `absorber`, `ideal_beam_splitter` (`reflectance_s`, `reflectance_p`), `ideal_polarizer` (`transmission_axis`, `extinction_ratio`), `ideal_retarder` (`fast_axis`, `retardance_waves`) oder `{"type": "coating", "name": "KATALOG:NAME"}`. Die Achsen von Polarisator und Retarder stehen in Elementkoordinaten. |

Zahlenwerte dürfen auch als `{"value": …, "variable": true, "min": …, "max": …}` stehen: Das
markiert sie als Optimierungsvariable mit Grenzen (beide optional, in der Einheit des Werts). Ein
Zahlenwert kann stattdessen an eine Zeile der Parametertabelle gebunden sein,
`{"param": "NAME"}`, ohne eigenen Wert (siehe unten). Zernike-Terme (`shape.terms`) kennt das
Schema, die Engine kompiliert sie in Version 0.4.0 aber noch nicht (Meldung: "not supported before
M8"). Gitter und Phasenflächen (`phases`) übernimmt die Kompilierung; der Tracer verfolgt seit #127
jede Beugungsordnung `order` nach der lokalen Gittergleichung (ADR 0025) und rechnet
`diffraction_efficiency` in das Gewicht ein; evaneszente Ordnungen enden mit dem Status
`Evanescent`. Kristalle mit `ordinary` und `extraordinary` verfolgt der Tracer seit #132
(ADR 0026).

## Relative Platzierung (ab Schema 0.4, ADR 0028)

Jede Pose kann einen anderen Bezug als den Elternknoten haben und die Reihenfolge von Drehung und
Verschiebung umkehren. Dann ist eine Luftdicke ein einziger Wert, und alles danach wandert mit.

| Feld der Pose | Werte |
| --- | --- |
| `reference` | `absolute` (Standard: der Elternknoten), `relative_to_preceding` (die letzte Fläche davor in Baumreihenfolge, auch über Element- und Baugruppengrenzen), `relative_to_sibling` (das vorangehende Geschwister in derselben Baugruppe bzw. demselben Element) |
| `order` | `translate_first` (Standard: erst verschieben, dann um den Pivot im eigenen KS drehen), `rotate_first` (erst um den Pivot im KS des Bezugs drehen, dann längs der gedrehten Achsen verschieben: ein Koordinatensprung, z. B. hinter einem Faltspiegel) |

```json
{"type": "lens", "name": "L2",
 "pose": {"reference": "relative_to_preceding", "position": [0, 0, 12.5]}, "...": "..."}
```

Die erste Fläche eines Elements bleibt absolut, sie legt das KS des Elements fest. Die Prüfung
meldet `pose.relative_first_surface`, `pose.no_preceding` und `pose.no_sibling`. Die Kompilierung
setzt die globale Lage in Baumreihenfolge zusammen: global(X) = global(Bezug) · Pose(X) (#163).

## Parametertabelle und Konfigurationen (ab Schema 0.4, ADR 0029)

Die Tabelle `parameters` steht vor `root`. Jede Zeile hat einen `name` (`[A-Za-z_][A-Za-z0-9_]*`,
eindeutig) und genau eine Form:

| Form | Bedeutung |
| --- | --- |
| `"value": 40.0` | derselbe Wert in allen Konfigurationen |
| `"values": [20.0, 5.0]` | ein Wert je Konfiguration, in der Reihenfolge von `configurations` |
| `"expression": "TOTAL - G"` | aus früheren Zeilen berechnet (`+ - * /`, Klammern, unäres Minus; nur Namen weiter oben) |

Dazu optional `variable` (Optimierungsvariable, nicht bei `expression`) und `min`/`max` (nicht bei
`expression`). Ein Zahlenwert im Baum bindet sich mit `{"param": "NAME"}` an eine Zeile:

```json
"configurations": [{"name": "wide"}, {"name": "tele"}],
"parameters": [
  {"name": "TOTAL", "value": 40.0},
  {"name": "G", "values": [20.0, 5.0], "variable": true, "min": 2.0, "max": 30.0},
  {"name": "B", "expression": "TOTAL - G"}
],
"root": {"type": "assembly", "name": "zoom", "children": [
  {"type": "lens", "name": "L2",
   "pose": {"reference": "relative_to_preceding", "position": [0.0, 0.0, {"param": "G"}]},
   "...": "..."}
]}
```

`pickup` gibt es nicht mehr: Beim Lesen einer Datei vor 0.4 fällt er weg, der Wert bleibt, und
`rtt validate`, `rtt format` und `rt.load` melden die Warnung `io.pickup_dropped` mit der Stelle.
In einer 0.4-Datei ist `pickup` ein Fehler. Bis #164/#165 wertet die Bibliothek die Ausdrücke
noch nicht aus; die Kompilierung meldet gebundene Werte als `param.unresolved`.

## Pfade

`"events": "auto"` besucht alle Flächen in Baumreihenfolge: Linsen und Platten brechen, Spiegel
reflektieren, dünne Elemente, Blenden und Detektoren lassen durch. Zwei Fälle lehnt "auto" mit einem
Kompilierfehler ab: eine Platte aus einem Material mit mehr als zwei Flächen (Prisma, Würfel) und
ein Spiegel mit Substrat und mehreren Flächen (Mangin-Spiegel). Dort und wo "auto" sonst nicht
reicht (Strahlteiler, Doppeldurchgang), steht eine explizite Liste:

```json
"events": [{"surface": "P.S1"}, {"surface": "P.S2", "kind": "reflect"}]
```

`kind` ist `refract` (Standard), `reflect`, `transmit`, `ordinary` oder `extraordinary`. Eine
Beugungsordnung steht als `"order": m` an einem beliebigen Ereignis einer Fläche mit `phases`
(ADR 0025), z. B. `{"surface": "G.S1", "kind": "reflect", "order": 1}` für ein
Reflexionsgitter; an Flächen ohne Phasenschicht ist `order` ≠ 0 ein Fehler. Die frühere Art
`diffract` (bis Schema 0.2) liest `rtt` als `transmit` mit derselben Ordnung.

## Beispiel 1: Plankonvex-Singlet

Das ist die kleinste sinnvolle Datei (f ≈ 100 mm; geprüft: EFL = 100,0527 mm in Luft, der
Standardumgebung). Die Blende liegt bei z = 0, die Linse bei z = 5 mm, ihre zweite Fläche 4 mm
dahinter, der Detektor im paraxialen Fokus bei z = 106,442 mm.

```json
{
  "schema_version": "0.4.0",
  "name": "Plankonvex-Singlet",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [{"um": 0.4861}, {"um": 0.5876, "reference": true}, {"um": 0.6563}],
  "aperture": {"type": "epd", "value": 20.0},
  "fields": {"type": "angle_deg", "points": [{}, {"y": 3.5}, {"y": 5.0}]},
  "root": {"type": "assembly", "name": "system", "children": [
    {"type": "stop", "name": "stop",
     "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]},
    {"type": "lens", "name": "L1", "pose": {"position": [0, 0, 5.0]}, "material": "CONST:1.5168",
     "surfaces": [
       {"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 51.68}},
        "aperture": {"type": "circular", "radius": 12.7}},
       {"id": "L1.S2", "pose": {"position": [0, 0, 4.0]},
        "aperture": {"type": "circular", "radius": 12.7}}]},
    {"type": "detector", "name": "image", "pose": {"position": [0, 0, 106.442]},
     "surfaces": [{"id": "IMG"}]}
  ]},
  "paths": [{"name": "main", "events": "auto"}]
}
```

## Beispiel 2: verkittetes Dublett mit Katalogglas

Eine Linse mit drei Flächen und zwei Segmenten. `material` ist eine Liste mit einem Eintrag
je Segment. Die Positionen der Flächen sind relativ zum Element (4,0 mm und 6,5 mm hinter der
ersten Fläche). Es braucht einen Katalog mit N-LAK9 und N-SF5 (EFL ≈ 60,59 mm).

```json
{"type": "lens", "name": "L1", "pose": {"position": [0, 0, 5.0]},
 "material": ["SCHOTT:N-LAK9", "SCHOTT:N-SF5"],
 "surfaces": [
   {"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 62.0}}},
   {"id": "L1.S2", "pose": {"position": [0, 0, 4.0]},
    "shape": {"base": {"type": "conic", "radius": -44.0}}},
   {"id": "L1.S3", "pose": {"position": [0, 0, 6.5]},
    "shape": {"base": {"type": "conic", "radius": -130.0}}}]}
```

## Beispiel 3: Hohlspiegel mit Blende

Der Spiegel liegt bei z = 0 und der Detektor im Fokus bei z = −100 mm. Die Blende davor ist nötig,
sonst lassen sich keine Strahlen über die Pupille erzeugen. Die Schnittweite (BFL) ist |R|/2 =
100 mm; die EFL ist 1/Φ und in Luft |R|/(2 n_Luft) ≈ 99,97 mm.

```json
"children": [
  {"type": "stop", "name": "stop", "pose": {"position": [0, 0, -30.0]},
   "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]},
  {"type": "mirror", "name": "M",
   "surfaces": [{"id": "M",
     "shape": {"base": {"type": "conic", "radius": -200.0, "conic": -1.0}},
     "aperture": {"type": "circular", "radius": 20.0},
     "interaction": {"type": "ideal_mirror"}}]},
  {"type": "detector", "name": "focus", "pose": {"position": [0, 0, -100.0]},
   "surfaces": [{"id": "F"}]}
]
```

## Dateien per Python-Skript erzeugen

Am robustesten baut man das Dictionary im Code und schreibt es mit `json.dump`. Dieses Skript
(geprüft, EFL = 97,631 mm für R1 = 100, R2 = −100) erzeugt ein Singlet aus Parametern:

```python
import json

def singlet(r1, r2, dicke, glas="CONST:1.5168", epd=20.0, bild_abstand=95.0, halbdurchmesser=12.7):
    """Eine Linse hinter einer Blende. Radien in mm (0 = plan), Dicke in mm."""
    ap = {"type": "circular", "radius": halbdurchmesser}

    def flaeche(id_, radius, z=None):
        f = {"id": id_, "aperture": ap}
        if z:
            f["pose"] = {"position": [0, 0, z]}
        if radius:
            f["shape"] = {"base": {"type": "conic", "radius": radius}}
        return f

    z_linse = 5.0
    return {
        "schema_version": "0.4.0", "name": f"Singlet R1={r1} R2={r2}",
        "units": {"length": "mm", "wavelength": "um"},
        "wavelengths": [{"um": 0.4861}, {"um": 0.5876, "reference": True}, {"um": 0.6563}],
        "aperture": {"type": "epd", "value": epd},
        "fields": {"type": "angle_deg", "points": [{}, {"y": 5.0}]},
        "root": {"type": "assembly", "name": "system", "children": [
            {"type": "stop", "name": "stop",
             "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": epd / 2}}]},
            {"type": "lens", "name": "L1", "pose": {"position": [0, 0, z_linse]}, "material": glas,
             "surfaces": [flaeche("L1.S1", r1), flaeche("L1.S2", r2, z=dicke)]},
            {"type": "detector", "name": "image",
             "pose": {"position": [0, 0, z_linse + dicke + bild_abstand]}, "surfaces": [{"id": "IMG"}]},
        ]},
        "paths": [{"name": "main", "events": "auto"}],
    }

with open("meine_linse.rtt.json", "w", encoding="utf-8") as f:
    json.dump(singlet(100.0, -100.0, 5.0), f, indent=2)
```

Alternativ: im Explorer unter **System bauen** eine Flächentabelle ausfüllen und über
**System-Datei → Als .rtt.json speichern** die fertige Datei herunterladen.

## Prüfen und laden

```powershell
rtt validate meine_linse.rtt.json                  # Format und Konsistenz
python -c "import raytatouille as rt; s = rt.load('meine_linse.rtt.json'); print(rt.validate(s))"
python -m streamlit run rtt_explorer.py            # dann links "Datei hochladen"
```

Eine leere Liste bei `rt.validate` heißt, es gibt keine Befunde. Fehler nennen den JSON-Pointer der
betroffenen Stelle, z. B. `/root/children/1/surfaces/1/id`. Für Glaskataloge und Coatings
(`KATALOG:NAME`) braucht `rt.compile` die passenden Bibliotheken, der Explorer lädt sie anhand der
Referenzen in der Datei aus `tests/catalogs` oder aus hochgeladenen Dateien.

## Häufige Fehler

| Meldung oder Symptom | Ursache |
| --- | --- |
| `exactly one wavelength must be the reference, found 0` | keine oder mehrere Wellenlängen mit `"reference": true` |
| `duplicate surface id` | zwei Flächen mit gleicher `id` |
| `unknown material reference 'SCHOTT:…'` | Glas fehlt im geladenen Katalog (oder falscher Katalogname) |
| `the path has no stop to aim at` | Es gibt keine Blende. Spot, Fans und OPD brauchen eine, im Explorer hilft ein freies Bündel. |
| `Zernike sag terms are not supported before M8` | Zernike-Terme sind im Schema, werden aber noch nicht kompiliert |
| `tilted against the z axis` | Das System ist gekippt, paraxiale Daten und Seidel gibt es nur für rotationssymmetrische Systeme. Spot und Trace gehen trotzdem. |
| Die Fläche erscheint im Layout riesig | Fläche ohne `aperture`: Die Bibliothek schneidet sie an ihrem Formbereich (halbe Kugel). `aperture` angeben. |
