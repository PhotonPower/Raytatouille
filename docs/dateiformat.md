# Systemdateien `.rtt.json` erzeugen

Kurzanleitung für das Dateiformat, das die Engine, die CLI (`rtt`) und der Raytatouille Explorer
laden. Alle Beispiele hier wurden gegen Schema `0.2.x` und Bibliotheksversion 0.4.0 geprüft.
Die maßgebliche Beschreibung sind `schema/raytatouille.schema.json` und `docs/architecture.md`.

> Hinweis: Der Raytatouille Explorer (Streamlit-Oberfläche `rtt_explorer.py`) ist nicht Teil dieses
> Repos. Verweise darauf gelten für Nutzer dieses Werkzeugs; alles andere funktioniert mit der
> Bibliothek und der CLI allein.

## Grundregeln

- Einheiten sind fest: Längen in **mm**, Wellenlängen in **µm**, Winkel in **Grad**.
- Koordinaten sind rechtshändig, die optische Achse ist **+z**.
- Das System ist ein Baum: **Baugruppe → Element → Fläche**. Jede Ebene hat eine optionale Pose.
  Die Position einer Fläche gilt relativ zu ihrem Element, die eines Elements relativ zur Baugruppe.
- Die Pose besteht aus `position` [x, y, z] in mm und `rotation_deg` [rx, ry, rz]. Die Rotationen
  laufen intrinsisch X → Y → Z, mit optionalem `pivot`.
- Radien haben das Vorzeichen der Standardkonvention: positiv, wenn der Krümmungsmittelpunkt
  bei +z liegt. Ein fehlendes `shape` bedeutet eine Planfläche.
- Jede Fläche braucht eine eindeutige `id`. Auf sie beziehen sich Pfade, Ergebnisse und
  Fehlermeldungen.

## Aufbau auf oberster Ebene

| Feld | Pflicht | Inhalt |
| --- | --- | --- |
| `schema_version` | ja | `"0.2.0"` (Muster `0.2.x`) |
| `name` | nein | Anzeigename |
| `units` | ja | genau `{"length": "mm", "wavelength": "um"}` |
| `environment` | nein | `temperature_c`, `pressure_atm`, `medium` (Standard: Luft nach Ciddor, 20 °C, 1 atm) |
| `object` | nein | `{"at_infinity": true}` ist der Standard; sonst `distance` |
| `wavelengths` | ja | Liste von `{"um": …, "weight": …, "reference": true}`. **Genau eine** Wellenlänge ist `reference`. |
| `aperture` | ja | `{"type": "epd" \| "image_fnumber" \| "object_na" \| "stop_size", "value": …}` |
| `fields` | ja | `{"type": "angle_deg" \| "object_height" \| "paraxial_image_height", "points": [{"x":…, "y":…, "weight":…}]}`. `{}` ist der Achspunkt. |
| `root` | ja | die oberste Baugruppe |
| `paths` | ja | mindestens ein Pfad, meist `{"name": "main", "events": "auto"}` |

## Elemente (`type`)

| Typ | Verwendung | Material |
| --- | --- | --- |
| `lens` | Linse; N Flächen ergeben N − 1 Glassegmente (Kittglieder) | `material` |
| `plate` | Platte oder Prisma | `material` |
| `mirror` | Spiegel, meist mit `"interaction": {"type": "ideal_mirror"}` | – |
| `thin_element` | dünnes Element ohne Dicke (Polarisator, Retarder, Strahlteiler) | – |
| `stop` | Aperturblende | – |
| `detector` | Bildebene | – |

**Material** ist ein Text `KATALOG:NAME`. Ein Text gilt für alle Segmente, eine Liste enthält
einen Eintrag pro Segment. Gültig sind `AIR`, `VACUUM`, `CONST:1.5168` und Katalogreferenzen
wie `SCHOTT:N-BK7`. Der Katalogname ist der Dateiname der AGF-Datei in Großbuchstaben (oder
ein Alias), die Datei lädt man beim Kompilieren oder im Explorer.

## Flächen

| Feld | Inhalt |
| --- | --- |
| `shape.base` | `{"type": "plane"}`, `{"type": "conic", "radius": R, "conic": k}` oder `{"type": "even_asphere", "radius": R, "conic": k, "coefficients": [A4, A6, …]}` |
| `aperture` | `circular` (`radius`, optional `inner_radius`), `rectangular` (`half_width_x`, `half_width_y`) oder `elliptical` (`semi_axis_x`, `semi_axis_y`) |
| `interaction` | `fresnel`, `ideal_mirror`, `ideal_anti_reflection`, `absorber`, `ideal_beam_splitter`, `ideal_polarizer`, `ideal_retarder` oder `{"type": "coating", "name": "KATALOG:NAME"}` |

Zahlenwerte dürfen auch als `{"value": …, "variable": true}` stehen, das markiert sie später als
Optimierungsvariable. Zernike-Terme (`shape.terms`) kennt das Schema, die Engine kompiliert
sie in Version 0.4.0 aber noch nicht (Meldung: "not supported before M8"). Gitter und Phasenflächen
(`phases`) lassen sich laden und kompilieren; ob ein Pfad sie auswertet, habe ich nicht geprüft.

## Pfade

`"events": "auto"` besucht alle Flächen in Baumreihenfolge: Linsen brechen, Spiegel reflektieren,
Blenden und Detektoren lassen durch. Wo das nicht reicht (Prismen mit Reflexion, Strahlteiler,
Doppeldurchgang), steht eine explizite Liste:

```json
"events": [{"surface": "P.S1"}, {"surface": "P.S2", "kind": "reflect"}]
```

`kind` ist `refract`, `reflect`, `transmit`, `ordinary`, `extraordinary` oder `diffract`.

## Beispiel 1: Plankonvex-Singlet

Das ist die kleinste sinnvolle Datei (f ≈ 100 mm; geprüft: EFL = 100,0527 mm). Die Blende liegt
bei z = 0, die Linse bei z = 5 mm, ihre zweite Fläche 4 mm dahinter.

```json
{
  "schema_version": "0.2.0",
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
    {"type": "detector", "name": "image", "pose": {"position": [0, 0, 106.363]},
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
sonst lassen sich keine Strahlen über die Pupille erzeugen (EFL ≈ 99,97 mm, also |R|/2).

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
        "schema_version": "0.2.0", "name": f"Singlet R1={r1} R2={r2}",
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
