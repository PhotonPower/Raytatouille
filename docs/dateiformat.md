# Systemdateien `.rtt.json` erzeugen

Kurzanleitung für das Dateiformat, das die Engine, die CLI (`rtt`) und der Raytatouille Explorer
laden. Alle Beispiele hier wurden gegen Schema `0.10.x` geprüft (Stand nach Bibliothek 0.9.0, Schema 0.10). Dateien mit
Schema 0.1 bis 0.9 liest `rtt` weiter und schreibt sie mit `rtt format` als 0.10; ein `pickup` aus
einer Datei vor 0.4 fällt dabei weg, mit der Warnung `io.pickup_dropped` (siehe
„Parametertabelle“).
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
| `schema_version` | ja | `"0.10.0"` (Muster `0.10.x`) |
| `name` | nein | Anzeigename |
| `units` | ja | genau `{"length": "mm", "wavelength": "um"}` |
| `environment` | nein | `temperature_c`, `pressure_atm`, `medium` (Standard: Luft nach Ciddor, 20 °C, 1 atm) |
| `object` | nein | `{"at_infinity": true}` ist der Standard; sonst `distance` (mm, Param): Abstand vom Objekt zum **globalen Ursprung** längs −z. Das Objekt liegt also bei z = −`distance` im globalen System, unabhängig davon, wo die erste Fläche liegt; positiv heißt vor dem Ursprung. Soll der Abstand vor einer Fläche gelten, nennt `reference` diese Fläche (ab Schema 0.8, ADR 0038): `"reference": {"surface": "L1.S1"}` legt das Objekt `distance` vor den **nominalen Scheitel** dieser Fläche, längs −z; es zählt nur die z-Koordinate des Scheitels (das Objekt bleibt auf der globalen z-Achse). `"origin"` ist der Standard und wird nicht geschrieben. Ein Bezug nur bei endlichem Objekt (sonst `object.reference_at_infinity`) und auf eine Fläche des Systems (sonst `object.reference_surface_unknown`). Bei der Optimierung folgt das Objekt der Fläche, bei der Toleranzierung bleibt es am nominalen Ort. Optional `medium` (ab Schema 0.7, ADR 0037): Material des Objektraums als `KATALOG:NAME`, nur bei endlichem Objekt (sonst Fehler `object.medium_infinite`). Fehlt es, ist der Objektraum die Umgebung. Eine Quelle in Glas (n > 1) strahlt dann auch jenseits des Grenzwinkels Glas/Luft; `object_na` wird mit dem Index dieses Mediums gelesen (NA darf größer als 1 sein, aber nicht größer als n) |
| `wavelengths` | ja | Liste von `{"um": …, "weight": …, "reference": true}`. **Genau eine** Wellenlänge ist `reference`. `weight` ≥ 0. |
| `aperture` | ja | `{"type": "epd" \| "image_fnumber" \| "object_na" \| "stop_size", "value": …}`. `epd`: Durchmesser der Eintrittspupille in mm. `image_fnumber`: F-Zahl bei unendlichen Konjugierten, EPD = |EFL| / F#, auch bei endlichem Objekt. `object_na`: objektseitige numerische Apertur, nur bei endlichem Objekt, **paraxial gelesen**: die paraxiale Randstrahlsteigung vom axialen Objektpunkt ist u = NA / n, n der Brechungsindex des Objektraums (Greivenkamp, OPTI-502, Abschn. 9, S. 9-34: NA = n sin U ≈ n u; ein realer Randstrahl mit n sin U = NA weicht bei großer NA ab). Bei paraxialer Zielung hat der Randstrahl tan U = NA / n; bei realer Zielung trifft er den paraxialen Blendenrand R_s, und sein tan U weicht um die Pupillenaberration davon ab. Optional `"reading": "real"` (ab 0.8, ADR 0040, nur bei `object_na`; `"paraxial"` ist der Standard und wird nicht geschrieben) liest die NA für die Zielung **real**: Der Randstrahl vom axialen Objektpunkt in der Meridionalebene hat bei der Referenzwellenlänge n_ref sin U = NA, n_ref der Index des Objektraums dort. Die Eintrittspupille der Zielung ist der Kegel mit diesem U, das Ziel auf der Blende der Blendentreffer dieses Strahls bei der Referenzwellenlänge, für alle Wellenlängen derselbe. Das Pupillenraster ist gleichmäßig in den Blendenkoordinaten (reale Zielung) bzw. in tan U (paraxiale Zielung), nicht in sin U. Die paraxialen Größen (`first_order`, Seidel, Bericht, Geister, `image_fnumber`) bleiben bei u = NA / n und beschreiben damit die paraxiale Lesart, bei großer NA also nicht dasselbe Bündel wie die Strahlen. Ist NA ≥ n_ref, wirft jede Zielung den Fehler `aperture.na_unreachable`; erreicht der reale Randstrahl die Blende nicht, wirft ihn die reale Zielung. Bei objektseitiger Telezentrie legen nur `object_na` und `stop_size` das Bündel fest. `stop_size` braucht kein `value` (die Blendengröße kommt aus der Apertur der Blende). Optional `"stop": {"surface": "<id>"}` (ab 0.8, ADR 0039): Die Apertur dieser Fläche ist die Systemblende, z. B. die Fassung einer Linsen- oder Spiegelfläche. Zulässig sind Flächen von `lens`, `mirror`, `plate` und `thin_element`; die Fläche braucht eine Apertur; ein System hat entweder ein `stop`-Element oder `aperture.stop`. Das Ereignis der Fläche behält seine Art (eine Linsenfläche bricht); die Blendenebene ist ihre Scheitelebene, die Austrittspupille das Bild der Blende durch die Fläche selbst und alle folgenden. Bei der Toleranzierung folgt die Blende der gestörten Fläche |
| `fields` | ja | `{"type": "angle_deg" \| "object_height" \| "paraxial_image_height", "points": [{"x":…, "y":…, "weight":…}]}`. `{}` ist der Achspunkt. `weight` ≥ 0. |
| `configurations` | nein | Konfigurationen (ADR 0029), z. B. `[{"name": "wide"}, {"name": "tele"}]`; fehlt der Abschnitt, gibt es nur die Nominalkonfiguration. Nie leer. |
| `parameters` | nein | Parametertabelle (ADR 0029), siehe unten |
| `root` | ja | die oberste Baugruppe |
| `paths` | ja | mindestens ein Pfad, meist `{"name": "main", "events": "auto"}` |
| `optimization` | nein | Merit-Funktion (ADR 0030), siehe unten; fehlt ohne Merit-Funktion |
| `tolerances` | nein | Toleranzsatz (ADR 0035), siehe unten; fehlt ohne Toleranzen |

## Elemente (`type`)

| Typ | Verwendung | Material |
| --- | --- | --- |
| `lens` | Linse; N Flächen ergeben N − 1 Glassegmente (Kittglieder) | `material` |
| `plate` | Platte oder Prisma (Prisma: expliziter Pfad nötig, siehe Pfade) | `material` |
| `mirror` | Spiegel; Standard-Interaktion `fresnel` (siehe unten) | optional `material` (Substrat) |
| `thin_element` | dünnes Element ohne Dicke (Polarisator, Retarder, Strahlteiler) | – |
| `stop` | freistehende Aperturblende; Blende auf einer Fläche siehe `aperture.stop` | – |
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
| `shape.base` | `{"type": "plane"}`, `{"type": "conic", "radius": R, "conic": k}` oder `{"type": "even_asphere", "radius": R, "conic": k, "coefficients": [A4, A6, …]}`; ab Schema 0.9 (ADR 0042) die anamorphen Formen `{"type": "biconic", "radius_x": R_x, "radius_y": R_y, "conic_x": k_x, "conic_y": k_y}`, `{"type": "toroid", "radius": R, "conic": k, "coefficients": [A4, A6, …], "radius_of_rotation": R_rot}` und `{"type": "cylinder", "radius": R, "conic": k, "coefficients": [A4, A6, …]}` sowie (ADR 0043) `{"type": "odd_asphere", "curvature": c, "conic": k, "coefficients": [A1, A2, A3, …]}`, `{"type": "axicon", "angle_deg": α}`, `{"type": "qcon", "curvature": c, "conic": k, "normalization_radius": ρ, "coefficients": [a0, a1, …]}` und `{"type": "qbfs", …}` mit denselben Schlüsseln (siehe unten) |
| `aperture` | `circular` (`radius`, optional `inner_radius`), `rectangular` (`half_width_x`, `half_width_y`) oder `elliptical` (`semi_axis_x`, `semi_axis_y`) |
| `phases` | Phasenschichten: `{"type": "linear_grating", "lines_per_mm": G, "orientation_deg": ψ}` oder `{"type": "radial_phase", "normalization_radius": R, "coefficients": [c1, c2, …]}` (ADR 0025). Einheiten: G in 1/mm, ψ in Grad von der lokalen x- zur lokalen y-Achse (ψ = 0: Rillen parallel zur lokalen y-Achse), φ = 2πG(x cos ψ + y sin ψ); R in mm, Koeffizienten c_k in **rad** (nicht in Wellen, kein Faktor 2π), c1 gehört zu ρ², φ = Σ c_k ρ^{2k} mit ρ = r/R. Die Werte gelten unverändert für jede Wellenlänge; das Vorzeichen von φ und der Ordnung nach ADR 0025, Punkte 1 bis 3 |
| `diffraction_efficiency` | nur an Flächen mit `phases`: Leistungsanteil je Beugungsordnung, z. B. `[{"order": 1, "efficiency": 0.8}]`; fehlt das Feld, hat jede Ordnung 1, sonst haben nicht aufgeführte Ordnungen 0 |
| `interaction` | `fresnel` (Standard), `ideal_mirror`, `ideal_anti_reflection`, `absorber`, `ideal_beam_splitter` (`reflectance_s`, `reflectance_p`), `ideal_polarizer` (`transmission_axis`, `extinction_ratio`), `ideal_retarder` (`fast_axis`, `retardance_waves`), `ideal_lens` (`focal_length`, `object_distance`), `ideal_cylinder_lens` (`focal_length`, `axis_deg`, `object_distance`; ab Schema 0.5, siehe unten) oder `{"type": "coating", "name": "KATALOG:NAME"}`. Die Achsen von Polarisator und Retarder stehen in Elementkoordinaten. |

Zahlenwerte dürfen auch als `{"value": …, "variable": true, "min": …, "max": …}` stehen: Das
markiert sie als Optimierungsvariable mit Grenzen (beide optional, in der Einheit des Werts). Ein
Zahlenwert kann stattdessen an eine Zeile der Parametertabelle gebunden sein,
`{"param": "NAME"}`, ohne eigenen Wert (siehe unten). Zernike-Terme (`shape.terms`) kennt das
Schema, die Engine kompiliert sie aber noch nicht, bis M8 (Meldung: "not supported before
M8"). Gitter und Phasenflächen (`phases`) übernimmt die Kompilierung; der Tracer verfolgt seit #127
jede Beugungsordnung `order` nach der lokalen Gittergleichung (ADR 0025) und rechnet
`diffraction_efficiency` in das Gewicht ein; evaneszente Ordnungen enden mit dem Status
`Evanescent`. Kristalle mit `ordinary` und `extraordinary` verfolgt der Tracer seit #132
(ADR 0026).

## Anamorphe Flächen (ab Schema 0.9, ADR 0042)

Biconic, Toroid und Zylinder haben in den Schnitten x-z und y-z verschiedene Krümmungen. Sie
gelten im lokalen System der Fläche; die Drehung der Fläche um z legt die Richtung der Schnitte
fest. Alle Radien sind in mm, endlich und ≠ 0 (sonst `shape.radius_invalid`), mit dem
Vorzeichen wie bei `conic`: positiv, wenn der Krümmungsmittelpunkt bzw. die Drehachse auf der
+z-Seite des Scheitels liegt. Alle Werte sind Params.

```json
{"type": "biconic", "radius_x": 40.0, "radius_y": -60.0, "conic_x": -0.5}
{"type": "toroid", "radius": 48.0, "coefficients": [2e-6], "radius_of_rotation": 55.0}
{"type": "cylinder", "radius": -75.0}
```

- `biconic`: z = (c_x x² + c_y y²) / (1 + √(1 − (1 + k_x) c_x² x² − (1 + k_y) c_y² y²)) mit
  c_x = 1/`radius_x`, c_y = 1/`radius_y`; `conic_x` und `conic_y` sind optional (Standard 0).
  Asphärenterme gibt es nicht.
- `toroid`: Das Profil z_y(y) in der y-z-Ebene ist eine gerade Asphäre wie `even_asphere`
  (`radius`, optional `conic` und `coefficients` A4, A6, …). Es dreht sich um die Achse parallel
  zu y durch (x, z) = (0, `radius_of_rotation`). Mit `radius_of_rotation` = `radius`, ohne Konik
  und Koeffizienten ist das eine Kugel.
- `cylinder`: das Profil des Toroids ohne Drehung, also flach in x (`radius`, optional `conic`
  und `coefficients`).

Die rotationssymmetrische Paraxialrechnung (`first_order`, Seidel, Bericht) lehnt eine Fläche
mit verschiedenen Schnitten ab; ein Biconic mit gleichen Schnitten und ein Toroid als Kugel
rechnen dort wie `conic`. Die Toleranzklassen erzeugen an diesen Typen weder
`curvature_relative` noch `conic` (Grund `ShapeType`); die Toleranzarten `radius`,
`curvature_relative` und `conic` sind an ihnen `tolerances.target_invalid`.

## Ungerade Asphäre und Axicon (ab Schema 0.9, ADR 0043)

Beide sind rotationssymmetrisch um die lokale z-Achse; alle Werte sind Params. Wie die
anamorphen Formen sind sie in einer Datei vor Schema 0.9 ein Lesefehler.

```json
{"type": "odd_asphere", "curvature": 0.02, "conic": -0.5, "coefficients": [0.0, 1e-4, 2e-6]}
{"type": "axicon", "angle_deg": -5.0}
```

- `odd_asphere` (OpticStudio „Odd Asphere“): z = c r²/(1 + √(1 − (1 + k) c² r²)) + A1 r +
  A2 r² + A3 r³ + … mit r² = x² + y², A_n in mm^(1−n).
  - Die Basis trägt ihre **Krümmung** `curvature` c in mm⁻¹, keinen Radius: c ist endlich
    (sonst `shape.curvature_invalid`), 0 ist eine ebene Basis, und c bleibt als Variable auch
    durch 0 hindurch stetig. Bei einem Import ist c = 1/R; ein ZMX-Radius 0 heißt c = 0.
  - `conic` (Standard 0) wirkt bei c = 0 nicht und gibt dann die Warnung
    `shape.conic_without_curvature`. `coefficients` beginnen bei A1 und schließen A2 ein; eine
    leere Liste gibt `shape.asphere_without_coefficients`.
  - Die Scheitelkrümmung ist c + 2 A2, nicht c. Mit A1 ≠ 0 ist der Scheitel eine Kegelspitze
    ohne paraxiales Modell (`first_order` wirft); mit A3 ≠ 0 gibt es die erste Ordnung, aber
    keine Seidel-Summen.
- `axicon`: der Kegel z = r tan α mit der Spitze im Scheitel. `angle_deg` ist der Basiswinkel α
  zwischen Mantel und der Ebene z = 0, mit dem Vorzeichen einer Krümmung: α > 0 heißt, der
  Mantel steigt nach +z; α = 0 ist die Ebene. |α| < 90° (sonst `shape.axicon_angle_invalid`).
  Ein Kegel mit Konik ist eine `odd_asphere` mit A1 = tan α.
  - Die Spitze hat kein paraxiales Modell: Zielung, Pupillen, OPD und die paraxialen Operanden
    werfen für einen Pfad mit Axicon; eigene Startstrahlen mit `trace` laufen.

Toleranzen (ADR 0035, Nachtrag #314): An `odd_asphere` stört `curvature_relative` die Krümmung,
C′ = C (1 + δ), und `conic` die Konik; `radius` ist dort `tolerances.target_invalid`. Die
Klassen überspringen beide bei c = 0 mit `Plane`. Eine Krümmung, die 0 erreicht oder das
Vorzeichen wechselt, macht den Prüfling nicht ungültig; deshalb ist an einer Krümmung auch eine
Grenze `curvature_relative` ≥ 1 erlaubt (ADR 0043, Nachtrag PR C), an einem Radius nicht. Am Axicon sind `radius`,
`curvature_relative` und `conic` `tolerances.target_invalid` (Klassen: Grund `ShapeType`); der
Winkel geht über die Art `param` auf `/shape/base/angle_deg`, und ein Prüfling, dessen Winkel
das Vorzeichen wechselt oder |α| ≥ 90° erreicht, ist ungültig.

## Forbes Q-con und Q-bfs (ab Schema 0.9, ADR 0043)

Rotationssymmetrische Asphären nach G. W. Forbes mit orthogonalen Polynomen in u² und
u = r/ρ (`docs/quellen.md`). Die Basis trägt wie bei `odd_asphere` ihre Krümmung `curvature`
c in mm⁻¹ (endlich, 0 ist eine ebene Basis), `conic` k ist optional, `normalization_radius` ρ
in mm ist Pflicht (endlich und > 0, sonst `shape.normalization_radius_invalid`), die
`coefficients` a₀, a₁, … sind in mm. Alle Werte sind Params.

```json
{"type": "qcon", "curvature": 0.015625, "conic": -0.8, "normalization_radius": 12.7,
 "coefficients": [-0.05, 0.01, -0.002]}
{"type": "qbfs", "curvature": 0.0166667, "conic": -0.3, "normalization_radius": 12.7,
 "coefficients": [0.1, -0.02]}
```

- `qcon` (Forbes 2007, Gl. (2), (3), (8)): z = c r²/(1 + φ) + u⁴ Σ a_m Q^con_m(u²) mit
  φ = √(1 − (1 + k) c² r²). Forbes 2007 schreibt die Konik als „conic parameter“ ε = 1 + k.
  Die Scheitelkrümmung ist c.
- `qbfs` (Forbes 2011, Gl. (2.4)): z = c r²/(1 + φ) + u²(1 − u²)/σ · Σ a_m Q^bfs_m(u²) mit
  1/σ = √(1 − k c² r²)/φ. Mit k = 0 ist c die Krümmung der Best-Fit-Sphäre, die die Fläche im
  Scheitel und bei u = 1 trifft (Forbes 2007, Gl. (10), (11)). Die Fläche ist nur für φ > 0
  definiert; (1 + k) c² ρ² ≥ 1 ist `shape.normalization_radius_invalid`. Die Scheitelkrümmung
  ist c + 2 S(0)/ρ² mit S(0) = Σ a_m Q_m(0), nicht c.
- ρ ist bei Forbes die halbe freie Apertur. Jenseits von ρ (u > 1) werden die Polynome
  extrapoliert; reicht die Apertur der Fläche weiter als ρ (Radius, halbe Diagonale bzw. große
  Halbachse, mehr als 10⁻⁹ relativ), warnt `shape.normalization_radius_exceeded`.
- Toleranzen wie bei `odd_asphere`: `curvature_relative` stört c, `conic` stört k, `radius` ist
  `tolerances.target_invalid`; ρ und die Koeffizienten gehen über `param`.

## Ideale Linse und ideale Zylinderlinse (ab Schema 0.5, ADR 0031)

Eine ideale Linse ohne Dicke ist eine Interaktion an der einzigen Planfläche eines
`thin_element` (ohne `shape.terms` und ohne `phases`, sonst `interaction.ideal_lens_not_allowed`):

```json
{"type": "thin_element", "name": "IL",
 "surfaces": [{"id": "IL",
               "aperture": {"type": "circular", "radius": 12.0},
               "interaction": {"type": "ideal_lens", "focal_length": 50.0,
                               "object_distance": 75.0}}]}
```

- `focal_length` (Pflicht, mm, Param): Abstand von der Linse zum hinteren Brennpunkt im
  umgebenden Medium, längs der Ausbreitung; f > 0 sammelt, f hängt nicht von der Wellenlänge ab.
  Variabel, mit Grenzen und an die Parametertabelle bindbar.
- `object_distance` (optional, mm, Param): Auslegungs-Konjugierte für den optischen Weg,
  positiv für ein reelles Objekt vor der Linse; das Vorzeichen ist dasselbe wie bei
  `object.distance`, gemessen wird aber von der Linse aus, nicht vom globalen Ursprung (ADR 0031,
  Nachtrag #279). Fehlt das Feld, liegt das Objekt im Unendlichen.
- `ideal_cylinder_lens` hat dazu `axis_deg` (Standard 0, im kanonischen Text bei 0 weggelassen):
  die Zylinderachse in lokalen Flächenkoordinaten von x nach y. Die Linse fokussiert nur in der
  Richtung (−sin ψ, cos ψ, 0); bei 0 also im yz-Schnitt, die Brennlinie liegt parallel zu x.
- An der Linse wirkt nur `transmit`; `refract` oder `reflect` dort meldet compile mit
  `paths.ideal_lens_event`, eine Linse in einem Kristall mit `crystal.unsupported`.
- **Grenze:** Die Richtung ist für jeden Strahl exakt (kollineare Abbildung). Der optische Weg,
  also die OPD, ist nur für Objekte in der Ebene von `object_distance` exakt; bei einer anderen
  Konjugierten bleibt ein Rest ∝ r⁴, im Beispiel f = 50 mm, Objekt bei 75 mm statt im
  Unendlichen etwa 11 λ am Rand. Seidel zeigt für die ideale Linse trotzdem 0. Auf einem
  Rückweg gilt die Konjugierte gegen dessen Ausbreitung.
- Paraxial ist die ideale Linse ein Schritt mit φ = |n|/f, n der Index des umgebenden Mediums:
  `rear_focal_length` = f, EFL = f/|n|, Seidel-Beiträge 0. Die Zylinderlinse ist paraxial ein
  `ParaxialError` (nicht rotationssymmetrisch); ihre Pfade brauchen eigene Startstrahlen.

## Bauteile aus der Rechnung nehmen: `active` (ab Schema 0.10, ADR 0049)

Jede Baugruppe und jedes Element kann `"active": false` tragen. Der Knoten bleibt im Modell und
behält seine Geometrie (Lage, relative Bezüge, Layout), wird aber aus der Rechnung genommen;
mit `false` an einer Baugruppe gilt das für alle ihre Kinder. Standard ist `true`, die kanonische
Form schreibt den Schlüssel nur bei `false`, die Bearbeitungsform immer (ein Patch
`{"op": "replace", "path": "/root/children/2/active", "value": false}` schaltet um). An der
Wurzel ist `false` der Fehler `node.root_inactive`. Eine Datei vor 0.10 mit dem Schlüssel ist ein
Lesefehler.

```json
{"type": "plate", "name": "Filter", "active": false, "material": "SCHOTT:N-BK7", "surfaces": [ … ]}
```

Wie Pfade, Blende, Optimierung, Toleranzierung und Export inaktive Knoten behandeln, legt
ADR 0049 fest; die Rechnung beachtet das Flag ab der nächsten Ausbaustufe von #380 (bis dahin
lehnt der Export Systeme mit inaktiven Knoten ab).

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
In einer 0.4-Datei ist `pickup` ein Fehler. compile wertet die Tabelle für die gewählte
Konfiguration aus (`compile(system, materials, coatings, configuration)`, #165); ein Index jenseits
der Konfigurationen ist `config.unknown`.

## Merit-Funktion `optimization` (ab Schema 0.4, ADR 0030)

Der Abschnitt steht nach `paths` und hat zwei Listen, `operands` und `generators`; eine leere Liste
und der ganze Abschnitt ohne Einträge fallen weg. Ein Operand liefert einen Wert v, sein Residuum ist
√w·(v − t). Ein Generator erzeugt Residuen je Strahl mit dem Ziel 0. Die Einstellungen des Lösers
(Toleranzen, Iterationen) stehen nicht in der Datei, sie sind Argumente des Laufs.

Gemeinsame Schlüssel:

- `target` (Pflicht bei Operanden, Einheit des Werts; Generatoren haben kein `target`);
- `weight` (Standard 1, ≥ 0);
- `configuration` (Name einer Konfiguration; fehlt: Konfiguration 0, also die erste bzw. die
  Nominalkonfiguration);
- `path` (Name eines Pfads; nicht bei `param_value`).
- **Indizes:** `field` ist ein Index in `fields.points` (Standard 0), `wavelength` ein Index in
  `wavelengths`; fehlt `wavelength`, gilt die Referenzwellenlänge.
- **Pupillenkoordinaten** `px`, `py` (Standard 0) sind normiert wie bei `make_rays`: Der
  Einheitskreis ist der Rand der paraxialen Eintrittspupille, +y ist meridional.

| `type` | Wert | Eigene Schlüssel (Standard) |
| --- | --- | --- |
| `efl`, `bfl` | paraxiale Brenn- bzw. Schnittweite, mm | `wavelength` |
| `image_fnumber`, `magnification` | bildseitige F-Zahl bzw. Abbildungsmaßstab, dimensionslos | `wavelength` |
| `ray_x`, `ray_y` | x bzw. y eines realen Strahls im lokalen KS einer Fläche, mm | `surface` (Pflicht), `occurrence`, `field`, `px`, `py`, `wavelength` |
| `spot_rms` | RMS-Spotradius eines Felds, mm | `field`, `wavelength`, `polychromatic` (false), `reference` (`centroid` oder `chief`), `rings` (6, hexapolar) |
| `opd_rms` | RMS der Wellenfront eines Felds, Wellen | `field`, `wavelength`, `grid` (33) |
| `param_value` | Wert einer Zeile der Parametertabelle | `parameter` (Pflicht) |

Generatoren `rms_spot` (zwei Residuen je Strahl, Bezug `reference`) und `rms_wavefront` (eines je
Strahl, in Wellen) tasten die Pupille mit Gauß-Quadratur ab: `rings` (3) Ringe in ρ², `arms` (6) Arme
je Ring. `fields` und `wavelengths` wählen Felder bzw. Wellenlängen per Index; fehlen sie, gelten alle
mit ihren Gewichten. Eine leere Liste ist ein Fehler.

`occurrence` (0-basiert) wählt bei `ray_x`/`ray_y` das Ereignis, wenn der Pfad die Fläche mehrmals
trifft (Doppeldurchgang, Geisterbild); ein automatischer Pfad trifft jede Fläche einmal.
`polychromatic: true` nimmt alle Wellenlängen mit ihren Gewichten und verträgt kein `wavelength`.
Die Namen, Indexbereiche, Gewichte und Abtastungen prüft `rtt validate` (Codes `merit.*`).

```json
"optimization": {
  "operands": [
    {"type": "efl", "path": "main", "target": 100.0},
    {"type": "ray_y", "path": "main", "surface": "IMG", "field": 1, "py": 1.0, "target": 0.0},
    {"type": "param_value", "parameter": "G", "configuration": "tele", "target": 5.0, "weight": 0.1}
  ],
  "generators": [{"type": "rms_spot", "path": "main", "fields": [0, 1]}]
}
```

## Toleranzsatz `tolerances` (ab Schema 0.6, ADR 0035)

Der Abschnitt steht nach `optimization`; alle Schlüssel sind optional, leere Listen und ein Satz
ohne Einträge fallen weg. Er beschreibt, was toleriert wird. Die Analyse (Sensitivität, Monte
Carlo) und ihre Laufparameter (Threads, Perzentile, LM-Optionen) stehen nicht in der Datei
(ADR 0036).

| Schlüssel | Inhalt |
| --- | --- |
| `class` | Toleranzklasse des Systems: `coarse`, `normal`, `precise` oder `none` (wie fehlend) |
| `overrides` | Überschreibungen je Knoten oder Fläche: Klasse, Einzelwert oder Pivot |
| `items` | explizite Toleranzen |
| `compensators` | Kompensatoren mit Verstellweg |
| `compensator_merit` | Merit-Funktion der Kompensatoren, aufgebaut wie `optimization`; fehlt sie, gilt `optimization` |
| `criteria` | Kriterien: Operanden mit Grenzen |
| `monte_carlo` | Vorgaben `seed` und `trials` |

**Ziele** werden über Anker adressiert, nie über Indizes: genau einer von `surface` (Flächen-`id`),
`node` (Name eines Elements oder einer Baugruppe, nicht der Wurzel) und `parameter` (Name einer
Zeile der Parametertabelle). Ein Param wird als Anker plus angehängter Pointer `param` adressiert,
relativ zum Objekt des Ankers in der Bearbeitungsform, z. B. `{"surface": "L1.S2", "param":
"/shape/base/coefficients/0"}` oder `{"node": "image", "param": "/pose/position/2"}`. So bleiben
Ziele gültig, wenn Knoten eingefügt oder verschoben werden.

**Items** haben `id`, ein Ziel, `kind`, `limit` (Grenze T > 0; die Störung liegt in [−T, T]) und
`distribution`: `normal` (Standard, σ = T/2, abgeschnitten bei ±T), `uniform` oder `endpoint`
(±T). Die `id` bestimmt die Zufallsfolge eines Items; sie darf nicht mit `class:` beginnen (das
ist der Namensraum der erzeugten Items).

| `kind` | Ziel | Einheit | Wirkung |
| --- | --- | --- | --- |
| `radius` | Fläche mit Konik oder gerader Asphäre | mm | R′ = R + ΔR |
| `curvature_relative` | wie `radius` | dimensionslos | C′ = C (1 + δ) |
| `conic` | wie `radius` | dimensionslos | κ′ = κ + Δκ |
| `thickness` | Fläche oder Knoten | mm | Offset auf `pose.position[2]`; relativ platzierte Folgeknoten wandern mit |
| `param` | Anker + `param` | Einheit des Params | Offset auf den Wert des Params (z. B. Asphärenkoeffizient) |
| `parameter` | Zeile | Einheit der Zeile | Offset auf die Zeile, in allen Konfigurationen derselbe |
| `decenter_x`, `decenter_y`, `shift_z` | Fläche oder Knoten | mm | starre Verschiebung im Nennrahmen des Ziels und seines Teilbaums |
| `tilt_x_deg`, `tilt_y_deg`, `tilt_z_deg` | Fläche oder Knoten | Grad | starre Drehung um die Achse des Nennrahmens durch den Pivot |
| `index`, `abbe` | Element mit isotropem Material (Klassen: nur Lens und Plate) | dimensionslos | Δn_d bzw. Δν_d; `segment` (0-basiert) bei einer Liste `material`; eine Liste aus lauter gleichen Einträgen ist ein Material, dort gilt nur `segment` 0 |

`thickness` ist ein Maß (Offset auf ein Param), `shift_z` eine Lage (starr): Eine Dickentoleranz
schiebt relativ platzierte Folgeknoten mit, eine starre Verschiebung nicht. Params auf
Systemebene (Objektabstand, Systemapertur) sind keine Toleranz- und keine Kompensatorziele,
ebenso eine Zeile der Parametertabelle, von der sie abhängen (direkt gebunden oder über Ausdrücke).

**Klassen und Überschreibungen.** Die Klasse erzeugt Items aus einer festen Tabelle der
Bibliothek (Flächen: `curvature_relative`, `conic` an Asphären, `thickness` ab der zweiten
Fläche; Elemente: `thickness`, `decenter_x/y`, `tilt_x/y_deg`, bei Lens und Plate `index` und
`abbe`; Baugruppen: `decenter_x/y`, `shift_z`, `tilt_x/y_deg`). Ein Eintrag in `overrides` hat
genau eine von drei Formen:

- `{"node": …, "class": …}`: andere Klasse für den Knoten und seinen Teilbaum; `none` schaltet die
  Erzeugung dort ab. Der nächste Vorfahr gewinnt.
- `{"node" | "surface": …, "kind": …, "limit": …}` (optional `distribution`, `segment`): ein
  Einzelwert für eine erzeugte Art an genau diesem Ziel.
- `{"node" | "surface": …, "pivot": [x, y, z]}`: der Pivot der Kippungen dieses Ziels (mm, in
  seinen Koordinaten), nicht vererbt; ohne Eintrag gilt der Pivot seiner Pose.

Vorrang: explizites Item > Einzelwert > Klasse des nächsten Vorfahren > `class` des Systems,
verglichen wird die gestörte Größe (ein `radius`-Item ersetzt das erzeugte `curvature_relative`
derselben Fläche).

**Kompensatoren** sind ein ungebundenes Param (Anker plus `param`) oder eine Zeile der Form
`value` (`parameter`), mit `travel_min` < 0 < `travel_max` relativ zum Wert des Prüflings vor der
Kompensation. **Kriterien** haben `id`, `operand` (ein Operand wie in `optimization`, aber ohne
`target` und `weight`) und mindestens eine Grenze: `min`, `max` (absolut) oder `max_change`
(|Wert − Nennwert|). **`monte_carlo`**: `seed` ganzzahlig in [0, 2⁵³ − 1], `trials` in [1, 2³²].

Was zusammenpasst (Art und Ziel, Segment, Grenzen, ids, Kompensatoren, Bereiche), prüft
`rtt validate` (Codes `tolerances.*`); der Leser prüft nur Schlüssel und Typen. Eine Datei vor
Schema 0.6 mit `tolerances` ist ein Lesefehler.

```json
"tolerances": {
  "class": "normal",
  "overrides": [
    {"node": "objective", "class": "precise"},
    {"surface": "L1.S2", "kind": "curvature_relative", "limit": 0.005},
    {"node": "objective", "pivot": [0.0, 0.0, 9.5]}
  ],
  "items": [
    {"id": "S3-wedge", "surface": "L1.S3", "kind": "tilt_x_deg", "limit": 0.02},
    {"id": "bfd", "parameter": "BFD", "kind": "parameter", "limit": 0.1, "distribution": "endpoint"}
  ],
  "compensators": [{"parameter": "BFD", "travel_min": -1.0, "travel_max": 1.0}],
  "criteria": [
    {"id": "spot-edge", "operand": {"type": "spot_rms", "path": "main", "field": 2}, "max": 0.03},
    {"id": "efl", "operand": {"type": "efl", "path": "main"}, "max_change": 0.5}
  ],
  "monte_carlo": {"seed": 20261012, "trials": 200}
}
```

Vollständiges Beispiel: `tests/reference/m7/achromat_tolerances.rtt.json`.

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
  "schema_version": "0.10.0",
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
        "schema_version": "0.10.0", "name": f"Singlet R1={r1} R2={r2}",
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

## Flächentabelle (CSV, ADR 0045, Punkt 7)

Eine Linsentabelle für Tabellenprogramme. Sie enthält nur die Flächen eines koaxialen Systems;
Wellenlängen, Felder, Systemapertur und Umgebung bleiben in der `.rtt.json`-Datei. Beim Einlesen
ersetzt die Tabelle Baum, Objekt und Blende eines Basissystems.

```powershell
rtt export system.rtt.json linsen.csv --catalog schott.agf          # schreiben
rtt import linsen.csv neu.rtt.json --base system.rtt.json --catalog schott.agf
```

C++: `rtt::io::write_surface_table`, `read_surface_table` in `rtt/io/surface_table.hpp`.

**Format** (RFC 4180, Abschnitt 2):
- UTF-8, Komma als Trenner, jede Zeile mit gleich vielen Feldern.
- Geschrieben wird CRLF, gelesen werden CRLF und LF; ein BOM wird übergangen.
- Felder mit Komma, Anführungszeichen oder Zeilenumbruch stehen in Anführungszeichen, ein
  Anführungszeichen darin wird verdoppelt.
- Leerzeichen gehören zum Feld: `" 50"` ist keine Zahl.
- Zahlen in der kürzesten Form, die beim Lesen denselben Wert gibt.

**Kopfzeile:** genau diese Spalten, gelesen in beliebiger Reihenfolge. Eine fehlende, unbekannte
oder doppelte Spalte ist der Fehler `table.column`.

| Spalte | Inhalt |
|---|---|
| `row` | `object` (erste Zeile nach der Kopfzeile) oder `surface` |
| `surface` | Flächen-ID |
| `element`, `kind` | Elementname und Art (`lens`, `plate`, `mirror`, `thin_element`, `stop`, `detector`); aufeinanderfolgende Zeilen mit demselben Namen bilden ein Element |
| `shape` | `plane`, `conic`, `even_asphere` |
| `radius`, `conic` | R in mm, Konik |
| `a4` … `a16` | Asphärenkoeffizienten in mm^(1−n) |
| `thickness` | Abstand zum Scheitel der nächsten Zeile längs z in mm, nach einem Spiegel negativ; beim Objekt der Abstand zur ersten Fläche oder `Infinity`; in der letzten Zeile leer |
| `material` | Medium hinter der Fläche (`KATALOG:NAME`, `CONST:…`, leer = Umgebung); bei Spiegeln, Blenden, dünnen Elementen, Detektoren und der letzten Fläche einer Linse leer; beim Objekt das Objektmedium |
| `aperture`, `aperture_1`, `aperture_2` | `circular` (Radius, Innenradius), `rectangular` bzw. `elliptical` (Halbbreiten x, y) oder leer |
| `stop` | `true` an der Blende, sonst leer |
| `focal_length`, `object_distance` | ideale Linse eines dünnen Elements, mm |

**Regeln:**
- Die Objektzeile enthält nur `thickness` und `material`.
- Ein leeres `conic`, `a4` … `a16` oder ein leerer Innenradius ist 0, und 0 wird dort leer
  geschrieben. Ein leerer `radius` gehört zu `plane`.
- `stop` = `true` an einer Zeile der Art `stop` macht das Element zur Blende, an jeder anderen
  Zeile setzt es `aperture.stop` auf ihre ID.
- Beim Einlesen entsteht die kanonische Form: flacher Baum, erstes Element bei z = 0, alle
  weiteren Elemente und Flächen `relative_to_preceding` mit (0, 0, t), ein automatischer Pfad.
  Spiegel werden `ideal_mirror`.
- Die Tabelle hält nur Werte: Optimierungsvariablen und Grenzen der Parameter, Pivot und
  Reihenfolge der Posen und ein ausdrücklich gesetztes Objektmedium gleich der Umgebung
  entfallen ohne Warnung.
- Eine leere Zeile am Ende der Datei wird übergangen.
- Was im Basissystem auf den alten Baum verweist, entfällt mit der Warnung `table.base_dropped`:
  weitere Pfade, `optimization`, `tolerances`, `parameters`, `configurations`. Gebundene Werte des
  Basissystems behalten den Wert der Konfiguration 0.
- **Fehler** tragen die Zeile, in der das Feld beginnt, und die Nummer des Feldes (`table.parse`,
  `table.column`, `table.value`, `table.duplicate`), siehe `docs/diagnostics.md`. Am Ende laufen
  `validate` und `compile`; ein Fehler dort wird `table.invalid_result`.
- **Schreiben:** nur koaxiale Systeme, mit denselben Fehlern `export.*` wie der ZMX-Export.
  Anders als dort sind `CONST:`, `VACUUM`, jede Umgebung und eine ideale Linse mit
  `object_distance` erlaubt. Eine `odd_asphere` und andere Formen sind `export.surface_unsupported`.
- **Excel mit deutscher Einstellung** schreibt `;` und Dezimalkomma; das gibt `table.parse` mit
  einem Hinweis. Als CSV mit Komma und Dezimalpunkt speichern. Tabellenprogramme speichern oft nur
  15 Stellen; das ist eine Grenze des Programms, nicht des Formats.

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
