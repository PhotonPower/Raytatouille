# Lehrbeispiele für Optikdesign

Sechzehn Systemdateien (`.rtt.json`) in acht Themen. Jede Datei zeigt genau einen Zusammenhang, den man mit
dem Raytatouille Explorer (oder der Python-API) nachmessen kann. Die Zahlen unten wurden mit
raytatouille 0.4.0 gemessen und lassen sich in den Ansichten des Explorers nachprüfen.

> Hinweis: Der Raytatouille Explorer (Streamlit-Oberfläche `rtt_explorer.py`) ist nicht Teil dieses Repos. Die
> Beispiele funktionieren auch mit der Bibliothek und der CLI allein (`rtt validate`, `rt.load`, `rt.compile`).

## So benutzt man die Beispiele

1. Die Dateien liegen in `examples/lehrbeispiele/`. Im Explorer: **System → Datei hochladen**.
2. `02a` und `02b` brauchen die Gläser N-BK7 und F2 (`tests/catalogs/schott.agf`). Der Explorer wählt den Katalog selbst,
   wenn der Repo-Pfad gesetzt ist; sonst einen eigenen Schott-AGF-Katalog in der Seitenleiste hochladen.
3. Alle Längen sind in mm, Wellenlängen in µm. Alle Systeme haben eine **Blende** (außer `08`), die Bildebene liegt
   im **paraxialen Fokus**.

**Messweise.** RMS-Werte sind die „RMS um Hauptstrahl“ der Spot-Ansicht (hexapolar, 14 Ringe, 631 Strahlen), in µm,
am paraxialen Fokus und nicht am besten Fokus. Das ist Absicht: So bleibt der Vergleich zwischen den Dateien fair, und die
Studierenden können die Bildebene selbst verschieben (`pose.position` des Detektors in der Datei).

| Thema | Dateien | Lernziel |
| --- | --- | --- |
| 01 | `01a` bis `01d` | Sphärische Aberration hängt von der **Form** der Linse ab |
| 02 | `02a`, `02b` | Farblängsfehler und Achromat |
| 03 | `03` | Eine Asphäre (Hyperbel) beseitigt die sphärische Aberration |
| 04 | `04a`, `04b` | Kugel- und Parabolspiegel, Koma |
| 05 | `05a`, `05b` | Cassegrain und Ritchey-Chrétien, Seidel-Koeffizienten |
| 06 | `06a`, `06b` | Teleobjektiv: kurze Baulänge bei langer Brennweite |
| 07 | `07a`, `07b` | Blendenlage bestimmt die Verzeichnung |
| 08 | `08` | Brewster-Winkel und Polarisation |

---

## 01 Die Form der Linse (sphärische Aberration)

**Dateien.** Vier Linsen mit gleicher Brennweite (f ≈ 100 mm), gleichem Glas (n = 1,5168, ohne Dispersion), gleicher
Öffnung (EPD 20 mm, f/5) und gleicher Dicke (4 mm). Nur die Radien unterscheiden sich („Biegung“).

| Datei | R₁ (mm) | R₂ (mm) | RMS auf der Achse (µm) | RMS bei 3° (µm) |
| --- | --- | --- | --- | --- |
| `01a_biegung_plan_zuerst` | ∞ (plan) | −51,68 | **248,5** | 282,4 |
| `01b_biegung_symmetrisch` | 100,00 | −105,50 | 87,3 | 117,1 |
| `01c_biegung_gewoelbt_zuerst` | 51,68 | ∞ (plan) | 61,6 | 90,0 |
| `01d_biegung_bestform` | 58,82 | −415,88 | **57,1** | 83,3 |

**Beobachten.** Ansicht *Spot* (Feld 0) und *Ray Fans*, dazu *Seidel* (S_I).

**Fragen.**
- Welche Seite einer Plankonvex-Linse muss zur Lichtquelle zeigen? (Vergleich `01a` mit `01c`, Faktor ≈ 4.)
- Wie unterscheidet sich die Bestform (`01d`) von der Plankonvex-Linse `01c`? Vergleiche die Radien und die Seidel-Beiträge je Fläche.
- Verschiebe in `01a` die Bildebene: Der kleinste Fleck liegt nicht im paraxialen Fokus. Gemessen: bei −3,0 mm (also 3,0 mm vor dem paraxialen Fokus) nur noch 82,7 µm statt 248,5 µm. Warum liegt der beste Fokus **vor** dem paraxialen? (Randstrahlen werden stärker gebrochen.)

## 02 Farbfehler und Achromat

**Dateien.** Beide f ≈ 100 mm, EPD 20 mm, Wellenlängen F, d, C (0,4861 / 0,5876 / 0,6563 µm). `02a` ist ein Singlet
aus N-BK7 in Bestform. `02b` ist ein verkitteter Achromat aus N-BK7 und F2, dessen Radien so berechnet sind, dass
der paraxiale Farblängsfehler zwischen F und C null wird.

| | `02a` Singlet | `02b` Achromat |
| --- | --- | --- |
| Radien (mm) | 58,82 / −413,41 | 44,39 / −42,00 / −830,99 |
| Brennweite (mm) | 99,97 | 100,01 |
| Farblängsfehler paraxial (F–C) | **1,53 mm** | **0,2 µm** (≈ 0 per Konstruktion) |
| Farblängsfehler real (Randzone) | 1,52 mm | 0,135 mm (Restfehler, „Sphärochromasie“) |
| RMS auf der Achse (µm) | 85,6 | **4,3** |
| RMS bei 2,5° (µm) | 101,2 | 39,7 |
| RMS bei 5° (µm) | 154,0 | 110,8 |

**Beobachten.** Ansicht *Farbfehler* (Paar F–C), *Layout* mit „alle (Farben)“, *Spot* polychromatisch.

**Fragen.**
- Die Faustregel „Farblängsfehler ≈ f / ν“ ergibt für N-BK7 (ν_d = 64,17) etwa 1,56 mm. Passt das zu `02a`?
- Der Achromat ist auf der Achse mit 4,3 µm RMS in der Größenordnung des Airy-Radius (3,6 µm bei f/5, nur ein grober Vergleichswert), bei 5° aber 26-mal schlechter. Sieh dir dazu die Seidel-Summen und die Ray Fans bei 5° an.
- Der Farblängsfehler zwischen F und C ist null. Probiere in der Ansicht *Farbfehler* andere Wellenlängenpaare und die Randzone (Pupillenzone 1,0): Was bleibt übrig?

## 03 Die Descartes-Linse (Asphäre)

**Datei.** `03_descartes_hyperbel` ist die Linse `01a` (plane Seite zuerst, gleiche Radien), nur mit einer
**Hyperbel** statt der Kugel auf der zweiten Fläche: Konik **k = −n² = −2,3007**.
Parallele Strahlen laufen im Glas ungebrochen auf die zweite Fläche zu. Die Hyperbel (Exzentrizität e = n) fokussiert sie
exakt.

| Konik k der Rückfläche | −1,5 | −2,0 | **−2,3007** | −2,6 | −3,0 |
| --- | --- | --- | --- | --- | --- |
| RMS auf der Achse (µm) | 83,2 | 30,8 | **0,13** | 30,5 | 70,4 |

Vergleich: dieselbe Linse mit Kugel (`01a`) hat 248,5 µm. Das ist rund 1900-mal besser.

**Beobachten.** Ändere `conic` der Fläche `L1.S2` in der Datei und sieh dir *Spot* und *Ray Fans* an. Das Minimum liegt scharf bei k = −n².

**Fragen.**
- Warum verschwindet der Fehler nicht ganz (0,13 µm)? Probiere k = −(n/1,00027)² = −2,2994: Dann sinkt der Rest auf etwa 0,001 µm, weil die Luft n = 1,00027 hat.
- Außeraxial ist die Linse trotzdem schlecht (108,2 µm bei 3°). Sieh dir den Ray Fan bei 3° an: Ist er symmetrisch oder unsymmetrisch? Was folgt daraus für den Fehler?

## 04 Kugelspiegel und Parabolspiegel

**Dateien.** Hohlspiegel mit R = −200 mm (f ≈ 100 mm), EPD 40 mm (f/2,5), Felder 0°, 0,5°, 1°. `04a` ist eine Kugel, `04b` eine Parabel (k = −1).

| RMS (µm) | 0° | 0,5° | 1° |
| --- | --- | --- | --- |
| `04a` Kugel | 56,3 | 57,5 | 61,0 |
| `04b` Parabel | **0,0** | 12,2 | 24,4 |

(Zum Vergleich: Airy-Radius bei f/2,5 ≈ 1,8 µm.)

**Beobachten.** *Layout* (Strahlen bündeln sich bei der Parabel in einem Punkt), *Spot* und *Ray Fans* bei 1°.

**Fragen.**
- Warum ist die Parabel auf der Achse perfekt, aber nicht außeraxial? Der Fehler wächst linear mit dem Feldwinkel (12,17 → 24,40 µm bei Verdopplung des Winkels). Das ist **Koma**.
- Wie ändert sich der Kugelspiegel-Fehler, wenn man die Öffnung halbiert? Gemessen (`04a`, auf der Achse, EPD und Blende geändert): EPD 40 mm: 56,3 µm, EPD 20 mm: 6,94 µm, EPD 10 mm: 0,86 µm. Das ist Faktor 8,1 bei jeder Halbierung, also ∝ Öffnung³.

## 05 Cassegrain und Ritchey-Chrétien

**Dateien.** Zwei Spiegelteleskope mit gleicher Geometrie: Primärspiegel R = −800 mm, Spiegelabstand 300 mm, Brennweite
≈ 1440 mm (f/14,4), Felder 0°, 0,25°, 0,5°, 1°. Der Fangspiegel hat R = −276,92 mm. Unterschied sind nur die Kegelschnitte:

| | `05a` Cassegrain | `05b` Ritchey-Chrétien |
| --- | --- | --- |
| Konik Primärspiegel | −1,0000 (Parabel) | −1,0514 |
| Konik Fangspiegel | −3,1302 | −3,6764 |
| Seidel S_I (sphärische Aberration) | 5,80·10⁻⁸ | 1,05·10⁻⁶ |
| Seidel S_II (**Koma**) | **−5,26·10⁻⁴** | **2,79·10⁻⁸** |
| Seidel S_III (Astigmatismus) | 1,67·10⁻³ | 1,90·10⁻³ |
| RMS (µm) bei 0° / 0,25° / 0,5° / 1° | 0,0 / 5,4 / 19,5 / 76,1 | 0,0 / 5,0 / 20,0 / 80,3 |

Die Konik von `05b` wurde so gelöst, dass die Summen S_I und S_II zu null werden (Seidel-Ansicht).

**Beobachten.** *Seidel* (Tabelle und Balken je Fläche), *Layout* (gefalteter Strahlengang).

**Fragen.**
- Warum ist die Fleckgröße in beiden Teleskopen fast gleich, obwohl `05b` keine Koma hat? Auf einer **ebenen** Bildebene überwiegen Astigmatismus und Bildfeldwölbung. Der Astigmatismus ist bei `05b` sogar etwas größer (1,90·10⁻³ gegen 1,67·10⁻³).
- Wozu dient dann die Koma-Freiheit? (Man kann einen Korrektor oder eine gewölbte Bildfläche ergänzen.)

*Einschränkung des Modells:* Das sequenzielle Modell kennt weder die Abschattung durch den Fangspiegel noch das Loch im Primärspiegel.

## 06 Teleobjektiv: kurze Baulänge

**Dateien.** `06a` besteht aus einer Sammellinse (f ≈ 80 mm) und einer Zerstreuungslinse (f ≈ −50 mm) im Abstand
48,67 mm. `06b` ist zum Vergleich eine einfache Linse mit gleicher Brennweite. EPD 20 mm (f/10), Felder 0°, 1°, 2°.

| | `06a` Teleobjektiv | `06b` einfache Linse |
| --- | --- | --- |
| Brennweite (mm) | 199,96 | 200,77 |
| Baulänge (erste Fläche bis Bild, mm) | **127,9** | 203,5 |
| Verhältnis Baulänge / Brennweite | **0,64** | 1,01 |
| RMS auf der Achse (µm) | 0,9 | 22,1 |
| RMS bei 1° / 2° (µm) | 77,6 / 156,5 | 25,9 / 36,4 |

**Beobachten.** *Layout* (der Fokus liegt nur 127,9 mm hinter der ersten Fläche, obwohl die Brennweite 200 mm beträgt) und die *Paraxialen Kenndaten* im Kopf der Seite.

**Fragen.**
- Wie groß ist der Telefoto-Effekt? (Das Baulängen-Verhältnis liegt bei 0,64 statt 1.)
- Auf der Achse ist `06a` besser als die einfache Linse, außeraxial deutlich schlechter. Die Biegungen beider Linsen sind so gewählt, dass sich die Aberrationen auf der Achse aufheben (Radien: 52,63, −187,79, −18,87, −73,72 mm). Das gilt außeraxial nicht mehr. Wie viele Linsen braucht ein echtes Objektiv?

## 07 Blendenlage und Verzeichnung

**Dateien.** Dieselbe Linse (f ≈ 100 mm, Bestform): `07a` mit der Blende 20 mm **vor** der Linse, `07b` mit der Blende 20 mm **dahinter**.
Felder 0°, 5°, 10°, 15°. Alle Strahlen erreichen in beiden Fällen die Bildebene.

| Feldwinkel | 5° | 10° | 15° |
| --- | --- | --- | --- |
| `07a` Verzeichnung | −0,29 % | −1,17 % | −2,74 % |
| `07b` Verzeichnung | +0,43 % | +1,82 % | +4,46 % |

Negatives Vorzeichen: **Tonnenverzeichnung**; positives: **Kissenverzeichnung**.

**Beobachten.** Ansicht *Verzeichnung & Feldkrümmung* (Tabelle unter dem Diagramm), dazu *Layout*: Hauptstrahlen (andere Felder) laufen durch die Mitte der Blende.

**Fragen.**
- Die Linse ist in beiden Dateien identisch. Warum ändert sich trotzdem das Vorzeichen der Verzeichnung? (Der Hauptstrahl läuft durch die Blendenmitte und trifft die Linse an anderer Stelle.)
- Wo liegt die Blende, wenn die Verzeichnung null ist? Gemessen bei 15° (Abstand der Blende vom ersten Linsenscheitel; `pose.position` der Blende ändern, die Linse ist 5 mm dick):

| Blendenlage | Verzeichnung bei 15° |
| --- | --- |
| 20 mm vor dem 1. Scheitel | −2,74 % |
| 10 mm vor dem 1. Scheitel | −1,41 % |
| 5 mm vor dem 1. Scheitel | −0,76 % |
| 2 mm vor dem 1. Scheitel | −0,37 % |
| 7 mm nach dem 1. Scheitel | +0,62 % |
| 10 mm nach dem 1. Scheitel | +1,09 % |
| 20 mm nach dem 1. Scheitel | +3,08 % |

  Die Verzeichnung wechselt das Vorzeichen, wo die Blende durch die Linse wandert. Dort ist sie klein.

Hinweis: Die Blende hat in beiden Dateien den Radius 6 mm. Weil der Eintrittspupillen-Durchmesser (12 mm) festgelegt ist, schneidet sie keine Strahlen ab. Die Blendenlage wirkt über den Hauptstrahl.

## 08 Brewster-Winkel und Polarisation

**Datei.** `08_brewster_platte` ist eine 10 mm dicke Glasplatte (n = 1,5168), um die x-Achse um **56,60°** gekippt
(Brewster-Winkel: arctan(n/n_Luft) = 56,60°). Es gibt keine Blende, der Explorer nimmt automatisch das **freie
kollimierte Bündel**. Die Einfallsebene ist die yz-Ebene: Licht mit E-Feld in x ist **s-polarisiert**, Licht mit E-Feld in y ist **p-polarisiert**.

Transmission der Platte (beide Flächen, ohne Mehrfachreflexion), gemessen für verschiedene Kippwinkel (`rotation_deg` in der Datei ändern):

| Kippwinkel | 0° | 30° | 45° | 56,6° | 70° |
| --- | --- | --- | --- | --- | --- |
| s-Licht | 0,9176 | 0,8825 | 0,8174 | **0,7138** | 0,4811 |
| p-Licht | 0,9176 | 0,9472 | 0,9817 | **1,0000** | 0,9183 |

**So geht's im Explorer.** Ansicht *Polarisation*, Quelle *Linear*, Winkel 0° (s) oder 90° (p), Strahlquelle *Kollimiertes
Bündel*. Die Karten zeigen die Transmission über dem Bündel. Ansicht *Layout* zeigt den Strahlversatz in der Platte.

**Fragen.**
- Beim senkrechten Einfall sind s und p gleich (0,9176). Warum? (Je Fläche gehen ca. 4 % verloren, 0,96² ≈ 0,92.)
- Bei 56,6° wird p-Licht verlustfrei durchgelassen, s-Licht nur zu 0,71. Die analytische Rechnung (Fresnel, zwei Flächen) gibt (1 − R_s)² = 0,7138, also stimmt der Tracer.
- Wozu nutzt man das? (Brewster-Fenster in Lasern, einfacher Polarisator aus mehreren gekippten Platten.)

---

## Weitere Beispiele im Repo

`tests/reference/m3/` enthält Systeme für Polarisation und Beschichtungen: Malus-Gesetz (`malus`), λ/4- und λ/2-Platten
(`quarter_wave`, `half_wave`), Fresnel-Rhombus (`fresnel_rhomb`), Totalreflexion im Prisma (`tir_prism`) und ein Singlet mit
MgF₂-Entspiegelung (`ar_singlet`, braucht den Coating-Katalog `DEMO`).

## Neue Beispiele erzeugen

Die Radien und Kegelschnitte wurden numerisch mit der Bibliothek gelöst (Formel der dicken Linse
Φ = (n−1)(c₁−c₂) + (n−1)² t c₁ c₂ / n für c₂ bei gegebenem c₁, Newton-Verfahren für den Achromaten, Seidel-Summen für den
Ritchey-Chrétien); das dazu benutzte Skript liegt nicht im Repo. Die Struktur der Dateien beschreiben
`schema/raytatouille.schema.json` und der Abschnitt „Dateiformat“ in `docs/architecture.md`.

Die Dateien liegen bewusst nicht unter `tests/reference/`: Dort liest jeder Test jede Datei (Rundreise in C++, bitgleicher
Python-Vergleich, Schema-Test). Die Messwerte oben sind von Hand geprüft, aber nicht als Tests hinterlegt.
