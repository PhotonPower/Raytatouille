# Komplexe Modelle

Systeme mit vielen Flächen, Kittgliedern und Katalogglas, die eher einer realen Optik entsprechen als die
Lehrbeispiele. Die Zahlen unten wurden mit raytatouille 0.4.0 gemessen.

| Datei | System | Flächen | Gläser |
| --- | --- | --- | --- |
| `mikroskopobjektiv_20x.rtt.json` | Mikroskopobjektiv 20x/0,45, unendlich-korrigiert | 20 + Blende | `NIKON-HIKARI` |

## Mikroskopobjektiv 20x/0,45

**Herkunft.** Das Objektiv ist nach einem Schnittbild nachgebaut, einer Illustration eines
Mikroskopobjektivs mit acht Gruppen. Für das Bild wurde ein Maßstab von 0,1 mm/px angenommen. Daraus
stammen die Lage der Gruppen, die Mittendicken, die Luftabstände und die Wölbungsrichtungen. Die
Radien sind auf Abbildungsqualität optimiert (gedämpfte kleinste Quadrate über den RMS-Spot bei F, d, C und
drei Feldern, Nebenbedingungen f = 10 mm, Randdicken ≥ 0,25 mm). **Es ist kein Herstellerdatensatz.**

**Erzeugen.** `erzeuge_mikroskopobjektiv_20x.py` enthält die Prescription vorwärts (Präparat links, wie im
Schnittbild) und schreibt daraus die Systemdatei. Radien, Abstände, NA oder Feld lassen sich dort ändern.

```bash
python examples/komplexe_modelle/erzeuge_mikroskopobjektiv_20x.py
```

### Warum rückwärts gerechnet

Die Datei rechnet von der Tubusseite zum Präparat: ein paralleles Bündel (EPD 9 mm) fällt auf die
Hinterlinse G8 und wird auf das Präparat fokussiert (Detektor `OBJ`). So rechnet man unendlich-korrigierte
Objektive üblicherweise, und hier ist es auch nötig.

Die Blende sitzt wie bei realen Objektiven in der hinteren Brennebene, im Luftspalt zwischen G7 und G8. In
Lichtrichtung Präparat → Tubus liegt sie damit hinter dem hinteren Brennpunkt der Gruppe G1–G7, und die
Eintrittspupille ist virtuell und liegt etwa 90 mm hinter dem Präparat. Mit Version 0.4.0 starten die
Strahlen in diesem Fall rückwärts (#93). Ein Vorwärtsmodell mit Tubuslinse lässt sich deshalb erst nach
der Behebung von #93 analysieren.

### Aufbau

Die Gruppen und Flächen sind wie im Schnittbild von links (Präparat) nach rechts (Tubus) nummeriert. In der
Datei erscheinen sie in Strahlreihenfolge, also von G8 bis G1. Radien gelten vorwärts, in der Datei haben
sie das umgekehrte Vorzeichen.

| Gruppe | Art | Radien vorwärts (mm) | Gläser (NIKON-HIKARI) |
| --- | --- | --- | --- |
| G1 | Frontlinse, Meniskus | −2,885 / −5,455 | E-LAKH1 |
| G2 | aplanatischer Meniskus | −72,42 / −6,735 | NICF-V |
| G3 | Kittglied, drei Linsen | 13,173 / 103,777 / 27,842 / −13,665 | NICF-V, J-SFH1, NICF-V |
| G4 | Plankonvex | plan / −15,97 | J-FK5 |
| G5 | Bikonkavlinse | −42,889 / 11,805 | E-KZFH1 |
| G6 | Kittglied, kugelförmig | 12,392 / −6,596 / −13,452 | J-FK5, J-SFH1 |
| G7 | Kittglied, Meniskus | −14,225 / −7,563 / −16,94 | J-SFH1, E-LAKH1 |
| G8 | Hinterlinse, Meniskus | −63,261 / −41,623 | E-LAKH1 |

NICF-V ist ein fluoritähnliches Glas (n_d 1,434, ν_d 95). Alle Gläser stehen in
`tests/catalogs/nikon/nikon-hikari.agf`.

### Kenndaten

| Größe | Wert |
| --- | --- |
| Brennweite | 9,997 mm (20x bei 200-mm-Tubuslinse) |
| NA im Präparat | 0,45 (EPD 9,0 mm) |
| Arbeitsabstand (BFL) | 0,509 mm |
| Objektfeld | ±0,55 mm (Sehfeldzahl 22), Felder 0 / 0,385 / 0,55 mm |
| Baulänge G1.S1 bis G8.S2 | 33,7 mm |

**Messweise.** RMS um den Hauptstrahl, hexapolar 14 Ringe, polychromatisch F, d, C (1893 Strahlen),
am paraxialen Fokus. OPD als Karte mit 31 × 31 Punkten.

| Feld | RMS-Spot (µm) | OPD RMS (Wellen) |
| --- | --- | --- |
| 0 | 0,308 | 0,014 |
| 0,385 mm | 0,351 | 0,065 |
| 0,55 mm | 0,583 | 0,146 |

Zum Vergleich: Der Airy-Radius bei NA 0,45 und 587,6 nm beträgt 0,80 µm. Bis 70 % des Feldes ist das
Objektiv beugungsbegrenzt (Maréchal-Grenze 0,07 Wellen RMS), am Feldrand nicht mehr ganz.

| Größe | Wert |
| --- | --- |
| Farblängsfehler F−C, paraxial / Randzone | −1,17 µm / −0,56 µm |
| Verzeichnung am Feldrand | −0,665 % |
| Astigmatismus am Feldrand | 5,90 µm |

### Laden

```python
import raytatouille as rt
from raytatouille import analysis as an

lib = rt.MaterialLibrary()
lib.add_catalog("tests/catalogs/nikon/nikon-hikari.agf")
system = rt.compile(rt.load("examples/komplexe_modelle/mikroskopobjektiv_20x.rtt.json"), lib)
spot = an.spot(system, path="main", field=2, rays="hexapolar:14")
print(f"RMS {1000 * spot.stats.rms_chief:.3f} µm")
```

Im Explorer: **System → Datei hochladen**. Den Katalog NIKON-HIKARI wählt der Explorer aus
`tests/catalogs`, sonst `nikon-hikari.agf` in der Seitenleiste hochladen.
