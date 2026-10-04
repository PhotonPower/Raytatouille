# Geprüfte Quellen

Regel 6 aus [AGENTS.md](../AGENTS.md) verlangt für jede Formel eine Literaturstelle mit
Gleichungsnummer. Diese Liste enthält Quellen, deren Gleichung jemand **selbst gelesen und mit
dem Code verglichen** hat. Agenten verwenden zuerst diese Einträge.

## Neue Quelle eintragen

1. Frei zugängliche Quelle mit Gleichungsnummer suchen: Open-Access-Artikel (z. B. Optics
   Express, Optica, arXiv), Hochschulskripte, Normen-Entwürfe mit freiem Zugang.
2. Die Gleichung im Volltext lesen und Symbol für Symbol mit dem Code vergleichen:
   Vorzeichen, Normalenrichtung, Einheiten, Indexzählung.
3. Eintrag in diese Liste im selben PR wie der Code, mit Link, Gleichungsnummer,
   Konvention und bekannten Fallstricken.
4. Erst wenn keine prüfbare Quelle auffindbar ist: Frage-Issue (Label `frage`).

Nicht erlaubt: Gleichungsnummern aus dem Gedächtnis oder aus Zitaten anderer Arbeiten
übernehmen, ohne den Volltext gesehen zu haben.

## Geprüft

| Thema | Quelle | Stelle | Konvention und Fallstricke | Geprüft |
| --- | --- | --- | --- | --- |
| Sag einer Konik | G. W. Forbes, „Manufacturability estimates for optical aspheres“, *Opt. Express* 19(10), 9923–9942 (2011), [Volltext](https://opg.optica.org/oe/fulltext.cfm?uri=oe-19-10-9923&id=213662) | Gl. (2.1) und die direkt folgende Definition von φ(ρ) | z = c ρ² / (1 + φ), φ = √(1 − (1 + κ) c² ρ²); c Scheitelkrümmung, κ konische Konstante. Diese Form ist für c → 0 und κ = −1 stabil, die Form (1 − φ)/((1 + κ) c) nicht. Die Nummer der φ-Definition ist nicht geprüft, nicht zitieren. | #13 |
| Reflexion in Vektorform | B. de Greve, „Reflections and Refractions in Ray Tracing“ (2006), [PDF](https://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf) | Gl. (13) | d′ = d − 2 (d·N) N; unabhängig von der Richtung von N. | #19 |
| Brechung in Vektorform | wie oben | Gl. (22) mit (23) und (28); Zusammenfassung Abschnitt 6 | d′ = (n/n′) d + ((n/n′) cos θ − cos θ′) N mit **N zum Einfallsmedium** (N·d < 0) und cos θ = −d·N. Unsere Geometrie liefert N in +z am Scheitel: vor der Formel N so drehen, dass N·d < 0. **Fehler in der Quelle:** nach Gl. (20a) steht cos θ = i·n, richtig ist −i·n (wie in Abschnitt 6). | #19 |
| Totalreflexion, Grenzwinkel | wie oben | Gl. (24) und (17); Grenzwinkel Gl. (25) | TIR, wenn sin² θ′ = (n/n′)² (1 − cos² θ) > 1. | #19 |
| Fresnel-Reflexionsgrade (Leistung, unpolarisiert gemittelt) | wie oben | Gl. (27a), (27b), (29a) | Nur Leistungen, keine Amplituden und keine Phasen. Für die Polarisation (M3) **nicht ausreichend**. | – |

## Im Code zitiert, aber nicht online geprüft

Buchquellen und kostenpflichtige Normen. Sie dürfen als zusätzliche Angabe stehen bleiben;
wer sie am Original prüft, verschiebt sie nach oben.

| Thema | Quelle | Verwendet in |
| --- | --- | --- |
| Newton-Verfahren, Abbruchkriterium | Press et al., *Numerical Recipes*, 3. Aufl., Abschn. 9.4 | `rtt-geom/intersect.hpp` |
| Stabile Lösung der quadratischen Gleichung | Press et al., *Numerical Recipes*, 3. Aufl., Abschn. 5.6, Gl. 5.6.4–5.6.5 | `rtt-geom/intersect.hpp` |
| Gerade Asphäre | ISO 10110-12 | `rtt-geom/asphere.hpp` |
| Konik, allgemein | W. T. Welford, *Aberrations of Optical Systems*, Kap. 2 | `rtt-geom` |
| Brechung, allgemein | M. Born, E. Wolf, *Principles of Optics*, 7. Aufl., §3.2.2 | `rtt-trace` (#6) |
| Linsenmacherformel | E. Hecht, *Optics*, Kap. 6 | `rtt-paraxial` (#7) |
| Polarisations-Raytracing | Chipman, Lam, Young, *Polarized Light and Optical Systems* (2018) | `rtt-trace/ray_batch.hpp` |

## Noch ohne geprüfte Quelle (bald gebraucht)

Wer an diesen Themen arbeitet, sucht zuerst eine Quelle nach dem Verfahren oben.

| Thema | Meilenstein |
| --- | --- |
| Paraxialer y-nu-Trace, Hauptebenen, Pupillen | M1 (#7) |
| Dispersionsformeln (Sellmeier, Schott, Herzberger, Conrady), AGF-Format | M2 |
| Brechzahl der Luft (Ciddor), dn/dT nach Schott | M2 |
| Seidel-Koeffizienten | M2 |
| Fresnel-Amplituden mit Phase, auch für komplexen Index | M3 |
| Transfermatrix-Methode für Schichtsysteme | M3 |
| 3D-Polarisations-Raytracing-Matrix | M3 |
