# Raytatouille – Architektur

> Maßgebliche Fassung. Änderungen nur per PR; Architekturentscheidungen zusätzlich als ADR unter `docs/adr/`.

## Ziel, Scope und Leitprinzipien

Raytatouille ist eine headless Raytracing-Engine für Optikdesign in C++20 mit Python-API, lizenziert unter MIT. Die GUI wird später in Python auf dieser API gebaut. Dieses Dokument ist die verbindliche Vorgabe für Programmieragenten: Es legt Module, Schnittstellen, Konventionen und Abnahmekriterien fest.

**Scope v1:** bildgebende und polarisationsoptische Systeme aus Linsen (sphärisch, asphärisch, Freiform), Spiegeln, Beschichtungen, Strahlteilern, Polarisatoren und Retardern. Tracing paraxial, sequenziell und multi-sequenziell (mehrere benannte Strahlpfade in einem Modell). Analysen, lokale Optimierung und Toleranzierung.

**Später (v2+):** als erste Priorität nicht-sequenzielles Tracing mit Streuung, danach Python-GUI, globale Optimierung, Beamlet-Wellenoptik, GPU-Backend.

**Nicht geplant:** rigorose EM-Solver (RCWA, FDTD), Beleuchtungsdesign, Struktur-FEA. Dafür sind Import-Schnittstellen (Jones-Tabellen, Coating-Dateien) vorgesehen.

**Leitprinzipien**

1. **Polarisation ist nativ.** Jeder reale Strahl trägt ab dem ersten Release eine 3×3-Polarisations-Raytracing-Matrix. Es gibt keinen separaten „Polarisationsmodus“.
2. **Flächen sind Stacks.** Eine Fläche = Form + Apertur + optionale Phase + Interaktion (Fresnel, Coating, Polarisator, Retarder, Strahlteiler). Neue Eigenschaften sind neue Layer, keine neuen Flächentypen.
3. **Pfade sind Daten.** Ein Strahlpfad ist eine explizite Liste von Flächen-Ereignissen (brechen, reflektieren, o-/e-Strahl, Beugungsordnung). Strahlteiler und Doppelbrechung erzeugen mehrere Pfade, keine Sonderfälle im Tracer.
4. **Hierarchie statt Flächenliste.** System → Baugruppe → Element → Fläche, jede Ebene mit eigenem Koordinatensystem und Pivot. Toleranzen wirken auf diese Hierarchie.
5. **Headless und deterministisch.** Keine UI- und keine Zufallsabhängigkeit ohne expliziten Seed. Gleiche Eingabe → bitgleiche Ausgabe auf derselben Plattform.
6. **Physik vor Features.** Kein Modul gilt als fertig ohne Referenztest gegen eine analytische Lösung oder ein publiziertes Ergebnis.
7. **Agentenfreundlich.** Kleine Bibliotheken mit schmalen Schnittstellen, damit Agenten parallel arbeiten können.

## Technologie-Entscheidungen

Details und Begründungen stehen in `docs/adr/`. Kurzfassung:

| ADR | Entscheidung |
| --- | --- |
| 0001 | Kern in C++20, CMake ≥ 3.25, Abhängigkeiten über vcpkg (Manifest `vcpkg.json`) |
| 0002 | Python-Bindings mit nanobind und scikit-build-core, NumPy-Views ohne Kopie |
| 0003 | Eigen (≥ 3.4) für 3-Vektoren und 3×3-Matrizen, `std::complex<double>`, durchgängig `double` |
| 0004 | Parallelität mit oneTBB, `static_partitioner`, deterministisch |
| 0005 | Strahlen als Structure-of-Arrays |
| 0006 | Physik-Bibliotheken als Templates über den Skalartyp (`rtt::math::Real`; nur `double` und später Dual-Zahlen, ADR 0015) |
| 0007 | Ableitungen in v1 per zentraler finiter Differenz |
| 0008 | Dateiformat JSON (nlohmann/json); strenger C++-Parser ist Referenz, JSON-Schema für externe Werkzeuge |
| 0009 | Exceptions nur an API-Grenzen, Strahlprobleme als Status-Flags |
| 0010 | Lizenz MIT, nur kompatible Abhängigkeiten |
| 0011 | Linux und Windows; getestet mit GCC 13, Clang 18, MSVC 2022; kein `-ffast-math` |
| 0012 | GUI später in Python (PySide6, pyvista/VTK, matplotlib) nur über die öffentliche API |
| 0013 | Keine GPU in v1 |
| 0014 | Laufzeit-Interface `Material` in `double`; Formelkerne als Templates |
| 0015 | `rtt::math::Real` umfasst nur `double` (später Dual-Zahlen) |
| 0016 | `CompiledSystem` in eigener Bibliothek `rtt-compile`; feste Reihenfolge in der Schicht Tracing |
| 0017 | Kittglieder: Material je Segment (`Lens`, `Plate`), Schema 0.2 |
| 0018 | Python-Build- und Testwerkzeuge (nanobind, scikit-build-core, numpy, pytest, mypy) über `pyproject.toml` statt vcpkg |
| 0019 | Coating-Kataloge als eigene JSON-Dateien (`CATALOG:NAME`), Auflösung in `compile()`, Substrat = Inneres des Elements |
| 0020 | Gemeinsamer strikter JSON-Parser `rtt-json` (header-only, Schicht Basis): doppelte Schlüssel und Zahlen-Overflow sind Fehler mit Pointer |
| 0021 | Interaktionen im Tracer (Fresnel, Coatings, ideale Elemente je Ereignisart, `EventImpossible` für den Rest); `prt` leistungsnormiert, `weight` = s·½‖P_T‖²_F (angenommen) |
| 0022 | Stabile Diagnosecodes in Punktnotation mit Registry `rtt-diagnostics` (header-only, Schicht Basis), Warnungen von `compile()` als Daten, Fehlerorte in `ParaxialError`/`AnalysisError`, eine Klasse `NoStopError` (angenommen) |
| 0023 | Ergebnisformat `raytatouille-result` (JSON mit dtype/shape für Arrays, NaN als Zeichenkette, Roundtrip = Daten) und Strahlverluste als Daten (angenommen) |
| 0024 | Modell lesen über unveränderliche typisierte Kopien, ändern über JSON Patch (RFC 6902) auf der Bearbeitungsform, Adressierung mit JSON-Pointern und verankerten Verweisen, Undo und Redo über inverse Patches (angenommen) |
| 0025 | Gitter und diffraktive Flächen als Phasenschicht: Ordnung am Ereignis, lokale Gittergleichung, OPL-Beitrag der Ordnung, Effizienz nur im Gewicht, Status `Evanescent` (angenommen) |
| 0026 | Einachsiger Kristall: Material mit n_O und n_E, optische Achse am Element, Strahlzustand mit Wellennormale und Modenindex, Projektionsmodell für P (angenommen) |
| 0027 | Ghost-Generator: Zweifachreflexions-Ghosts eines Pfads als erzeugte explizite Pfade, Medien durch compile (angenommen) |
| 0028 | Relative Platzierung: `Pose.reference` (absolut, relativ zur vorangehenden Fläche oder zum vorangehenden Geschwister) und `Pose.order` (`translate_first` wie bisher, `rotate_first` als Koordinatensprung) (angenommen) |
| 0029 | Parametertabelle mit Ausdrücken nur über frühere Zeilen, Konfigurationen als Spalten (nur Parameter), Modell-Param verweist mit `{"param": …}` auf eine Zeile, `pickup` entfällt mit Warnung (angenommen) |
| 0030 | Optimierung: Merit-Funktion als kleinste Quadrate im Abschnitt `optimization`, Operanden und Gauß-Quadratur-Generatoren, Levenberg-Marquardt (Madsen/Nielsen/Tingleff, Dämpfung nach Nielsen, feste Skalierung nach MINPACK-1), zentrale Differenz parallel, Grenzen per MINUIT-Transformation, Abbruch mit bestem Stand, Ergebnis als ein Patch (angenommen) |

**Konventionen (verbindlich für alle Bibliotheken)**

- **Einheiten:** Längen in mm, Wellenlängen in µm, Winkel intern in rad (API-Parameter und Dateifelder in Grad tragen das Suffix `_deg`), Temperatur in °C, Druck in atm. Hilfsfunktionen in `rtt/math/units.hpp`.
- **Koordinaten:** rechtshändig, optische Achse +z, y ist die meridionale Richtung.
- **Radius:** positiv, wenn der Krümmungsmittelpunkt auf der +z-Seite des Flächenscheitels liegt (übliche Konvention in der Optik). Im Modell steht der Radius; der Tracer rechnet mit c = 1/R. Ebene Flächen sind `Plane`, nie Radius 0 oder unendlich.
- **Transformationen:** Translation, dann Rotation um den Pivot; Rotationen intrinsisch X → Y → Z, also R = Rx·Ry·Rz. p_parent = position + pivot + R·(p_child − pivot). Implementiert in `rtt::math::Isometry3::from_pose`. Mit `Pose::order` = `rotate_first` (ADR 0028) gilt stattdessen p_ref = Q + R·(position + p − Q) mit dem Pivot Q im KS des Bezugs (Koordinatensprung); `Pose::reference` legt den Bezug fest (Elternknoten, vorangehende Fläche in Baumreihenfolge oder vorangehendes Geschwister). Relative Bezüge wertet compile ab #163 aus.
- **Felder und Phase:** Ebene Wellen ~ exp(i(k·r − ωt)). Komplexer Index ñ = n + iκ mit κ ≥ 0 für Absorption. Positive Retardance heißt: die langsame Achse ist phasenverzögert.
- **Polarisation und Fresnel** (festgelegt für #56 und #58; Quelle S. J. Byrnes, arXiv:1603.02720v5, siehe `docs/quellen.md`):
  - **Tangentialinvariante:** ξ = n sin θ ist reell und an jeder Grenzfläche erhalten; bei absorbierendem Einfallsmedium ist ξ = Re(ñ_i) sin θ_i, wie bei der Brechung im Tracer. ξ ist die kanonische Eingabe von `rtt-polar` und `rtt-coating`; Varianten mit cos θ_i (`fresnel_at_angle`) gelten nur für nicht absorbierende Einfallsmedien.
  - **Normalkomponente:** q_j = ñ_j cos θ_j = √(ñ_j² − ξ²) mit Im q ≥ 0 und, falls Im q = 0, Re q ≥ 0. Die Welle klingt ab bzw. trägt Leistung vorwärts (Byrnes, Anhang D). Das gilt auch im Einfallsmedium. Totalreflexion: q_t = i√(ξ² − n_t²), |r_s| = |r_p| = 1.
  - **Vorzeichen (Convention A, Byrnes Gl. (6) und Anhang A; Anhang A, Fußnote 6, nennt Lehrbücher mit dieser Konvention):** r_s = (ñ_i cos θ_i − ñ_t cos θ_t)/(ñ_i cos θ_i + ñ_t cos θ_t), r_p = (ñ_t cos θ_i − ñ_i cos θ_t)/(ñ_t cos θ_i + ñ_i cos θ_t), t_s = 2ñ_i cos θ_i/(ñ_i cos θ_i + ñ_t cos θ_t), t_p = 2ñ_i cos θ_i/(ñ_t cos θ_i + ñ_i cos θ_t). Bei senkrechtem Einfall ist r_p = −r_s.
  - **Basis:** Mit der Flächennormale N ins Medium t sind die p-Einheitsvektoren in Byrnes Gl. (5) für die einfallende, die reflektierte und die transmittierte Welle genau p = k × s mit demselben s = k_in × N/|k_in × N|. Kehrt man N um, wechseln s und alle p gemeinsam das Vorzeichen, r und t bleiben gleich. r und t sind damit direkt a_s, a_p in J = diag(a_s, a_p, 1) der PRT-Matrix (Abschnitt Polarisation).
  - **Begründung über den idealen Leiter:** Er reflektiert E_r = −E_tan + E_norm, also E_r = (−I + 2NNᵀ)E_i. Daraus folgt (−I + 2NNᵀ)s = −s und (−I + 2NNᵀ)p_in = p_out für jeden Einfallswinkel, also a_s = −1 und a_p = +1. Auf transversalen Feldern ist P = −R mit der geometrischen Spiegelung R = I − 2NNᵀ: Beide Eigenpolarisationen bekommen dieselbe Phase π, es gibt weder Diattenuation noch Retardance. Fresnel mit ñ = 1 + iκ, κ → ∞, ergibt in Convention A genau r_s → −1 und r_p → +1 (Test in #56). Convention B gäbe r_p → −1, also scheinbar die Retardance π.
  - **Energie:** R = |r|², T_s = |t_s|² Re(q_t)/Re(q_i), T_p = |t_p|² Re(ñ_t cos θ_t*)/Re(ñ_i cos θ_i*) (Byrnes Gl. (21)–(23), mit Konjugation). Bei verlustfreiem Einfallsmedium gilt R + T = 1, auch wenn das zweite Medium absorbiert. Bei absorbierendem Einfallsmedium kann R + T ≠ 1 sein; das ist kein Fehler (Byrnes Abschnitt 4.4.1, Anhang B).
  - **Phasendifferenz (#56):** bei Reflexion Δ = arg(−r_p/r_s) (Ellipsometrie, Byrnes Gl. (16); 0 bei senkrechtem Einfall), bei Transmission Δ = arg(t_p/t_s), jeweils in (−π, π]. Δ > 0 heißt p verzögert; bei exp(−iωt) bedeutet ein größeres arg eine Verzögerung. Die Retardance mit langsamer Achse aus P folgt in #57.
  - **Fehlerbehandlung:** Die Kerne sind `noexcept`, mit dokumentierten Vorbedingungen und `assert` im Debug-Build. An der API-Grenze gibt es geprüfte Varianten (`fresnel_checked`), die `std::invalid_argument` werfen.
  - **Händigkeit und Stokes** (festgelegt für #60; Quelle Lam, Tab. 2.1, S. 69, und Abb. 1.1, S. 48): In der Transversalbasis (e1, e2 = k × e1, k) ist E = (i·e1 + e2)/√2 **rechtszirkular**, „evolving clockwise in time as one looks into the beam“ (Lam). Mit e^{−iωt} gilt dann E(t) ∝ (sin ωt, cos ωt): Gegen die Ausbreitung gesehen dreht sich E im Uhrzeigersinn. (e1 + i·e2)/√2 ist linkszirkular. Stokes-Parameter (Projektdefinition, Lam gibt keine Gleichung; das Vorzeichen von S3 ist so gewählt, dass Lams rechtszirkular (i, 1, 0) S3 = +S0 ergibt): S0 = |E1|² + |E2|², S1 = |E1|² − |E2|², S2 = 2 Re(E1* E2), S3 = 2 Im(E1 E2*); **S3 > 0 heißt rechtszirkular**. Der lineare Polarisationswinkel ist ½·atan2(S2, S1), von e1 nach e2. In der Literatur gibt es die entgegengesetzte Händigkeitskonvention (Blick in Ausbreitungsrichtung) und das entgegengesetzte Vorzeichen von S3; Vergleichswerte aus anderen Quellen sind deshalb vorher umzurechnen. Probe nach Lam Abb. 1.1: Eine λ/4-Platte mit langsamer Achse x macht 45°-lineares Licht rechtszirkular (Test in #60).
- **Größen:** OPL in mm, OPD und Wellenfront in Wellen bei der Referenzwellenlänge, Strahlgewicht als Leistung (Quelle normiert auf 1).
- **OPD-Vorzeichen** (festgelegt in #29): W = OPL_ref − OPL_Strahl, beide bis zur Referenzsphäre gemessen; **W > 0 heißt voreilend**, die Welle krümmt sich stärker als die Referenzsphäre (Wyant & Creath, Abschnitt I, siehe `docs/quellen.md`). Ein zur Austrittspupille hin verschobener Detektor sieht W < 0; unterkorrigierte sphärische Aberration mit Referenz im paraxialen Fokus ergibt W > 0 am Rand.
- **Feldwinkel und Pupille** (festgelegt in #8): Ein Feldwinkel (θx, θy) in Grad gibt die Richtung des Hauptstrahls im Objektraum d ∝ (tan θx, tan θy, 1) an; θy > 0 heißt, der Strahl steigt in +y (der Objektpunkt liegt bei −y, bei virtueller EP hinter dem Objekt (z_EP < z_Objekt) bei +y). Eine Objekthöhe (x, y) in mm ist der Objektpunkt (x, y, −Objektabstand), positiv = +y. Eine paraxiale Bildhöhe wird über den paraxialen Hauptstrahl in Winkel bzw. Objekthöhe umgerechnet. Feldwerte werden bei der Referenzwellenlänge in Richtung bzw. Objektpunkt umgerechnet (paraxiale Daten der Referenzwellenlänge); ein Feldpunkt ist damit für alle Wellenlängen derselbe physikalische Punkt, während Pupille, Blendenziel und Startebene zur Wellenlänge des Strahls gehören (#31). Normierte Pupillenkoordinaten (px, py): der Einheitskreis ist der Rand der paraxialen Eintrittspupille, +y meridional; Ray Aiming zielt auf (px, py)·R_s auf der Blendenfläche (R_s paraxialer Blendenradius zur Eintrittspupille, vorzeichenbehaftet bei invertierter Abbildung). Ein Feldwinkel bei endlichem Objekt legt den Objektpunkt auf den Hauptstrahl durch die EP-Mitte. Bei Objekt im Unendlichen starten die Strahlen eines Feldes auf einer Ebene senkrecht zur Feldrichtung (ebene Welle, OPL 0), so platziert, dass das Bündel bis 3 EP-Radien um den Hauptstrahl mindestens 1 mm vor der EP und vor jeder erreichbaren Fläche des Pfades startet; bei endlichem Objekt im Objektpunkt. Die EP darf im System hinter Flächen liegen (paraxiales Bild der Blende). Bei endlichem Objekt ist der Strahl zum Pupillenpunkt Q die Gerade durch den Objektpunkt P und Q; er läuft in der Ausbreitungsrichtung des Objektraums, nach der Paraxial-Konvention also in +z (bei parabasalen, gekippten Systemen später entlang der Objektraumachse). Liegt die EP virtuell auf der anderen Seite des Objekts (z_EP < z_Objekt, etwa bei einer Blende hinter dem hinteren Brennpunkt der Gruppe davor), ist d ∝ P − Q statt Q − P (#93). Weil die Konjugation Geraden betrifft, bleibt die Zuordnung (px, py) → (px, py)·R_s auf der Blende gleich; da P und alle Q auf den Ebenen z_Objekt und z_EP liegen und die Pfade rotationssymmetrisch sind, gilt die Entscheidung für alle Strahlen aller Felder. **Objektseitig telezentrisch** (endliches Objekt, paraxiale EP im Unendlichen, also a = 0 in der Frontmatrix bis zur Blende; #96): Der Blendenpunkt hängt dann nur von der Objektraum-Steigung ab (y_Blende = b n₁ u). Pupillenkoordinaten sind deshalb Steigungen: Der Strahl von (px, py) startet in P mit d ∝ (px·u_m, py·u_m, 1), der paraxiale Hauptstrahl läuft parallel zur Achse (der real gezielte weicht um die Pupillenaberration ab, ∝ h³), und py > 0 heißt Steigung nach oben. u_m ist die paraxiale Randstrahlsteigung: NA/n₁ bei `object_na`, r_Blende/|s| bei `stop_size` (s = Blendenhöhe des Strahls aus dem axialen Objektpunkt mit Steigung 1, R_s = u_m·s). `entrance_pupil_diameter` und `image_fnumber` legen bei EP im Unendlichen kein Bündel fest und ergeben einen Fehler, ebenso ein Feldwinkel bei endlichem Objekt (alle Hauptstrahlen sind achsparallel); eine paraxiale Bildhöhe wird über den achsparallelen Hauptstrahl umgerechnet. Telezentrisch gilt nur bei a = 0 exakt, ohne eigene Schwelle. Der Übergang ist stetig: Liegt die Blende knapp neben der Brennebene, ist die EP endlich, aber sehr weit weg, und das Bündel geht stetig in das telezentrische über. Die Beschriftung der Pupillenpunkte ist dabei stetig von der Seite z_EP → +∞; von der Seite z_EP → −∞ (virtuelle EP hinter dem Objekt, #93) ist sie gespiegelt, (px, py) → (−px, −py). Das ist die bekannte Umkehr des Pupillenbilds beim Durchgang durch Unendlich. Verschiebt man die Blende durch die Brennebene (etwa beim Optimieren oder in einem Parameter-Sweep), kippen deshalb Ray Fans und OPD-Karten in (px, py), obwohl sich das Strahlbündel selbst stetig ändert. **Pupillenkoordinaten sind an der Referenzwellenlänge orientiert** (#96, präzisiert #93): Bei Dispersion kann die EP einer Wellenlänge auf der anderen Seite des Objekts liegen als die der Referenzwellenlänge (typisch bei objektseitiger Telezentrie: blau mit Blende hinter dem Fokus, rot davor). Dann wird σ = sign(z_EP,ref − z_Objekt)·sign(z_EP,λ − z_Objekt) (EP im Unendlichen zählt als „+“), und der Pupillenpunkt (px, py) ist der EP-Punkt σ·(px, py)·r_EP bzw. die Steigung σ·(px, py)·u_m mit dem Ziel σ·(px, py)·R_s. So meint ein Pupillenpunkt für jede Wellenlänge dieselbe Blendenseite, und Fächer und OPD-Karten bleiben über λ vergleichbar. Bei Objekt im Unendlichen gibt es keine EP-Seite, σ = +1; ohne Seitenwechsel über λ ist σ = +1 und alles bitgleich wie zuvor. Liegt die Blende so nah an der Brennebene, dass die Rundung über das Vorzeichen von a entscheidet, ist auch diese Beschriftung zufällig; das Bündel bleibt richtig. Damit die reale Zielung dort konvergiert, hat die Schrittweite ihrer Jacobi-Matrix eine Untergrenze in Pupilleneinheiten (Nachtrag #96 zu ADR 0007). Zur OPD bei Austrittspupille im Unendlichen (bild- bzw. beidseitig telezentrisch) siehe den OPD-Absatz unter Analyse (#102). Die erreichbare Zone einer Fläche wird iterativ bestimmt (Bündelradius am tiefsten gefundenen z); konvergiert das nur langsam (linear mit Faktor tan θ·|dsag/dr| nahe 1, z. B. L2.S1 des Cooke-Triplets bei 14°), gilt nach 50 Runden das Minimum über die Kreisscheibe mit dem Außenradius der Apertur bzw. dem Formgebiet der Fläche als sichere Schranke. Ein Fehler entsteht nur bei einer unbegrenzten Fläche ohne Apertur oder mit einer Apertur ab 1 km (Startebene sonst ohne nutzbare Genauigkeit), die sich gegen ein steiles Feld zurückkrümmt (#72).

## Systemübersicht

17 CMake-Bibliotheken bzw. -Programme in sieben Schichten. Jede darf nur Bibliotheken aus tieferen Schichten verwenden. Ausnahmen: die Schicht Tracing mit der festen Reihenfolge `rtt-compile` < `rtt-paraxial` < `rtt-trace`, dort darf eine Bibliothek zusätzlich die in dieser Reihenfolge vor ihr stehenden verwenden (ADR 0016); in der Schicht Schnittstellen dürfen `rtt-py` und `apps/rtt-cli` `rtt-io` verwenden (Dateien lesen und schreiben).

| Schicht | Bibliotheken | Status |
| --- | --- | --- |
| Schnittstellen | `rtt-py` (Python-API), `apps/rtt-cli` (Kommandozeile), `rtt-io` (Dateien) | `rtt-cli`, `rtt-io`: M0; `rtt-py`: M2 |
| Workflows | `rtt-optim`, `rtt-tolerance` | M5, M7 |
| Auswertung | `rtt-analysis` | M2, M6 |
| Tracing | `rtt-compile` < `rtt-paraxial` < `rtt-trace` | M1 |
| Modell | `rtt-model` | M0 |
| Physik | `rtt-geom`, `rtt-material`, `rtt-coating`, `rtt-polar` | M1–M4 |
| Basis | `rtt-math`, `rtt-json` (strikter JSON-Parser für `rtt-io` und `rtt-coating`, ADR 0020), `rtt-diagnostics` (Registry der Diagnosecodes, ADR 0022) | `rtt-math`: M0; `rtt-json`: M3; `rtt-diagnostics`: #86 |

Die vier Physik-Bibliotheken kennen weder Modell noch Tracer und können daher parallel gebaut werden. `rtt-model` beschreibt nur, `rtt-compile` macht daraus ein unveränderliches `CompiledSystem`, `rtt-paraxial` und `rtt-trace` rechnen, `rtt-analysis` wertet aus. CMake-Targets heißen `rtt_<name>` mit Alias `rtt::<name>`.

## Datenmodell (`rtt-model`, umgesetzt in M0)

Das Modell ist reine Datenstruktur ohne Tracing-Logik. Header: `libs/rtt-model/include/rtt/model/`.

**Hierarchie**

- **`System`:** Name, Umgebung (Temperatur, Druck, Umgebungsmedium), Objektraum (unendlich oder Abstand), Wellenlängen (genau eine Referenz), Systemapertur, Feldpunkte, Wurzel-Baugruppe, Pfade.
- **`Assembly`:** Name, `Pose`, Kinder (`Node` = Assembly oder Element).
- **`Element`:** physischer Körper mit `ElementKind`: `Lens` (≥ 2 Flächen, Material je Segment), `Mirror` (≥ 1 Fläche, optional ein Substratmaterial), `Plate` (≥ 2 plane Flächen, Material je Segment; Platten, Würfel, Prismen), `ThinElement` (1 Fläche ohne Dicke), `Stop` (1 Fläche mit Apertur, höchstens einer pro System), `Detector` (1 Fläche).
- **Segmente (ADR 0017):** Ein Element mit N Flächen hat N − 1 Segmente; Segment i liegt zwischen Fläche i und i + 1 (nullbasiert). `Lens` und `Plate` geben entweder ein Material für alle Segmente an (`material`, Kurzform) oder eine Liste mit genau N − 1 Einträgen (`segment_materials`, Kittglieder); `Element::segment_material(i)` liefert das Material von Segment i in beiden Formen.
- **`Surface`:** eindeutige `SurfaceId` (lesbarer String wie `"L1.S1"`), `Pose` im Element, Layer-Stack.

**Layer-Stack einer Fläche**

| Layer | Anzahl | Umgesetzt (Schema 0.3) | Geplant |
| --- | --- | --- | --- |
| `shape.base` | genau 1 | `Plane`, `Conic`, `EvenAsphere` | Q-con, Q-bfs, Biconic, Toroid, Axicon, XY-Polynom, ungerade Asphäre |
| `shape.terms` | 0..n, additiv | `ZernikeSag` (Noll) | Grid-Sag, Zernike Fringe |
| `aperture` | 0..1 | kreisförmig (auch Ring), rechteckig, elliptisch | Polygon, Obskurationen |
| `phases` | 0..n, additiv | `LinearGrating`, `RadialPhase`; dazu `diffraction_efficiency` je Ordnung (ADR 0025) | Zernike-Phase, Metaflächen-Tabellen |
| `interaction` | genau 1 | Fresnel (Default), ideal: Spiegel, AR, Strahlteiler, Polarisator, Retarder; Absorber; Coating-Referenz | Jones-/Mueller-Tabellen |
| `scatter` | 0..1 | – | Lambert, Gauß, ABg (M9) |

**Medien:** innerhalb eines Elements das Material des jeweiligen Segments (bei der Kurzform überall dasselbe), außerhalb das Umgebungsmedium. Material ist ein Katalogverweis `"KATALOG:NAME"` oder ein einachsiger Kristall `CrystalMaterial { ordinary, extraordinary }` mit `optic_axis` am Element (ADR 0026; nur `Lens` und `Plate`, ein Kristall für alle Segmente). Die Kittfläche i (0 < i < N − 1) ist genau ein Übergang von Segment i − 1 zu Segment i.

**Pfade:** `Path` = Name + Liste von `Event { SurfaceId surface; EventKind kind; int order; }` oder `automatic` (alle Flächen in Baumreihenfolge). `EventKind`: `Refract`, `Reflect`, `Transmit`, `Ordinary`, `Extraordinary`. Jedes Ereignis an einer Fläche mit Phasenschicht darf eine Beugungsordnung `order` tragen, Ordnung 0 ist das Ereignis ohne Beugung (ADR 0025; bis Schema 0.2 gab es dafür die Art `Diffract`). Dieselbe Fläche darf mehrfach vorkommen (Doppeldurchgang, Interferometer). Der Ghost-Generator (M4, #123, ADR 0027) erzeugt aus einem Pfad alle Zweifach-Reflexions-Ghosts als explizite Pfade (`rtt/compile/ghosts.hpp`: `ghost_paths`, `compile_with_ghosts`). Ghost-Flächen sind die Refract-Ereignisse; je Paar i < j wird an j zurück- und an i wieder vorwärts reflektiert, N Refract-Ereignisse ergeben N(N − 1)/2 Ghosts; Pfade mit Beugungsordnungen oder Kristallmoden sind in M4 ausgeschlossen. Die Medien bestimmt compile mit den Regeln für explizite Pfade. Ghosts sind abgeleitete Daten und stehen nur in einer Datei, wenn der Nutzer sie speichert. Das Ghost-Ranking (#124, `rtt/analysis/ghosts.hpp`: `ghost_ranking`, Nachtrag zu ADR 0027) verfolgt alle Ghosts mit den Startstrahlen des Basispfads und sortiert sie nach der Bestrahlungsstärke am Bild relativ zum Nutzbild, ρ = (P_g/P_b)(r_b² + r₀²)/(r_g² + r₀²) mit dem geometrischen RMS-Radius r und der Auflösungsgrenze r₀ des Detektors; paraxiale Fokuslage und Unschärfe sind Diagnosen.

**Parameter:** Jeder optimierbare Wert ist ein `Param { double value; bool variable; std::optional<double> min, max; std::optional<std::string> param; }` (ADR 0029, Schema 0.4). Ein ungebundenes Param trägt seinen Wert, `variable` und optional Grenzen; ein gebundenes (`param` gesetzt) verweist auf eine Zeile der Parametertabelle `System::parameters` und hat keinen eigenen Wert. Die Tabelle hat Zeilen mit genau einer Form (`value` für alle Konfigurationen, `values` je Konfiguration, `expression` aus früheren Zeilen); `System::configurations` sind ihre Spalten. Ausgewertet wird sie ab #164/#165; bis dahin meldet compile gebundene Params als `param.unresolved`. `pickup` gibt es seit Schema 0.4 nicht mehr (Migration mit `io.pickup_dropped`).

**Validierung:** `rtt::model::validate(system)` liefert `Diagnostic`s mit JSON-Pointer auf die betroffene Stelle (Wellenlängen, Apertur, Element-Regeln inkl. Länge der Materialliste, eindeutige IDs und Namen, gültige Radien und Aperturen, Pfadverweise). Jede Diagnose trägt einen stabilen Code aus der Registry `rtt-diagnostics`, z. B. `material.unknown` (ADR 0022, Liste in [diagnostics.md](diagnostics.md)). `compile()` behält die Warnungen in `CompiledSystem::diagnostics()`, kennt den JSON-Pointer jeder Fläche (`CompiledSurface::location`) und meldet Fehler mit Code in `CompileError`. `ParaxialError` und `AnalysisError` tragen den Ort, wo er bekannt ist; „keine Blende“ ist überall `rtt::compile::NoStopError` aus einer einzigen Prüfung `require_stop`.

## Dateiformat (`rtt-io`, umgesetzt in M0)

Ein System ist eine Datei `*.rtt.json`. Die vollständige Struktur steht in `schema/raytatouille.schema.json`; Referenz ist aber der C++-Parser.

- **Strikt:** unbekannte und doppelte Schlüssel, falsche Typen, Zahlen jenseits des `double`-Bereichs, andere Einheiten und inkompatible `schema_version` sind Fehler mit JSON-Pointer (doppelte Schlüssel und Overflow über `rtt-json`, ADR 0020).
- **Versionierung:** `schema_version` SemVer; aktuell 0.3.0. Vor 1.0 müssen Major und Minor exakt passen oder zu einer älteren Version gehören, die `rtt-io` migriert (derzeit 0.1 und 0.2 → 0.3: `diffract` wird `transmit` mit derselben Ordnung; geschrieben wird immer die aktuelle Version). Jede Formatänderung erhöht die Version und bringt eine getestete Migration mit. Das JSON-Schema beschreibt nur die aktuelle Version.
- **Material:** `"material"` ist ein String (ein Material für alle Segmente) oder ein Array mit einem Eintrag je Segment, z. B. `"material": ["SCHOTT:N-BK7", "SCHOTT:F2"]` (ADR 0017), oder ab 0.3 ein Kristallobjekt `{"ordinary": …, "extraordinary": …}` mit `"optic_axis"` am Element (ADR 0026). Die gelesene Form wird unverändert geschrieben.
- **Bearbeitungsform** (ADR 0024): `rtt::io::to_edit_json` schreibt jeden Wert, auch Standardwerte, und jedes `Param` als Objekt; optionale Mitglieder ohne Wert fehlen. Jeder Pointer in die kanonische Datei gilt auch hier. `rtt::io::apply_patch` ändert ein System mit JSON Patch (RFC 6902) auf dieser Form, ganz oder gar nicht, und liefert den inversen Patch für Undo; abgelehnt wird ein Patch, der Fehler von `validate()` hinzufügt (`rtt::io::EditError`, Codes `edit.*`).
- **Kanonisch:** `rtt::io::to_json` schreibt immer dieselben Bytes: 2 Leerzeichen Einzug, LF, abschließender Zeilenumbruch, Standardwerte weggelassen, kleine Objekte aus Skalaren auf einer Zeile. Für jede kanonische Datei gilt `to_json(parse(text)) == text`. `rtt format` bringt Dateien in diese Form.
- **Parameter:** als Zahl (`51.68`) oder als Objekt (`{"value": 30.0, "variable": true}`).
- **Import v1:** AGF-Glaskataloge der Hersteller (M2, `rtt-material`); ZMX-Dateien (M8) als reines Austauschformat (sequenzieller Teil), nicht unterstützte Flächentypen ergeben eine klare Fehlermeldung. **Export v1:** CSV für Analysedaten. **Später:** STEP-Export der Flächen, ISO-10110-Zeichnung.

Beispiel (gekürzt aus `tests/reference/m0/singlet.rtt.json`):

```json
{
  "schema_version": "0.3.0",
  "name": "Plano-convex singlet f = 100 mm",
  "units": {"length": "mm", "wavelength": "um"},
  "wavelengths": [
    {"um": 0.4861},
    {"um": 0.5876, "reference": true},
    {"um": 0.6563}
  ],
  "aperture": {"type": "epd", "value": 20.0},
  "fields": {
    "points": [
      {},
      {"y": 5.0}
    ]
  },
  "root": {
    "type": "assembly",
    "name": "system",
    "children": [
      {
        "type": "lens",
        "name": "L1",
        "pose": {"position": [0.0, 0.0, 5.0]},
        "material": "SCHOTT:N-BK7",
        "surfaces": [
          {
            "id": "L1.S1",
            "shape": {
              "base": {"type": "conic", "radius": 51.68}
            },
            "aperture": {"type": "circular", "radius": 12.7}
          },
          {
            "id": "L1.S2",
            "pose": {"position": [0.0, 0.0, 4.0]},
            "aperture": {"type": "circular", "radius": 12.7}
          }
        ]
      }
    ]
  },
  "paths": [
    {"name": "main", "events": "auto"}
  ]
}
```

Alle Modellelemente zeigt `tests/reference/m0/feature_tour.rtt.json`, Mehrfachpfade `michelson.rtt.json`, ein Kittglied mit Material je Segment `tests/reference/m2/achromat.rtt.json`.

## Strahlen und Trace-Engines (ab M1)

Tracing läuft immer über ein unveränderliches `CompiledSystem`: Das Modell wird einmal kompiliert, dann tracen beliebig viele Threads lesend darauf.

**Kompilierung (`rtt-compile`, ADR 0016)**

1. Baum flach machen: globale Transformation jeder Fläche (mit den Bezügen nach ADR 0028), Parametertabelle und Konfiguration auflösen (ADR 0029).
2. Materialien bei allen Wellenlängen, Temperatur und Druck auswerten (komplexer Index; bei Doppelbrechung n_o, n_e und optische Achse).
3. Beschichtungen vorberechnen: wahlweise Tabellen über Einfallswinkel je Wellenlänge oder exakt pro Strahl.
4. Flächentypen in `std::variant` auflösen, damit der Hot Path ohne virtuelle Aufrufe auskommt.
5. Ergebnis ist unveränderlich, thread-sicher lesbar und hält keine Zeiger ins Modell. Ein Hash für Caching ist vorgesehen, aber noch nicht umgesetzt (Stand M4; ADR 0016, Folgen).

**Medien entlang eines Pfads** (festgelegt in #5, für Segmente erweitert in #27): Der Strahl startet im Umgebungsmedium. `Refract`, `Ordinary` und `Extraordinary` an einer Fläche eines Elements mit Material wechseln zwischen dem Inneren dieses Elements und der Umgebung (bei mehreren Segmenten siehe unten); `Reflect` und `Transmit` (mit jeder Ordnung) und alle Events an Elementen ohne Material behalten das Medium. Der automatische Pfad besucht alle Flächen in Baumreihenfolge: `Lens`/`Plate` → `Refract`, `Mirror` → `Reflect`, `Stop`/`Detector`/`ThinElement` → `Transmit`. Kittglieder (ADR 0017, festgelegt in #27): Eine `Lens`/`Plate` mit N Flächen hat N − 1 Segmente; Segment i liegt zwischen Fläche i und i + 1, die Kurzform `material` gilt für alle Segmente. Ein `Mirror` mit Substratmaterial und mehr als einer Fläche (Mangin-Spiegel) braucht einen expliziten Pfad (`Refract`, `Reflect`, `Refract`); auf dem automatischen Pfad ist er ein Kompilierfehler (#6). Stand M1: Formen `Plane`, `Conic` und `EvenAsphere`; Pickups werden noch nicht ausgewertet; der Hash folgt, wenn Caching gebraucht wird.

Für `Refract`, `Ordinary` und `Extraordinary` hängt die Regel nur von der Art des Elements und seinen Materialien ab, nie von der Art des Pfads; dieselben Events ergeben also auf automatischen und expliziten Pfaden dieselben Medien. „Einheitlich“ heißt Kurzform oder Liste mit lauter gleichen Einträgen.

| Element | automatischer Pfad | expliziter Pfad |
| --- | --- | --- |
| `Lens`, einheitlich oder gemischt | Segmentregel | Segmentregel |
| `Plate`, einheitlich (Prisma, Würfel) | 2 Flächen: Umschaltregel; mehr als 2: Kompilierfehler | Umschaltregel |
| `Plate`, gemischt (z. B. gekitteter Würfel) | Segmentregel | Segmentregel |
| `Mirror` mit Substrat | 1 Fläche: Umschaltregel; mehr: Kompilierfehler (#6) | Umschaltregel |

*Umschaltregel* (Regel aus #5): An jeder Fläche wechselt der Strahl zwischen dem Inneren und der Umgebung; damit lässt sich ein Prisma durch jede Fläche betreten und verlassen. Bei 2 Flächen ist sie gleich der Segmentregel. Eine einheitliche `Plate` mit mehr als 2 Flächen braucht einen expliziten Pfad, weil die Baumreihenfolge ihrer Flächen keinen sinnvollen Strahlweg ergibt.

*Segmentregel*: An Fläche i führt das Event aus Segment i − 1 ins Segment i (an der letzten Fläche ins Umgebungsmedium) und aus Segment i ins Segment i − 1 (an der ersten Fläche ins Umgebungsmedium); aus einem anderen Segment desselben Elements, das die Fläche nicht begrenzt, ins Umgebungsmedium. Von außen führt die erste Fläche ins Segment 0, die letzte ins Segment N − 2 und eine innere Fläche i in das Material der Segmente i − 1 und i, wenn beide gleich sind; sonst ist das Event ein Kompilierfehler (Seite mehrdeutig, z. B. die Kittfläche eines Achromaten von außen). Eine einheitliche `Lens` mit 3 Flächen ergibt damit Luft → Glas → Glas → Luft, ein Achromat Luft → Glas A → Glas B → Luft.

*Mehrere Elemente:* Elemente schachteln nicht. Ein brechendes Event an Element B, während der Strahl in Element A ist, verlässt A und betritt B wie aus der Umgebung (Medium vor dem Event: das von A); das Verlassen von B führt ins Umgebungsmedium, nicht zurück nach A. Folgen Austrittsfläche von A und Eintrittsfläche von B direkt aufeinander, läuft der Strahl über das Umgebungsmedium, auch bei Nullabstand; Kittglieder werden ausschließlich als ein Element mit Segmenten beschrieben.

**Strahl-Batch (Structure-of-Arrays)**

| Feld | Typ je Strahl | Bedeutung |
| --- | --- | --- |
| `pos`, `dir` | 2 × 3 × `double` | Ort und Einheitsrichtung, global; im Kristall ist `dir` die Energierichtung Ŝ |
| `wave`, `mode_index` | 3 × `double`, `double` | Wellennormale k̂ (global) und n der Kristallmode; isotrop `wave` = `dir`, `mode_index` = 0 (ADR 0026, Punkt 3) |
| `wl` | `uint16_t` | Index der Systemwellenlänge |
| `opl` | `double` | akkumulierter optischer Weg in mm |
| `prt` | 9 × `std::complex<double>` | akkumulierte 3×3-Polarisations-Raytracing-Matrix P, leistungsnormiert |
| `weight` | `double` | Leistung für eine unpolarisierte Quelle (Start 1) |
| `status` | `uint8_t` (enum) | `Alive`, `Missed`, `Vignetted`, `Tir`, `NoConvergence`, `Absorbed`, `EventImpossible`, `Evanescent` (ADR 0025, angehängt; `kRayStatusCount` = 8) |
| `field`, `pupil` | `uint16_t`, 2 × `double` | Herkunft für Analysen |
| `last_surface` | `uint32_t` | letzte getroffene Fläche |

Der Eingangs-Jones-Vektor wird nicht gespeichert: P wirkt auf jeden Eingangszustand (|P E_in|² ist der Leistungsanteil, die Phasen sind die der Feld-PRT; die Feldamplitude unterscheidet sich um die Faktoren √c der Grenzflächen, siehe unten), damit sind polarisiert und unpolarisiert mit demselben Trace auswertbar.

**Semantik von `prt` und `weight` (ADR 0021, umgesetzt in #61):** P ist **leistungsnormiert**: Jede Grenzfläche trägt prt_matrix mit a_s√c_s, a_p√c_p bei (c nach Byrnes Gl. (21), (22) bei Transmission, c = 1 bei Reflexion), sodass |P E|² für |E| = 1 der Leistungsanteil ist und die Phasen unverändert bleiben. P ist damit nicht die Feldamplituden-PRT-Matrix nach Chipman/Lam; die Jones-Pupille (M6) muss die Faktoren √c berücksichtigen, wenn sie Feldamplituden braucht. `weight` = s·½‖P_T‖²_F mit P_T = P − k k₀ᵀ und s den polarisationsunabhängigen Faktoren (Effizienz η einer Beugungsordnung, ADR 0025; Volumenabsorption exp(−4πκt/λ_vac), Absorber 0); polarisierte Leistung = weight·|P E|²/(½‖P_T‖²_F). Interaktionen je Event: Fresnel bei Refract (before → after) und Reflect (gegen `medium_beyond`); ein Mirror ohne Material reflektiert als idealer Leiter (−1, +1), jede andere Fläche ohne Medienwechsel mit r = 0 (weight 0, Status Alive). IdealMirror nur Reflect, IdealAntiReflection Refract (1, 1) und Reflect (0, 0) (verschwindende Reflexion wie bei Fresnel ohne Medienwechsel: weight 0, Status Alive), IdealBeamSplitter Reflect (−√R_s, √R_p) sowie Refract/Transmit (√(1−R_s), √(1−R_p)), CoatingRef Refract/Reflect mit dem Stapel aus `rtt-coating` (aus dem Substrat umgekehrt, ADR 0019), Polarisator und Retarder nur Transmit (Achse von compile global gedreht, `CompiledSurface::ideal_axis`). Transmit an Fresnel-, AR- und Coating-Flächen ist ein Dummy-Durchgang ohne Wirkung. Jede andere Kombination ergibt `EventImpossible`.

**Layout, Geometrie-Export (umgesetzt in #81, ohne ADR):** `rtt/compile/layout.hpp` und `CompiledSystem::elements()` liefern die Geometrie zum Zeichnen eines Linsenschnitts.
- **Elemente und Flächen:**
  - `CompiledElement` enthält Name, Art, die Flächen [first_surface, + surface_count), die Medien des Körpers und das Flag `segmented`.
  - `CompiledSurface` hat zusätzlich `element`, `medium_front` und `medium_back`.
- **Medien vor und hinter einer Fläche:** in der Flächenreihenfolge des Elements, so wie eine Brechung durch die Flächen in dieser Reihenfolge, von der Umgebung kommend, sie sähe (Regeln von ADR 0017).
  - Linsen und segmentierte Platten: Segmente.
  - Platten aus einem Material, Prismen und Spiegel mit Substrat: Wechsel zwischen innen und Umgebung.
  - Elemente ohne Material: Umgebung auf beiden Seiten.
  - Die Werte sind unabhängig von den Pfaden. Die Medien eines Pfads stehen in seinen Events; Reflexionen und rückwärts laufende Pfade können eine Fläche anders sehen.
  - Bei Platten und Prismen aus einem Material mit mehr als zwei Flächen ist das nur die Konvention der Flächenreihenfolge (jede Brechung wechselt innen ↔ Umgebung). Welche Seite Glas ist, ergibt sich aus `CompiledElement::media` und `segmented`.
- **sag und Normale:** in lokalen Koordinaten, vektorisiert. Die Normale (−∂z/∂x, −∂z/∂y, 1)/|…| zeigt am Scheitel nach +z wie in `rtt-geom`; sie gibt es lokal oder global. Außerhalb des Formbereichs kommt NaN; die Apertur wird nicht angewendet.
- **Profile:** nur in Schnittebenen, die zur lokalen z-Achse der Fläche parallel sind (|n·z_lokal| ≤ 1e−12). Das umfasst Meridionalschnitte, Kippungen in der Ebene, Faltspiegel und gegen die Achse versetzte Ebenen.
  - Dann ist der Schnitt exakt die Kurve z = sag entlang der Geraden „Ebene ∩ lokale x-y-Ebene“, ohne Nullstellensuche, in Richtung t = z_lokal × n.
  - Begrenzung: Die Apertur wird analytisch geschnitten (Kreis, Ring mit zwei Stücken, Rechteck, Ellipse), dazu der Formbereich `max_radius`. Eine Fläche ohne beides ist ein Fehler.
  - Am Formbereich endet das Profil um 1e−12·`max_radius` innerhalb, weil die Sag dahinter NaN ist und am Rand eine unendliche Steigung hat. Eine Fläche ohne Apertur wird bis zum Formbereich (z. B. Halbkugel) gezeichnet; für Zeichnungen gibt man Linsenflächen Aperturen.
  - Andere Ebenen ergeben einen Fehler mit Hinweis auf M10; allgemeine Schnitte kommen mit den 3D-Meshes.
- **Umrisse:** ein geschlossenes Polygon je Glassegment j, zusammengesetzt aus dem Profil von Fläche j, der Randkante, dem Profil von Fläche j + 1 rückwärts und der zweiten Randkante.
  - **Randkante als Stufe:** gerechnet in den Schnittkoordinaten von Fläche j. Enden die Ränder an unterschiedlicher Stelle der Schnittgeraden, läuft die Kante parallel zur Achse von Fläche j beim äußersten der beiden Ränder (vom Glas aus gesehen) und dann entlang der Schnittgeraden zum anderen.
    - Am Außenrand ist das der Zylinderrand beim größeren Radius mit Schulter, wie in Linsenzeichnungen.
    - An einer Bohrung ist es die Bohrung beim kleineren Radius.
    - Eine umgedrehte Fläche j + 1 wird in die Richtung von Fläche j gedreht.
  - Ein Ring auf beiden Flächen ergibt je Seite der Bohrung ein Polygon.
  - Nur Linsen und segmentierte Platten haben Umrisse. Alle Platten aus einem Material (auch Fenster mit zwei Flächen, Prismen, Rhomben), Spiegel, dünne Elemente, Blenden und Detektoren liefern eine leere Liste; für sie zeichnet man die Profile. Umrisse für Platten aus einem Material sind ein Folgeauftrag.

**Engine 1, paraxial (`rtt-paraxial`):** y-nu-Trace für rotationssymmetrische Systeme: EFL, BFL, Hauptebenen, Pupillen, Vergrößerung, Seidel-Koeffizienten pro Fläche (M2, umgesetzt in #30: S_I–S_V je Event und als Summe aus Randstrahl (axialer Objektpunkt, Rand der paraxialen Eintrittspupille bei +y) und Hauptstrahl des maximalen Feldes (radialer Feldwert, Richtung nach der Feldkonvention aus #8 (Winkel: Hauptstrahl steigt in +y, bei endlichem Objekt liegt der Objektpunkt dann bei −y; Höhen: Punkt bei +y), durch die EP-Mitte; paraxiale Bildhöhen und Feldwinkel bei endlichem Objekt werden wie in `rtt-trace` bei der Referenzwellenlänge umgerechnet, auch wenn Seidel bei einer anderen Wellenlänge gerechnet wird (#35); S_II und S_V wechseln mit dieser Richtung das Vorzeichen), Beiträge von Konik und A4, chromatische Terme C_L und C_T für ein wählbares Wellenlängenpaar; Normierung W040 = S_I/8 usw. in mm, Quellen in docs/quellen.md). v2: parabasal für dezentrierte Systeme. Konventionen (festgelegt in #7): u = dy/dz global; der Index trägt ein Vorzeichen (positiv bei Ausbreitung in +z) und wechselt es an jedem Spiegel, sodass n·u der optische Richtungskosinus ist und Brechung n′u′ = nu − yφ mit φ = c(n′ − n) sowie Transfer mit globalem Δz auch nach Spiegeln gelten. Brennweiten sind positiv für sammelnde Systeme (Hohlspiegel: EFL = |R|/2); alle Lagen als globale z-Koordinaten, FFL/BFL ab dem ersten bzw. letzten strahlverändernden Scheitel. Nicht rotationssymmetrische Pfade (dezentriert, gekippt, Phasenschichten, Beugung, Doppelbrechung) ergeben einen `ParaxialError`. Da Brechzahlen absolut sind (`AIR` = Ciddor-Luft, n ≈ 1,00027), ist EFL = 1/Φ im absoluten Sinn und `rear_focal_length` = |n′|/Φ der physikalische Abstand; Programme, die relativ zu Luft rechnen (z. B. OpticStudio), geben für ein System in Luft eine um den Faktor n_Luft größere EFL an. Dasselbe gilt für die Systemapertur `image_fnumber`: EPD = EFL/F# mit der absoluten EFL, also bei gleichem F# eine um den Faktor n_Luft kleinere EPD als in OpticStudio. Literaturwerte mit Luft = 1 vergleicht man daher mit `VACUUM` als Umgebung oder rechnet um (#25).

Prescription-Daten (#84, `rtt/paraxial/prescription.hpp`, Python `rt.paraxial.prescription`):

- Je Event des Pfads Scheitel-z, vorzeichenbehafteter Index nach dem Event, Rand- und Hauptstrahl (y, u nach dem Event, paraxialer Einfallswinkel i = u + y c mit u vor dem Event, Sasian L4 S. 24) und die Lagrange-Invariante n′(ū′y − u′ȳ) nach dem Event.
- Es sind dieselben Strahlen wie in `seidel` (eine interne Konstruktion).
- `total_track` = Σ|Δz| zwischen aufeinanderfolgenden Event-Scheiteln, bei Spiegeln also die abgewickelte Länge; `object_distance` = erster Scheitel − Objekt-z bei endlichem Objekt.
- `paraxial_working_f_number` = 1/(2|n′u′|) und `paraxial_image_na` = |n′u′| des Randstrahls nach dem letzten Event (Greivenkamp OPTI-502 S. 9-34 bis 9-37). Die Namen sind ausdrücklich paraxial: Ein reales Arbeits-F/# aus einem realen Randstrahl kann später dazukommen. Das Arbeits-F/# gilt für jeden Lichtkegel (Greivenkamp S. 9-36), also auch afokal mit endlichem Objekt. Afokal mit Objekt im Unendlichen ist der Randstrahl parallel: NA = 0 exakt und kein F/#; das entscheidet das Afokal-Kriterium von `first_order`, nicht die Rundung von u′.
- H = n(ūy − uȳ) im Objektraum, wie `Seidel::lagrange`. Lateral- und Winkelvergrößerung kommen aus `first_order` (Feld `first_order` des Ergebnisses), ohne zweite Definition.
- Anders als `seidel` ist eine fehlende Blende oder Pupille kein Fehler. Werte ohne Strahl sind leer (Python `None`, in den Arrays NaN), wie bei `first_order`; der Randstrahl existiert bei einem Objekt im Unendlichen auch ohne Blende, wenn die Systemapertur den EPD festlegt. Das bleibt bewusst so, wenn #86 `NoStopError` einführt.

**Engine 2, sequenziell/multi-sequenziell (`rtt-trace`):** pro Pfad und Event: lokal transformieren → Schnitt → Apertur → Normale → Phasengradient → Event (Brechung, Reflexion, Beugung, o-/e-Strahl) → Interaktion auf P und `weight` → global transformieren. Ist ein Event physikalisch unmöglich (z. B. TIR statt Brechung), endet der Strahl mit Status. Aperturränder sind inklusiv mit der Toleranz `kApertureTolerance` = 1e−9 mm (gleich der Aiming-Toleranz), am Außen- wie am Innenrand; ein auf den Blendenrand gezielter Strahl wird damit nicht durch Rundung vignettiert (#50). Ray Aiming per Newton auf die Blendenmitte; Aperturtypen EPD, Bildraum-F/#, objektseitige NA, Blendengröße; Feldtypen Winkel, Objekthöhe, paraxiale Bildhöhe; Verteilungen hexapolar, Gitter, Zufall mit Seed, Fans, Einzelstrahl.

**Strahlpfad-Aufzeichnung (umgesetzt in #80, ohne ADR):** `SequentialTracer::trace(system, path, rays, paths, record_rays, max_recorded_rays)` zeichnet zusätzlich in `RayPaths` (`rtt/trace/ray_paths.hpp`) auf. `RayBatch` und die Tracer-Schnittstelle bleiben unverändert.
- **Recorder-Policy:** Die Strahlschleife ist ein Template über einen Recorder. Das normale `trace()` nutzt `NoRecord`, dessen Aufrufe leere Inline-Funktionen sind. Ohne Aufzeichnung ist der Hot Path damit derselbe, und die Ergebnisse sind bitgleich. Mit Aufzeichnung endet der Batch ebenfalls bitgleich.
- **Layout:** Structure of Arrays in C-Reihenfolge, damit Python Views ohne Kopie bekommt. Es gibt S = Ereignisse + 1 Slots: Slot 0 ist der Startzustand, Slot k der Zustand nach Ereignis k − 1.
  - Je Strahl und Slot: Position, Richtung (global), OPL, `weight` und Status. Kein P; P je Ereignis folgt mit M6.
  - Je Strahl: `count` (gültige Slots) und `lost_at`. Je Pfad: `event_surfaces`.
  - Ein bei Ereignis j verlorener Strahl hat `count` = j + 2 und `lost_at` = j. Slot j + 1 ist der Zustand beim Verlust mit dessen Status, die späteren Slots sind NaN mit dem Verluststatus.
  - Ein durchgekommener Strahl hat `count` = S und `lost_at` = −1. Ein Strahl, der schon beim Start nicht Alive war, hat `count` = 1 und `lost_at` = −1; seine späteren Slots sind NaN mit dem Startstatus.
  - Einen Verlust erkennt man an `lost_at` ≥ 0, nicht an `count` < S: Ein Verlust beim letzten Ereignis hat `count` = S. Die Fläche ist `event_surfaces[lost_at]`.
  - Slot `count` − 1 ist bitgleich mit dem Endzustand im `RayBatch`.
  - Bei einer Exception bleibt `paths` unverändert (starke Garantie).
- **Auswahl statt Abschneiden:** `record_rays` wählt die Strahlen explizit aus, in der angegebenen Reihenfolge und ohne Duplikate; `ray_indices` gibt die Zuordnung zum `RayBatch`. „Alle Strahlen“ heißt in C++ `std::nullopt` und in Python `record_rays=None`; eine leere Auswahl ist in beiden ein Fehler (entschieden in der Prüfung von #88). Ohne Auswahl werden alle Strahlen aufgezeichnet, höchstens `max_recorded_rays` (Standard 10 000, etwa 65 Byte je Strahl und Slot). Sonst gibt es einen `invalid_argument` mit dem Speicherbedarf. Die ersten Strahlen eines Hexapolar-Bündels wären nur die inneren Ringe; stilles Abschneiden würde eine GUI falsch zeichnen lassen.
- **Determinismus:** Jeder Strahl schreibt nur seine eigenen Slots; das Ergebnis ist unabhängig von der Thread-Anzahl (ADR 0004).

**Abbruch und Fortschritt (#83, `rtt/trace/run_control.hpp`):** Lange Rechnungen lassen sich aus einer GUI abbrechen und melden Fortschritt.
- **Typen:** `CancelToken` (geteiltes Flag, `request_cancel()` aus jedem Thread), `RunControl` (Token, Callback `progress(Progress{done, total, stage})`, Mindestabstand der Meldungen, `block_size`) und die Exception `Cancelled`.
- **Wo:** Überladungen mit `const RunControl&` als letztem Parameter: `SequentialTracer::trace`, `make_rays` und alle Bündel- und Verlaufsanalysen (`spot`, `ray_fan`, `opd_map`, `opd_fan`, `distortion`, `field_curvature`, `longitudinal_colour`, `lateral_colour`) sowie die Pfadauswertung (`path_transmission`, `opl_difference`, #122). Die alten Signaturen bleiben und rufen die neuen mit leerem Control.
- **Abbruch und Fortschritt sind Laufzeitsteuerung, nie Teil der Options** (sie gehören nicht zur Rechnung). Ein leerer Control rechnet exakt wie ohne; ein aktiver ändert kein Ergebnis, bitgleich und unabhängig von Blockgröße und Thread-Anzahl (ADR 0004, Nachtrag #83).
- **Ablauf:** Mit aktivem Control wird in Blöcken von `block_size` Strahlen gerechnet (Standard 256).
  - Vor jedem Block wird das Token geprüft. Nach einer Abbruchanfrage beginnt kein neuer Block, und nach dem parallelen Teil kommt `Cancelled`, nie aus einem Worker (Regel 3). Die Rückkehr erfolgt nach höchstens einem Block je Worker. **Dokumentierte Abbruchzeit je Stufe:** `trace` deutlich unter 1 ms (256 Strahlen typischer Linsen), `aim` bei realer Zielung etwa 5 ms (256 × etwa 20 µs je Worker, Cooke-Triplett; seit #119 zielt `make_rays` parallel, ADR 0004), `field` und `wavelength` ein Feldpunkt bzw. eine Wellenlänge (eine Handvoll Strahlen). Der Test prüft mit Reserve < 1 s.
  - Der Callback wird serialisiert, gedrosselt und aus einem beliebigen Thread gerufen. Seine Exception beendet den Lauf und wird danach weitergereicht.
  - Stufen: `aim` und `trace` je Strahlbündel, `field` für Feldverläufe, `wavelength` für die Farbanalysen. Jede vollständige Stufe endet mit done = total.
- **Ohne Abbruch** bleiben Seidel, `first_order`, `prescription` und die Einzelpunkt-Funktionen (`*_at`); sie dauern Mikro- bis Millisekunden.
- **Python:** `rt.CancelToken` (`cancel()`, `cancelled`) und die Schlüsselwortargumente `cancel=` und `progress=` an `rt.trace.trace`, `rt.trace.make_rays` und den acht Analysen; `progress(done, total, stage)` hat den Typ `rt.trace.ProgressCallback`. Ein Abbruch kommt als `rt.errors.Cancelled`, eine Exception des Callbacks unverändert.
  - GIL: Das Control wird mit GIL gebaut, die Rechnung läuft ohne GIL. Der Callback holt sich das GIL im rufenden Worker und ruft das Python-Objekt über einen Zeiger, damit kein Worker ohne GIL dessen Referenzzähler ändert. Argumente und Control werden mit GIL zerstört.
  - `trace(..., record_path=True)` nimmt weder `cancel` noch `progress` (ValueError): Die Aufzeichnung ist für die wenigen Strahlen einer Zeichnung gedacht.
- **Benchmark:** `apps/rtt-bench` misst den Spot des Cooke-Tripletts mit 1027 Strahlen (Zielung, Trace, Spot; 1 und alle Threads). Die CI schreibt die Zeiten in die Job-Summary und warnt erst über dem 3-Fachen des Richtwerts von 10 ms. Der Richtwert ist laut Maintainer nur ein grober Anhaltspunkt; die CI schlägt deshalb nie fehl.

**Engine 3, nicht-sequenziell (erste Priorität nach v1.0):** BVH, Ray Splitting, Strahlbaum mit Abbruch nach Tiefe und Energie, Quellen, Detektoren. Gleiche Physik wie Engine 2; deshalb ist die Event-Ausführung als `apply_event(ray, hit, kind)` gekapselt und jedes Review ab M1 prüft, dass Schnitt, Interaktion und Event-Ausführung pfadunabhängig bleiben.

```cpp
class Tracer {
 public:
  virtual ~Tracer() = default;
  virtual TraceStats trace(const CompiledSystem& sys, PathId path, RayBatch& rays) const = 0;
};
```

## Physik-Module (ab M1)

Jede Physik-Bibliothek ist eine reine Funktionssammlung ohne Kenntnis von Modell oder Tracer und wird isoliert gegen Referenzwerte getestet.

**Geometrie (`rtt-geom`)**

```cpp
template <rtt::math::Real T>
class Shape {
 public:
  virtual ~Shape() = default;
  virtual T sag(T x, T y) const = 0;
  virtual std::pair<T, T> grad(T x, T y) const = 0;  // dz/dx, dz/dy
  virtual std::pair<T, T> base_conic() const = 0;    // (c, k) als Startwert
  virtual std::optional<T> max_radius() const = 0;
};
```

- Schnitt: analytisch mit der Basis-Konik, dann Newton auf F(t) = z(t) − sag(x(t), y(t)); konvergiert bei |F| ≤ 1e−12 mm oder |Δt| ≤ 8·ε_double·(1 + |t| + ‖o‖) (Rundungsgrenze bei fernem Strahlursprung), höchstens 30 Iterationen, sonst `NoConvergence`. Die nächste Lösung mit t > ε in Ausbreitungsrichtung gilt; ε = 1e−9 mm (`kDefaultTMin` in `rtt-geom`), überall im Schnitt einheitlich zu verwenden.
- Normale aus dem Gradienten: n = (−∂z/∂x, −∂z/∂y, 1) normiert.
- Phasenschichten (#126, ADR 0025, `rtt/geom/phase.hpp`): `LinearGratingPhase` (φ = 2πG(x cos ψ + y sin ψ)) und `RadialPhasePolynomial` (φ = Σ c_k ρ^{2k}) mit φ in rad und ∇φ in rad/mm über den lokalen lateralen Koordinaten; die Phase einer Fläche ist die Summe ihrer Schichten (`phase`, `phase_grad`). `tangential_gradient` projiziert ∇φ mit der Einheitsnormalen in die Tangentialebene, (I − N Nᵀ)∇φ, wie es die lokale Gittergleichung braucht (Mansuripur Gl. (11) für Rotationsflächen mit rotationssymmetrischem φ, `docs/quellen.md`). Die Konstruktoren prüfen die Werte (`std::invalid_argument`), die Auswertung wirft nicht.
- Flächentypen v1: Ebene, Sphäre/Konik, gerade und ungerade Asphäre, Forbes Q-con und Q-bfs, Biconic, Toroid/Zylinder, Axicon, Zernike-Sag (Standard, Fringe), XY-Polynom, Grid-Sag (bikubisch).

**Materialien (`rtt-material`):** Sellmeier (1–5), Schott, Herzberger, Conrady, Extended 2 und 3, Cauchy, tabellierte n/k. AGF-Import; Herstellerkataloge werden nicht mitgeliefert. Unterstützte AGF-Formelnummern sind nur die mit belegter Reihenfolge der CD-Koeffizienten (`docs/quellen.md`): 1, 2 (#24), 6, 12 (a₇ = 0) und 13 (#42); CD- und TD-Datensätze dürfen auf Folgezeilen ohne Kennung umbrechen (#42). Für eine GUI (#85): `MaterialLibrary::catalog(name)` liefert jeden geladenen Katalog so, wie er in der Datei steht (alle Gläser mit NM-Extras, GC, ED, TD, MD, OD, IT; N(d)/ν(d) aus dem NM-Datensatz, berechnete Werte über `index_many`, bitgleich zum skalaren `index`); `add_catalog_text` lädt aus dem Speicher. **Katalog-Alias (Entscheidung #85):** `add_catalog(datei, name)` lädt eine einzelne Datei unter `name` statt des Dateistamms (z. B. zwei `schott.agf` als `SCHOTT` und `SCHOTT_M2`); Referenzen bleiben `KATALOG:GLAS`, der Alias ersetzt nur KATALOG. Ohne Alias ändert sich nichts. Namen werden wie alle Referenzen mit Groß- und Kleinschreibung verglichen (`schott` ≠ `SCHOTT`); ein Alias ist nicht leer, ohne `:` und Leerraum und nicht `CONST`; ein Name für ein Verzeichnis ist ein Fehler. Doppelte Glasnamen (#71, Nachtrag ADR 0008): ein Glas mit identischen Daten wird einmal gelesen, ein Glas mit verschiedenen Blöcken ist mehrdeutig (`UnknownMaterial` mit beiden Zeilen und Alias-Hinweis, der Rest des Katalogs lädt); Textzeilen vor dem ersten Datensatz und ein einzelnes Wort vor dem nächsten NM werden übersprungen, jeweils mit `LoadWarning` (`MaterialLibrary::load_warnings()`). Luft nach Ciddor, Glasdaten relativ zu Luft, dn/dT nach Schott-Modell, CTE. Umgesetzt in #25: `AIR` ist trockene Luft (0 % Feuchte, 450 ppm CO₂) nach Ciddor bei Systemtemperatur und -druck; Katalogglas wird mit λ_rel = λ_vac/n_Luft(λ_vac, T_ref, 1 atm) ausgewertet und als n_abs = n_rel·n_Luft + Δn_abs(T − T_ref) (SCHOTT TIE-19) absolut geliefert; `VACUUM` bleibt 1, `CONST:` und eigene Formeln bleiben absolut und temperaturunabhängig. Metalle als n + iκ. Einachsig doppelbrechend in v1, zweiachsig in v2. Umgesetzt in #129 (ADR 0014 Punkt 4, ADR 0026 Punkt 1): `UniaxialMaterial` (`rtt/material/uniaxial.hpp`) hält zwei gewöhnliche Materialien für n_O und n_E, `MaterialLibrary::resolve_uniaxial(ordinary, extraordinary)` löst beide über `resolve` auf (Fehler mit Präfix `ordinary: ` bzw. `extraordinary: `); der Wellenlängenbereich ist die Schnittmenge, κ wird unverändert durchgereicht (κ ≠ 0 meldet compile ab #131 als `crystal.absorbing`). Die optische Achse gehört zum Element, nicht zum Material.

**Beschichtungen (`rtt-coating`):** Transfermatrix-Methode (Abelès) für s und p mit komplexen Indizes, Ausgabe r_s, r_p, t_s, t_p. Dicke physikalisch oder als QWOT. Ideale und tabellierte Coatings. Invariante R + T + A = 1. Umgesetzt in #58: Transfermatrix in der Grenzflächen-Formulierung nach Byrnes (Gl. 8–15, nicht die Admittanzen nach Macleod, weil diese für p meist die andere Vorzeichenkonvention haben); Schichten vom Einfallsmedium zum Substrat, Dicke in µm (`_um`) oder als QWOT d = q·λ₀/(4 Re n(λ₀)); Kerne als Templates über `Real`, noexcept, Prüfung an der API-Grenze mit `check_stack`. Ideale Coatings mit Reflexionsgrad R ∈ [0, 1] ohne Retardance (r_s = −√R, r_p = +√R, Transmission verlustfrei). Bewusste Modellannahme (festgelegt in #58): Unter Totalreflexion verhält sich ein ideales Coating mit R < 1 wie die unbeschichtete Grenzfläche, einschließlich der TIR-Retardance; der ideale Spiegel (R = 1) bleibt auch dann bei r_s = −1, r_p = +1. Tabellierte Coatings über Einfallswinkel (Gitter muss [0, π/2] abdecken) und Wellenlänge, bilinear in Real- und Imaginärteil, Auswertung noexcept; der Wellenlängenbereich wird bei `compile()` geprüft (#61). Kataloge (ADR 0019, #59): eigene JSON-Dateien (`schema/raytatouille-coatings.schema.json`), geladen mit `CoatingLibrary::add_catalog` wie Glaskataloge und referenziert als `CATALOG:NAME` aus `CoatingRef`; je Design Schichten von der Umgebungsseite zum Substrat mit Materialverweis und `thickness_um` oder `qwot` sowie Designwellenlänge. `compile(system, materials, coatings)` löst die Designs auf, wertet die Schichten je Systemwellenlänge bei Umgebungstemperatur und -druck aus (`CompiledSystem::coatings()`) und legt je Fläche das Substratmedium fest: das Innere des Elements (angrenzendes Segment, Material der Platte bzw. des Spiegelsubstrats); Licht von innen sieht den umgekehrten Stapel. Innere Flächen zwischen zwei Segmenten, Spiegel ohne Material, Stop, Detector und ThinElement mit Coating sind Kompilierfehler.

**Polarisation (`rtt-polar`):** 3D-Polarisations-Raytracing nach Chipman: P_q = O_out · J · O_inᵀ, O_in = [s, p_in, k_in], O_out = [s, p_out, k_out], s = k_in × n normiert, J = diag(a_s, a_p, 1). Bei senkrechtem Einfall s aus stabiler Referenzrichtung. P_gesamt = P_q · P_gesamt. Polarisator (Achse projiziert senkrecht zu k, Extinktionsverhältnis), Retarder (linear, zirkular, allgemeine Jones-Matrix). Depolarisierende Mueller-Elemente erst v2. Auswertung: Diattenuation, Retardance, Jones-Pupille, Stokes. **Umgesetzt in #57** (Quelle: W.-S. T. Lam, Dissertation, Arizona, Gl. (3.1)–(3.9), (4.3)–(4.6), Abschnitt 4.5.2): `prt_basis`, `prt_matrix` und `geometric_transform` sind Templates über `Real`; `diattenuation`, `retardance` und `physical_retardance` gibt es nur für double (Auswertung, nicht im Ableitungspfad). Für einfallende und ausfallende Welle gilt dasselbe s (bei isotroper Fläche ist s′ = s); p_out = k_out × s. Bei senkrechtem Einfall (|k_in × N| < 1e−12) ist s = k_in × e normiert, mit der globalen Achse e mit kleinstem |k_in·e|, bei Gleichstand in der Reihenfolge x, y, z; das ist zustandslos und deterministisch. P hängt dann nicht von s ab, solange J in der Transversalebene isotrop ist (Fresnel, Coatings): Reflexion mit a_p = −a_s und p_out = −p_in ergibt P = a_s(I − kkᵀ) − kkᵀ, Transmission mit a_p = a_s ergibt P = a_s(I − kkᵀ) + kkᵀ. Q = s sᵀ ± p_out p_inᵀ + k_out k_inᵀ (+ bei Brechung, − bei Reflexion; das Element (3,3) ist 1, Lam druckt in Gl. (4.6) 0). Diattenuation D = (Λ₁² − Λ₂²)/(Λ₁² + Λ₂²) aus den Singulärwerten des Transversalanteils P − k_out k_inᵀ. Sie beruht auf Feldamplituden: Bei Reflexion und bei Transmission zwischen verlustfreien Medien ist sie gleich der Leistungs-Diattenuation, nur bei Transmission mit Absorption weicht sie ab. Retardance δ ∈ [0, π] über die Polarzerlegung (M_R = U V† aus der SVD), gerechnet im 2×2-Block der Transversalebene, damit die schnelle Achse auch bei entarteten Eigenwerten transversal bleibt. Die schnelle Achse ist der Eigenvektor mit der kleineren Phase; bei δ = 0 oder π ist sie nicht eindeutig und folgt der Reihenfolge des Lösers. Die physikalische Retardance ist die von Q⁻¹P und entspricht bei einer einzelnen Reflexion |arg(−r_p/r_s)|. **Retardance einzelner Reflexionen (Entscheidung zu #57):** Wegen der Referenz Q mit diag(1, −1) bei Reflexion ist die physikalische Retardance einer einzelnen Reflexion um π gegenüber Hechts relativer Phase δ_p − δ_s verschoben. Beispiele: Totalreflexion am Grenzwinkel (r_s = r_p = 1) hat Retardance π; eine Totalreflexion mit Hechts relativer Phase 45° (Fresnel-Rhombus) hat 135°. Über zwei Reflexionen hebt sich der π-Anteil auf; der Rhombus wirkt insgesamt als λ/4. Abnahmetests (#60, #63) prüfen deshalb beobachtbare, konventionsfreie Größen: Ausgangszustand bzw. Stokes-Größen, Malus-Transmission, Zirkularität nach λ/4, Drehwinkel nach λ/2, bei Rhombus und Spiegelfolgen den Gesamtdurchgang. Die Retardance einzelner Reflexionen wird nur zusammen mit dieser Konvention geprüft. **Ideale Elemente (umgesetzt in #60, `ideal.hpp`):** dünne Elemente ohne Ablenkung (P k = k); eine allgemeine Jones-Matrix wird in der Basis (e1 = Achse senkrecht zu k projiziert, e2 = k × e1) eingebettet. Polarisator: Durchlassachse mit Amplitude 1 ohne Phase, gekreuzte Achse mit √ε, ε = T_min/T_max als Leistung (D = (1 − ε)/(1 + ε)). Linearer Retarder: schnelle Achse e^{−iδ/2}, langsame e^{+iδ/2} mit δ = 2π·Retardance in Wellen; die symmetrische Phase fügt der Wellenfront keinen Piston hinzu. Zirkularer Retarder: rechtszirkular um δ verzögert, also Drehung linearer Polarisation um +δ/2 von e1 nach e2. Die λ/2-Platte entspricht Lams Matrix (S. 123) bis auf die globale Phase −i. Geprüfte Varianten werfen `std::invalid_argument` (Achse ∥ k oder null, ε ∉ [0, 1], |k| ≠ 1, nicht endlich).

**Doppelbrechung:** o- und e-Wellenvektor an der Grenzfläche; der e-Strahl läuft entlang des Poynting-Vektors, der OPL über den Wellenvektor. Amplituden aus den Eigenpolarisationen (McClain/Chipman). Der Pfad wählt `Ordinary` oder `Extraordinary`. Umgesetzt in #130 (`rtt/polar/birefringence.hpp`, Lam Gl. (2.11), (2.39), (2.41)): `uniaxial_mode` löst k = κ + βN aus dem Tangentialanteil κ (mit Beugungsterm, ADR 0025) für die o-Mode (|k| = n_O, S = k) und die e-Mode (K-Fläche, S als ihre Normale mit S·N > 0), `isotropic_from_tangential` den Austritt; dazu die Eigenpolarisationen und die PRT-Matrizen der Projektion (ADR 0026, Punkt 5). Keine reelle Lösung und streifender Austritt ergeben `std::nullopt`.

**Kristalle im Pfad (ADR 0026 Punkt 4, umgesetzt in compile mit #131):** Jedes Kristallelement hat ein eigenes Medium (`CompiledMedium` mit n_O, n_E und globaler, normierter optischer Achse); jedes Ereignis im Kristall trägt die Mode seines Eintritts (`CompiledEvent::crystal_mode`). Die Medienregeln bleiben wie oben. B ist das Medium vor dem Ereignis, A danach, „jenseits“ die andere Seite der Fläche; Fehler stehen am Pointer des Ereignisses.

| Fall | erlaubt | sonst |
| --- | --- | --- |
| Eintritt (B isotrop, A Kristall) | `ordinary` oder `extraordinary`, `order` nach ADR 0025 | `refract`: `paths.crystal_mode_required` |
| Austritt (B Kristall, A isotrop) | `refract`, `order` nach ADR 0025 | `ordinary`/`extraordinary`: `paths.mode_without_crystal` |
| `ordinary`/`extraordinary` an einem Element ohne Kristall | – | `paths.mode_without_crystal` |
| Kristall → Kristall | – | `crystal.unsupported` |
| `reflect` mit B Kristall | – | `crystal.unsupported` |
| `reflect` von außen (jenseits Kristall) | nur mit `ideal_anti_reflection` (r = 0) | `crystal.unsupported` |
| `transmit` mit B isotrop an einer Kristallfläche | ja (Durchgang ohne Wirkung) | – |
| `transmit` mit B Kristall | nur `order` = 0, Mode bleibt | `crystal.unsupported` |
| automatischer Pfad durch ein Kristallelement | – | `crystal.unsupported` (am Pfad) |

Interaktionen an Kristallflächen in M4: `fresnel` und `ideal_anti_reflection`, sonst `crystal.interaction_unsupported`; ein Teil mit κ ≠ 0 bei einer Systemwellenlänge ist `crystal.absorbing`. Haben die Wellenlängenbereiche von n_O und n_E keinen gemeinsamen Bereich, meldet `material.wavelength_out_of_range` das am Material. **Tracer (umgesetzt in #132):** Der Strahl läuft im Kristall entlang `dir` = Ŝ, die Phase entlang `wave` = k̂: OPL += `mode_index` · l · (k̂·Ŝ) (Lam Gl. (2.17)); isotrope Medien behalten n · l und sind bitgleich. Leseregel: `wave` und `mode_index` gelten nur bei `mode_index` > 0, sonst liest jedes Ereignis `dir`, und `trace` setzt beim Start `wave := dir`; ein Batch, der nur `dir` setzt, ist damit gültig. `apply_event` löst beim Eintritt (`ordinary`/`extraordinary`) die Mode mit `rtt/polar/birefringence.hpp` aus dem Tangentialimpuls (`incident_tangential`, beim Austritt `mode_index` · `wave`, plus `order_momentum`) in lokalen Flächenkoordinaten; P ist die Projektion von ADR 0026, Punkt 5 (`prt_crystal_entry`/`prt_crystal_exit`), die Eigenpolarisation aus der globalen Achse. Keine reelle Lösung: `Tir` bei Ordnung 0, `Evanescent` bei einer Ordnung m ≠ 0, deren Ordnung 0 existiert. `transmit` im Kristall lässt Mode, `wave`, `dir` und P unverändert. Die Kristalldaten (n_O, n_E, Achse, Mode) füllt `event_media.hpp` für Tracer und Ray Aiming (`stop_hit`) an einer Stelle (`EventMedia::crystal_before`, `crystal_after`, `crystal_mode`). Die M1-Überladungen von `apply_event` haben keine Kristalldaten; `ordinary`/`extraordinary` geben dort `EventImpossible`. **Grenze in M4:** `first_order` lehnt Pfade in einen Kristall ab (`ParaxialError`), und damit alles, was darauf aufbaut: `make_rays` mit Ray Aiming, OPD-Karte und OPD-Fächer sowie die Komfortformen von `spot` und `path_transmission`. Solche Pfade brauchen eigene Startstrahlen (Hauptform von `path_transmission`, #122). Eine OPD im Kristall ist ein Folgepunkt (#35).

**Phasen und Beugung (ADR 0025, umgesetzt in #126 und #127):** lokale Gittergleichung n′t′_∥ = n t_∥ + m λ₀ g_∥/(2π) mit g_∥ = (I − N Nᵀ)∇φ (`rtt/geom/phase.hpp`, `apply_event`: `incident_tangential`, `order_momentum`, `order_direction`); die Phasenfunktionen legt compile in `CompiledSurface::phase_functions` ab. Ordnung 0 ist das Ereignis ohne Phasenschicht, bitgleich. Bei m ≠ 0 zuerst Tir der Ordnung 0, dann `Evanescent` für |τ| ≥ n′ (Gleichheit eingeschlossen); OPL + m φ λ₀/(2π); P = R(k₀ → k_m)·P₀ mit der kleinsten Drehung (`rotation_between`, Diebel); Effizienz je Ordnung (`diffraction_efficiency`) nur in `weight`. Ray Aiming (`stop_hit`) nutzt dieselben Ereignisdaten wie der Tracer (privat `event_media.hpp`) und geht durch Phasenflächen mit Ordnung 0 vor der Blende; `first_order` nimmt Phasenflächen mit Ordnung 0 an. Für Ordnungen ≠ 0 braucht ein Pfad eigene Startstrahlen (ADR 0025, Punkt 4; Hauptform von `path_transmission`, #122).

**Strahlteiler:** kein eigenes Modul, sondern Kombination aus Element und Interaktion (Platte, Würfel, Pellicle, PBS). Der Pfad wählt `Reflect` oder `Refract`/`Transmit`.

**Literatur** (vor Nutzung prüfen): Chipman, Lam, Young: *Polarized Light and Optical Systems* (2018); Yun, McClain, Chipman: *Three-dimensional polarization ray-tracing calculus I/II*, Applied Optics 2011; Forbes: *Shape specification for axially symmetric optical surfaces*, Optics Express 2007; Macleod: *Thin-Film Optical Filters*.

## Analyse, Optimierung, Toleranzierung

Analysen liefern Datenobjekte, niemals Plots. Optimierung und Toleranzierung arbeiten nur über Parameter-Pfade und Analyse-Operanden.

**Spot und Ray Fans (`rtt-analysis`, festgelegt in #28):** Bildfläche ist die Fläche des letzten Events im Pfad; Koordinaten x, y in mm im lokalen KS dieser Fläche (Achsen nach ihrer Pose). Bezug ist der Hauptstrahl (Pupillenmitte) bei der Referenzwellenlänge, für alle Spots und Fans, damit der Farbquerfehler sichtbar bleibt. Polychromatisch zählen die Wellenlängengewichte des Modells, normiert auf Summe 1, mal dem Leistungsgewicht des Strahls; Feldgewichte spielen erst bei feldübergreifenden Größen (Merit-Funktion) eine Rolle. Kennzahlen: Schwerpunkt, RMS und GEO jeweils um Schwerpunkt und Hauptstrahl, GEO nur über Punkte mit Gewicht > 0; Anteil vignettierter Strahlen ungewichtet (alles, was nicht mit Status `Alive` auf der Bildfläche ankommt). Ray Fans: tangential ε(p_y) bei p_x = 0, sagittal ε(p_x) bei p_y = 0, je Punkt ε_x, ε_y und Status; nicht angekommene Punkte haben ε = 0 und ihren Status. Kein angekommener Hauptstrahl → `AnalysisError`.

**OPD (`rtt-analysis`, festgelegt in #29):** Referenzsphäre um den realen Durchstoßpunkt des Hauptstrahls der Referenzwellenlänge auf der Bildfläche, Radius bis zur Mitte der paraxialen Austrittspupille; die OPL eines Strahls bis zur Sphäre ist seine OPL auf der Bildfläche minus |n′|·s (s vorzeichenbehafteter Weg von der Sphäre zur Bildfläche entlang des Strahls). Bezugsstrahl ist der Hauptstrahl derselben Wellenlänge (W(0, 0) = 0 je Wellenlänge); Piston und Tilt werden nicht entfernt, der Farbquerfehler erscheint als Tilt. Einheit: Wellen bei der Referenzwellenlänge. OPD-Karte auf einem gleichmäßigen Gitter im Einheitskreis und OPD-Fans; RMS = Standardabweichung über die angekommenen Punkte (Piston entfernt, **Tilt nicht**), PV = max − min; vignettierte Punkte haben W = 0, ihren Status und zählen nicht mit. **Austrittspupille im Unendlichen** (bild- bzw. beidseitig telezentrisch, #102): Die Referenz ist der Grenzfall R → ∞ der Sphäre um denselben Mittelpunkt C; `ReferenceSphere::radius` ist dann +∞. Für jeden Strahl gilt W∞ = OPL_Haupt − OPL_Strahl + |n′|·d·(p − C), also die OPL bis zum Fußpunkt des Lots von C auf den Strahl; das ist eine Taylorentwicklung des Schnittparameters mit der Sphäre (Herleitung in `opd.cpp`, `docs/quellen.md`). Das ist bewusst **keine Ebene senkrecht zum Hauptstrahl** (Issue #102 sprach von „Referenzebene“): Gegen eine solche Ebene erschiene schon eine perfekte, auf C konvergierende Kugelwelle als Defokus. Der Grenzwert der Sphäre um C dagegen gibt für sie W ≡ 0, denn jeder ihrer Strahlen läuft durch C (p − C ∥ d). Die OPL bis zur Sphäre wird für jedes R in der auslöschungsfreien Form OPL − |n′|·(s − R) = OPL − |n′|·(b + q/(√(q + R²) + R)) gerechnet (b = d·(p − C), q = b² − |p − C|²; der gemeinsame Term |n′|R fällt in W heraus). So ist der Übergang von endlicher zu unendlicher Austrittspupille auch numerisch stetig; mit OPL − |n′|·s bei s ≈ R ≈ 1e14 mm wäre allein die Rundung bis zu Dutzende Wellen. Seit #102 rechnet auch der endliche Fall so; gegenüber vorher ändern sich die OPD-Werte bestehender Systeme nur in der Rundung. Die Pupillenkoordinaten der OPD-Karte sind die Startkoordinaten der Strahlquelle (rtt-trace) und hängen nicht von der Austrittspupille ab.

**Feldabhängige Analysen (`rtt-analysis`, festgelegt in #31):** Feldverläufe laufen über das relative Feld 0 … 1 entlang +y bis zum größten Feldpunkt (linear im Feldwert, Grad bzw. mm). Höhen sind vorzeichenbehaftet: Abstand vom Bildflächenscheitel, projiziert auf die Feldrichtung in der globalen x-y-Ebene. *Verzeichnung* D = (h_real − h_par)/h_par · 100 % mit dem realen Hauptstrahl auf der Bildfläche und dem paraxialen Hauptstrahl **in der Ebene des Bildflächenscheitels** (nicht in der paraxialen Bildebene); D = 0 auf der Achse. Realer und paraxialer Hauptstrahl starten am selben Feldpunkt: Der Feldwert wird wie in `rtt-trace` bei der Referenzwellenlänge umgerechnet (#50), beide Strahlen laufen bei der ausgewerteten Wellenlänge. *Bildfeldwölbung* numerisch: je Feld zwei Nachbarstrahlen bei ±δ (δ = 1e−3 normierte Pupille) tangential (zum Feldpunkt hin) und sagittal; Fokus = Mitte der kürzesten Verbindung der beiden Geraden im Bildraum, angegeben als z-Komponente des Abstands vom Bildflächenscheitel (parallel zur Achse) mit dem Vorzeichen der Ausbreitung im Bildraum; Fehler O(δ²), ungerade Terme (Koma) heben sich auf. *Farblängsfehler* für ein Wellenlängenpaar (Standard erste minus letzte Systemwellenlänge): paraxial aus den Brennpunkten von `rtt-paraxial`, real über den Achsenschnitt des Zonenstrahls (0, zone) des Achsfeldes, beide längs der Ausbreitung. *Farbquerfehler*: realer Hauptstrahl jeder Wellenlänge auf der Bildfläche und sein Versatz gegen den der Referenzwellenlänge.

| Gruppe | v1 | v2 |
| --- | --- | --- |
| Erste Ordnung | Systemdaten, Pupillen, Vergrößerung, Seidel pro Fläche | parabasal |
| Geometrisch | Spot (RMS, GEO, Schwerpunkt), Ray Fans, Footprint, Verzeichnung, Bildfeldwölbung, Farbfehler, Vignettierung | Bildsimulation |
| Wellenfront | OPD gegen Referenzsphäre, Zernike-Fit, Strehl-Näherung | – |
| Beugung | FFT-PSF, FFT-MTF, Through-Focus, Encircled Energy | Huygens-PSF, Beamlets, Faserkopplung |
| Polarisation | Jones-Pupille, Diattenuations- und Retardance-Karten, Transmission je Pfad | Stokes-Detektoren |
| Pfade | Transmission je Pfad, Ghost-Ranking | Interferogramm zweier Pfade |

**Pfadauswertung (`rtt/analysis/paths.hpp`, #122):** Transmission je Pfad und OPL-Differenz zweier Pfade für dieselben Startstrahlen; Grundlage der Michelson-Abnahme, kein Interferogramm.
- **Hauptform:** Der Aufrufer gibt die Startstrahlen vor (`RayBatch`); sie werden je Pfad kopiert und verfolgt. Das geht auch für gefaltete und gekippte Pfade (Strahlteiler, Interferometer), die die paraxiale Zielung von `make_rays` nicht annimmt. Die Komfortform zielt über `make_rays` und gilt deshalb nur für rotationssymmetrische Pfade (sonst `ParaxialError`).
- **Startstrahlen:**
  - Jeder Strahl zählt als gestartet. Ist er schon zu Beginn nicht `Alive`, gilt er mit seinem Status als verloren (weight 0).
  - Jede Systemwellenlänge ist erlaubt, auch gemischt.
  - Startgewicht und Start-OPL werden mitgeführt.
- **Transmission:** `mean` ist der Mittelwert der Endgewichte über alle gestarteten Strahlen, verlorene zählen 0. Das ist der übertragene Leistungsanteil bei gleichmäßig ausgeleuchteter Pupille und unpolarisierter Quelle (ADR 0021); bei Startgewichten ≠ 1 die apodisierte Leistung, kein Verhältnis zum Start. Dazu `min`, `max` und das Gewicht je Strahl, Verluste und Warnungen nach ADR 0023.
- **OPL-Differenz** OPL_b − OPL_a je Startstrahl; die beiden Strahlen dürfen an verschiedenen Punkten der Bildfläche landen (im Michelson fallen sie zusammen):
  - **In mm, nicht in Wellen:** Es ist eine Weglänge auf der Bildfläche, keine Wellenfront; Wellen gehören zum Interferogramm.
  - **Nur für zwei Pfade mit derselben Bildfläche** (Fläche des letzten Events), sonst `invalid_argument`.
- **Abbruch und Fortschritt:** Überladungen mit `RunControl` (#83), Stufen `aim` (Komfortform) und `trace` je Pfad.
- **Referenz:** `tests/reference/m4/michelson_offset.rtt.json` mit Stop, Testarm 7,5 mm länger und beiden Ausgängen (Kamera und Rückweg durch den Stop). Je Strahl summieren sich die vier Pfade zu (R + T)² = 1, die Arme zu R·T, und die OPL-Differenz ist 2Δ = 15 mm.

**Optimierung (`rtt-optim`, M5):** Variablen sind alle `Param` mit `variable = true`. Merit-Funktion als gewichtete Summe von Operanden über beliebige Pfade und Konfigurationen; Generatoren für RMS-Spot und RMS-Wellenfront. Levenberg-Marquardt mit paralleler zentraler Differenz; Grenzen per Variablentransformation. v2: globale Suche, Glassubstitution, exakte Gradienten.

**Toleranzierung (`rtt-tolerance`, M7):** Toleranzen auf Fläche, Element und Baugruppe (Kippung um den jeweiligen Pivot), Kompensatoren, Sensitivität, inverse Sensitivität, Monte Carlo mit Seed. Eine Baugruppen-Toleranz bewegt alle Kinder starr; Doppeldurchgänge sehen dieselbe Störung, weil Pfade auf `SurfaceId` verweisen.

## Python-API und CLI

**Python (`rtt-py`, Paket `raytatouille`, ab M2):** dünne Spiegelung der C++-Typen, Ergebnisse als NumPy-Arrays, vollständige Typ-Stubs (`.pyi`). Plots in `raytatouille.plot`. Der ganze Modellbaum ist als unveränderliche typisierte Kopien lesbar (`raytatouille.model`, `System.root`, `paths`, `fields`, `aperture`, `object_space`; ADR 0024). Geändert wird mit JSON Patch auf der Bearbeitungsform (`System.to_dict`, `json_at`, `locate_surface`/`locate_node`/`locate_path`): `rt.apply_patch` als Skriptbefehl, `rt.Editor` mit Undo und Redo und einer als JSON exportierbaren, wieder abspielbaren History (ADR 0024). Pfadauswertung und Ghosts (#133): `rt.analysis.path_transmission` und `opl_difference` nehmen entweder Startstrahlen (`start=`, die Hauptform, auch für gefaltete Pfade) oder Feld, Wellenlänge, Abtastung (`rays`) und Zielung für `make_rays`; beides zugleich ist ein `ValueError`. `rt.compile_with_ghosts` und `rt.ghost_paths` liegen wie `rt.compile` auf oberster Ebene, `rt.analysis.ghost_ranking` bewertet einen `GhostSystem`. Listen in den Ergebnissen sind Spalten (NumPy-Kopien); ein fehlender Diagnosewert ist in einer Spalte NaN, ein fehlender Einzelwert None (wie bei `Prescription`, #108). Beugungsordnungen (#134) brauchen keine eigene Funktion: Ein Pfad mit `order` an einer Fläche mit Phasenschicht wird wie jeder Pfad mit `rt.trace.trace` oder `path_transmission` gerechnet; die Strahlen tragen den Status `EVANESCENT`, wenn die Ordnung keine reelle Richtung hat. Kristalle (#134): `RayBatch` trägt die Wellennormale (`wave_x/y/z`) und den Modenindex (`mode_index`, 0 ohne Mode) neben der Energierichtung `dir`; ein Batch, der nur `dir` setzt, ist gültig (Leseregel, ADR 0026). `CompiledMedium` zeigt `is_crystal`, `index_extraordinary` und die globale `optic_axis`. Die spätere GUI (`raytatouille.gui`, PySide6) nutzt ausschließlich die öffentliche API; jede GUI-Aktion ist als Skriptbefehl reproduzierbar.

```python
import raytatouille as rt

sys = rt.load("doublet.rtt.json")
spot = rt.analysis.spot(sys, path="main", field=1, rays="hexapolar:12")
pupil = rt.analysis.jones_pupil(sys, path="main", field=0, grid=64)
mf = rt.optim.MeritFunction.rms_wavefront(sys, rings=6, arms=8)
mf.add("EFL", target=100.0, weight=1.0)
rt.optim.local(sys, mf, max_iter=200)
sys.save("doublet_opt.rtt.json")
```

**Stand M2 (#32):** Paket `raytatouille` unter `libs/rtt-py` (C++ in `src/`, Paket in `python/raytatouille/`, Tests in `tests/`), Erweiterung `raytatouille._core` mit nanobind, Build mit scikit-build-core (`pyproject.toml` im Wurzelverzeichnis, `pip install .`) oder in der CMake-Baumstruktur mit `RTT_BUILD_PYTHON=ON` (Presets `ci-linux-python`, `ci-windows-python`; Werkzeuge vorher mit `pip install -r libs/rtt-py/tests/requirements.txt`, ADR 0018). `pip install .` braucht die C++-Abhängigkeiten aus vcpkg über `CMAKE_ARGS="-DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows"` (Linux: Triplet `x64-linux` oder Systempakete). Unter Windows nicht aus einer Visual-Studio-Entwicklerumgebung für x86 bauen: scikit-build-core übernimmt `Platform=x86` bzw. `VSCMD_ARG_TGT_ARCH=x86` und konfiguriert dann Win32, vcpkg baut alle Abhängigkeiten für x86 neu, und das Linken gegen die 64-Bit-Python scheitert; vorher `Platform` und `VSCMD_ARG_TGT_ARCH` auf `x64` setzen oder eine normale Shell nehmen. Vorhanden: `rt.load`, `rt.save`, `rt.validate`, `rt.System` (Name und Umgebung änderbar; ganzes Modell lesbar und mit `rt.Editor` bzw. `rt.apply_patch` änderbar, ADR 0024), `rt.MaterialLibrary` (`add_catalog`, `index`), `rt.compile`, `rt.paraxial.first_order`, `rt.trace.make_rays` und `rt.trace.trace` (Pfad als Name oder Index, Wellenlänge `None` = Referenz, `threads` ohne Einfluss auf das Ergebnis). `RayBatch`-Spalten sind beschreibbare NumPy-Views ohne Kopie; jede View hält den Batch am Leben, die Größe ist aus Python fest (kein `resize`), `prt_matrices()` ist eine dokumentierte Kopie; während `trace` ist das GIL freigegeben. C++-Fehler kommen als Klassen aus `raytatouille.errors` (`RaytatouilleError` als Basis; `ParseError` mit `pointer`, `CompileError` mit `diagnostics`, `ParaxialError`, `UnknownMaterial`, `AgfError` mit `file` und `line`) mit der C++-Meldung; Dateifehler als `OSError`, sonst ValueError/IndexError wie in nanobind. Die Typ-Stubs `_core.pyi` erzeugt `nanobind.stubgen` beim Build. Tests: Roundtrip Datei → Python → Datei bitgleich für alle Referenzdateien; Ergebnisse aus Python bitgleich zu C++ über das Testprogramm `rtt_py_reference` aus demselben CMake-Baum (1 und 4 Threads); mypy --strict auf die Tests.

**Stand M2 (#33):** `rt.analysis` mit `spot`, `ray_fan`, `opd_map`, `opd_fan`, `longitudinal_colour`, `lateral_colour`, `distortion`, `distortion_at`, `field_curvature`, `field_curvature_at` und `seidel` (dasselbe Objekt wie `rt.paraxial.seidel`). Erstes Argument ist ein `System` (bei jedem Aufruf mit `materials=` kompiliert) oder ein `CompiledSystem`; für mehrere Analysen einmal `rt.compile()` aufrufen. Wellenlänge `None` heißt Referenzwellenlänge, nur bei `spot` polychromatisch mit den Wellenlängengewichten des Modells (wie in C++). Pupillenverteilungen als Objekt oder Kurzform `"hexapolar:N"`, `"grid:N"`, `"fan_x:N"`, `"fan_y:N"`, `"random:N[:SEED]"`, `"single:PX,PY"` (ValueError mit der Kurzform bei Fehlern). Ergebnisse sind die gebundenen C++-Structs; Listen darin (Spotpunkte, Fans, OPD-Punkte, Fokuslagen, Feldverläufe, Seidel-Beiträge je Fläche) kommen als read-only NumPy-Kopien, je Größe ein Array. `AnalysisError` gehört zu `raytatouille.errors`. `raytatouille.plot` zeichnet die Ergebnisse mit matplotlib ohne eigene Auswertung; matplotlib ist optional (`raytatouille[plot]`, ADR 0018 Nachtrag). Jede Analyse ist gegen C++ bitgleich getestet (`rtt_py_reference`, 1 und 4 Threads, auch mit Vignettierung).

**Layout in Python (#81):** `rt.layout` enthält:
- `surfaces(cs)`: `SurfaceLayout` mit Id, Element, Art, `rotation`/`translation` bzw. `to_global` (4 × 4), Form (Typ, c, k, A4 …, `max_radius`), Apertur (Typ und Parameter) und `medium_front`/`medium_back`.
- `elements(cs)`.
- `sag` und `normal` mit Broadcasting über NumPy.
- `profile` und `outlines`. Die Fläche wird über Index oder Id gewählt, die Schnittebene über `"yz"` (Richtung +y), `"xz"` (Richtung +x) oder (Punkt, Normale), global.
- Beispiel `examples/python/layout.py` (Linsenschnitt des Achromaten mit Kittglied).

**Strahlpfade in Python (#80):** `rt.trace.trace(..., record_path=True, record_rays=None, max_recorded_rays=10000)` gibt `(TraceStats, RayPaths)` zurück. `RayPaths` hat schreibgeschützte NumPy-Views ohne Kopie: `position` und `direction` (N, S, 3), `opl`, `weight` und `status` (N, S), `count`, `lost_at` und `ray_indices` (N,), `event_surfaces` (S − 1,). Bitgleich gegen `rtt_py_reference` in allen Fällen von `test_bitwise.py`, die dabei mit Aufzeichnung tracen.

**Stand M3 (#62):** Coatings und Polarisation in Python.
- `rt.CoatingLibrary` mit `add_catalog` (Datei oder Verzeichnis) und `in`; `rt.compile(system, materials=None, coatings=None)`. Ohne `coatings` ist eine vergütete Fläche ein `CompileError` (ADR 0019). Fehlerhafte Kataloge ergeben `rt.CoatingCatalogError` mit `file` und `pointer`.
- `RayBatch.prt`/`prt_matrices()`/`weight` sind nach ADR 0021 dokumentiert: P leistungsnormiert, `weight` = Leistung bei unpolarisierter Quelle. `prt_matrices()` bleibt eine Kopie, weil P als neun Spalten gespeichert ist.
- `rt.polar` mit Batch-Funktionen auf einer verfolgten `RayBatch` (Schleifen in C++, `libs/rtt-py/src/polar_batch.cpp`, auch vom Referenzprogramm benutzt):
  - `initial_directions` (k₀ = Re(Pᵀk), weil P = P_T + k k₀ᵀ mit kᵀP_T = 0)
  - `transverse_polarization` (Zustand senkrecht zu k₀ projiziert und normiert)
  - `transmission` (unpolarisiert = `weight`, polarisiert = weight·|P E|²/(½‖P_T‖²))
  - `diattenuation` (bei leistungsnormiertem P die Leistungs-Diattenuation)
  - `retardance` (nur für Strahlen mit k·k₀ ≥ 1 − 1e−12, sonst NaN)
  - `stokes`
- Einzelmatrix-Funktionen: `prt_matrix`, `geometric_transform`, `diattenuation(p, k_in, k_out)`, `retardance(m, k)`, `physical_retardance(p, q, k_in)`, `stokes(e, axis, k)`.
- Eingaben werden an der API-Grenze geprüft: Form, endlich, Einheitsvektoren und Transversalität auf 1e−12; sonst `ValueError`.
- Die physikalische Retardance eines ganzen Pfads braucht das akkumulierte Q und folgt mit der Jones-Pupille (M6).
- Tests:
  - Bitgleich gegen `rtt_py_reference` mit P und den `rt.polar`-Größen je Fall. Neue Fälle: `ar_singlet` mit Coatings, absorbierende AR-Platte unter verschiedenen Winkeln, beide Michelson-Arme, Polarisator mit λ/4 und Analysator.
  - Analytische Prüfungen: Malus, (1 + ε)/4, |S3| = S0, π/2.
  - Beispiel `examples/python/polarization.py`.

**Ergebnisse als Daten (#86, ADR 0023):** Jedes Ergebnisobjekt der Python-API hat `to_dict()` und `to_json()` im Format `raytatouille-result` (Schema `schema/raytatouille-result.schema.json`, Beispiele unter `tests/reference/results/`); `rt.results.load_json` liest die Daten bitgleich zurück, ohne die C++-Objekte zu rekonstruieren. Spot, Fächer und OPD zählen in `losses` die Strahlen je Status und die Fläche mit den meisten Verlusten und melden in `warnings` die Codes `rays.lost` (über `lost_warning_fraction`, Standard 50 %) und `stop.clips_beam`; Python gibt sie zusätzlich als `RaytatouilleWarning` aus. In den rtt-py-Tests ist eine unerwartete `RaytatouilleWarning` ein Fehler.

**CLI (`rtt`):** vorhanden: `rtt validate` (Ausgabe `error [code] /pointer: Meldung`, ADR 0022), `rtt format [--check]`, `rtt --version`. Geplant: `rtt trace`, `rtt analyze`, `rtt optimize`, `rtt import <zmx> <rtt.json>`.

## Validierung und Tests

Ein Feature ist erst fertig, wenn es mindestens einen Referenztest gegen eine unabhängige Lösung besteht. Pflichtfälle (Testsysteme unter `tests/reference/`):

| Fall | Erwartung | Toleranz |
| --- | --- | --- |
| Paraboloid, Objekt im Unendlichen, auf Achse | perfekter Fokus bei R/2 | RMS-Spot < 1e−9 mm, OPD < 1e−6 λ |
| Dicke Einzellinse | EFL, BFL, Hauptebenen nach Linsenmacherformel | relativ 1e−12 (paraxial) |
| Kleine Apertur, real vs. paraxial | Konvergenz gegen paraxial | Fehler ∝ NA² |
| Sphärische Aberration, Plankonvexlinse | Seidel = Fit an reale Strahlen bei kleiner NA | relativ 1e−3 |
| Fresnel an N-BK7 (n = 1,5168) | senkrecht R ≈ 4,22 %; Brewster arctan(n) mit R_p = 0 | 1e−12 |
| Totalreflexion | Grenzwinkel arcsin(1/n), Phasensprung s/p | 1e−12 |
| λ/4-AR MgF₂ (1,38) auf 1,52 | R ≈ 1,26 % bei Designwellenlänge | 1e−10 |
| Verlustfreier Schichtstapel | R + T = 1 für s und p | 1e−12 |
| Zwei Polarisatoren | Malus: T = cos²θ | 1e−12 |
| λ/4-Platte unter 45° | linear → zirkular | 1e−12 |
| λ/2-Platte unter θ | Drehung um 2θ | 1e−12 |
| Metallspiegel unter 45° | Retardance, Diattenuation aus Fresnel mit ñ | 1e−10 |
| Calcit-Platte, Achse 45° | Walk-off nach analytischer Formel | 1e−10 rad |
| Lineares Gitter | Gittergleichung, evaneszent → Status | 1e−12 |
| Beugungsbegrenztes System | PSF-Nullstelle 1,22 λ F/#, analytische MTF | relativ 1e−3 |
| Michelson | R + T = 1, OPD-Differenz = 2 × Armdifferenz | 1e−10 |

**Abnahme M4 (#135):** `libs/rtt-analysis/tests/test_m4_acceptance.cpp` (Tag `[m4]`) prüft jedes Kriterium der Roadmap-Zeile durch den Tracer, mit Referenzsystem unter `tests/reference/m4/` und analytischem Sollwert: Michelson-Bilanz und OPL-Differenz 2Δ (`michelson_offset`), Calcit-Walk-off t·tan ρ (`calcite_walkoff`), Gittergleichung mit evaneszenter Ordnung (`grating_transmission`), Ghost-Ranking an zwei Platten und Fokuslage eines Singlet-Ghosts (`ghost_plates`, `ghost_singlet`), dazu die erste Ordnung der Feature-Tour (`m0/feature_tour`, `refract` mit Ordnung) gegen Mansuripur (7b). Python ist für diese Fälle bitgleich zu C++: `test_bitwise_paths.py` (#133, #135), `test_bitwise_diffraction.py` (#134 A) und `test_bitwise_crystal.py` (#134 B).

**Teststufen:** Unit-Tests mit Catch2 v3 in jeder Bibliothek; Property-Tests (RapidCheck, ab M1) für Energieerhaltung, |dir| = 1, P·k_ein = k_aus, Umkehrbarkeit, Invarianz bei starrer Bewegung; Golden-Tests mit Standardsystemen aus der Literatur; NSC-Tests (M9) gegen analytische Fälle (Lambert-Strahler, Etendue, Integrationskugel, Leistungsbilanz); Benchmarks mit Google Benchmark (Ziel ≥ 5 Mio. Strahl-Flächen-Schnitte pro Sekunde und Kern bei Sphären inkl. Polarisation, Regression > 10 % blockiert); Schema-Tests mit pytest.

**CI (`.github/workflows/ci.yml`):** clang-format, JSON-Schema, Linux GCC (Release), Linux Clang mit clang-tidy (Release), Linux Clang mit ASan/UBSan (Debug), Windows MSVC (Release); Warnungen sind Fehler. Job Python (ADR 0018): Paket `raytatouille` in der CMake-Baumstruktur (Presets `ci-linux-python` mit Python 3.10 und 3.13, `ci-windows-python` mit 3.13), pytest einschließlich Bitvergleich mit C++ und mypy --strict über ctest; je OS ein Smoke-Test mit `pip install .`.

## Roadmap und Abnahmekriterien

Nach M1 können bis zu drei Agenten parallel arbeiten. M4 und M6 brauchen M2 und M3; M5 braucht nur M2. Ein Meilenstein ist erreicht, wenn alle Kriterien in der CI grün sind.

| Meilenstein | Umfang | Abnahme |
| --- | --- | --- |
| **M0 Fundament** ✅ | CMake mit vcpkg, `rtt-math`, `rtt-model`, `rtt-io` (JSON mit Schema), `rtt` CLI, CI Linux und Windows | Roundtrip Datei → Modell → Datei bitgleich; CI grün |
| **M1 Sequenzieller Kern** ✅ | Ebene, Sphäre, Konik, gerade Asphäre; Schnitt, Brechung, Reflexion; `CompiledSystem`; Ray Aiming; paraxialer Trace | Paraboloid-Fokus, Linsenmacherformel, Konvergenz real → paraxial |
| **M2 Materialien, Analyse** ✅ | Dispersionsformeln, AGF-Import, Luft, dn/dT; Spot, Ray Fans, OPD, Verzeichnung, Seidel; Python-Bindings dafür | Seidel-Fall, Farblängsfehler eines Achromaten; Cooke-Triplet als Golden-Test |
| **M3 Polarisation** ✅ | P-Matrix, Fresnel mit ñ, Transfermatrix-Coatings, ideale Coatings, Polarisator, Retarder | Fresnel, Brewster, TIR, AR, Malus, λ/4, λ/2, Metallspiegel |
| **GUI-Grundlagen** ✅ (0.5.0) | G1–G5, G9, G10 aus `docs/feedback/gui-anforderungen.md`: Strahlpfade, Geometrie-Export, Modell lesen und ändern mit Undo (ADR 0024), Abbruch und Fortschritt, Prescription-Daten, Materialbibliothek, Diagnosecodes und Ergebnisformat (ADR 0022, 0023); telezentrische Systeme, paralleles Ray Aiming | Abnahmekriterien der Issues #80–#86, #93, #96, #102, #119 |
| **M4 Multi-Path** ✅ (0.6.0) | explizite Pfade, Strahlteiler, einachsige Kristalle, Gitter, Ghost-Generator | Michelson-Bilanz, Calcit-Walk-off, Gittergleichung, Ghost-Ranking |
| M5 Optimierung | Parameter-Pfade, Operanden, Merit-Generatoren, LM, Grenzen, Pickups, Konfigurationen | Singlet auf EFL + minimalen RMS-Spot; reproduzierbar |
| M6 Beugung | Zernike-Fit, FFT-PSF/MTF, Encircled Energy, Jones-Pupille, Retardance-Karten | Airy-Nullstelle, analytische MTF, Jones-Pupille eines Faltspiegels |
| M7 Toleranzierung | Toleranzen auf allen Ebenen, Kompensatoren, Sensitivität, Monte Carlo | Baugruppen-Kippung bewegt alle Kinder; Monte Carlo reproduzierbar |
| M8 Release v1.0 | Forbes-, Zernike-, XY-, Grid-Flächen, ZMX-Import, CLI komplett, Doku, Beispiele | alle Referenzfälle grün; Benchmark-Ziel; ZMX-Testdateien korrekt importiert |
| M9 Nicht-sequenziell | BVH, Quellen, Detektoren, Ray Splitting, Streumodelle, Strahlpfad-Reports | analytische NSC-Testfälle |
| M10 Python-GUI | PySide6-Anwendung auf der öffentlichen API | – |

**Erweiterungskandidaten (Entwurf, nicht verbindlich):** Aus dem Vergleich mit einem kommerziellen Optikdesign-Programm sind Ergänzungen zu M5 bis M9 und weiterführende Pakete gesammelt in [roadmap-erweiterungen.md](roadmap-erweiterungen.md). Erst nach Freigabe durch den Maintainer wandern Einträge in die Tabelle oben.
