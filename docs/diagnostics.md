# Diagnosecodes

Jede Diagnose (`rtt::model::Diagnostic`, in Python `rt.Diagnostic`) trägt einen stabilen Code in Punktnotation `<gruppe>.<was>` (ADR 0022). Dazu kommen Schwere, JSON-Pointer und Meldung. Der Code ist die maschinenlesbare Ursache; die Meldung ist für Menschen und darf sich ändern.

**Stabilitätsversprechen:** Ein Code wird nie umbenannt, wiederverwendet oder umgedeutet. Neue Codes kommen nur hinzu. Ein Code, der nicht mehr erzeugt wird, bleibt in dieser Liste und wird als „nicht mehr erzeugt“ markiert. Ausnahme: Ein Zwischencode zwischen zwei Releases (in der Bedeutung als „Zwischenstand“ markiert) darf vor dem nächsten Release wieder entfallen, weil er nie stabil veröffentlicht war (ADR 0022).

**Registry:** `libs/rtt-diagnostics/include/rtt/diagnostics/codes.hpp` (Schicht Basis, header-only), in Python `rt.diagnostics.CODES`.
- Ein Erzeuger kann nur registrierte Codes verwenden. Der Typ `DiagnosticCode` prüft das beim Kompilieren, ein Tippfehler ist ein Compile-Fehler.
- Die Schwere einer Diagnose ist immer die ihres Registry-Eintrags.
- Ein pytest gleicht diese Tabelle mit der Registry ab (Code, Schwere und Erzeuger).
- Catch2-Tests erzeugen jeden Code an seinem Ort: `libs/rtt-compile/tests/test_diagnostic_codes.cpp` die Codes von validate und compile (bei mehreren Erzeugungsorten jeden), `libs/rtt-analysis/tests/test_warnings.cpp` die der Analysen, `libs/rtt-optim/tests/test_optimize.cpp` die von `optimize` (`optim.rays_lost` und `merit.operand_unsupported` an Generatoren in `libs/rtt-optim/tests/test_generators.cpp`).

**Ausgabe:** `to_string` und `rtt validate` schreiben `error [code] /pointer: Meldung`, z. B. `error [material.unknown] /root/children/1/material: …`.

**Wer erzeugt:**
- `validate`: `rtt::model::validate`, auch `rt.validate` und `rtt validate`.
- `compile`: `rtt::compile::compile`, als `CompileError` bzw. Warnung in `CompiledSystem::diagnostics()`.
- `analysis`: Spot, Strahlfächer, OPD-Karte und OPD-Fächer (`rtt-analysis`), als Warnung im Feld `warnings` des Ergebnisses (ADR 0023).
- `edit`: `rtt::io::apply_patch` (JSON Patch auf der Bearbeitungsform, ADR 0024), als `rtt::io::EditError` mit `code()`, `location()` (Pointer in die Bearbeitungsform) und `op_index()`. Fügt ein Patch Fehler von `validate` hinzu, trägt `EditError` deren Code und Ort (z. B. `surface.id_duplicate`).
- `io`: der Leser von `rtt-io` (`parse_system`, `load_system` mit Warnliste), als Diagnose mit Schwere Warnung; die Überladungen ohne Warnliste verwerfen sie. `rtt validate` und `rtt format` schreiben sie nach stderr, Python (`rt.load`, `System.from_json`) gibt sie als `RaytatouilleWarning` aus.
- `optim`: `rtt::optim::optimize` (ADR 0030). Startfehler kommen als `rtt::optim::OptimError` mit `diagnostics()`; Meldungen des Laufs stehen in `OptimResult::diagnostics`, zusammen mit den Warnungen der Analysen am Endstand.
- `agf`: der AGF-Leser von `rtt-material` (`parse_agf`, `MaterialLibrary::add_catalog`, `add_catalog_text`), als `rtt::material::LoadWarning{code, file, line, message}` in `AgfCatalog::warnings` und `MaterialLibrary::load_warnings()` (#71); Ort ist Datei und Zeile.
- Warnungen von `compile()`, den Analysen und dem AGF-Leser gibt Python zusätzlich als `RaytatouilleWarning` mit `code` und `location` aus (beim AGF-Leser `datei:zeile`).

Pointer-Platzhalter: `…/el` steht für ein Element, z. B. `/root/children/1`; `…/s` für eine Fläche, z. B. `/root/children/1/surfaces/0`; `i`, `k` für Listenindizes.

| Code | Schwere | Bedeutung | Erzeuger | Ort |
| --- | --- | --- | --- | --- |
| `agf.duplicate_glass` | Warnung | Glas mit gleichen Werten aller gelesenen Datensätze wiederholt; die Wiederholung wird übergangen | agf | Datei:Zeile der Wiederholung |
| `agf.duplicate_glass_conflict` | Warnung | Glas zweimal mit verschiedenen Daten; es ist mehrdeutig und nicht auflösbar (`UnknownMaterial` mit beiden Zeilen und Alias-Hinweis), beide Blöcke stehen im Listing | agf | Datei:Zeile des zweiten Blocks |
| `agf.preamble_skipped` | Warnung | Textzeile vor dem ersten Datensatz (CC/NM) übersprungen, z. B. die RadiantZemax-Kopfzeile | agf | Datei:Zeile |
| `agf.stray_line` | Warnung | einzelnes Wort vor dem nächsten NM-Datensatz übersprungen, nur Leerzeilen dazwischen (NIKON-HIKARI_201911 Z. 10102) | agf | Datei:Zeile |
| `aperture.na_not_physical` | Warnung | objektseitige NA ≥ 1 in Luft | validate | `/aperture/value` |
| `aperture.stop_missing` | Fehler | Aperturtyp `stop_size` ohne Blendenelement | validate | `/aperture/type` |
| `aperture.value_invalid` | Fehler | Wert der Systemapertur nicht endlich oder ≤ 0 | validate | `/aperture/value` |
| `bounds.invalid` | Fehler | `min` nicht kleiner als `max`, oder Grenzen an einer abgeleiteten Zeile (nur über die API) | validate | `…/min` (bei einer Zeile ohne `min`: `…/max`) |
| `bounds.value_outside` | Warnung | Wert eines ungebundenen Params oder einer Zeile außerhalb seiner Grenzen (ein Startwert darf außerhalb liegen) | validate | `…/value` bzw. `/parameters/i/values/k` |
| `coating.design_wavelength_out_of_range` | Fehler | QWOT-Designwellenlänge außerhalb des Bereichs des Schichtmaterials | compile | `…/s/interaction/name` |
| `coating.layer_material_unknown` | Fehler | Material einer Schicht nicht auflösbar | compile | `…/s/interaction/name` |
| `coating.no_substrate` | Fehler | Beschichtung auf Linse, Platte oder Spiegel ohne Material | compile | `…/s/interaction` |
| `coating.not_allowed` | Fehler | Beschichtung auf dünnem Element, Blende oder Detektor (kein Substrat) | compile | `…/s/interaction` |
| `coating.substrate_ambiguous` | Fehler | Beschichtung auf einer inneren Fläche zwischen zwei Segmenten | compile | `…/s/interaction` |
| `coating.thickness_invalid` | Fehler | Schichtdicke nicht berechenbar (z. B. QWOT mit Re n ≤ 0) | compile | `…/s/interaction/name` |
| `coating.unknown` | Fehler | Beschichtungsreferenz nicht auflösbar | compile | `…/s/interaction/name` |
| `coating.wavelength_out_of_range` | Fehler | Systemwellenlänge außerhalb des Bereichs eines Schichtmaterials | compile | `…/s/interaction/name` |
| `config.unknown` | Fehler | Konfigurationsindex von compile größer oder gleich der Zahl der Konfigurationen (ADR 0029, Punkt 5) | compile | `/configurations`, ohne Abschnitt leer |
| `configurations.name_duplicate` | Fehler | Name einer Konfiguration doppelt | validate | `/configurations/k/name` (zweites Vorkommen) |
| `configurations.name_invalid` | Fehler | Name einer Konfiguration leer oder nur Leerraum | validate | `/configurations/k/name` |
| `crystal.absorbing` | Fehler | Teil eines Kristalls mit κ ≠ 0 bei einer Systemwellenlänge (in M4 nicht unterstützt, ADR 0026) | compile | `…/el/material/ordinary` bzw. `…/extraordinary` |
| `crystal.interaction_unsupported` | Fehler | andere Interaktion als `fresnel` oder `ideal_anti_reflection` an einer Kristallfläche | compile | `…/s/interaction` |
| `crystal.kind_not_allowed` | Fehler | Kristallmaterial an einem Element, das keine Linse und keine Platte ist (ADR 0026) | validate | `…/el/material` |
| `crystal.material_conflict` | Fehler | Kristall und isotropes Material zugleich gesetzt (nur über die API möglich) | validate | `…/el/material` |
| `crystal.optic_axis_invalid` | Fehler | optische Achse null oder nicht endlich | validate | `…/el/optic_axis` |
| `crystal.optic_axis_missing` | Fehler | Kristall ohne optische Achse | validate | `…/el/optic_axis` |
| `crystal.optic_axis_not_allowed` | Fehler | optische Achse an einem Element ohne Kristall | validate | `…/el/optic_axis` |
| `crystal.unsupported` | Fehler | Kristallfall, der noch nicht verfolgt werden kann (ADR 0026, Punkt 4): Kristall → Kristall, Reflexion im Kristall oder von außen (außer `ideal_anti_reflection`), Beugungsordnung im Kristall, automatischer Pfad durch einen Kristall | compile | `/paths/i/events/k` bzw. `/paths/i/events` |
| `edit.base_not_representable` | Fehler | das zu ändernde System ist in der Bearbeitungsform nicht darstellbar (nicht endliche Zahl, Material und Segmentliste oder Kristall zugleich, kein gültiges UTF-8) | edit | leer (ganzes Dokument) |
| `edit.invalid_value` | Fehler | die geänderte Bearbeitungsform ist keine gültige Systemdatei (Typ, Schlüssel, Aufzählungswert); Ort vom strengen Leser | edit | z. B. `/name`, `…/s/shape/base/conic` |
| `edit.patch_invalid` | Fehler | kein gültiges JSON-Patch-Dokument oder keine gültige Operation (RFC 6902): kein Array, unbekannte Operation, fehlendes Mitglied, ungültiger Pointer oder Array-Index (keine Ziffern, führende Null, `-` außerhalb eines Ziels), Verschieben in ein eigenes Kind, doppelte Schlüssel | edit | `path` bzw. `from` der Operation oder leer |
| `edit.path_not_found` | Fehler | Ziel, Quelle oder Eltern-Wert einer Operation fehlt, ein Array-Index liegt außerhalb, oder der Eltern-Wert ist kein Objekt bzw. Array | edit | `path` bzw. `from` der Operation |
| `edit.read_only` | Fehler | eine Operation außer `test` schreibt auf oder unter `/schema_version` oder `/units` oder auf das ganze Dokument | edit | `path` bzw. `from` der Operation |
| `edit.test_failed` | Fehler | eine `test`-Operation des Patches ist fehlgeschlagen | edit | `path` der Operation |
| `element.material_both` | Fehler | ein Material für alle Segmente und eine Segmentliste zugleich | validate | `…/el/material` |
| `element.material_empty` | Fehler | leere Materialreferenz | validate | `…/el/material`, `…/el/material/i` |
| `element.material_list_not_allowed` | Fehler | Spiegel mit Liste von Segmentmaterialien | validate | `…/el/material` |
| `element.material_missing` | Fehler | Linse oder Platte ohne Material | validate | `…/el/material` |
| `element.material_not_allowed` | Fehler | dünnes Element, Blende oder Detektor mit Material | validate | `…/el/material` |
| `element.plate_surface_not_plane` | Fehler | Fläche einer Platte ist keine Ebene | validate | `…/el/surfaces/i/shape` |
| `element.segment_count` | Fehler | Anzahl der Segmentmaterialien passt nicht zu den Flächen | validate | `…/el/material` |
| `element.surface_count` | Fehler | falsche Flächenzahl für die Elementart | validate | `…/el/surfaces` |
| `environment.medium_empty` | Fehler | leeres Umgebungsmedium | validate | `/environment/medium` |
| `environment.pressure_invalid` | Fehler | Druck nicht endlich oder < 0 atm | validate | `/environment/pressure_atm` |
| `environment.temperature_invalid` | Fehler | Temperatur nicht endlich oder nicht über dem absoluten Nullpunkt | validate | `/environment/temperature_c` |
| `fields.coordinate_invalid` | Fehler | Feldkoordinate nicht endlich | validate | `/fields/points/i` |
| `fields.empty` | Fehler | kein Feldpunkt | validate | `/fields/points` |
| `fields.weight_invalid` | Fehler | Feldgewicht nicht endlich oder < 0 | validate | `/fields/points/i/weight` |
| `interaction.axis_invalid` | Fehler | Achse von Polarisator oder Verzögerer null oder nicht endlich | validate | `…/s/interaction/transmission_axis`, `…/fast_axis` |
| `interaction.coating_name_empty` | Fehler | leere Beschichtungsreferenz | validate | `…/s/interaction/name` |
| `interaction.extinction_ratio_invalid` | Fehler | Extinktionsverhältnis des Polarisators außerhalb [0, 1] | validate | `…/s/interaction/extinction_ratio` |
| `interaction.reflectance_invalid` | Fehler | Reflexionsgrad des Strahlteilers außerhalb [0, 1] | validate | `…/s/interaction` |
| `interaction.retardance_invalid` | Fehler | Verzögerung nicht endlich | validate | `…/s/interaction/retardance_waves` |
| `io.pickup_dropped` | Warnung | `pickup` eines Params aus einer Datei vor Schema 0.4 beim Lesen verworfen, der Wert bleibt (ADR 0029, Punkt 6); der Text steht in der Meldung | io | `…/pickup` in der gelesenen Datei |
| `material.unknown` | Fehler | Materialreferenz nicht auflösbar | compile | `/environment/medium`, `…/el/material`, `…/el/material/i`, `…/el/material/ordinary`, `…/el/material/extraordinary` |
| `material.wavelength_out_of_range` | Fehler | Systemwellenlänge außerhalb des Bereichs eines Materials auf einem Pfad | compile | `/environment/medium`, `…/el/material`, `…/el/material/i` |
| `merit.index_out_of_range` | Fehler | Feld- oder Wellenlängenindex der Merit-Funktion außerhalb des Systems | validate | `…/field`, `…/wavelength`, `…/fields/k`, `…/wavelengths/k` |
| `merit.operand_unsupported` | Fehler | Operand oder Generator, den die Optimierung beim Start nicht auswerten kann: `magnification` bei Objekt im Unendlichen, ein Operand außer `param_value` oder ein Generator auf einem Pfad, den `first_order` ablehnt (Kristall, Ordnung ≠ 0, nicht rotationssymmetrisch; ADR 0030, Punkt 3), ein Generator, dessen gewählte Feld- oder Wellenlängengewichte sich zu 0 summieren (ADR 0030, Nachtrag #168) | optim | `/optimization/operands/i`, `/optimization/generators/i` |
| `merit.polychromatic_wavelength` | Fehler | polychromatischer `spot_rms` mit `wavelength` (in einer Datei ein Lesefehler) | validate | `…/wavelength` |
| `merit.sampling_invalid` | Fehler | `rings`, `arms` oder `grid` < 1 | validate | `…/rings`, `…/arms`, `…/grid` |
| `merit.selection_empty` | Fehler | leere Liste `fields` oder `wavelengths` eines Generators (in einer Datei ein Lesefehler) | validate | `…/fields`, `…/wavelengths` |
| `merit.surface_ambiguous` | Fehler | `ray_x`/`ray_y` an einer Fläche, die der Pfad mehrmals trifft, ohne `occurrence` | validate | `…/surface` |
| `merit.surface_not_on_path` | Fehler | `ray_x`/`ray_y` an einer Fläche, die der Pfad nicht trifft, oder `occurrence` zu groß | validate | `…/surface`, `…/occurrence` |
| `merit.unknown_configuration` | Fehler | Merit-Funktion nennt eine unbekannte Konfiguration | validate | `…/configuration` |
| `merit.unknown_parameter` | Fehler | `param_value` nennt eine unbekannte Zeile der Parametertabelle | validate | `…/parameter` |
| `merit.unknown_path` | Fehler | Merit-Funktion nennt einen unbekannten Pfad | validate | `…/path` |
| `merit.weight_invalid` | Fehler | Gewicht der Merit-Funktion nicht endlich oder < 0 | validate | `…/weight` |
| `node.name_duplicate` | Fehler | Name einer Baugruppe oder eines Elements doppelt | validate | `…/name` |
| `node.name_empty` | Fehler | leerer Name einer Baugruppe oder eines Elements | validate | `…/name` |
| `object.distance_invalid` | Fehler | endlicher Objektabstand nicht endlich oder ≤ 0 mm | validate | `/object/distance` |
| `optim.evaluation_failed` | Warnung | Auswertungen der Merit-Funktion sind im Lauf gescheitert (Analysefehler oder nicht definierter Wert); diese Versuche galten als abgelehnt (ADR 0030, Punkt 10) | optim | `/optimization` |
| `optim.jacobian_failed` | Fehler | eine Spalte der Jacobi-Matrix hat kein gültiges Auswertungspaar; der Lauf endet mit Status `failed` und dem letzten angenommenen Stand | optim | Pointer der Variablen (kleinster Index) |
| `optim.no_operands` | Fehler | Optimierung ohne Operanden und ohne Generatoren (ADR 0030, Nachtrag #167) | optim | leer |
| `optim.no_variables` | Fehler | Optimierung ohne variables Param und ohne variable Tabellenzeile | optim | leer |
| `optim.parameter_at_bound` | Warnung | Ergebniswert einer Variablen liegt an ihrer Grenze | optim | Pointer der Variablen |
| `optim.rays_lost` | Warnung | Strahlen eines Generators gingen am Endstand verloren; sie haben die Residuen 0 und senken die Merit-Funktion (ADR 0030, Nachtrag #168) | optim | `/optimization/generators/i` |
| `param.bound_conflict` | Fehler | gebundenes Param (`param`) mit `variable` oder Grenzen (nur über die API) | validate | Pointer des Params |
| `param.unknown_parameter` | Fehler | Param an eine Zeile gebunden, die es nicht gibt | validate | `…/param` |
| `parameters.expression_syntax` | Fehler | Ausdruck einer Zeile nicht nach der Grammatik von ADR 0029 (Zeichenposition in der Meldung), länger als 1000 Zeichen, tiefer als 64 geschachtelt oder mit einem Literal, das weder 0 noch ein endliches normales double ist (Überlauf, Unterlauf, subnormal) | validate | `/parameters/i/expression` |
| `parameters.forward_reference` | Fehler | Ausdruck nutzt die eigene oder eine spätere Zeile | validate | `/parameters/i/expression` |
| `parameters.name_duplicate` | Fehler | Name einer Zeile der Parametertabelle doppelt | validate | `/parameters/i/name` (zweites Vorkommen) |
| `parameters.name_invalid` | Fehler | Name einer Zeile nicht von der Form `[A-Za-z_][A-Za-z0-9_]*` | validate | `/parameters/i/name` |
| `parameters.not_finite` | Fehler | Ergebnis eines Ausdrucks nicht endlich, einmal je betroffener Konfiguration (Name in der Meldung); Zeilen, die eine solche Zeile nutzen, melden nichts | validate | `/parameters/i/expression` |
| `parameters.unknown_name` | Fehler | Ausdruck nutzt einen Namen, der keine Zeile ist | validate | `/parameters/i/expression` |
| `parameters.values_count` | Fehler | `values` mit einer anderen Anzahl als Konfigurationen | validate | `/parameters/i/values` |
| `parameters.variable_expression` | Fehler | Zeile mit `expression` als `variable` markiert | validate | `/parameters/i/variable` |
| `paths.crystal_mode_required` | Fehler | Eintritt in einen Kristall mit `refract` statt `ordinary` oder `extraordinary` (ADR 0026) | compile | `/paths/i/events/k` |
| `paths.empty` | Fehler | kein Pfad | validate | `/paths` |
| `paths.events_empty` | Fehler | expliziter Pfad ohne Ereignisse | validate | `/paths/i/events` |
| `paths.inner_surface_ambiguous` | Fehler | innere Fläche zwischen verschiedenen Materialien von außen betreten | compile | `/paths/i/events/k` |
| `paths.mangin_mirror_automatic` | Fehler | Spiegel mit Substrat und mehreren Flächen auf einem automatischen Pfad | compile | `…/el/surfaces` |
| `paths.mode_without_crystal` | Fehler | `ordinary` oder `extraordinary` an einem Ereignis, dessen Medium danach kein Kristall ist (Element ohne Kristall oder Austritt) | compile | `/paths/i/events/k` |
| `paths.name_duplicate` | Fehler | Pfadname doppelt | validate | `/paths/i/name` |
| `paths.name_empty` | Fehler | leerer Pfadname | validate | `/paths/i/name` |
| `paths.order_not_allowed` | Fehler | Beugungsordnung an einem Ereignis, das nicht beugt: Fläche ohne Phasenschicht (ADR 0025) | validate | `/paths/i/events/k/order` |
| `paths.uniform_plate_automatic` | Fehler | Platte aus einem Material mit mehr als 2 Flächen auf einem automatischen Pfad | compile | `…/el/surfaces` |
| `paths.unknown_surface` | Fehler | Ereignis an einer unbekannten Flächen-ID | validate | `/paths/i/events/k/surface` |
| `phase.lines_per_mm_invalid` | Fehler | Liniendichte des Gitters nicht endlich oder ≤ 0 | validate | `…/s/phases/i/lines_per_mm` |
| `phase.radius_invalid` | Fehler | Normierungsradius der Phase nicht endlich oder ≤ 0 mm | validate | `…/s/phases/i/normalization_radius` |
| `pose.no_preceding` | Fehler | `relative_to_preceding`, aber in Baumreihenfolge steht keine Fläche davor | validate | `…/pose/reference` |
| `pose.no_sibling` | Fehler | `relative_to_sibling` am ersten Kind oder an der Wurzel | validate | `…/pose/reference` |
| `pose.relative_first_surface` | Fehler | erste Fläche eines Elements relativ platziert (beide Arten; einzige Diagnose an diesem Ort) | validate | `…/s/pose/reference` |
| `rays.lost` | Warnung | mehr Strahlen verloren als die Schwelle `lost_warning_fraction` der Analyse (Standard 50 %; Vignettierung am Feldrand ist gewollt) | analysis | Fläche, an der die meisten verlorenen Strahlen enden (`…/s`), sonst leer |
| `report.paraxial_unavailable` | Warnung | Systemdaten-Report ohne paraxiale Daten, weil `prescription` den Pfad ablehnt (`ParaxialError`, dessen Meldung der Grund ist: z. B. nicht rotationssymmetrisch, Kristall, Ordnung ≠ 0, nicht kreisförmige Stop-Apertur, Feld ohne Hauptstrahl); die übrigen Angaben bleiben | analysis | `/paths/i` |
| `shape.asphere_without_coefficients` | Warnung | gerade Asphäre ohne Koeffizienten | validate | `…/s/shape/base/coefficients` |
| `shape.radius_invalid` | Fehler | Radius null oder nicht endlich | validate | `…/s/shape/base/radius` |
| `shape.zernike_radius_invalid` | Fehler | Zernike-Normierungsradius nicht endlich oder ≤ 0 mm | validate | `…/s/shape/terms/i/normalization_radius` |
| `shape.zernike_unsupported` | Fehler | Zernike-Pfeilhöhenterme noch nicht unterstützt (M8) | compile | `…/s/shape/terms/i` |
| `shape.zernike_without_coefficients` | Warnung | Zernike-Term ohne Koeffizienten | validate | `…/s/shape/terms/i/coefficients` |
| `stop.aperture_missing` | Fehler | Blende ohne Apertur | validate | `…/el/surfaces/0/aperture` |
| `stop.clips_beam` | Warnung | Strahlen der Abtastung enden an der Blende mit VIGNETTED: die Blendenöffnung beschneidet das Bündel, das die Systemapertur festlegt | analysis | Blendenfläche `…/s` |
| `stop.multiple` | Fehler | mehr als ein Blendenelement | validate | `…/el` (zweite Blende) |
| `stop.not_on_path` | Warnung | ein Pfad besucht das Blendenelement nicht; Zielen, Pupillen, Seidel-Summen und OPD werfen auf ihm `NoStopError` | compile | `/paths/i` |
| `surface.efficiency_invalid` | Fehler | Beugungseffizienz ungültig: leere Liste, Wert außerhalb [0, 1], Ordnung doppelt oder Fläche ohne Phasenschicht (ADR 0025) | validate | `…/s/diffraction_efficiency`, `…/k/efficiency`, `…/k/order` |
| `surface.id_duplicate` | Fehler | Flächen-ID doppelt | validate | `…/s/id` |
| `surface.id_empty` | Fehler | leere Flächen-ID | validate | `…/s/id` |
| `surface_aperture.half_width_invalid` | Fehler | halbe Breite der Rechteckapertur nicht endlich oder ≤ 0 mm | validate | `…/s/aperture` |
| `surface_aperture.inner_radius_invalid` | Fehler | Innenradius nicht endlich, < 0 oder ≥ Radius | validate | `…/s/aperture/inner_radius` |
| `surface_aperture.radius_invalid` | Fehler | Radius der Kreisapertur nicht endlich oder ≤ 0 mm | validate | `…/s/aperture/radius` |
| `surface_aperture.semi_axis_invalid` | Fehler | Halbachse der elliptischen Apertur nicht endlich oder ≤ 0 mm | validate | `…/s/aperture` |
| `value.not_finite` | Fehler | Zahl des Modells ohne eigene Prüfung nicht endlich (Param-Werte, Grenzen, Pivot, `orientation_deg`, Werte der Parametertabelle, `target`, `px` und `py` der Merit-Funktion; nur über die API, JSON kennt keine nicht endlichen Zahlen) | validate | Pointer der Zahl, bei einem Param `…/value` |
| `wavelengths.empty` | Fehler | keine Wellenlänge | validate | `/wavelengths` |
| `wavelengths.reference_count` | Fehler | nicht genau eine Referenzwellenlänge | validate | `/wavelengths` |
| `wavelengths.too_many` | Fehler | mehr als 65535 Wellenlängen | compile | `/wavelengths` |
| `wavelengths.value_invalid` | Fehler | Wellenlänge nicht endlich oder ≤ 0 µm | validate | `/wavelengths/i/um` |
| `wavelengths.weight_invalid` | Fehler | Wellenlängengewicht nicht endlich oder < 0 | validate | `/wavelengths/i/weight` |

## Ausnahmen mit Ortsangabe

Nicht jede Ursache ist eine Diagnose. Diese Ausnahmen tragen den Ort als Daten (ADR 0022):

| Ausnahme | Daten | Bedeutung |
| --- | --- | --- |
| `rtt::compile::NoStopError` (`std::invalid_argument`); Python `NoStopError(ParaxialError, AnalysisError, ValueError)` | `path_name`, `location` (`/paths/i`) | Der Pfad hat keine Blende, gebraucht für Zielen, Pupillen, Seidel-Summen und OPD. Eine Prüfstelle: `rtt::compile::require_stop`, vor dem ersten Strahl; eigene Argumentprüfungen der Analyse (Pfad, Feld, Wellenlänge, Rotationssymmetrie, Feldpunkt außerhalb der Achse) dürfen vorher melden. |
| `rtt::paraxial::ParaxialError` | `surface`, `location` (Fläche oder z. B. `/fields/points/i`, `/fields/type`) | Paraxiale Daten nicht definiert (nicht rotationssymmetrisch, Blende nicht kreisförmig, Feldangaben). |
| `rtt::analysis::AnalysisError` | `surface`, `location`, `ray_status`, `field`, `wavelength` | Analyse ohne Ergebnis. Die Daten beschreiben den verlorenen Haupt- oder Zonenstrahl. `surface` ist die letzte erreichte Fläche: bei VIGNETTED, ABSORBED, TIR, EVENT_IMPOSSIBLE, EVANESCENT die Fläche, an der er endet; bei MISSED, NO_CONVERGENCE die davor. |
| `rtt::io::ParseError`, `rtt::material::AgfError`, `rtt::coating::CoatingCatalogError` | Pointer bzw. Datei und Zeile | Lesefehler; ohne Code. |

## Gruppe `agf.*`

Lade-Warnungen der Glaskataloge (#71) tragen einen eigenen Typ `rtt::material::LoadWarning{code, file, line, message}`, weil `rtt-material` unter `rtt-model` liegt; die Codes stehen in derselben Registry (Erzeuger `agf`). Die Regeln (Nachtrag ADR 0008) und Belege stehen in `docs/quellen.md`.
