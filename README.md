# Raytatouille

Raytracing-Engine für Optikdesign: Linsen (sphärisch, asphärisch, Freiform), Spiegel,
Beschichtungen, Strahlteiler, Polarisatoren und Retarder; Polarisation von Anfang an.
C++20 mit Python-API, Linux und Windows, MIT-Lizenz.

**Stand:** Version 0.9.0, Toleranzierung erreicht (M7): Toleranzsatz in der Systemdatei mit
Toleranzklassen, Einzelwerten und expliziten Toleranzen auf Fläche, Element und Baugruppe (Kippung
um den Pivot, starre Bewegung aller Kinder), auf Brechzahl, Abbe-Zahl und Konfigurationszeilen;
Kompensatoren mit Levenberg-Marquardt; Sensitivität mit RSS und inverser Sensitivität; Monte Carlo
mit Seed, bitgleich für jede Threadzahl, Statistik, Perzentile und Ausbeute mit Wilson-Intervall;
alles auch aus Python. Davor: Beugung (M6) mit Pupillenfunktion, Zernike-Fit (Noll) und Strehl,
FFT-PSF, MTF und Encircled Energy, Jones-Pupille mit Diattenuations-, Retardance- und Stokes-Karten,
Interferogramm zweier Pfade und die ideale Linse; Optimierung mit Levenberg-Marquardt,
Merit-Funktion und Parametertabelle mit Konfigurationen (M5; exakte Bedingungen wie eine feste
Brennweite als Bindung in der Parametertabelle, #196); explizite Pfade durch Strahlteiler, Ghosts,
Gitter und einachsige Kristalle (M4); sequenzielles Raytracing mit Ray Aiming und paraxialer Analyse
(M1); Glaskataloge (AGF), Luft und dn/dT, Spot, Ray Fans, OPD, Verzeichnung, Seidel und Python-API
(M2); Polarisations-Raytracing mit Fresnel (komplexer Index), Coatings (Transfermatrix, Kataloge),
idealen Polarisatoren und Retardern (M3). Für eine GUI: Strahlpfade und Geometrie-Export, Modell
lesen und ändern mit Undo (JSON Patch), Abbruch und Fortschritt, Prescription-Daten, durchsuchbare
Materialbibliothek, stabile Diagnosecodes und ein JSON-Ergebnisformat, telezentrische Systeme
(GUI-Grundlagen, 0.5.0). Als Nächstes: Release v1.0 (M8). Der Plan steht in
[docs/architecture.md](docs/architecture.md).

## Bauen

```bash
# Ubuntu 24.04 mit Systempaketen
sudo apt-get install cmake ninja-build g++ libeigen3-dev nlohmann-json3-dev catch2
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
```

Mit vcpkg (Linux oder Windows, wie die CI): `VCPKG_ROOT` setzen, dann Preset
`ci-linux-gcc`, `ci-linux-clang` oder `ci-windows-msvc` verwenden.

## Benutzen

```bash
rtt validate tests/reference/m0/singlet.rtt.json   # Struktur und Semantik prüfen
rtt format meine_optik.rtt.json                     # in kanonische Form bringen
```

Beispielsysteme liegen unter `tests/reference/`, mit einem Unterordner je Meilenstein (`m0/` bis
`m3/`).

## Dokumentation

- [Architektur, Konventionen, Roadmap](docs/architecture.md)
- [Roadmap-Erweiterungen](docs/roadmap-erweiterungen.md)
- [Systemdateien `.rtt.json` schreiben (Kurzanleitung)](docs/dateiformat.md)
- [Lehrbeispiele: 16 Systeme mit Messwerten](examples/lehrbeispiele/README.md)
- [Architekturentscheidungen](docs/adr/README.md)
- [Regeln für Programmieragenten](AGENTS.md)
- [Agenten einsetzen (Anleitung)](docs/agenten.md)
- [Geprüfte Literaturquellen](docs/quellen.md)
- [Dateiformat (JSON-Schema)](schema/raytatouille.schema.json)
- [Komplexe Modelle: Mikroskopobjektiv 20x/0,45](examples/komplexe_modelle/README.md)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
