# Raytatouille

Raytracing-Engine für Optikdesign: Linsen (sphärisch, asphärisch, Freiform), Spiegel,
Beschichtungen, Strahlteiler, Polarisatoren und Retarder; Polarisation von Anfang an.
C++20 mit Python-API, Linux und Windows, MIT-Lizenz.

**Stand:** Meilenstein M3 (Polarisation) erreicht. Sequenzielles Raytracing mit Ray Aiming und
paraxialer Analyse (M1); Glaskataloge (AGF), Luft und dn/dT, Spot, Ray Fans, OPD, Verzeichnung,
Seidel und Python-API (M2); Polarisations-Raytracing mit Fresnel (komplexer Index), Coatings
(Transfermatrix, Kataloge), idealen Polarisatoren und Retardern (M3). Als Nächstes: Multi-Path
(M4). Der Plan steht in [docs/architecture.md](docs/architecture.md).

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
- [Architekturentscheidungen](docs/adr/README.md)
- [Regeln für Programmieragenten](AGENTS.md)
- [Agenten einsetzen (Anleitung)](docs/agenten.md)
- [Geprüfte Literaturquellen](docs/quellen.md)
- [Dateiformat (JSON-Schema)](schema/raytatouille.schema.json)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
