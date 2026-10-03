# Raytatouille

Raytracing-Engine für Optikdesign: Linsen (sphärisch, asphärisch, Freiform), Spiegel,
Beschichtungen, Strahlteiler, Polarisatoren und Retarder; Polarisation von Anfang an.
C++20 mit Python-API, Linux und Windows, MIT-Lizenz.

**Stand:** Meilenstein M0 (Fundament) – Datenmodell, Dateiformat `.rtt.json`,
Validierung und Kommandozeile `rtt`. Noch kein Raytracing; das folgt mit M1.

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

Beispielsysteme liegen unter `tests/reference/m0/`.

## Dokumentation

- [Architektur, Konventionen, Roadmap](docs/architecture.md)
- [Architekturentscheidungen](docs/adr/README.md)
- [Regeln für Programmieragenten](AGENTS.md)
- [Dateiformat (JSON-Schema)](schema/raytatouille.schema.json)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
