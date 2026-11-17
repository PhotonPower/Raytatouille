# Lizenzen der eingebetteten Bestandteile

Raytatouille selbst steht unter der MIT-Lizenz (`LICENSE` im Wurzelverzeichnis). Die Wheels auf
PyPI enthalten zusätzlich übersetzten Code der folgenden Bibliotheken; ihre Lizenztexte liegen in
diesem Verzeichnis und stehen in jedem Wheel unter `raytatouille-<version>.dist-info/licenses/`
(ADR 0018, Nachtrag #404).

| Bestandteil | Version (Pin) | Lizenz | Datei | Im Wheel |
|---|---|---|---|---|
| oneTBB | 2023.1.0 (vcpkg) | Apache-2.0 | `oneTBB-Apache-2.0.txt` | Linux, macOS statisch gelinkt; Windows `tbb12.dll` |
| Eigen | 5.0.1 (vcpkg) | MPL-2.0 | `Eigen-MPL-2.0.txt` | Header, in `_core` übersetzt |
| nlohmann-json | 3.12.0 (vcpkg) | MIT | `nlohmann-json-MIT.txt` | Header, in `_core` übersetzt |
| darin Grisu2 (F. Loitsch) und UTF-8-Decoder (B. Hoehrmann) | wie nlohmann-json | MIT | `nlohmann-json-third-party-MIT.txt` | in `_core` übersetzt (`dump()`) |
| pocketfft | `0fa0ef5` (vcpkg, `cmake/RttPocketfft.cmake`) | BSD-3-Clause | `pocketfft-BSD-3-Clause.txt` | Header, in `_core` übersetzt |
| nanobind | 2.13.0 (`libs/rtt-py/tests/requirements.txt`) | BSD-3-Clause | `nanobind-BSD-3-Clause.txt` | in `_core` übersetzt |
| robin-map (in nanobind) | wie nanobind | MIT | `robin-map-MIT.txt` | in `_core` übersetzt |

**Microsoft Visual C++ Runtime (nur Windows):** Das Windows-Wheel enthält `msvcp140.dll`
(vendored von delvewheel), weil `_core` und `tbb12.dll` gegen die DLL-Laufzeit gebaut sind und
CPython nur `vcruntime140*.dll` mitbringt. Die DLL ist eine weitergebbare Datei des Visual
Studio Redistributable und wird nach dessen Lizenzbedingungen unverändert weitergegeben.
Die übrigen Systemlaufzeiten (libstdc++ aus gcc-toolset unter der GCC Runtime Library Exception,
libc++ unter Apache-2.0 WITH LLVM-exception) verlangen in Binärform keinen Hinweis.

**SPDX-Ausdruck:** `license` in `pyproject.toml` beschreibt die Wheels. Das sdist enthält
keinen Code dieser Bibliotheken (vcpkg und nanobind holen ihn beim Bau); PEP 639 kennt nur einen
Ausdruck für beide.

**Eigen (MPL-2.0, Abschnitt 3.2):** Eigen wird unverändert verwendet. Der Quelltext der
verwendeten Version steht unter <https://gitlab.com/libeigen/eigen> (Tag `5.0.1`).

Die Texte stammen aus den gepinnten Versionen. Wer einen Pin ändert, prüft, ob sich der
Lizenztext geändert hat, und passt Tabelle und Datei an.
