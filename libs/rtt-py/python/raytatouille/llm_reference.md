# raytatouille – compact reference for language models and agents

Version 0.10.0. Read this first; then `help(obj)` for details (all docstrings are English) and
`raytatouille/_core.pyi` for exact types. Everything below works with the installed wheel alone.
The code blocks run in order as one session (a test checks that).

## What it is

Sequential raytracing engine for optical design (C++20 core, NumPy results): lenses, mirrors,
prisms, gratings, crystals, coatings, polarization, paraxial optics, aberrations, wavefront,
PSF/MTF, optimization (Levenberg-Marquardt) and tolerancing (sensitivity, Monte Carlo).
Headless, deterministic, thread-safe. `import raytatouille as rt`.

## Lifecycle (always the same)

1. **Describe** an `rt.System`: `rt.load(path)`, or `rt.System.from_json(text)` from a dict.
2. **Check**: `rt.validate(system)` → list of `Diagnostic` (`severity`, `code`, `location` =
   JSON pointer, `message`). Never raises.
3. **Compile**: `cs = rt.compile(system, materials, coatings, configuration=0)` →
   immutable `CompiledSystem` (raises `rt.CompileError` with `.diagnostics`).
4. **Analyse** with `rt.paraxial`, `rt.analysis`, `rt.trace`, `rt.polar` on `cs` (most also
   accept a `System` plus `materials=` and compile it themselves; compile once for many calls).
5. **Change** with JSON Patch on the edit form: `rt.apply_patch(system, patch)` or `rt.Editor`
   (undo/redo); results have `to_dict()` / `to_json()`.

## Conventions

| Topic | Rule |
| --- | --- |
| Units | mm, µm (vacuum wavelength), degrees in files and API (`*_deg`), °C, atm |
| Axes | right-handed, optical axis +z, y meridional |
| Rotation | intrinsic X→Y→Z, R = Rx·Ry·Rz; pose default: translate first |
| Radius sign | R > 0 if the centre of curvature lies at +z of the vertex; `plane` has no radius |
| Indices | absolute; `AIR` = Ciddor air at 20 °C, 1 atm (environment) |
| Fields | `angle_deg` (direction ∝ (tan θx, tan θy, 1)), `object_height` mm, `paraxial_image_height` mm |
| Pupil | unit circle = rim of the paraxial entrance pupil, +py meridional |
| OPD | W = OPL_ref − OPL, positive = leading, in waves at the reference wavelength |
| Ray weight | power for an unpolarized source, starts at 1 |
| Indices of lists | `field`, `wavelength` are 0-based indices; `path` is an index or a name |
| Determinism | bitwise equal for any thread count; randomness only with a seed |

## Minimal system file (schema 0.10.0)

```python
import json
import raytatouille as rt

SINGLET = {
    "schema_version": "0.10.0",
    "name": "Singlet",
    "units": {"length": "mm", "wavelength": "um"},
    "wavelengths": [{"um": 0.4861}, {"um": 0.5876, "reference": True}, {"um": 0.6563}],
    "aperture": {"type": "epd", "value": 20.0},
    "fields": {"type": "angle_deg", "points": [{}, {"y": 3.5}, {"y": 5.0}]},
    "root": {"type": "assembly", "name": "system", "children": [
        {"type": "stop", "name": "stop",
         "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": 10.0}}]},
        {"type": "lens", "name": "L1", "pose": {"position": [0.0, 0.0, 5.0]},
         "material": "CONST:1.5168",
         "surfaces": [
             {"id": "L1.S1", "shape": {"base": {"type": "conic", "radius": 51.68}},
              "aperture": {"type": "circular", "radius": 12.7}},
             {"id": "L1.S2", "pose": {"position": [0.0, 0.0, 4.0]},
              "aperture": {"type": "circular", "radius": 12.7}}]},
        {"type": "detector", "name": "image", "pose": {"position": [0.0, 0.0, 106.442]},
         "surfaces": [{"id": "IMG"}]}]},
    "paths": [{"name": "main", "events": "auto"}],
}
s = rt.System.from_json(json.dumps(SINGLET))
assert not [d for d in rt.validate(s) if d.severity == rt.Severity.ERROR]
cs = rt.compile(s)                      # CONST:, AIR, VACUUM need no catalog
```

Format essentials (`rt.schema("system")` returns the full JSON schema as a dict):
- Top level: `schema_version`, `name`, `units`, `environment` (`temperature_c`, `pressure_atm`,
  `medium`), `object` (`{"at_infinity": true}` default, or `{"at_infinity": false, "distance": d,
  "reference": {"surface": id}, "medium": m}`), `wavelengths` (exactly one `reference`),
  `aperture` (`type`: `epd` | `image_fnumber` | `object_na` | `stop_size`; optional
  `"stop": {"surface": id}`, `"reading": "real"` for `object_na`), `fields`, `configurations`,
  `parameters`, `root`, `paths`, `optimization`, `tolerances`.
- Tree: `assembly` (`children`) and elements `lens` (N surfaces, N−1 materials: string or list),
  `plate`, `mirror`, `thin_element` (one surface: polarizer, retarder, beam splitter, grating,
  ideal lens), `stop`, `detector`. Names and surface `id`s are unique. `"active": false` takes a
  node out of the computation. `pose`: `position` [x, y, z] mm, `rotation_deg`, `reference`
  (`absolute` | `relative_to_preceding` | `relative_to_sibling`), `order`.
- Surface `shape.base`: `plane`, `conic` (`radius`, `conic`), `even_asphere` (+ `coefficients`
  A4, A6, …), `odd_asphere` (`curvature`, `conic`, A1, A2, …), `axicon` (`angle_deg`), `qcon` /
  `qbfs` (`curvature`, `conic`, `normalization_radius`, `coefficients`), `biconic`, `toroid`,
  `cylinder`. `aperture`: `circular` (`radius`, `inner_radius`), `rectangular`, `elliptical`.
  `interaction`: `fresnel` (default), `ideal_mirror`, `absorber`, `ideal_beam_splitter`,
  `ideal_polarizer`, `ideal_retarder`, `ideal_lens`, `coating` (`"name": "CAT:NAME"`), … .
- Materials: `CONST:n` or `CONST:n,k`, `AIR`, `VACUUM`, `CATALOG:NAME` from a loaded AGF file
  (catalog name = file stem upper-case). No glass catalogs ship with the package.
- Numbers that may vary are Params: `51.68`, `{"value": 51.68, "variable": true, "min": 40}`
  or `{"param": "ROW"}` (row of the parameter table).
- Paths: `{"name": "main", "events": "auto"}` (all active surfaces in tree order) or explicit
  events `{"surface": id, "kind": "refract" | "reflect" | "transmit" | "ordinary" |
  "extraordinary", "order": m}`.

## Recipes

### Paraxial data, spot, wavefront

```python
fo = rt.paraxial.first_order(cs, path="main")          # FirstOrder
print(round(fo.efl, 4), round(fo.bfl, 4), fo.entrance_pupil.diameter)
spot = rt.analysis.spot(cs, field=2, rays="hexapolar:8")  # all wavelengths, weighted
print(spot.rays_arrived, spot.stats.rms_centroid)        # mm
opd = rt.analysis.opd_map(cs, field=0, grid=33)          # waves
print(round(opd.rms, 3), round(opd.pv, 3))
seidel = rt.paraxial.seidel(cs)                          # .sum.s1 … .s5, .surfaces
fan = rt.analysis.ray_fan(cs, field=1)                   # .tangential, .sagittal
pres = rt.paraxial.prescription(cs)                      # surface table, total track, F#
```

Other analyses in `rt.analysis`: `opd_fan`, `opd_points`, `distortion`, `field_curvature`,
`longitudinal_colour`, `lateral_colour`, `pupil_function`, `zernike_fit` (Noll), `psf`,
`psf_at`, `through_focus`, `mtf`, `encircled_energy`, `jones_pupil`, `stokes_pupil`,
`path_transmission`, `opl_difference`, `interferogram`, `ghost_ranking`, `irradiance_map`,
reports `system_report`, `raytrace_report`, `dimension_report`. Rays are given as a string
(`"hexapolar:N"` = 1+3N(N+1) rays, `"grid:N"`, `"fan_x:N"`, `"fan_y:N"`, `"random:N:SEED"`,
`"single:PX,PY"`) or a sampling object from `rt.trace`.

### Diffraction

```python
psf = rt.analysis.psf(cs, field=0, grid=64)              # .intensity (L, L), .strehl
mtf = rt.analysis.mtf(cs, field=0, grid=64)              # .tangential, .sagittal (MtfCurve)
ee = rt.analysis.encircled_energy(cs, [0.005, 0.01, 0.02], field=0)   # radii in mm
print(round(psf.strehl, 3), ee.energy.shape)
```

This singlet (1.4 waves RMS on axis) warns `psf.pupil_undersampled` and `psf.aliasing` (raise
`grid`) and `pupil.mapping_nonlinear` (a fast system; the pupil-to-image mapping is not linear).
Warnings say how far to trust a result; they never stop a run.

### Rays

```python
rays = rt.trace.make_rays(cs, rt.trace.HexapolarPupil(rings=6), fields=[0, 1])  # RayBatch
rt.trace.trace(cs, rays)                    # in place; columns are NumPy views
ok = rays.status == rt.trace.RayStatus.ALIVE
print(int(ok.sum()), rays.size, float(rays.pos_y[ok].max()))
stats, paths = rt.trace.trace(cs, rt.trace.make_rays(cs, rt.trace.FanYPupil(n=5)),
                              record_path=True)  # RayPaths: .position (event, ray, 3), .status
```

`RayBatch` columns: `pos_x/y/z`, `dir_x/y/z`, `wl` (index), `opl`, `weight`, `field`,
`pupil_x/y`, `status`, `last_surface`. Ray problems are status values (`MISSED`, `VIGNETTED`,
`TIR`, `NO_CONVERGENCE`, `ABSORBED`, `EVENT_IMPOSSIBLE`, `EVANESCENT`), never exceptions. Own
start rays: fill a `rt.trace.RayBatch(n)` and call `trace`; needed for tilted, decentred or
folded paths and for an axicon, where automatic aiming raises `rt.ParaxialError`.

### Read and change the model

```python
ptr = s.locate_surface("L1.S1") + "/shape/base/radius"   # a Param: {"value": …, "variable": …}
print(s.root.children[1].name, s.json_at(ptr + "/value"))
s2 = rt.apply_patch(s, [{"op": "replace", "path": ptr, "value": 50.0}])   # new System
ed = rt.Editor(s)
ed.apply([{"op": "replace", "path": ptr, "value": 52.0}])
ed.undo()
assert ed.system.to_json() == s.to_json()
rt.save(s2, "singlet.rtt.json")            # canonical text; rt.load reads it back
assert rt.load("singlet.rtt.json").to_json() == s2.to_json()
```

The edit form (`s.to_dict()`, pointers from `locate_surface`, `locate_node`, `locate_path`,
`locate_parameter`) writes every Param as an object, so values sit at `…/value`. Read-only typed
views: `s.root`, `s.paths`, `s.fields`, `s.aperture`, `s.object_space` (classes in `rt.model`).

### Materials

```python
lib = rt.MaterialLibrary()
lib.add_catalog_text(
    "NM N-BK7 2 517642.251 1.5168 64.17 0 1\n"
    "CD 1.03961212 6.00069867E-03 0.231792344 2.00179144E-02 1.01046945 103.560653 0 0\n"
    "LD 0.3 2.5\n", "SCHOTT")              # or lib.add_catalog("schott.agf")
n = lib.index("SCHOTT:N-BK7", 0.5875618)    # complex absolute index (array in -> array out)
print(round(n.real, 5))
cs_glass = rt.compile(rt.apply_patch(s, [{"op": "replace",
    "path": s.locate_node("L1") + "/material", "value": "SCHOTT:N-BK7"}]), lib)
```

### Optimization

Mark Params `variable` and add an `optimization` section (operands `efl`, `bfl`,
`image_fnumber`, `magnification`, `ray_x`, `ray_y`, `spot_rms`, `opd_rms`, `param_value`;
generators `rms_spot`, `rms_wavefront`):

```python
d = json.loads(s.to_json())
lens = d["root"]["children"][1]
lens["surfaces"][0]["shape"]["base"]["radius"] = {"value": 51.68, "variable": True}
d["root"]["children"][2]["pose"]["position"][2] = {"value": 106.442, "variable": True}
d["optimization"] = {
    "operands": [{"type": "efl", "path": "main", "target": 100.0, "weight": 100.0}],
    "generators": [{"type": "rms_spot", "path": "main"}]}
res = rt.optim.optimize(rt.System.from_json(json.dumps(d)),
                        options=rt.optim.OptimizeOptions(max_iterations=30))
print(res.status.name, [round(v.end, 3) for v in res.variables])
# res.system is the optimized System, res.patch the JSON Patch (text) from the input
```

`optimize(..., cancel=rt.CancelToken(), progress=f, on_iteration=g)`; a cancel ends with status
`CANCELLED` and the last accepted state, no exception. `rt.optim.MeritFunction` evaluates the
merit at arbitrary values.

### Tolerancing

```python
t = json.loads(s.to_json())
t["tolerances"] = {
    "items": [{"id": "R1", "surface": "L1.S1", "kind": "radius", "limit": 0.2},
              {"id": "tilt", "node": "L1", "kind": "tilt_x_deg", "limit": 0.05}],
    "compensators": [{"node": "image", "param": "/pose/position/2",
                      "travel_min": -1.0, "travel_max": 1.0}],
    "compensator_merit": {"generators": [{"type": "rms_spot", "path": "main", "fields": [0]}]},
    "criteria": [{"id": "efl", "operand": {"type": "efl", "path": "main"}, "max_change": 0.3}],
    "monte_carlo": {"seed": 1, "trials": 50}}
ts = rt.System.from_json(json.dumps(t))
sens = rt.tolerance.sensitivity(ts)        # .criteria[i].rss, .items[k] (ItemContribution)
mc = rt.tolerance.monte_carlo(ts)          # .yield_.fraction, .criteria[i].percentiles
print(round(sens.criteria[0].rss, 4), mc.yield_.good, "/", mc.yield_.total)
```

Item `kind`: `radius`, `curvature_relative`, `conic`, `thickness`, `param`, `parameter`,
`decenter_x/y`, `shift_z`, `tilt_x/y/z_deg`, `index`, `abbe`. Targets are anchors
(`surface`, `node`, `parameter`) plus an optional `param` pointer, never indices. `"class":
"coarse" | "normal" | "precise"` generates items from a table (`rt.tolerance.expand` lists them).

### Configurations, layout, results, plots

```python
lay = rt.layout.surfaces(cs)               # SurfaceLayout: id, rotation, translation, shape, …
prof = rt.layout.profile(cs, "L1.S1")      # list of (M, 3) polylines in the y-z plane, global
text = spot.to_json()                      # result format "raytatouille-result"
back = rt.results.load_json(text)          # Result(type, schema_version, data); data as to_dict()
print(back.type, len(lay), prof[0].shape)
```

Configurations: a `configurations` list plus `parameters` rows with `values` per configuration;
`rt.compile(s, configuration="tele")` (index or name). `rt.plot` (needs matplotlib) has one
function per result type, `rt.plot.spot(diagram, ax=None) -> Axes` and likewise `ray_fan`,
`opd_fan`, `opd_map`, `psf`, `mtf`, `through_focus`, `encircled_energy`, `distortion`,
`field_curvature`, `longitudinal_colour`. Layout drawings come from `rt.layout.outlines`.

## Errors and warnings

- Exceptions (all `rt.RaytatouilleError`): `ParseError` (file text, with pointer), `CompileError`
  and `OptimError` (`.diagnostics`, `.codes`), `ParaxialError`, `AnalysisError`, `NoStopError`,
  `ApertureError`, `EditError`, `AgfError`, `CoatingCatalogError`, `UnknownMaterial`.
  `Cancelled` only from analyses with a `cancel` token.
- Warnings are `rt.RaytatouilleWarning` (Python `warnings`), e.g. lost rays or extrapolated
  material data; results also keep them in `.warnings`.
- Diagnostic codes are stable strings (`material.unknown`, `shape.radius_invalid`, …); the
  registry is `rt.diagnostics.CODES` (code → severity, producer, summary).

## Limits of 0.10

No freeform surfaces (Zernike sag, XY, grid) yet; automatic aiming only for rotationally
symmetric and anamorphic paths; sequential tracing only; scalar diffraction; ZMX import and
export only in the command-line tool `rtt` (not in the wheel); only `System` pickles.

## API index (generated)

One line per public function. `…` stands for the keywords many functions share:
`materials`, `coatings`, `threads`, `cancel`, `progress` and the warning limits
`lost_warning_fraction`, `undersampled_waves`, `mapping_warning`, `mapping_limit`
(see `help(function)`). `SystemLike` = `System | CompiledSystem`.

<!-- BEGIN GENERATED: python tools/llm_reference.py --write -->
### raytatouille
```text
apply_patch(system: System, patch: Patch, *, check: Check = 'no_new_errors') -> System
apply_patch_with_inverse(system: System, patch: Patch, *, check: Check = 'no_new_errors') -> PatchResult
compile(system: System, materials: MaterialLibrary | None = None, coatings: CoatingLibrary | None = None, *, configuration: int | str = 0) -> CompiledSystem
compile_with_ghosts(system: System, base: int | str = 0, materials: MaterialLibrary | None = None, coatings: CoatingLibrary | None = None, *, max_paths: int = 10000, configuration: int | str = 0) -> GhostSystem
ghost_paths(system: CompiledSystem, base: int | str = 0, *, max_paths: int = 10000) -> list[Path]
llm_reference() -> str
load(file: str | os.PathLike) -> System
save(system: System, file: str | os.PathLike) -> None
schema(name: str = 'system') -> dict[str, Any]
validate(system: System) -> list[Diagnostic]
```
Classes: `CancelToken`, `CoatingLibrary`, `CompiledMedium`, `CompiledSystem`, `Diagnostic`, `Editor`, `Environment`, `Field`, `GhostPaths`, `GhostSystem`, `LoadWarning`, `MaterialLibrary`, `PatchResult`, `System`, `Wavelength`.

### rt.paraxial
```text
first_order(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, …) -> FirstOrder
first_order_section(system: SystemLike, section: Section, path: int | str = 0, wavelength: int | None = None, *, …) -> FirstOrder
is_anamorphic(system: CompiledSystem, path: int | str = 0) -> bool
prescription(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, …) -> Prescription
seidel(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, pair: ChromaticPair | tuple[int, int] | None = None, …) -> Seidel
trace_ray_section(system: CompiledSystem, section: Section, z0: float, h0: float, u0: float, path: int | str = 0, wavelength: int | None = None) -> list[RayAtEvent]
```
Classes: `ChromaticPair`, `FirstOrder`, `Prescription`, `PrescriptionSurfaces`, `Pupil`, `RayAtEvent`, `RayStart`, `Seidel`, `SeidelSurfaces`, `SeidelTerms`.

### rt.trace
```text
gauss_pupil_weights(gauss: GaussPupil) -> list[float]
make_rays(system: CompiledSystem, sampling: PupilSampling, *, path: int | str = 0, fields: Sequence[int] | None = None, wavelength: int | None = None, aiming: Aiming = Aiming.REAL, …) -> RayBatch
trace(system: CompiledSystem, rays: RayBatch, *, path: int | str = 0, record_path: bool = False, record_rays: npt.ArrayLike | None = None, polarization: bool = False, …) -> TraceStats | tuple[TraceStats, RayPaths]
```
Classes: `CancelToken`, `FanXPupil`, `FanYPupil`, `GaussPupil`, `GridPupil`, `HexapolarPupil`, `PointListPupil`, `RandomPupil`, `RayBatch`, `RayPaths`, `SinglePupilPoint`, `TraceStats`.

### rt.analysis
```text
dimension_report(system: SystemLike, *, …) -> DimensionReport
distortion(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, samples: int = 11, aiming: Aiming = Aiming.REAL, …) -> DistortionSweep
distortion_at(system: SystemLike, field: Field | tuple[float, float], path: int | str = 0, wavelength: int | None = None, *, aiming: Aiming = Aiming.REAL, …) -> DistortionPoint
encircled_energy(system: SystemLike | PupilFunction, radii: Sequence[float] | int | np.integer[Any], path: int | str | None = None, field: int | None = None, *, relative: bool = False, centre: EeCentre = EeCentre.REFERENCE, wavelengths: Sequence[int] | None = None, grid: int | None = None, apodization: float | None = None, focus_shift: float | None = None, aiming: Aiming | None = None, …) -> EncircledEnergy
field_curvature(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, samples: int = 11, delta: float = 0.001, aiming: Aiming = Aiming.REAL, …) -> FieldCurvatureSweep
field_curvature_at(system: SystemLike, field: Field | tuple[float, float], path: int | str = 0, wavelength: int | None = None, *, delta: float = 0.001, aiming: Aiming = Aiming.REAL, …) -> FieldCurvaturePoint
ghost_ranking(ghosts: GhostSystem, field: int = 0, wavelength: int | None = None, *, rays: str | PupilSampling = 'hexapolar:6', aiming: Aiming = Aiming.REAL, resolution_radius: float = 0.005, …) -> GhostRanking
interferogram(system: SystemLike, path_a: int | str, path_b: int | str, *, start: RayBatch | None = None, field: int | None = None, wavelength: int | None = None, rays: str | PupilSampling | None = None, aiming: Aiming | None = None, …) -> Interferogram
irradiance_map(system: SystemLike, rays: RayPaths | RayBatch, *, surface: str | None = None, path: int | str = 0, nx: int = 64, ny: int = 64, extent: tuple[float, float, float, float] | None = None, event: int | None = None, wavelength: int | None = None, …) -> IrradianceMap
jones_pupil(system: SystemLike, path: int | str = 0, *, start: RayBatch | None = None, chief: int | None = None, field: int | None = None, wavelength: int | None = None, rays: str | PupilSampling | None = None, aiming: Aiming | None = None, basis: JonesBasis = JonesBasis.DOUBLE_POLE, remove_global_phase: bool = False, field_amplitude: bool = False, …) -> JonesPupil
lateral_colour(system: SystemLike, path: int | str = 0, field: int = 0, *, aiming: Aiming = Aiming.REAL, …) -> LateralColour
longitudinal_colour(system: SystemLike, path: int | str = 0, *, pair: ChromaticPair | tuple[int, int] | None = None, zone: float = 1.0, aiming: Aiming = Aiming.REAL, …) -> LongitudinalColour
mtf(system: SystemLike | PupilFunction, path: int | str | None = None, field: int | None = None, *, wavelengths: Sequence[int] | None = None, grid: int | None = None, map: bool = False, complex_otf: bool = False, apodization: float | None = None, focus_shift: float | None = None, aiming: Aiming | None = None, …) -> Mtf
noll_from_nm(n: int, m: int, sine: bool = False) -> int
noll_to_nm(j: int) -> NollIndex
opd_fan(system: SystemLike, path: int | str = 0, field: int = 0, wavelength: int | None = None, *, points: int = 21, aiming: Aiming = Aiming.REAL, focus_shift: float = 0.0, …) -> OpdFan
opd_map(system: SystemLike, path: int | str = 0, field: int = 0, wavelength: int | None = None, *, grid: int = 33, aiming: Aiming = Aiming.REAL, focus_shift: float = 0.0, …) -> OpdMap
opd_points(system: SystemLike, path: int | str = 0, field: int = 0, wavelength: int | None = None, *, rays: str | PupilSampling = 'hexapolar:6', aiming: Aiming = Aiming.REAL, focus_shift: float = 0.0, …) -> OpdPupilPoints
opl_difference(system: SystemLike, path_a: int | str, path_b: int | str, *, start: RayBatch | None = None, field: int | None = None, wavelength: int | None = None, rays: str | PupilSampling | None = None, aiming: Aiming | None = None, …) -> PathOplDifference
path_transmission(system: SystemLike, path: int | str = 0, *, start: RayBatch | None = None, field: int | None = None, wavelength: int | None = None, rays: str | PupilSampling | None = None, aiming: Aiming | None = None, …) -> PathTransmission
psf(system: SystemLike, path: int | str = 0, field: int = 0, *, wavelengths: Sequence[int] | None = None, grid: int = 64, padding: int = 4, apodization: float = 0.0, focus_shift: float = 0.0, aiming: Aiming = Aiming.REAL, aliasing_factor: float = 3.0, …) -> Psf
psf_at(pupil: PupilFunction, points: npt.ArrayLike) -> npt.NDArray[np.float64]
psf_positions(result: Psf) -> npt.NDArray[np.float64]
pupil_cells(grid: int, obscuration: float = 0.0) -> PupilCells
pupil_function(system: SystemLike, path: int | str = 0, field: int = 0, wavelength: int | None = None, *, grid: int = 64, apodization: float = 0.0, focus_shift: float = 0.0, aiming: Aiming = Aiming.REAL, with_prt: bool = False, …) -> PupilFunction
ray_fan(system: SystemLike, path: int | str = 0, field: int = 0, wavelength: int | None = None, *, points: int = 21, aiming: Aiming = Aiming.REAL, …) -> RayFan
raytrace_report(system: SystemLike, path: int | str = 0, *, start: RayBatch, …) -> RaytraceReport
refocus(system: SystemLike, pupil: PupilFunction, focus_shift: float, *, …) -> PupilFunction
sampling(rays: str | PupilSampling) -> PupilSampling
seidel(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, pair: ChromaticPair | tuple[int, int] | None = None, …) -> Seidel
spot(system: SystemLike, path: int | str = 0, field: int = 0, wavelength: int | None = None, *, rays: str | PupilSampling = 'hexapolar:6', aiming: Aiming = Aiming.REAL, …) -> SpotDiagram
stokes_pupil(jones: JonesPupil, e1: complex, e2: complex) -> StokesPupil
system_report(system: SystemLike, path: int | str = 0, wavelength: int | None = None, *, …) -> SystemReport
zernike(j: int, rho: float, theta: float) -> float
zernike_fit(pupil: PupilFunction, *, terms: int = 37, incomplete_warning_fraction: float = 0.01) -> ZernikeFit
zernike_radial(n: int, m: int, rho: float) -> float
```
Classes: `DimensionReport`, `DimensionSegments`, `DistortionPoint`, `DistortionSweep`, `EncircledEnergy`, `FanPoints`, `FieldCurvaturePoint`, `FieldCurvatureSweep`, `Foci`, `GhostEntries`, `GhostRanking`, `Interferogram`, `InterferogramPoints`, `IrradianceLosses`, `IrradianceMap`, `JonesPupil`, `LateralColour`, `LongitudinalColour`, `Mtf`, `MtfCurve`, `NollIndex`, `OpdFan`, `OpdMap`, `OpdPoints`, `OpdPupilPoints`, `OplDifferencePoints`, `PathOplDifference`, `PathRays`, `PathTransmission`, `Point2`, `Points2`, `Psf`, `PsfPeak`, `PsfWavelengths`, `PupilCells`, `PupilFunction`, `PupilMapping`, `RayFan`, `RayLosses`, `RaytraceReport`, `RaytraceRows`, `ReferenceSphere`, `SpotDiagram`, `SpotStatistics`, `StokesPupil`, `SystemReport`, `ThroughFocus`, `ZernikeFit`.

### rt.polar
```text
diattenuation(rays: RayBatch | npt.ArrayLike, k_in: npt.ArrayLike | None = None, k_out: npt.ArrayLike | None = None) -> Diattenuations | Diattenuation
geometric_transform(k_in: npt.ArrayLike, k_out: npt.ArrayLike, normal: npt.ArrayLike, reflection: bool) -> FloatArray
initial_directions(rays: RayBatch) -> FloatArray
physical_retardance(p: npt.ArrayLike, q: npt.ArrayLike, k_in: npt.ArrayLike) -> Retardance
prt_matrix(k_in: npt.ArrayLike, k_out: npt.ArrayLike, normal: npt.ArrayLike, a_s: complex, a_p: complex) -> ComplexArray
retardance(rays: RayBatch | npt.ArrayLike, k: npt.ArrayLike | None = None) -> Retardances | Retardance
stokes(rays: RayBatch | npt.ArrayLike, polarization: npt.ArrayLike, axis: npt.ArrayLike) -> FloatArray
transmission(rays: RayBatch, polarization: npt.ArrayLike | None = None) -> FloatArray
transverse_polarization(rays: RayBatch, polarization: npt.ArrayLike) -> ComplexArray
```
Classes: `Diattenuation`, `Diattenuations`, `Retardance`, `Retardances`.

### rt.optim
```text
optimize(system: System, *, options: OptimizeOptions | None = None, on_iteration: IterationCallback | None = None, …) -> OptimResult
variables(system: System) -> list[Variable]
```
Classes: `GeneratorStats`, `GeneratorValue`, `IterationState`, `MeritEvaluation`, `MeritFunction`, `OperandValue`, `OptimIterations`, `OptimResult`, `OptimizeOptions`, `Variable`, `VariableValue`.

### rt.tolerance
```text
abbe_share(tolerance_class: ToleranceClass) -> float
class_limit(tolerance_class: ToleranceClass, kind: ToleranceKind, level: ToleranceLevel) -> float
expand(system: System, *, …) -> ToleranceItems
monte_carlo(system: System, *, trials: int | None = None, seed: int | None = None, percentiles: Sequence[float] = (0.5, 0.9, 0.95, 0.99), z: float = 1.959963984540054, keep_trials: bool = True, compensators: OptimizeOptions | None = None, …) -> MonteCarloResult
sensitivity(system: System, *, confirm: bool = True, compensators: OptimizeOptions | None = None, …) -> ToleranceSensitivity
sigma_ratio(distribution: ToleranceDistribution) -> float
```
Classes: `ClassTolerance`, `ClassTolerances`, `CompensatorValue`, `Confirmation`, `CriterionSensitivity`, `CriterionStatistics`, `CriterionValue`, `InverseSensitivity`, `ItemContribution`, `ItemOrigin`, `ItemTrials`, `MonteCarloResult`, `OverrideStatus`, `Percentile`, `ToleranceItems`, `ToleranceSensitivity`, `TrialResult`, `TrialTable`, `Yield`.

### rt.layout
```text
elements(system: CompiledSystem) -> list[CompiledElement]
node_frames(system: CompiledSystem) -> list[NodeFrame]
normal(system: CompiledSystem, surface: int | str, x: npt.ArrayLike, y: npt.ArrayLike, frame: str = 'local') -> FloatArray
outlines(system: CompiledSystem, element: int, plane: SectionPlane = 'yz', samples: int = 101) -> list[FloatArray]
profile(system: CompiledSystem, surface: int | str, plane: SectionPlane = 'yz', samples: int = 101) -> list[FloatArray]
reference_frame(system: CompiledSystem, node: str) -> FloatArray
sag(system: CompiledSystem, surface: int | str, x: npt.ArrayLike, y: npt.ArrayLike) -> FloatArray
surfaces(system: CompiledSystem) -> list[SurfaceLayout]
```
Classes: `CompiledElement`, `NodeFrame`, `SurfaceLayout`.

### rt.reports
```text
to_csv(report: Report) -> str
```

### rt.results
```text
encode_json(type_name: str, data: dict[str, Any], indent: int | None = None) -> str
load_json(text: str) -> Result
to_dict(result: Any) -> dict[str, Any]
to_json(result: Any, indent: int | None = None) -> str
```
Classes: `Result`.

### rt.plot (matplotlib)
```text
distortion(sweep: DistortionSweep, ax: Axes | None = None) -> Axes
encircled_energy(result: EncircledEnergy, ax: Axes | None = None) -> Axes
field_curvature(sweep: FieldCurvatureSweep, ax: Axes | None = None) -> Axes
longitudinal_colour(colour: LongitudinalColour, ax: Axes | None = None) -> Axes
mtf(result: Mtf, ax: Axes | None = None) -> Axes
opd_fan(fan: OpdFan, ax: Axes | None = None) -> Axes
opd_map(opd: OpdMap, ax: Axes | None = None) -> Axes
psf(result: Psf, ax: Axes | None = None, *, log: bool = True, floor: float = 1e-06) -> Axes
ray_fan(fan: RayFan, ax: Axes | None = None) -> Axes
spot(diagram: SpotDiagram, ax: Axes | None = None) -> Axes
through_focus(result: ThroughFocus, ax: Axes | None = None) -> Axes
```

### rt.System methods
```text
json_at(pointer: str) -> Any
locate_surface(id: str) -> str | None
locate_node(name: str) -> str | None
locate_path(name: str) -> str | None
locate_parameter(name: str) -> str | None
resolved(configuration: int | str = 0) -> System
from_json(text: str) -> System
save(file: str | os.PathLike) -> None
```

### rt.Editor methods
```text
apply(patch: Patch) -> None
set(pointer: str, value: Any) -> None
insert(pointer: str, value: Any) -> None
remove(pointer: str) -> None
move(from_: str, to: str) -> None
undo() -> None
redo() -> None
export_history() -> dict[str, Any]
from_history(data: Mapping[str, Any]) -> Editor
```

### rt.MaterialLibrary methods
```text
glasses(catalog: str) -> list[GlassInfo]
glass(reference: str) -> GlassInfo
glass_map(catalogs: Sequence[str] | None = None) -> GlassMap
```

Model classes (`rt.model`, read-only views of `System`): `Absorber`, `ApertureReading`, `Assembly`, `Axicon`, `Biconic`, `CircularAperture`, `Compensator`, `Configuration`, `CoatingRef`, `Conic`, `Criterion`, `CrystalMaterial`, `Cylinder`, `DiffractionEfficiency`, `Element`, `ElementKind`, `EllipticalAperture`, `EvenAsphere`, `Event`, `EventKind`, `FieldSet`, `FieldType`, `FirstOrderOperand`, `FirstOrderQuantity`, `Fresnel`, `IdealAntiReflection`, `IdealBeamSplitter`, `IdealCylinderLens`, `IdealLens`, `IdealMirror`, `IdealPolarizer`, `IdealRetarder`, `LinearGrating`, `MonteCarloDefaults`, `ObjectSpace`, `OddAsphere`, `OpdRmsOperand`, `Optimization`, `Param`, `ParamValueOperand`, `ParameterRow`, `Path`, `Plane`, `Pose`, `PoseOrder`, `PoseReference`, `QBfs`, `QCon`, `RadialPhase`, `RayCoordinate`, `RayOperand`, `RectangularAperture`, `ShapeStack`, `SpotGenerator`, `SpotReference`, `SpotRmsOperand`, `Surface`, `SurfaceReference`, `SystemAperture`, `SystemApertureType`, `ToleranceClass`, `ToleranceClassOverride`, `ToleranceDistribution`, `ToleranceItem`, `ToleranceKind`, `TolerancePivotOverride`, `ToleranceSet`, `ToleranceValueOverride`, `Toroid`, `WavefrontGenerator`, `ZernikeSag`.
<!-- END GENERATED -->
