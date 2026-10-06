"""Material library with catalogue listing (#85): glasses of AGF catalogues and glass map.

``MaterialLibrary`` extends the C++ library (``add_catalog``, ``add_catalog_text``, ``catalogs``,
``index`` for a float or a NumPy array of wavelengths) with the listing of the loaded AGF
catalogues as read from the files. Record meanings after the Ansys OpticStudio User Guide 2025
R1, "The AGF & BGF File Formats"; units of ED, MD and IT checked against the SCHOTT N-BK7 data
sheet; placeholders ("_", "-") and "-1" mean "not available" (docs/quellen.md).

Example::

    import numpy as np
    import raytatouille as rt

    lib = rt.MaterialLibrary()
    lib.add_catalog("schott.agf")
    lib.add_catalog("old/schott.agf", name="SCHOTT_OLD")   # alias for a second SCHOTT file
    for g in lib.glasses("SCHOTT"):
        print(g.reference, g.nd, g.vd, g.status)
    n = lib.index("SCHOTT:N-BK7", np.linspace(0.4, 0.7, 301))   # complex128 array
    chart = lib.glass_map()                                      # n_d over v_d
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, NamedTuple

import numpy as np
import numpy.typing as npt

from . import _core
from .errors import UnknownMaterial

__all__ = [
    "ClassRange",
    "GlassInfo",
    "GlassMap",
    "MaterialLibrary",
    "MechanicalData",
    "ThermalData",
]

STATUS_NAMES = ("standard", "preferred", "obsolete", "special", "melt")


class ClassRange(NamedTuple):
    """A resistance class of the OD record: low == high for one class, low < high for a range
    such as "1-2" (SCHOTT)."""

    low: float
    high: float


class ThermalData(NamedTuple):
    """TD record: coefficients of the SCHOTT dn/dT model (rtt/material/thermal.hpp)."""

    d0: float
    d1: float
    d2: float
    e0: float
    e1: float
    ltk_um: float
    reference_temperature_c: float


class MechanicalData(NamedTuple):
    """MD record; None where the catalogue gives no value."""

    youngs_modulus_gpa: float | None
    poisson_ratio: float | None
    knoop_hardness: float | None
    """Knoop hardness HK in kgf/mm^2."""
    specific_heat_j_per_kg_k: float | None
    thermal_conductivity_w_per_m_k: float | None


@dataclass(frozen=True, eq=False)
class GlassInfo:
    """One glass of an AGF catalogue as written in the file (read-only; compared by identity,
    since ``transmission`` is an array).

    ``nd`` and ``vd`` are N(d) and V(d) of the NM record ("for reference only" in the format);
    values computed from the dispersion formula come from ``MaterialLibrary.index``.
    ``supported`` is False for a dispersion formula whose coefficient order is not verified
    (#42); ``unsupported_reason`` then says why and ``index`` raises UnknownMaterial.
    """

    reference: str
    catalog: str
    name: str
    line: int
    formula: int
    supported: bool
    unsupported_reason: str | None
    coefficients: tuple[float, ...]
    nd: float
    vd: float
    status: int | None
    """0 standard, 1 preferred, 2 obsolete, 3 special, 4 melt (see STATUS_NAMES)."""
    exclude_substitution: bool | None
    melt_frequency: int | None
    """Relative melt frequency 1..5; None if not given (the files also write -1 and 0)."""
    comment: str
    wavelength_range_um: tuple[float, float] | None
    thermal: ThermalData | None
    tce_m30_70: float | None
    """Thermal expansion -30..70 degC in 1e-6/K (ED)."""
    tce_100_300: float | None
    """Thermal expansion 100..300 degC in 1e-6/K (ED; SCHOTT writes alpha(+20/+300 degC) of its
    data sheet here, the format description calls it "currently not used")."""
    density_g_per_cm3: float | None
    dpgf: float | None
    ignore_thermal_expansion: bool | None
    mechanical: MechanicalData | None
    relative_cost: float | None
    climate_resistance: ClassRange | None
    stain_resistance: ClassRange | None
    acid_resistance: ClassRange | None
    alkali_resistance: ClassRange | None
    phosphate_resistance: ClassRange | None
    transmission: npt.NDArray[np.float64]
    """Internal transmittance, shape (n, 3): wavelength in um, tau_i, thickness in mm (read-only)."""


class GlassMap(NamedTuple):
    """Data for a glass map (n_d over v_d): one entry per glass, from the NM records."""

    reference: list[str]
    nd: npt.NDArray[np.float64]
    vd: npt.NDArray[np.float64]
    status: npt.NDArray[np.int64]
    """Status as in GlassInfo, -1 if not given."""


def _not_available(value: float | None) -> float | None:
    # The format writes -1 for "not available" in OD; the parser keeps it as written.
    return None if value is None or value == -1.0 else value


def _class_range(value: tuple[float, float] | None) -> ClassRange | None:
    if value is None or (value[0] == -1.0 and value[1] == -1.0):
        return None
    return ClassRange(value[0], value[1])


def _read_only(a: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    a.setflags(write=False)
    return a


def _item(values: Sequence[float] | None, i: int) -> float | None:
    return values[i] if values is not None and i < len(values) else None


def _glass_info(r: dict[str, Any]) -> GlassInfo:
    thermal = r["thermal"]
    extra = r["extra"]
    mechanical = r["mechanical"]
    other = r["other"] or {}
    exclude = r["exclude_substitution"]
    ignore = _item(extra, 4)
    return GlassInfo(
        reference=f"{r['catalog']}:{r['name']}",
        catalog=r["catalog"],
        name=r["name"],
        line=r["line"],
        formula=r["formula"],
        supported=r["unsupported_reason"] is None,
        unsupported_reason=r["unsupported_reason"],
        coefficients=tuple(r["coefficients"]),
        nd=r["nd"],
        vd=r["vd"],
        status=r["status"],
        exclude_substitution=None if exclude is None else bool(exclude),
        melt_frequency=r["melt_frequency"],
        comment=r["comment"],
        wavelength_range_um=r["wavelength_range_um"],
        thermal=ThermalData(*thermal) if thermal is not None and len(thermal) == 7 else None,
        tce_m30_70=_item(extra, 0),
        tce_100_300=_item(extra, 1),
        density_g_per_cm3=_item(extra, 2),
        dpgf=_item(extra, 3),
        ignore_thermal_expansion=None if ignore is None else bool(ignore),
        mechanical=MechanicalData(*mechanical) if mechanical is not None else None,
        relative_cost=_not_available(other.get("relative_cost")),
        climate_resistance=_class_range(other.get("cr")),
        stain_resistance=_class_range(other.get("fr")),
        acid_resistance=_class_range(other.get("sr")),
        alkali_resistance=_class_range(other.get("ar")),
        phosphate_resistance=_class_range(other.get("pr")),
        transmission=_read_only(np.asarray(r["transmission"], dtype=np.float64)),
    )


class MaterialLibrary(_core.MaterialLibrary):
    """Resolves material references ("VACUUM", "AIR", "CONST:<n>", catalogue glasses) and lists
    the loaded AGF catalogues. Pass it to ``compile``."""

    def glasses(self, catalog: str) -> list[GlassInfo]:
        """All glasses of a loaded catalogue in file order, as written in the AGF file.

        Raises UnknownMaterial if no catalogue of that name is loaded.
        """
        # Raw records from the C++ binding; the stub generator leaves out private methods.
        records: list[dict[str, Any]] = getattr(self, "_glass_records")(catalog)
        return [_glass_info(r) for r in records]

    def glass(self, reference: str) -> GlassInfo:
        """The glass ``CATALOG:NAME`` of a loaded catalogue.

        Raises UnknownMaterial if the catalogue is not loaded or has no such glass.
        """
        catalog, sep, name = reference.partition(":")
        if not sep:
            raise UnknownMaterial(f"'{reference}' is not a catalogue reference CATALOG:NAME")
        for info in self.glasses(catalog):
            if info.name == name:
                return info
        raise UnknownMaterial(f"glass '{name}' not found in catalog {catalog}")

    def glass_map(self, catalogs: Sequence[str] | None = None) -> GlassMap:
        """n_d and v_d of the NM records for a glass map, over all loaded catalogues or the given
        ones (in the given order, glasses in file order).

        The values are as written in the catalogues; some manufacturers write placeholders there
        (e.g. n_d = 1, v_d = 0 for infrared glasses), which a plot should leave out.
        """
        names = self.catalogs() if catalogs is None else list(catalogs)
        infos = [g for c in names for g in self.glasses(c)]
        return GlassMap(
            reference=[g.reference for g in infos],
            nd=np.array([g.nd for g in infos], dtype=np.float64),
            vd=np.array([g.vd for g in infos], dtype=np.float64),
            status=np.array([-1 if g.status is None else g.status for g in infos],
                            dtype=np.int64),
        )

