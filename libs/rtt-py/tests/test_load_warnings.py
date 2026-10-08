"""Load warnings of the glass catalogues in Python (ADR 0022, #71)."""

from __future__ import annotations

import warnings

import pytest

import raytatouille as rt

HEADER = "Reproduced here by permission of RadiantZemax (www.radiantzemax.com)."


def block(name: str, vd: str = "64.17") -> str:
    return (f"NM {name} 2 517642 1.5168 {vd} 0 1 1\n"
            "CD 1.03961212 0.00600069867 0.231792344 0.0200179144 1.01046945 103.560653\n"
            "LD 0.3 2.5\n")


def test_preamble_and_stray_line_are_python_warnings() -> None:
    lib = rt.MaterialLibrary()
    with pytest.warns(rt.RaytatouilleWarning) as caught:
        lib.add_catalog_text(HEADER + "\nCC excerpt\n" + block("A") + "E\n" + block("E-LAKH1"),
                             "N", "n.agf")
    assert [(w.message.code, w.message.location) for w in caught  # type: ignore[union-attr]
            ] == [("agf.preamble_skipped", "n.agf:1"), ("agf.stray_line", "n.agf:6")]
    assert all(w.filename == __file__ for w in caught)  # the line of the call
    data = lib.load_warnings
    assert [(w.code, w.file, w.line) for w in data] == [("agf.preamble_skipped", "n.agf", 1),
                                                       ("agf.stray_line", "n.agf", 6)]
    assert data[0].to_dict() == {"code": "agf.preamble_skipped", "file": "n.agf", "line": 1,
                                 "message": data[0].message}
    assert rt.diagnostics.CODES["agf.stray_line"].producer == "agf"
    assert lib.index("N:E-LAKH1", 0.5876) == lib.index("N:A", 0.5876)


def test_duplicates() -> None:
    lib = rt.MaterialLibrary()
    text = ("CC excerpt\n" + block("E-BAK1") + block("E-BAK1") + block("E-F2", "36.258938")
            + block("E-F2", "36.258932"))
    with pytest.warns(rt.RaytatouilleWarning) as caught:
        lib.add_catalog_text(text, "N", "n.agf")
    assert [w.message.code for w in caught] == [  # type: ignore[union-attr]
        "agf.duplicate_glass", "agf.duplicate_glass_conflict"]
    # Identical blocks: one glass. Different blocks: both listed, ambiguous, not resolvable.
    listing = lib.glasses("N")
    assert [g.name for g in listing] == ["E-BAK1", "E-F2", "E-F2"]
    assert listing[0].supported
    assert not listing[1].supported and not listing[2].supported
    assert listing[1].unsupported_reason is not None
    assert "ambiguous" in listing[1].unsupported_reason
    with pytest.raises(rt.UnknownMaterial, match="ambiguous") as info:
        lib.index("N:E-F2", 0.5876)
    assert "lines 8 and 11" in str(info.value) and "alias" in str(info.value)
    # The intended block, loaded on its own under an alias, can be used.
    lib.add_catalog_text("CC E-F2 of n.agf line 11\n" + block("E-F2", "36.258932"), "N_E_F2")
    assert lib.index("N_E_F2:E-F2", 0.5876).real > 1.5


def test_a_failed_catalogue_warns_about_nothing() -> None:
    lib = rt.MaterialLibrary()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(rt.AgfError):
            lib.add_catalog_text(HEADER + "\nCC c\n" + block("A") + "XX 1\n", "Z", "z.agf")
    assert lib.load_warnings == []
