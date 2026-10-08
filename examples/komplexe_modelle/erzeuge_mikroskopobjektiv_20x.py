"""Erzeugt mikroskopobjektiv_20x.rtt.json, nachgebaut nach einem Schnittbild eines Mikroskopobjektivs.

Schnittbild: R. Chandler, "Die Anatomie einer Objektivlinse", Olympus LS (heute Evident), 2020,
https://evidentscientific.com/de/insights/the-anatomy-of-an-objective-lens (nur Maße abgelesen).

Unendlich-korrigiertes Objektiv 20x/0,45, f = 10 mm, Arbeitsabstand ca. 0,51 mm, Objektfeld
Ø 1,1 mm (Sehfeldzahl 22). Gerechnet wird RÜCKWÄRTS: Das parallele Bündel kommt von der
Tubusseite (rechts im Bild) und fokussiert auf das Präparat (links im Bild). So liegt das Objekt
im Unendlichen und die Blende in der hinteren Brennebene kann direkt angezielt werden.

Die Gruppen G1 (Frontlinse) bis G8 (Hinterlinse) und die Flächen-ids sind wie im Bild
von links nach rechts nummeriert. Maßstab des Bildes: 0,1 mm/px. Abstände und Dicken stammen
aus dem Bild, die Radien sind auf Abbildungsqualität optimiert (kein Herstellerdatensatz).
"""
import json
from pathlib import Path

NK = "NIKON-HIKARI:"
F_OBJ, NA, OBJ_FELD = 10.0, 0.45, 0.55
ARBEITSABSTAND = 0.5095

# Vorwärts-Prescription (Objekt links, wie im Bild):
# (Gruppe, [(id, Radius | None = plan, Dicke danach)], [Material je Segment], Halbdurchmesser, Luft danach)
GRUPPEN = [
    ("G1", [("G1.S1", -2.885, 3.2), ("G1.S2", -5.455, None)], [NK + "E-LAKH1"], [1.2, 2.7], 0.3),
    ("G2", [("G2.S1", -72.42, 2.3), ("G2.S2", -6.735, None)], [NK + "NICF-V"], 3.4, 0.6),
    ("G3", [("G3.S1", 13.173, 4.0), ("G3.S2", 103.777, 1.2), ("G3.S3", 27.842, 3.2),
            ("G3.S4", -13.665, None)], [NK + "NICF-V", NK + "J-SFH1", NK + "NICF-V"], 5.0, 0.5),
    ("G4", [("G4.S1", None, 2.5), ("G4.S2", -15.97, None)], [NK + "J-FK5"], 4.6, 0.7),
    ("G5", [("G5.S1", -42.889, 2.3), ("G5.S2", 11.805, None)], [NK + "E-KZFH1"], 4.5, 0.3),
    ("G6", [("G6.S1", 12.392, 3.6), ("G6.S2", -6.596, 1.6), ("G6.S3", -13.452, None)],
     [NK + "J-FK5", NK + "J-SFH1"], 4.6, 0.3),
    ("G7", [("G7.S1", -14.225, 1.2), ("G7.S2", -7.563, 2.6), ("G7.S3", -16.94, None)],
     [NK + "J-SFH1", NK + "E-LAKH1"], 4.9, 0.3),
    ("G8", [("G8.S1", -63.261, 3.0), ("G8.S2", -41.623, None)], [NK + "E-LAKH1"], 5.0, None),
]
BLENDE_VOR_G8 = 0.23     # Blende im Luftspalt G7/G8 (hintere Brennebene), Abstand zum Scheitel G8.S1
BLENDE_RADIUS = 4.7


def vorwaerts_z():
    """Scheitel-z jeder Fläche in Vorwärtsrichtung, G1.S1 bei z = 0."""
    z, out = 0.0, []
    for _, flaechen, _, semi, luft in GRUPPEN:
        zs = []
        for sid, radius, dicke in flaechen:
            zs.append(z)
            if dicke:
                z += dicke
        out.append(zs)
        z += luft or 0.0
    return out


def rueckwaerts():
    zf = vorwaerts_z()
    z0 = zf[-1][-1]                        # G8.S2 wird z = 0 der Rückwärtsrechnung
    kinder = []
    for (name, flaechen, mat, semi, _), zs in reversed(list(zip(GRUPPEN, zf))):
        z_el = z0 - zs[-1]
        surfs = []
        for k in range(len(flaechen) - 1, -1, -1):
            sid, radius, _ = flaechen[k]
            f = {"id": sid}
            dz = zs[-1] - zs[k]
            if dz:
                f["pose"] = {"position": [0, 0, round(dz, 6)]}
            if radius is not None:
                f["shape"] = {"base": {"type": "conic", "radius": -radius}}
            f["aperture"] = {"type": "circular", "radius": semi[k] if isinstance(semi, list) else semi}
            surfs.append(f)
        m = list(reversed(mat))
        kinder.append({"type": "lens", "name": name, "pose": {"position": [0, 0, round(z_el, 6)]},
                       "material": m if len(m) > 1 else m[0], "surfaces": surfs})
        if name == "G8":
            kinder.append({"type": "stop", "name": "Blende (hintere Brennebene)",
                           "pose": {"position": [0, 0, round(z0 - zs[0] + BLENDE_VOR_G8, 6)]},
                           "surfaces": [{"id": "STO", "aperture": {"type": "circular", "radius": BLENDE_RADIUS}}]})
    kinder.append({"type": "detector", "name": "Präparat",
                   "pose": {"position": [0, 0, round(z0 + ARBEITSABSTAND, 6)]},
                   "surfaces": [{"id": "OBJ", "aperture": {"type": "circular", "radius": 1.0}}]})
    return kinder


system = {
    "schema_version": "0.3.0",
    "name": "Mikroskopobjektiv 20x/0,45 nach Schnittbild (rückwärts: Tubusseite → Präparat)",
    "units": {"length": "mm", "wavelength": "um"},
    "object": {"at_infinity": True},
    "wavelengths": [{"um": 0.4861, "weight": 1.0}, {"um": 0.5876, "weight": 1.0, "reference": True},
                    {"um": 0.6563, "weight": 1.0}],
    "aperture": {"type": "epd", "value": round(2 * F_OBJ * NA, 4)},
    "fields": {"type": "paraxial_image_height",
               "points": [{}, {"y": round(0.7 * OBJ_FELD, 4)}, {"y": OBJ_FELD}]},
    "root": {"type": "assembly", "name": "Objektiv 20x/0.45", "children": rueckwaerts()},
    "paths": [{"name": "main", "events": "auto"}],
}

if __name__ == "__main__":
    ziel = Path(__file__).with_name("mikroskopobjektiv_20x.rtt.json")
    with open(ziel, "w", encoding="utf-8") as f:
        json.dump(system, f, indent=2, ensure_ascii=False)
    try:  # kanonische Form wie `rtt format`, falls das Paket installiert ist
        import raytatouille as rt
        rt.save(rt.load(ziel), ziel)
    except ImportError:
        pass
    print(f"{ziel.name} geschrieben")
