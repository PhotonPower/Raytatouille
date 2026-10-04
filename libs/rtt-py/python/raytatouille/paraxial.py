"""First-order (paraxial) data of a compiled system (rtt-paraxial).

Indices are absolute: in AIR (Ciddor air, n about 1.00027) EFL = 1/power is smaller by the
factor n_air than in programs that compute relative to air.
"""

from ._core import FirstOrder, Pupil, first_order

__all__ = ["FirstOrder", "Pupil", "first_order"]
