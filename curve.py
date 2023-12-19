from py_ecc.fields.field_elements import FQ as Field
import py_ecc.bn128 as b
from typing import NewType

primitive_root = 5
G1Point = NewType("G1Point", tuple[b.FQ, b.FQ])
G2Point = NewType("G2Point", tuple[b.FQ2, b.FQ2])


class Scalar(Field):
    field_modulus = b.curve_order

    # Gets the first root of unity of a given group order
    @classmethod
    def root_of_unity(cls, group_order: int):
        return Scalar(primitive_root) ** ((cls.field_modulus - 1) // group_order)

    # Gets the full list of roots of unity of a given group order
    @classmethod
    def roots_of_unity(cls, group_order: int):
        o = [Scalar(1), cls.root_of_unity(group_order)]
        while len(o) < group_order:
            o.append(o[-1] * o[1])
        return o


Base = NewType("Base", b.FQ)


def ec_mul(pt, coeff):
    if hasattr(coeff, "n"):
        coeff = coeff.n
    return b.multiply(pt, coeff % b.curve_order)


def ec_lincomb(pairs):
    o = b.Z1
    for pt, coeff in pairs:
        o = b.add(o, ec_mul(pt, coeff))
    return o
