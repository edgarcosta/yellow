# Copyright 2017 Edgar Costa
# See LICENSE file for license details.

from sage.all import PolynomialRing, ZZ

def bad_factors_ec(E):
    """
    INPUT:
        ``E`` -- an elliptic cuver over a number field

    OUTPUT:
        A list of pairs (f, bad_factor) where bad_factor are the coefficients of the associated to the prime f

    EXAMPLEs:
        sage: bad_factors_ec(EllipticCurve("35a2"))
            [(5, [1, 1]), (7, [1, -1])]

        sage: K.<a> = NumberField(x^2 - x - 7)
        sage: bad_factors_ec( EllipticCurve(K, [a + 1, a - 1, 1, 275*a - 886, -4395*a + 13961]))
            [(Fractional ideal (2), [1, 1]), (Fractional ideal (a - 1), [1, 1]), (Fractional ideal (a + 5), [1, -1])]

        sage: K.<phi> = NumberField(x^2 - x - 1)
        sage: bad_factors_ec(EllipticCurve(K, [1, phi + 1, phi, phi, 0]))
            [(Fractional ideal (5*phi - 2), [1, 1])]
    """
    ZZT = PolynomialRing(ZZ, "T")
    T = ZZT.gen()
    F = E.conductor().factor()
    output = [None] * len(F);
    for i, (f, e) in enumerate(F):
        output[i] = (f, list(1 - E.local_data(f).bad_reduction_type() * T))

    return output
