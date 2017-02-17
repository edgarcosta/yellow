# Copyright 2017 Edgar Costa
# See LICENSE file for license details.

from sage.all import ZZ, magma, var

__magma_DiscriminantSurfaceP3_init__ = False

def DiscriminantSurfaceP3_init():
    global __magma_DiscriminantSurfaceP3_init__
    if not __magma_DiscriminantSurfaceP3_init__:
        magma.eval("""
R4<sage0, sage1, sage2, sage3> := PolynomialRing(IntegerRing(), 4);
DiscriminantSurfaceP3 := function(f)
    P<x, y, z, w> := PolynomialRing(IntegerRing(), 4);
    D := 1;
    for j := 1 to 4 do
        B := ideal<P | [ Derivative(P!f, i) : i in [1..4] ] cat [P.j - 1]>;
        G := GroebnerBasis(B);
        D := LCM(G[#G],D);
    end for;
    return D;
end function;
""");
        __magma_DiscriminantSurfaceP3_init__ = True;

def DiscriminantSurfaceP3(F):
    """
        INPUT:
            - ``F`` -- a polynomial in 4 variables
        OUTPUT:
            an integer D such that F mod p is smooth <=> p \ndivides D
        EXAMPLES:
            sage: R4.<x,y,z,w> = ZZ[];
            sage: F = x^4 - 3*x^2*y^2 - 3*x^2*z^2 - 3*x^2*w^2 - 2*x*y*z*w + y^4 - 3*y^2*z^2 - 3*y^2*w^2 + z^4 - 3*z^2*w^2 + w^4;
            sage: DiscriminantSurfaceP3(F)
                188160
    """
    assert len(F.variables()) == 4;
    assert F.base_ring() == ZZ;
    DiscriminantSurfaceP3_init();
    x, y, z, w= var('sage0 sage1 sage2 sage3');
    D = ZZ( magma.eval( "DiscriminantSurfaceP3(%s);" %  F(x, y, z, w) ) );
    return D;

    
