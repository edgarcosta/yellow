# Copyright 2017 Edgar Costa
# See LICENSE file for license details.

from sage.all import PolynomialRing, PowerSeriesRing, ZZ
from sage.all import ceil, log, next_prime

def DirichletCoefficients(euler_factors, bound = None):
    """
    INPUT:
        - ``euler_factors`` -- a dictionary p --> Lp.list() where Lp \in 1 + ZZ[t]

    OUTPUT:
        A list of length bound + 1 -- [0, a1, a2, ... , a_bound]
            
    EXAMPLES:
        sage: E17a2_euler_factors = {2: [1, 1, 2], 3: [1, 0, 3], 5: [1, 2, 5], 7: [1, -4, 7], 11: [1, 0, 11], 13: [1, 2, 13], 17: [1, -1], 19: [1, 4, 19], 23: [1, -4, 23], 29: [1, -6, 29], 31: [1, -4, 31], 37: [1, 2, 37], 41: [1, 6, 41], 43: [1, -4, 43], 47: [1, 0, 47], 53: [1, -6, 53], 59: [1, 12, 59], 61: [1, 10, 61], 67: [1, -4, 67], 71: [1, 4, 71], 73: [1, 6, 73], 79: [1, -12, 79], 83: [1, 4, 83], 89: [1, -10, 89], 97: [1, -2, 97]}
        sage: A = DirichletCoefficients(E17a2_euler_factors, bound  = 99)
        sage: A == [0, 1, -1, 0, -1, -2, 0, 4, 3, -3, 2, 0, 0, -2, -4, 0, -1, 1, 3, -4, 2, 0, 0, 4, 0, -1, 2, 0, -4, 6, 0, 4, -5, 0, -1, -8, 3, -2, 4, 0, -6, -6, 0, 4, 0, 6, -4, 0, 0, 9, 1, 0, 2, 6, 0, 0, 12, 0, -6, -12, 0, -10, -4, -12, 7, 4, 0, 4, -1, 0, 8, -4, -9, -6, 2, 0, 4, 0, 0, 12, 2, 9, 6, -4, 0, -2, -4, 0, 0, 10, -6, -8, -4, 0, 0, 8, 0, 2, -9, 0]
            True

        sage: E33a3_euler_factors = {2: [1, -1, 2], 3: [1, 1], 5: [1, 2, 5], 7: [1, -4, 7], 11: [1, -1], 13: [1, 2, 13], 17: [1, 2, 17], 19: [1, 0, 19], 23: [1, -8, 23], 29: [1, 6, 29], 31: [1, 8, 31], 37: [1, -6, 37], 41: [1, 2, 41], 43: [1, 0, 43], 47: [1, -8, 47], 53: [1, -6, 53], 59: [1, 4, 59], 61: [1, -6, 61], 67: [1, 4, 67], 71: [1, 0, 71], 73: [1, 14, 73], 79: [1, 4, 79], 83: [1, -12, 83], 89: [1, 6, 89], 97: [1, -2, 97]}
        sage: A = DirichletCoefficients(E33a3_euler_factors, bound  = 99)
        sage: A == [0, 1, 1, -1, -1, -2, -1, 4, -3, 1, -2, 1, 1, -2, 4, 2, -1, -2, 1, 0, 2, -4, 1, 8, 3, -1, -2, -1, -4, -6, 2, -8, 5, -1, -2, -8, -1, 6, 0, 2, 6, -2, -4, 0, -1, -2, 8, 8, 1, 9, -1, 2, 2, 6, -1, -2, -12, 0, -6, -4, -2, 6, -8, 4, 7, 4, -1, -4, 2, -8, -8, 0, -3, -14, 6, 1, 0, 4, 2, -4, 2, 1, -2, 12, 4, 4, 0, 6, -3, -6, -2, -8, -8, 8, 8, 0, -5, 2, 9, 1]
            True
    """
    if bound is None:
        bound = next_prime( max(euler_factors.keys() ) ) - 1;
    degree = 0;
    for p, L in euler_factors.iteritems():
        if len(L) > degree + 1:
            degree = len(L) - 1;
    R = PolynomialRing(ZZ, degree, "a")
    S = PowerSeriesRing(R, default_prec = ceil(log(bound)/log(2))); 
    P = [1] + list(R.gens());
    recursion = (1/S(P)).list();
    
    A = [None]*(bound + 1);
    A[0] = 0;
    A[1] = 1;
    i = 1;
    # sieve through [1, bound]
    while i <= bound:
        if A[i] is None:
            #i is a prime
            if i in euler_factors:
                euler_i = euler_factors[i][1:];
            else:
                euler_i = [];
                
            euler_i += [0] * (degree - len(euler_i)); # fill with zeros if necessary            
            # deal with its powers
            r = 1;
            ipower = i # i^r
            while ipower <= bound:
                A[ipower] = recursion[r](euler_i)
                ipower *= i;
                r += 1;
            # deal with multiples of its powers by smaller numbers 
            for m in range(2, bound):
                if A[m] is not None and m%i != 0:
                    ipower = i
                    while m*ipower <= bound:
                        assert A[m*ipower] is None
                        A[m*ipower] = A[m] * A[ipower]
                        ipower *= i
        i += 1
    return A;
