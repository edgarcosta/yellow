# Copyright 2017 Edgar Costa
# See LICENSE file for license details.

from sage.all import Matrix, identity_matrix, vector

# Outputs a matrix with a given characteristic polynomial
def MatrixCharpoly(coefficients):
    """
    INPUT:
        - ``coefficients`` -- an array with n + 1, where the first entry is 1
    OUTPUT:
        - a matrix A such that det(1 - t A) = R[](coefficients)
    EXAMPLES:
        sage: list(reversed(yellow.MatrixCharpoly([1,2,3]).characteristic_polynomial().list())) == [1,2,3]
            True
    TESTS: 
        sage: test_pass = True
        sage: for i in range(10):
        ....:     cf = [1] + [ZZ.random_element() for _ in range(randint(1,10))];
        ....:     test_pass = test_pass and (list(reversed(MatrixCharpoly(cf).characteristic_polynomial().list())) == cf)
        sage: print test_pass
            True
    """
    assert coefficients[0] == 1
    n = len(coefficients) - 1;
    M = Matrix(1, n - 1).stack(identity_matrix(n-1,n-1));
    M = M.augment(-vector(reversed(coefficients[1:])))
    return M
