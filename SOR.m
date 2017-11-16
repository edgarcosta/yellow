// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

SOR := function(A, b, omega, phi0, epsilon, maxiterations)
    n := Ncols(A);
    phi := Vector(phi0);
    for k := 1 to maxiterations do
        for j := 1 to n do
        sigma := 0;
            for i := 1 to n do
                if j ne i then
                    sigma := sigma + phi[i] * A[i,j];
                end if;
            end for;
            phi[j] := (1 - omega) * phi[j] + omega * (b[j] - sigma)/A[j,j];
        end for;
        err := Norm(b - phi * A);
        if err le epsilon then
            print k;
            return phi;
        end if;
    print err;
    end for;
    print "failed to converge :(";
    return phi;
end function;
