package positiveDefinite

import (
    "fmt"
    "gonum.org/v1/gonum/mat"
    "math/cmplx"
)

type EigenComputer interface {
    Factorize(a mat.Matrix, left, right bool) bool
    Values(dst []complex128) []complex128
}

func isSymmetric(a mat.Matrix) bool {
    return mat.Equal(a, a.T())
}

func IsPositiveDefinite(eigen EigenComputer, candidate *mat.Dense) (isPositiveDefinite bool, err error) {
    m, n := candidate.Dims()
    c := mat.NewDense(m, n, nil)
    c.Clone(candidate)
    err = nil
    isPositiveDefinite = true

    if !isSymmetric(candidate) {
        isPositiveDefinite = false
        return
    }

    ok := eigen.Factorize(c, false, false)

    if ok {
        for _, val := range eigen.Values(nil) {
            r, theta := cmplx.Polar(val)
            if theta != 0 || r == 0 {
                isPositiveDefinite = false

                return
            }
        }

        return
    }

    return false, fmt.Errorf("eigen: Factorize unsuccessful %v", mat.Formatted(candidate, mat.Prefix("    "), mat.Squeeze()))
}
