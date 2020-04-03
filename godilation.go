// This project is licensed under the terms of the MIT license.
// It provides a function to calculate a unitary n-dilation for a given matrix contraction and a degree using the gonum library.
package godilation

import (
    "github.com/acra5y/go-dilation/internal/blockMatrix"
    "github.com/acra5y/go-dilation/internal/dilation"
    "github.com/acra5y/go-dilation/internal/positiveDefinite"
    "github.com/acra5y/go-dilation/internal/squareRoot"

    "gonum.org/v1/gonum/mat"
)

// returns a unitary n-dilation for the given square matrix contraction t or an error, if t is not a contraction or not a square matrix
func UnitaryNDilation(t *mat.Dense, n int) (*mat.Dense, error) {
    return dilation.UnitaryNDilation(positiveDefinite.IsPositiveDefinite, squareRoot.Calculate, blockMatrix.NewBlockMatrixFromSquares, t, n)
}
