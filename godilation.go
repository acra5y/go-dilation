// This project is licensed under the terms of the MIT license.
// It provides a function to calculate a unitary n-dilation for a given matrix contraction and a degree using the gonum library.
// The recipe to calculate the dilation follows the theory mentioned by Béla Szőkefalvi-Nagy in "Analyse harmonique des opérateurs de l'espace de Hilbert" (1967).
// The matrix square root needed for the dilation is calculated using the Exponential Method for Matrices.
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
