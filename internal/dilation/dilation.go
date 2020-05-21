package dilation

import (
    "fmt"
    "github.com/acra5y/go-dilation/internal/eye"
    "github.com/acra5y/go-dilation/internal/positiveDefinite"
    "gonum.org/v1/gonum/mat"
)

type isPositiveDefinite func(positiveDefinite.EigenComputer, *mat.Dense) (bool, error)

type squareRoot func(*mat.Dense) (*mat.Dense, error)

type newBlockMatrixFromSquares func([][]*mat.Dense) (*mat.Dense, error)

func defectOperatorSquared (t mat.Matrix) *mat.Dense {
    n, _ := t.Dims()
    eye := eye.OfDimension(n)

    defectSquared := mat.NewDense(n, n, nil)

    defectSquared.Product(t, t.T())

    defectSquared.Sub(eye, defectSquared)
    return defectSquared
}

func negativeTranspose(t *mat.Dense) *mat.Dense {
    m, n := t.Dims()
    data := make([]float64, m * n)

    for i := 0; i < m; i++ {
        for j := 0; j < n; j++ {
            data[m * i + j] = (-1) * t.At(j, i)
        }
    }
    return mat.NewDense(m, n, data)
}

// See E. Levy und O. M. Shalit: Dilation theory in finite dimensions: the possible, the impossible and the unknown. Rocky Mountain J. Math., 44(1):203-221, 2014

func UnitaryNDilation(isPD isPositiveDefinite, sqrt squareRoot, newBlockMatrix newBlockMatrixFromSquares, t *mat.Dense, degree int) (*mat.Dense, error) {
    m, n := t.Dims()

    if m != n {
        return nil, fmt.Errorf("Matrix does not have square dimension")
    }

    defectSquared := defectOperatorSquared(t)

    if pd, _ := isPD(&mat.Eigen{}, defectSquared); !pd {
        return nil, fmt.Errorf("Input is not a contraction")
    }

    defectSquaredOfTranspose := defectOperatorSquared(t.T())
    /*
        We can calculate the square root as it must exist as we assume T is a positive definite contraction:
        Let T be a positiv definite complex matrix of dimension n times n with ||T|| < 1 (where ||T|| denotes the operator norm).
        Let T^* denote the conjugate transpose of T.
        The equivalency (T^*T)^* = (T^*)(T^*)^* = (T^*)T shows T^*T is hermitian.
        Let I by the eye Matrix of the same dimension as T.
        It follows that for a given vector v and the euclidean norm |v|: v^*(I − T^*T)v = v^*Iv − v^*T^*Tv = v^*v − (Tv)^*Tv = |v|^2 − |Tv|^2 > 0.
        The last step ist based on the requirement that ||T|| < 1.
        For a positive definite matrix we then know that a square root must exist.
        See also "Harmonic Analysis of Operators on Hilbert Space" by  B. Sz.-Nagy, chapter I, 1. in section 3.
        (Please note this hint does not have the ambition to be a mathematical proof on its own).
    */
    defect, _ := sqrt(defectSquared)
    defectOfTransposed, _ := sqrt(defectSquaredOfTranspose)

    rows := make([][]*mat.Dense, degree + 1)

    firstRow := make([]*mat.Dense, degree + 1)
    secondRow := make([]*mat.Dense, degree + 1)

    blockDim := degree + 1
    firstRow[0] = t
    firstRow[blockDim - 1] = defectOfTransposed
    secondRow[0] = defect
    secondRow[blockDim - 1] = negativeTranspose(t)

    if degree > 1 {
        for i := 1; i < blockDim - 1; i++ {
            firstRow[i] = mat.NewDense(m, n, nil)
            secondRow[i] = mat.NewDense(m, n, nil)
        }

        for i := 2; i < blockDim; i++ {
            row := make([]*mat.Dense, blockDim)
            for j := 0; j < blockDim; j++ {
                if j == i - 1 {
                    row[j] = eye.OfDimension(m)
                } else {
                    row[j] = mat.NewDense(m, n, nil)
                }
            }
            rows[i] = row
        }
    }

    rows[0] = firstRow
    rows[1] = secondRow

    unitary, err := newBlockMatrix(rows)

    if err != nil {
        return nil, err
    }

    return unitary, nil
}
