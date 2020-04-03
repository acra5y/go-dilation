package dilation

import (
    "fmt"
    "github.com/acra5y/go-dilation/internal/positiveDefinite"
    "gonum.org/v1/gonum/mat"
    "reflect"
    "testing"
)

func testIsPositiveDefinite(t *testing.T, expected []*mat.Dense, isPD bool) isPositiveDefinite {
    calls := 0
    return func(a positiveDefinite.EigenComputer, candidate *mat.Dense) (bool, error) {
        if !mat.Equal(expected[calls], candidate) {
            t.Errorf("Unexpected argument in call to testIsPositiveDefinite. Got %v: ,want: %v", candidate, expected[calls])
        }
        calls++
        return isPD, nil
    }
}

func testSquareRoot(t *testing.T, expected []*mat.Dense) squareRoot {
    calls := 0
    return func(a *mat.Dense) (*mat.Dense, error) {
        if !mat.Equal(expected[calls], a) {
            t.Errorf("Unexpected argument in call %d to squareRoot. Got: %v, want: %v", calls + 1, a, expected[calls])
        }
        calls++
        return mat.NewDense(2, 2, nil), nil
    }
}

func testNewBlockMatrixFromSquares(t *testing.T, expected [][]*mat.Dense, errToReturn error) newBlockMatrixFromSquares {
    return func(rows [][]*mat.Dense) (*mat.Dense, error) {
        if len(rows) != len(expected) {
            t.Errorf("Unexpected argument in call to newBlockMatrixFromSquares. Wron length of rows, got: %d, want: %d", len(rows), len(expected))
        }
        for i, row := range rows {
            expectedRow := expected[i]
            if len(row) != len(expectedRow) {
            t.Errorf("Unexpected argument in call to newBlockMatrixFromSquares. Wron length of row %d, got: %d, want: %d", i, len(row), len(expectedRow))
            }
            for j, m := range row {
                if !mat.Equal(expectedRow[j], m) {
                    t.Errorf("Unexpected argument in call to newBlockMatrixFromSquares. Wrong block at position (%d, %d), got : %v, want: %v", i, j, m, expected[i][j])
                }
            }
        }
        return mat.NewDense(2, 2, nil), errToReturn
    }
}

func TestUnitaryNDilationSuccess(t *testing.T) {
    tables := []struct {
        desc string
        value *mat.Dense
        degree int
        expectedInSqrt []*mat.Dense
        expectedRows [][]*mat.Dense
    }{
        {
            desc: "Works for degree 1",
            value: mat.NewDense(2, 2, nil),
            degree: 1,
            expectedInSqrt: []*mat.Dense{mat.NewDense(2, 2, []float64{1,0,0,1,}),mat.NewDense(2, 2, []float64{1,0,0,1,}),},
            expectedRows: [][]*mat.Dense{
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),},
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),},
            },
        },
        {
            desc: "Uses negative transpose",
            value: mat.NewDense(2, 2, []float64{0.5,0.5,0,0.5,}),
            degree: 1,
            expectedInSqrt: []*mat.Dense{mat.NewDense(2, 2, []float64{0.5,-0.25,-0.25,0.75,}),mat.NewDense(2, 2, []float64{0.75,-0.25,-0.25,0.5,}),},
            expectedRows: [][]*mat.Dense{
                []*mat.Dense{mat.NewDense(2, 2, []float64{0.5,0.5,0,0.5,}),mat.NewDense(2, 2, nil),},
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, []float64{-0.5,0,-0.5,-0.5}),},
            },
        },
        {
            desc: "Works for degree > 1",
            value: mat.NewDense(2, 2, []float64{0.5,0.5,0,0.5,}),
            degree: 4,
            expectedInSqrt: []*mat.Dense{mat.NewDense(2, 2, []float64{0.5,-0.25,-0.25,0.75,}),mat.NewDense(2, 2, []float64{0.75,-0.25,-0.25,0.5,}),},
            expectedRows: [][]*mat.Dense{
                []*mat.Dense{mat.NewDense(2, 2, []float64{0.5,0.5,0,0.5,}),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),},
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, []float64{-0.5,0,-0.5,-0.5}),},
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, []float64{1,0,0,1,}),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2,nil),},
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, []float64{1,0,0,1,}),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),},
                []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),mat.NewDense(2, 2, []float64{1,0,0,1,}),mat.NewDense(2, 2, nil),},
            },
        },
    }

    for _, table := range tables {
        t.Run(table.desc, func(t *testing.T) {
            unitary, err := UnitaryNDilation(
                testIsPositiveDefinite(t, table.expectedInSqrt, true),
                testSquareRoot(t, table.expectedInSqrt),
                testNewBlockMatrixFromSquares(t, table.expectedRows, nil),
                table.value,
                table.degree,
            )

            if err != nil {
                t.Errorf("Unexpected err, want: %v, got: %v", nil, err)
            }

            if !mat.Equal(unitary, mat.NewDense(2, 2, nil)) {
                t.Errorf("Wrong matrix returned, want: %v, got: %v", mat.NewDense(2, 2, nil), unitary)
            }
        })
    }
}


func TestUnitaryNDilationErrors(t *testing.T) {
    expectedIsPDAndSQArgs := []*mat.Dense{mat.NewDense(2, 2, []float64{1,0,0,1,}),mat.NewDense(2, 2, []float64{1,0,0,1,}),}
    expectedBlockMatrixArgs := [][]*mat.Dense{
        []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),},
        []*mat.Dense{mat.NewDense(2, 2, nil),mat.NewDense(2, 2, nil),},
    }
    tables := []struct {
        desc string
        value *mat.Dense
        isPD bool
        blockMatrixErr error
        expectedError error
    }{
        {
            desc: "returns error when matrix is not square",
            value: mat.NewDense(2, 3, nil),
            isPD: true,
            blockMatrixErr: nil,
            expectedError: fmt.Errorf("Matrix does not have square dimension"),
        },
        {
            desc: "returns error when matrix not positive definite",
            value: mat.NewDense(2, 2, nil),
            isPD: false,
            blockMatrixErr: nil,
            expectedError: fmt.Errorf("Input is not a contraction"),
        },
        {
            desc: "returns error when block matrix can not be built",
            value: mat.NewDense(2, 2, nil),
            isPD: true,
            blockMatrixErr: fmt.Errorf("Some Error"),
            expectedError: fmt.Errorf("Some Error"),
        },
    }

    for _, table := range tables {
        t.Run(table.desc, func(t *testing.T) {
            _, err := UnitaryNDilation(
                testIsPositiveDefinite(t, expectedIsPDAndSQArgs, table.isPD),
                testSquareRoot(t, expectedIsPDAndSQArgs),
                testNewBlockMatrixFromSquares(t, expectedBlockMatrixArgs, table.blockMatrixErr),
                table.value,
                1,
            )

            if !reflect.DeepEqual(err, table.expectedError) {
                t.Errorf("Unexpected err, want: %v, got: %v", table.expectedError, err)
            }
        })
    }
}
