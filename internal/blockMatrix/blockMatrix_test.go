package blockMatrix

import (
    "gonum.org/v1/gonum/mat"
    "fmt"
    "math"
    "reflect"
    "testing"
)

func createRows(n, dim int) (rows [][]*mat.Dense) {
    rows = make([][]*mat.Dense, n)

    for i := range rows {
        row := make([]*mat.Dense, n)

        for j := range row {
        value := float64(n * i + j)
            data := make([]float64, int(math.Pow(float64(dim), 2)))
            for k, _ := range data {
                data[k] = value
            }
            row[j] = mat.NewDense(dim, dim, data)
        }

        rows[i] = row
    }
    return
}

func TestNewBlockMatrixFromSquares(t *testing.T) {
    tables := []struct {
        desc string
        rows [][]*mat.Dense
        expected *mat.Dense
        err error
    }{
        {
            rows: createRows(2, 2),
            expected: mat.NewDense(4, 4, []float64{0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3,}),
            desc: "returns correct matrix for 4 2x2 blocks",
            err: nil,
        },
        {
            rows: createRows(3, 2),
            expected: mat.NewDense(6, 6, []float64{0,0,1,1,2,2,0,0,1,1,2,2,3,3,4,4,5,5,3,3,4,4,5,5,6,6,7,7,8,8,6,6,7,7,8,8,}),
            desc: "returns correct matrix for 9 2x2 blocks",
            err: nil,
        },
        {
            rows: createRows(2, 3),
            expected: mat.NewDense(6, 6, []float64{0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,2,2,2,3,3,3,2,2,2,3,3,3,2,2,2,3,3,3,}),
            desc: "returns correct matrix for 4 3x3 blocks",
            err: nil,
        },
        {
            rows: createRows(5, 1),
            expected: mat.NewDense(5, 5, []float64{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,}),
            desc: "returns correct matrix for 25 1x1 blocks",
            err: nil,
        },
        {
            rows: createRows(1, 5),
            expected: mat.NewDense(5, 5, nil),
            desc: "returns correct matrix for 1 5x5 blocks",
            err: nil,
        },
        {
            rows: [][]*mat.Dense{
                []*mat.Dense{mat.NewDense(1, 1, nil), mat.NewDense(1, 2, nil),},
                []*mat.Dense{mat.NewDense(1, 1, nil), mat.NewDense(1, 1, nil),},
            },
            expected: mat.NewDense(0, 0, nil),
            desc: "validates all matrices have the same amount of columns",
            err: fmt.Errorf("Unexpected dimension: (1, 2) in row 0, col 1 (Expecting (1, 1))"),
        },
        {
            rows: [][]*mat.Dense{
                []*mat.Dense{mat.NewDense(1, 1, nil), mat.NewDense(2, 1, nil),},
                []*mat.Dense{mat.NewDense(1, 1, nil), mat.NewDense(1, 1, nil),},
            },
            expected: mat.NewDense(0, 0, nil),
            desc: "validates all matrices have the same amount of rows",
            err: fmt.Errorf("Unexpected dimension: (2, 1) in row 0, col 1 (Expecting (1, 1))"),
        },
        {
            rows: [][]*mat.Dense{
                []*mat.Dense{mat.NewDense(1, 1, nil), mat.NewDense(1, 1, nil),},
                []*mat.Dense{mat.NewDense(1, 1, nil),mat.NewDense(1, 1, nil),mat.NewDense(1, 1, nil),},
            },
            expected: mat.NewDense(0, 0, nil),
            desc: "validates length of each row is the same",
            err: fmt.Errorf("Unexpected length of row: 1 has length 3 (Expecting 2)"),
        },
    }

    for _, table := range tables {
        table := table
        t.Run(table.desc, func(t *testing.T) {
            t.Parallel()

            blockMatrix, err := NewBlockMatrixFromSquares(table.rows)

            if !reflect.DeepEqual(err, table.err) {
                t.Errorf("NewBlockMatrix returned wrong value for err, got: %v, want: %v.", err, table.err)
            }

            if !mat.Equal(blockMatrix, table.expected) {
                t.Errorf("NewBlockMatrix returned wrong value, got: %v, want: %v.", blockMatrix, table.expected)
            }
        })
    }
}
