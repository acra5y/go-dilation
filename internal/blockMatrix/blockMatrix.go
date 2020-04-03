package blockMatrix

import (
    "fmt"
    "math"
    "gonum.org/v1/gonum/mat"
)

func validateDims(rows [][]*mat.Dense) (bool, error) {
    d0, _ := rows[0][0].Dims()
    n0 := len(rows)

    for i, row := range rows {
        n := len(row)

        if n != n0 {
            return false, fmt.Errorf("Unexpected length of row: %d has length %d (Expecting %d)", i, n, n0)
        }

        for j, matrix := range row {
            d1, d2 := matrix.Dims()

            if d1 != d0 || d2 != d0 {
                return false, fmt.Errorf("Unexpected dimension: (%d, %d) in row %d, col %d (Expecting (%d, %d))", d1, d2, i, j, d0, d0)
            }
        }
    }

    return true, nil
}

func NewBlockMatrixFromSquares(rows [][]*mat.Dense) (*mat.Dense, error) {
    ok, err := validateDims(rows)

    if !ok {
       return mat.NewDense(0, 0, nil), err
    }

    var d0, d, entriesPerBlock, entriesPerBlockRow, blockDim int
    d0, _ = rows[0][0].Dims()
    entriesPerBlock = int(math.Pow(float64(d0), 2))
    blockDim = len(rows)
    entriesPerBlockRow = entriesPerBlock * blockDim

    d = d0 * blockDim

    var data []float64
    data = make([]float64, int(math.Pow(float64(d), 2)))

    for i, _ := range rows {
        offsetHandledBlockRows := i * entriesPerBlockRow
        for j := 0; j < d0; j++ {
            offsetCurrentBlocks := j * d
            for k, matrix := range rows[i] {
                raw := matrix.RawRowView(j)
                offsetHandledBlocks := k * d0
                index := offsetHandledBlockRows + offsetCurrentBlocks + offsetHandledBlocks

                data = append(data[:index], append(raw, data[(index + d0):]...)...)
            }
        }
    }

    return mat.NewDense(d, d, data), nil
}
