package squareRoot

import (
    "gonum.org/v1/gonum/mat"
    "testing"
)

var dummyMatrix = mat.NewDense(2, 2, nil)

func TestCalculate(t *testing.T) {
    tables := []struct {
        desc string
        value *mat.Dense
    }{
        {value: mat.NewDense(2, 2, []float64{1,0,0,1,}), desc: "square root"},
    }

    for _, table := range tables {
        table := table
        t.Run(table.desc, func(t *testing.T) {
            t.Parallel()
            res, err := Calculate(table.value)

            if err != nil {
                t.Errorf("Error: %v.", err)
            }

            if !mat.Equal(res, table.value) {
                t.Errorf("Wrong result, got: %v, want: %v", res, table.value)
            }
        })
    }
}
