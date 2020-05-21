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
        {value: mat.NewDense(2, 2, []float64{1,0,0,1,}), desc: "for eye matrix"},
        {value: mat.NewDense(3,3, []float64{0.9879, 0.0011, 0.0132, 0.0011, 0.9598, 0, 0.0132, 0, 0.9712}), desc:"for a random non diagonal p.d. matrix"},
        {value: mat.NewDense(3,3, []float64{0.1,0,0.3,0, 0.12, 0.1412,0, 0, 0.12}), desc:"for a random non diagonal p.d. matrix 2"},
        {value: mat.NewDense(3,3, []float64{2, -1, 0, -1, 2, -1,0, -1, 2}), desc:"for a random non diagonal p.d. matrix 3"},
    }

    for _, table := range tables {
        table := table
        t.Run(table.desc, func(t *testing.T) {
            t.Parallel()
            res, err := Calculate(table.value)

            if err != nil {
                t.Errorf("Error: %v.", err)
            }

            n, _ := res.Dims()
            value := mat.NewDense(n, n, nil)
            value.Pow(res, 2)
            diff := mat.NewDense(n, n, nil)
            diff.Sub(table.value, value)
            norm := mat.Norm(diff, 2)

            if norm >= 1e-3 {
                t.Errorf("Wrong result, got: %v, want: %v, norm difference %e", value, table.value, norm)
            }
        })
    }
}
