package godilation

import (
    "fmt"
    "gonum.org/v1/gonum/mat"
    "reflect"
    "testing"
)

func TestUnitaryNDilation(t *testing.T) {
    tables := []struct {
        desc string
        value *mat.Dense
        expectedValue *mat.Dense
        expectedErr error
    }{
        {
            value: mat.NewDense(2, 2, []float64{0.5,0,0,0.2,}),
            desc: "return a value for a contraction",
            expectedValue: mat.NewDense(6, 6, []float64{
                0.5, 0, 0, 0, 0.8660254037844386, 0,
                0, 0.2, 0, 0, 0, 0.9797958971132712,
                0.8660254037844386, 0, 0, 0, -0.5, 0,
                0, 0.9797958971132712, 0, 0, 0, -0.2,
                0, 0, 1, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
            }),
            expectedErr: nil,
        },
        {
            value: mat.NewDense(2, 2, []float64{0.5,0,0,2,}),
            desc: "return a value if matrix is not a contraction",
            expectedValue: nil,
            expectedErr: fmt.Errorf("Input is not a contraction"),
        },
    }
    for _, table := range tables {
        table := table
        t.Run(table.desc, func(t *testing.T) {
            t.Parallel()
            value, err := UnitaryNDilation(table.value, 2)

            if !reflect.DeepEqual(err, table.expectedErr) {
                t.Errorf("Wrong error, got: %v, want: %v", err, table.expectedErr)
            }

            if table.expectedValue == nil && value != nil {
                t.Errorf("Wrong result, got: %v, want: %v", value, table.expectedValue)
            }

            if table.expectedValue != nil && !mat.Equal(value, table.expectedValue) {
                t.Errorf("Wrong result, got: %v, want: %v", value, table.expectedValue)
            }
        })
    }
}
