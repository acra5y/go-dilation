package eye

import (
    "testing"
)

func TestOfDimensionReturnsCorrectDimension(t *testing.T) {
    d := 5

    eye := OfDimension(d)

    m, n := eye.Dims()

    if m != d {
        t.Errorf("Wrong dimension: Too many rows. Got %d, want %d", m, d)
    }

    if n != d {
        t.Errorf("Wrong dimension: Too many cols. Got %d, want %d", n, d)
    }
}

func TestOfDimensionOnesOnDiag(t *testing.T) {
    d := 5

    eye := OfDimension(d)
    m, _ := eye.Dims()

    for i := 0; i < m; i++ {
        x := eye.At(i, i)
        if x != 1 {
            t.Errorf("Result is not eye: Wrong entry at: (%d, %d). Got: %f, want: 1", i, i, x)
        }
    }
}

func TestOfDimensionZeros(t *testing.T) {
    d := 5

    eye := OfDimension(d)
    m, n := eye.Dims()

    for i := 0; i < m; i++ {
        for j := 0; j < n; j++ {
            if x := eye.At(i, j); i != j && x != 0 {
                t.Errorf("Result is not eye: Wrong entry at: (%d, %d). Got: %f, want: 0", i, i, x)
            }
        }
    }
}
