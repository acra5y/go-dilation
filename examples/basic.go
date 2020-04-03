package main

import (
	"fmt"
	"github.com/acra5y/go-dilation"
	"gonum.org/v1/gonum/mat"
)

// Print unitary n dilation of matrix T where T is a 2x2 square matrix with real values:
// T = | 0.5 0.5 |
//     | 0   0.5 |
func main() {
	dilation, err := godilation.UnitaryNDilation(mat.NewDense(2, 2, []float64{0.5,0.5,0,0.5,}), 2)

	if (err != nil) {
		fmt.Println(err)
	} else {
		fmt.Println(dilation)
	}
}
