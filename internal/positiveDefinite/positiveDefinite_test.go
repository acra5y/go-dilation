package positiveDefinite

import (
    "fmt"
    "gonum.org/v1/gonum/mat"
    "reflect"
    "testing"
)

type EigenMock struct {
    factorizeOk bool
    mockData []complex128
}

func (eigen EigenMock) Factorize(a mat.Matrix, left, right bool) bool {
    return eigen.factorizeOk
}

func (eigen EigenMock) Values(dst []complex128) []complex128 {
    return eigen.mockData
}

func createEigenMock(factorizeOk bool, values []complex128) EigenMock {
    eigen := EigenMock{ mockData: values, factorizeOk: factorizeOk }
    return eigen
}

var dummyMatrix = mat.NewDense(2, 2, nil)

func TestPdForMatrix(t *testing.T) {
    tables := []struct {
        desc string
        candidate *mat.Dense
        values []complex128
        isPd bool
        factorizeOk bool
    }{
        {values: []complex128{1+0i,5+0i}, isPd: true, factorizeOk: true, candidate: dummyMatrix, desc: "returns true"},
        {values: []complex128{1+0i,0-5i}, isPd: false, factorizeOk: true, candidate: dummyMatrix, desc: "returns false"},
        {values: []complex128{0+0i,5+7i}, isPd: false, factorizeOk: true, candidate: dummyMatrix, desc: "checks for factorize error"},
        {values: []complex128{}, isPd: false, factorizeOk: true, candidate: mat.NewDense(2, 2, []float64{0,1,0,0}), desc: "checks is symmetric"},
        {values: []complex128{}, isPd: false, factorizeOk: true, candidate: mat.NewDense(2, 3, nil), desc: "checks is square matrix"},
    }

    for _, table := range tables {
        table := table
        t.Run(table.desc, func(t *testing.T) {
            t.Parallel()
            isPd, err := IsPositiveDefinite(createEigenMock(table.factorizeOk, table.values), table.candidate)

            if err != nil {
                t.Errorf("IsPositiveDefinite returned unexpected error: %v", err)
            }
            if isPd != table.isPd {
                t.Errorf("IsPositiveDefinite was incorrect, got: %t, want: %t.", isPd, table.isPd)
            }
        })
    }
}

func TestPdFactorizeNotOk(t *testing.T) {
    isPd, err := IsPositiveDefinite(createEigenMock(false, []complex128{}), dummyMatrix)
    expectedError := fmt.Errorf("eigen: Factorize unsuccessful %v", mat.Formatted(dummyMatrix, mat.Prefix("    "), mat.Squeeze()))
    if !reflect.DeepEqual(err, expectedError) {
        t.Errorf("Wrong error returned, got: %t, want: %t.", err, expectedError)
    }
    if isPd {
        t.Errorf("IsPositiveDefinite was incorrect, got: %t, want: %t.", isPd, false)
    }
}