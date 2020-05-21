
# go-dilation

_Library to calculate a n-dilation of a square matrix contraction with real values_

This project is licensed under the terms of the MIT license.

If you stumbled across this and are interested, found a bug or have another idea how to contribute, feel free to open an issue or a pull request.

The recipe to calculate the dilation follows the theory mentioned by Béla Szőkefalvi-Nagy in "Analyse harmonique des opérateurs de l'espace de Hilbert" (1967). The matrix square root needed for the dilation is calculated using the Exponential Method for Matrices.

## Development

Run any common `go` tasks such as `go test ./...`.
This project was built on `go version go1.13.9`.

## Usage

Build n dilations by calling `UnitaryNDilation(m, n)`, where `t`  is of type `*mat.Dense` (see gonum) - the contraction that will be dilated, and `n` is of type `int` - the degree that the dilation will have.
Note that the defect operator of `t` must have real eigenvalues.
Import the library as github.com/acra5y/go-dilation.

For ill-conditionend matrices (e.g. with high eigenvalues or eigenvalues close to 0) the result might be imprecise.
