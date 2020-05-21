package squareRoot

import (
    "github.com/acra5y/go-dilation/internal/eye"
    "gonum.org/v1/gonum/mat"
    "math"
)

/*
    The algorithm to calculate the square root of a positive definite matrix is taken from
    "A New Algorithm for Computing the Square Rootof a Matrix"
    (https://scholarworks.rit.edu/cgi/viewcontent.cgi?article=10419&context=theses, chapter three):
    1. Declare some nonsingular matrix C with dimensions (n, n).
    2. Initialize i for number of iterations, S_0 = I and S_1 = C.
    3. Initialize Z = C − I.
    4. For i iterations or until S_i becomes too ill-conditioned, do S_{i+1} = 2S_i + (Z)(S_{i−1}),
    5. After iteration steps stop, find S_{i}^{−1}
    .
    6. Set n × n matrix Q = S_{i+1}(S_{i}^{−1}) − I.

    The step 6 is implemented by solving a linear equation to achieve a better numerical stability:
    For a matrix A we denote t(A) as the transposed Matrix
    Using the notation from above, we define Q' := Q - I, R := S_{i}, S := S_{i+1}
    We have the following equivalency: Q' = S*R^{-1} <=> Q'*R = S <=> t(R)*t(Q') = t(S).
    This is a linear equation that we can solve without computing a matrix.
    We assume that the matrix c fulfilles all necessary preconditions.
*/

func nextGuess(c, z, prePredecessor, predecessor *mat.Dense) (guess *mat.Dense) {
    n, _ := c.Dims()
    var p *mat.Dense
    guess = mat.NewDense(n, n, nil)
    p = mat.NewDense(n, n, nil)
    guess.Scale(2, predecessor)
    p.Product(z, prePredecessor)
    guess.Add(guess, p)
    return
}

func isIllConditioned(m* mat.Dense, iteration int) bool {
    n, _ := m.Dims()
    negative := mat.NewDense(n, n, nil)
    negative.Scale(-1, m)
    max := math.Max(mat.Max(m), mat.Max(negative))
    det := mat.Det(m)

    return math.Pow(max, float64(n)) / det > 1e15
}

func Calculate(c *mat.Dense) (sq *mat.Dense, err error) {
    err = nil
    n, _ := c.Dims()
    var m2, m3, eyeN, z *mat.Dense
    eyeN = eye.OfDimension(n)
    sq = mat.NewDense(n, n, nil)
    m2 = mat.NewDense(n, n, nil)
    m3 = mat.NewDense(n, n, nil)
    sq.CloneFrom(eyeN)
    m2.CloneFrom(c)
    z = mat.NewDense(n, n, nil)
    z.Sub(c, eyeN)

    for i := 1; i <= 100; i++ {
        m3 = nextGuess(c, z, sq, m2)
        sq.CloneFrom(m2)
        m2.CloneFrom(m3)

        if (isIllConditioned(m3, i)) {
            break;
        }
    }

    sq.Solve(sq.T(), m2.T())
    sq.Sub(sq.T(), eyeN)

    return
}
