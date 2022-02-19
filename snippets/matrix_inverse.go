package main

import (
	"fmt"
	"log"

	"gonum.org/v1/gonum/mat"
)

func main() {
	// Initialize a matrix A.
	a := mat.NewDense(2, 2, []float64{
		2, 1,
		6, 4,
	})

	// Compute the inverse of A.
	var aInv mat.Dense
	err := aInv.Inverse(a)
	if err != nil {
		log.Fatalf("A is not invertible: %v", err)
	}

	// Print the result using the formatter.
	fa := mat.Formatted(&aInv, mat.Prefix("       "), mat.Squeeze())
	fmt.Printf("aInv = %.2g\n\n", fa)

	// Confirm that A * A^-1 = I.
	var I mat.Dense
	I.Mul(a, &aInv)
	fi := mat.Formatted(&I, mat.Prefix("    "), mat.Squeeze())
	fmt.Printf("I = %v\n\n", fi)

	// The Inverse operation, however, should typically be avoided. If the
	// goal is to solve a linear system
	//  A * X = B,
	// then the inverse is not needed and computing the solution as
	// X = A^{-1} * B is slower and has worse stability properties than
	// solving the original problem. In this case, the SolveVec method of
	// VecDense (if B is a vector) or Solve method of Dense (if B is a
	// matrix) should be used instead of computing the Inverse of A.
	b := mat.NewDense(2, 2, []float64{
		2, 3,
		1, 2,
	})
	var x mat.Dense
	err = x.Solve(a, b)
	if err != nil {
		log.Fatalf("no solution: %v", err)
	}

	// Print the result using the formatter.
	fx := mat.Formatted(&x, mat.Prefix("    "), mat.Squeeze())
	fmt.Printf("x = %.1f", fx)

}
