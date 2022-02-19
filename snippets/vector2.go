package main

import (
	"fmt"

	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
)

func main() {
	var myvector []float64
	myvector = append(myvector, 11.0)
	myvector = append(myvector, 5.2)
	fmt.Println(myvector)

	myvec := mat.NewVecDense(2, []float64{11.0, 5.2})
	fmt.Println(*myvec)

	// === using gonum/floats
	// === lightweight, op on slices of floats
	vectorA := []float64{11.0, 5.2, -1.3}
	vectorB := []float64{-7.2, 4.2, 5.1}
	// dot product of A, B
	dotProduct := floats.Dot(vectorA, vectorB)
	fmt.Printf("DotP: %0.2f\n", dotProduct)

	// scale by 1.5
	floats.Scale(1.5, vectorA)
	fmt.Printf("Scaling A by 1.5: %v\n", vectorA)

	// norm/length of B
	normB := floats.Norm(vectorB, 2)
	fmt.Printf("norm of B: %0.2f\n", normB)

	// === using gonum/mat
	vectorC := mat.NewVecDense(3, []float64{11.0, 5.2, -1.3})
	vectorD := mat.NewVecDense(3, []float64{-7.2, 4.2, 5.1})

	dotProduct = mat.Dot(vectorC, vectorD)
	fmt.Printf("DotP: %0.2f\n", dotProduct)

	vectorC.ScaleVec(1.5, vectorC)
	fmt.Printf("Scaling C by 1.5: %v\n", vectorC)

}
