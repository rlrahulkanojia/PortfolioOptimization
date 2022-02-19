package main

import (
	"fmt"
	"math/rand"

	"gonum.org/v1/gonum/mat"
)

func myfunc(p, q int) (int, int, int) {
	return p - q, p * q, p + q
}

func main() {

	// Covariance

	// fmt.Println("Covariance computes the degree to which datasets move together")
	// fmt.Println("about their mean.")
	// x := []float64{8, -3, 7, 8, -4}
	// y := []float64{10, 2, 2, 4, 1}
	// cov := stat.Covariance(x, y, nil)
	// fmt.Printf("Cov = %.4f\n", cov)
	// fmt.Println("If datasets move perfectly together, the variance equals the covariance")
	// y2 := []float64{12, 1, 11, 12, 0}
	// cov2 := stat.Covariance(x, y2, nil)
	// varX := stat.Variance(x, nil)
	// fmt.Printf("Cov2 is %.4f, VarX is %.4f", cov2, varX)

	// Correlation

	// x := []float64{8, -3, 7, 8, -4}
	// y := []float64{10, 5, 6, 3, -1}
	// w := []float64{2, 1.5, 3, 3, 2}

	// fmt.Println("Correlation computes the degree to which two datasets move together")
	// fmt.Println("about their mean. For example, x and y above move similarly.")

	// c := stat.Correlation(x, y, w)
	// fmt.Printf("Correlation is %.5f\n", c)

	data := make([]float64, 18)
	for i := range data {
		data[i] = rand.NormFloat64()
	}
	a := mat.NewDense(3, 6, data)
	fmt.Println(a.Dims())

	fmt.Println(a.T().Dims())

	var myvar1, myvar2, myvar3 = myfunc(4, 2)

	// Display the values
	fmt.Printf("Result is: %d", myvar1)
	fmt.Printf("\nResult is: %d", myvar2)
	fmt.Printf("\nResult is: %d", myvar3)

}
