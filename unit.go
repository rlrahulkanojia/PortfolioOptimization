package main

import "fmt"

// Imports

// add edge case of multplile largest values
func findMaxElement(arr []float64) (float64, int) {
	min_num := arr[0]
	min_ind := 0

	for i := 0; i < len(arr); i++ {
		if arr[i] > min_num {
			min_num = arr[i]
			min_ind = i
		}
	}
	return min_num, min_ind
}

func getMatrices(covar []float64, f []int, num_assets int, expected_returns []float64) ([]float64, []float64, []float64) {

	covarF := reduceMatrix(covar, f, f, num_assets)
	temp := []int{0}
	meanF := reduceMatrix(expected_returns, f, temp, num_assets)
	b := getB(len(expected_returns), f)

	covarFB := reduceMatrix(covar, f, b, num_assets)
	// wB=reduceMatrix(w[-1],b,[0])
	return covarF, covarFB, meanF //, wB

}

func reduceMatrix(matrix []float64, listX []int, listY []int, num_assets int) []float64 {
	if len(listX) == 0 || len(listY) == 0 {
		return []float64{}
	}

	temp1 := []float64{}
	for ind_i := 0; ind_i < num_assets; ind_i++ {
		temp1 = append(temp1, matrix[(num_assets*ind_i)+listY[0]])
	}

	for ind_i := 0; ind_i < len(listY)-1; ind_i++ {
		for ind_j := 0; ind_j < num_assets; ind_j++ {
			temp1 = append(temp1, matrix[(num_assets*ind_j)+listY[1+ind_i]])
		}
	}

	shrink_param := len(listY)
	temp2 := []float64{}
	for ind_i := 0; ind_i < shrink_param; ind_i++ {
		temp2 = append(temp2, temp1[num_assets*ind_i+listX[0]])
	}

	for ind_i := 0; ind_i < len(listX)-1; ind_i++ {
		for ind_j := 0; ind_j < shrink_param; ind_j++ {
			temp2 = append(temp2, temp1[num_assets*ind_j+listX[ind_i+1]])
		}

	}
	return temp2
}

func contains(s []float64, e float64) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

func getB(length int, f []int) []int {
	ret := []float64{}
	ranged_array := []float64{}
	for ind_i := 0; ind_i < length; ind_i++ {
		ranged_array = append(ranged_array, float64(ind_i))
	}

	for ind_i := 0; ind_i < length; ind_i++ {
		if contains(f, arr[ind_i]) == false {
			ret = append(ret, arr[ind_i])
		}
	}
	return ret
}

func main() {
	expected_returns := []float64{0.00110055, 0.00133973, 0.00119665, 0.00071629}
	_, max_ind := findMaxElement(expected_returns)
	weights := []float64{0, 1, 2, 3}
	weights_ := []float64{3, 2}
	weights[max_ind] = 1
	// matrix := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}
	// listx := []int{3, 2}
	// listy := []int{3, 2}
	// fmt.Println(reduceMatrix(matrix, listx, listy, 4))

	fmt.Println(getB(weights, weights_))

}
