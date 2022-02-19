package main

import (
	"math/rand"

	"github.com/go-gota/gota/dataframe"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

func sum_float32(array []float32) float32 {
	result := float32(0)
	for _, v := range array {
		result += v
	}
	return result
}

func matrixReshape(nums [][]float64, r int, c int) [][]float64 {
	if canReshape(nums, r, c) {
		return reshape(nums, r, c)
	}
	return nums
}

func canReshape(nums [][]float64, r, c int) bool {
	row := len(nums)
	colume := len(nums[0])
	if row*colume == r*c {
		return true
	}
	return false
}

func reshape(nums [][]float64, r, c int) [][]float64 {
	newShape := make([][]float64, r)
	for index := range newShape {
		newShape[index] = make([]float64, c)
	}
	rowIndex, colIndex := 0, 0
	for _, row := range nums {
		for _, col := range row {
			if colIndex == c {
				colIndex = 0
				rowIndex++
			}
			newShape[rowIndex][colIndex] = col
			colIndex++
		}
	}
	return newShape
}

type matrix struct {
	dataframe.DataFrame
}

func (m matrix) At(i, j int) float64 {
	return m.Elem(i, j).Float()
}

func (m matrix) T() mat.Matrix {
	return mat.Transpose{m}
}

func cov(DF dataframe.DataFrame) []float64 { //
	// cov := stat.Covariance(DF.Col("AAPL").Float(), DF.Col("AAPL").Float(), nil)
	COV_DATA := make([]float64, DF.Ncol()*DF.Ncol())
	COLS := DF.Names()
	cov_ind := 0
	COV_DATA[cov_ind] = rand.NormFloat64()
	for col_ind_i := 0; col_ind_i < 4; col_ind_i++ {
		for col_ind_j := 0; col_ind_j < 4; col_ind_j++ {
			// fmt.Println(COLS[col_ind_i], COLS[col_ind_j])
			COV_DATA[cov_ind] = stat.Covariance(DF.Col(COLS[col_ind_i]).Float(), DF.Col(COLS[col_ind_j]).Float(), nil)
			cov_ind += 1
		}
	}
	// fmt.Println(COV_DATA)
	return COV_DATA
}
