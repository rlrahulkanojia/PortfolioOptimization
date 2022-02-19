package main

// Imports
import (
	"fmt"
	"log"
	"math/rand"
	"os"

	"github.com/go-gota/gota/dataframe"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

// General Function declaration

func read_csv(PATH string) dataframe.DataFrame {
	file, err := os.Open(PATH)
	if err != nil {
		log.Fatal(err)
	}
	df := dataframe.ReadCSV(file)
	return df
}

func df_mean(DF dataframe.DataFrame) []float64 {
	MEAN_DATA := make([]float64, DF.Ncol())
	COLS := DF.Names()
	mean_ind := 0
	// MEAN_DATA[mean_ind] = rand.NormFloat64()
	for col_ind_i := 0; col_ind_i < 4; col_ind_i++ {
		MEAN_DATA[mean_ind] = DF.Col(COLS[col_ind_i]).Mean()
		mean_ind += 1
	}

	return MEAN_DATA
}

func sum_float64(array []float64) float64 {
	result := 0.0
	for _, v := range array {
		result += v
	}
	return result
}

func randFloats(min, max float64, n int) []float64 {
	res := make([]float64, n)
	for i := range res {
		res[i] = min + rand.Float64()*(max-min)
	}
	return res
}

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

// Utility Function Declaration

func cov_matrix(DF dataframe.DataFrame) [][]float64 { //
	// cov := stat.Covariance(DF.Col("AAPL").Float(), DF.Col("AAPL").Float(), nil)
	// COV_DATA := make([][]float64, 4, 4) // DF.Ncol()*DF.Ncol()
	COV_DATA := make([][]float64, 4)
	for i := range COV_DATA {
		COV_DATA[i] = make([]float64, 4)
	}

	COLS := DF.Names()
	cov_ind := 0
	for col_ind_i := 0; col_ind_i < 4; col_ind_i++ {
		for col_ind_j := 0; col_ind_j < 4; col_ind_j++ {
			// fmt.Println(COLS[col_ind_i], COLS[col_ind_j])
			COV_DATA[col_ind_i][col_ind_j] = stat.Covariance(DF.Col(COLS[col_ind_i]).Float(), DF.Col(COLS[col_ind_j]).Float(), nil)
			cov_ind += 1
		}
	}
	return COV_DATA
}

// Given is a universe of n assets with
// Σ: an (n × n) positive definite covariance matrix,
// μ: n vector with the assets’ expected returns,
// w: n vector with the assets’ weights.
// l: an n vector containing the asset weights’ lower bounds (wi ≥ li, ∀i),
// u: an n vector containing the asset weights’ upper bounds (wi ≤ ui, ∀i),
// F: a subset of N = {1,2,... ,n} containing all assets’ indices where weights are within their bounds (li < wi < ui). We shall call the corresponding assets free assets.
// B: the subset of all asset indices where the weights lie on one of their bounds (wi = li or wi = ui). Thus, B = {1,2,... ,n} \ F. The sets of assets on their upper and lower bounds will be called U and L, respectively, with B = U ∪ L.
// k ≡ |F|: number of elements in F.

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

func CLA(cov_matrix []float64, expected_returns []float64) [][]float64 {

	// for each asset i
	// j ← argmaxi μi
	max_num, max_ind := findMaxElement(expected_returns)

	// w (1) j ← 1; ∀i 6= j : w (1) i ← 0
	weights := [4]float64{0} // make length constant

	fmt.Println("Max Num ", max_num)
	weights[max_ind] = 1

	// F ← {j}
	F := []float64{}
	F = append(F, float64(max_ind))

	// 	Cov Inverse method

	// cov_matrix := mat.NewDense(4, 4, cov_list)
	// var cov_matrix_inverse mat.Dense
	// _ = cov_matrix_inverse.Inverse(cov_matrix)
	// fmt.Println("cov_matrix_inverse ", cov_matrix_inverse)

	// λcurrent ← ∞
	// λcurrent = math.Inf(1)
	// t ← 1
	iteration := 0
	i_inside := []float64{}
	i_outside := []float64{}
	flag := true
	// until i inside = nil and i outside = nil
	for len(i_inside) != 0 && len(i_outside) != 0 || flag == true {
		flag = false
		iteration += 1
		fmt.Println("Iteration ", iteration)

		I_F := []float64{}
		for ind_i := 0; ind_i < len(F); ind_i++ {
			I_F = append(I_F, 1)
		}

		I_F_MAT := mat.NewDense(1, len(F), I_F)
		fmt.Println(I_F_MAT.Dims())

		COV := mat.NewDense(4, 4, cov_matrix)
		// var cov_matrix_inverse mat.Dense
		// _ = cov_matrix_inverse.Inverse(cov_matrix)
		// fmt.Println("cov_matrix_inverse ", cov_matrix_inverse)

		if len(F) > 0 {
			for ind_i := 0; ind_i < len(F); ind_i++ {
				// do Ci ← −(1 ′ F Σ −1 F 1F)(Σ −1 F μF)i + (1 ′ F Σ −1 F μF)(Σ −1 F 1F)i

				// C_i = -A + B
				// A = 1 ′ Σ −1  1)(Σ −1  μ)i
				// A := A1 * A2 * ( A3)
				// B := B1 * B2 * ( B3)

				// A1 > Transpose of 1
				fmt.Println("Calculating Transpose of Identity Matrix")
				A1 := I_F_MAT.T()

				// Creating inverse covariance Matrix
				fmt.Println("Calculating Inverse of Covariance Matrix")
				var COV_INV mat.Dense
				_ = COV_INV.Inverse(COV)
				rows, cols := COV_INV.Dims()
				inv_list := []float64{}
				fmt.Println(rows, cols)
				for ind_r := 0; ind_r < rows; ind_r++ {
					for ind_c := 0; ind_c < cols; ind_c++ {
						println("    :", ind_r, ind_c)
						fmt.Println("cov_matrix_inverse ", COV_INV.At(ind_r, ind_c))
						inv_list = append(inv_list, COV_INV.At(ind_r, ind_c))
					}
				}
				fmt.Println("DEBUG 1 ")

				COV_INV_matrix := mat.NewDense(rows, cols, inv_list)

				var A2 mat.Dense
				A2.Mul(A1, COV_INV_matrix)

				var A3 mat.Dense
				fmt.Println("Converting A2 Dense To matrix")
				rows, cols = A2.Dims()
				inv_list = []float64{}
				for ind_r := 0; ind_r < rows; ind_r++ {
					for ind_c := 0; ind_c < cols; ind_c++ {
						fmt.Println("cov_matrix_inverse ", COV_INV.At(ind_r, ind_c))
						inv_list = append(inv_list, A1.At(ind_r, ind_c))
					}
				}
				A3_matrix := mat.NewDense(rows, cols, inv_list)

				fmt.Println("Multiplying A3 with Identity Matrix")
				A3.Mul(A3_matrix, I_F_MAT)

				var A4 mat.Dense

				// Σ −1  μ

				required_return := []float64{}
				for ind_ret := 0; ind_ret < len(expected_returns); ind_ret++ {
					required_return = append(required_return, expected_returns[ind_ret])
				}

				// expected_returns

				// A4.Mul(COV_INV_matrix, )

				// ----------------------------------------------------------------------------------------

				part_1_A_list := []float64{}
				part_1_B_list := []float64{}
				rows, cols = A3.Dims()
				inv_list = []float64{}
				for ind_r := 0; ind_r < rows; ind_r++ {
					for ind_c := 0; ind_c < cols; ind_c++ {
						fmt.Println("cov_matrix_inverse ", COV_INV.At(ind_r, ind_c))
						part_1_A_list = append(part_1_A_list, A1.At(ind_r, ind_c))
					}
				}
				part_1_A := mat.NewDense(rows, cols, part_1_A_list)

				rows, cols = A4.Dims()
				inv_list = []float64{}
				for ind_r := 0; ind_r < rows; ind_r++ {
					for ind_c := 0; ind_c < cols; ind_c++ {
						fmt.Println("cov_matrix_inverse ", COV_INV.At(ind_r, ind_c))
						part_1_B_list = append(part_1_B_list, A1.At(ind_r, ind_c))
					}
				}
				part_1_B := mat.NewDense(rows, cols, part_1_B_list)

				var part_1 mat.Dense

				part_1.Mul(part_1_A, part_1_B)

				// ----------------------------------------------------------------------------------------

				// λi ← −(Σ −1 F 1F)i/Ci

			}
			// i inside ← argmaxi∈F{λi|λi < λcurrent}
		}
	}

	// ------------------------------------------------------------------------------------

	turning_points_weights := [][]float64{}
	turning_point_weights := []float64{}
	for ind_i := 0; ind_i < 1; ind_i++ {
		turning_point_weights = append(turning_point_weights, 1)
	}
	turning_points_weights = append(turning_points_weights, turning_point_weights)

	fmt.Println("TEST")
	return turning_points_weights

}

// Main Code
func main() {
	fmt.Println("Reading CSV .......")
	returns := read_csv("../returns.csv")
	returns_mean := df_mean(returns)

	cov_list := cov(returns)

	turning_points_weights := CLA(cov_list, returns_mean)
	fmt.Println(turning_points_weights)

}
