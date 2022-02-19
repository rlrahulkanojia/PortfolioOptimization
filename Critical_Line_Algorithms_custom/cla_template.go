package main

// Imports
import (
	"fmt"
	"log"
	"math"
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

func Identity_F(F int) mat.Dense {
	I_F := []float64{}
	for ind_i := 0; ind_i < F; ind_i++ {
		I_F = append(I_F, 1)
	}
	I_F_MAT := mat.NewDense(1, F, I_F)
	return *I_F_MAT
}

func Multiply_mat_dense(a mat.Matrix, b mat.Dense) []float64 {
	mat_rows, mat_cols := a.Dims()
	dense_rows, dense_cols := b.Dims()
	fmt.Println(mat_rows, mat_cols)
	fmt.Println(dense_rows, dense_cols)
	weights := []float64{}
	for row := 0; row < mat_rows; row++ {
		for col := 0; col < dense_cols; col++ {
			product := 0.0
			for key := 0; key < mat_cols; key++ {
				product += a.At(row, key) * b.At(key, col)
			}
			weights = append(weights, product)
		}
	}
	return weights
}
func Multiply_mat_mat(a mat.Matrix, b mat.Matrix) []float64 {
	mat_rows, mat_cols := a.Dims()
	dense_rows, dense_cols := b.Dims()
	fmt.Println(mat_rows, mat_cols)
	fmt.Println(dense_rows, dense_cols)
	weights := []float64{}
	for row := 0; row < mat_rows; row++ {
		for col := 0; col < dense_cols; col++ {
			product := 0.0
			for key := 0; key < mat_cols; key++ {
				product += a.At(row, key) * b.At(key, col)
			}
			weights = append(weights, product)
		}
	}
	return weights
}
func Multiply_dense_mat(a mat.Dense, b mat.Matrix) []float64 {
	mat_rows, mat_cols := a.Dims()
	dense_rows, dense_cols := b.Dims()
	fmt.Println(mat_rows, mat_cols)
	fmt.Println(dense_rows, dense_cols)
	weights := []float64{}
	for row := 0; row < mat_rows; row++ {
		for col := 0; col < dense_cols; col++ {
			product := 0.0
			for key := 0; key < mat_cols; key++ {
				product += a.At(row, key) * b.At(key, col)
			}
			weights = append(weights, product)
		}
	}
	return weights
}

func Multiply_dense_dense(a mat.Dense, b mat.Dense) []float64 {
	mat_rows, mat_cols := a.Dims()
	dense_rows, dense_cols := b.Dims()
	fmt.Println(mat_rows, mat_cols)
	fmt.Println(dense_rows, dense_cols)
	weights := []float64{}
	for row := 0; row < mat_rows; row++ {
		for col := 0; col < dense_cols; col++ {
			product := 0.0
			for key := 0; key < mat_cols; key++ {
				product += a.At(row, key) * b.At(key, col)
			}
			weights = append(weights, product)
		}
	}
	return weights
}

func Multiply_dense_list(a mat.Dense, b []float64) []float64 {
	mat_rows, mat_cols := a.Dims()
	dense_cols := len(b)

	fmt.Println(mat_rows, mat_cols)
	fmt.Println(dense_cols)
	weights := []float64{}
	for row := 0; row < mat_rows; row++ {
		for col := 0; col < dense_cols; col++ {
			product := 0.0
			for key := 0; key < mat_cols; key++ {
				product += a.At(row, key) * b[key]
			}
			weights = append(weights, product)
		}
	}
	return weights

}

func divide_list_float(l []float64, den float64) []float64 {
	weights := []float64{}
	for ind_i := 0; ind_i < len(l); ind_i++ {
		weights = append(weights, l[ind_i]/den)
	}
	return weights
}

func CLA(cov_matrix []float64, expected_returns []float64) [][]float64 {

	turning_point_weights := []float64{}
	turning_points_weights := [][]float64{}

	// for each asset i
	// j ← argmaxi μi
	λ_all := []float64{}
	_, max_ind := findMaxElement(expected_returns)

	// w (1) j ← 1; ∀i 6= j : w (1) i ← 0
	weights := []float64{0, 0, 0, 0} // make length constant
	// fmt.Println("Max Num ", max_num)
	weights[max_ind] = 1

	// F ← {j}
	F := []float64{}
	F = append(F, float64(max_ind))

	// λcurrent ← ∞
	λcurrent := math.Inf(1)
	// t ← 1
	iteration := 0
	i_inside := math.NaN()
	i_outside := math.NaN()
	flag := true
	// until i inside = nil and i outside = nil
	for i_inside != math.NaN() && i_outside != math.NaN() || flag == true {
		flag = false
		iteration += 1
		I_F_MAT := Identity_F(len(F))

		COV := mat.NewDense(4, 4, cov_matrix)
		var COV_INV mat.Dense
		_ = COV_INV.Inverse(COV)
		// fmt.Println("cov_matrix_inverse ", cov_matrix_inverse)

		if len(F) > 1 {
			for ind_i := 0; ind_i < len(F); ind_i++ {
				// do Ci ← −(1 ′ F Σ −1 F 1F)(Σ −1 F μF)i + (1 ′ F Σ −1 F μF)(Σ −1 F 1F)i
				// Breaking down to -> C_i = - A + B
				// A = −(1 ′  Σ −1  1)(Σ −1  μ)i
				// A1 = −(1 ′  Σ −1  1) ; A2 = (Σ −1  μ)i
				_, row := I_F_MAT.Dims()
				col, _ := COV_INV.Dims()

				A1_matrix := mat.NewDense(row, col, Multiply_mat_dense(I_F_MAT.T(), COV_INV)) // To be verified
				// fmt.Println("DEBUG A1 ", A1_matrix)
				A1 := mat.NewDense(row, col, Multiply_dense_dense(*A1_matrix, I_F_MAT)) // To be verified - scalar value required
				fmt.Println(A1)

				A2 := mat.NewDense(row, col, Multiply_dense_list(COV_INV, expected_returns))
				fmt.Println(A2)

				A := Multiply_mat_mat(A1, A2)
				fmt.Println(A)
				// add negative sign

				// B =  (1 ′  Σ −1  μ)(Σ −1  1)i
				// B1 = (1 ′  Σ −1  μ)
				B1_matrix := mat.NewDense(row, col, Multiply_mat_dense(I_F_MAT.T(), COV_INV)) // To be verified
				B1 := mat.NewDense(row, col, Multiply_dense_list(*B1_matrix, expected_returns))
				fmt.Println(B1)
				// B2 = (Σ −1  1)i
				B2 := mat.NewDense(row, col, Multiply_dense_dense(COV_INV, I_F_MAT))
				fmt.Println(B2)
				B := Multiply_mat_mat(B1, B2)
				fmt.Println(B)

				// C = -A + B
				λ_all = []float64{}
				// lambda_NUM := -1 * B1 / C
				λ_selected := []float64{}
				// i inside ← argmaxi∈F{λi|λi < λcurrent}
				for ind_j := 0; ind_i < len(λ_all); ind_j++ {
					if λ_all[ind_j] < λcurrent {
						λ_selected = append(λ_selected)
					}
				}
				_, argmax := findMaxElement(λ_selected)
				i_inside = float64(argmax)
			}
		}

		if len(F) < len(expected_returns) {
			for ind_i := 0; ind_i < len(expected_returns); ind_i++ {
				if ind_i != max_ind {
					F = append(F, float64(ind_i))
					I_F_MAT := Identity_F(len(F))
					// Ci ← −(1 ′ Fi Σ −1 Fi 1Fi)(Σ −1 Fi μFi)i + (1 ′ Fi Σ −1 Fi μFi)(Σ −1 Fi 1Fi)i
					// C = -A + B
					// A = − (1 ′ Fi Σ −1 Fi 1Fi) * (Σ −1 Fi μFi)i =>  -1 *  A1 * A2
					// do Ci ← −(1 ′ F Σ −1 F 1F)(Σ −1 F μF)i + (1 ′ F Σ −1 F μF)(Σ −1 F 1F)i
					_, row := I_F_MAT.Dims()
					col, _ := COV_INV.Dims()

					A1_matrix := mat.NewDense(row, col, Multiply_mat_dense(I_F_MAT.T(), COV_INV)) // To be verified
					// fmt.Println("DEBUG A1 ", A1_matrix)
					A1 := mat.NewDense(row, col, Multiply_dense_dense(*A1_matrix, I_F_MAT)) // To be verified - scalar value required
					fmt.Println(A1)

					A2 := mat.NewDense(row, col, Multiply_dense_list(COV_INV, expected_returns))
					fmt.Println(A2)

					A := Multiply_mat_mat(A1, A2)
					fmt.Println(A)
					// add negative sign

					// B =  (1 ′  Σ −1  μ)(Σ −1  1)i
					// B1 = (1 ′  Σ −1  μ)
					B1_matrix := mat.NewDense(row, col, Multiply_mat_dense(I_F_MAT.T(), COV_INV)) // To be verified
					B1 := mat.NewDense(row, col, Multiply_dense_list(*B1_matrix, expected_returns))
					fmt.Println(B1)
					// B2 = (Σ −1  1)i
					B2 := mat.NewDense(row, col, Multiply_dense_dense(COV_INV, I_F_MAT))
					fmt.Println(B2)
					B := Multiply_mat_mat(B1, B2)
					fmt.Println(B)

					// C = -A + B
					λ_all = []float64{}
					// lambda_NUM := -1 * B1 / C
					λ_selected := []float64{}
					// i inside ← argmaxi∈F{λi|λi < λcurrent}
					for ind_j := 0; ind_i < len(λ_all); ind_j++ {
						if λ_all[ind_j] < λcurrent {
							λ_selected = append(λ_selected)
						}
					}
					_, argmax := findMaxElement(λ_selected)
					i_outside = float64(argmax)

				}
			}
		}

		if i_inside != math.NaN() || i_outside != math.NaN() {
			iteration++
		}
		λcurrent = math.Max(λ_all[0], λ_all[0]) // FIX
		// if λi inside = max{λi inside,λi outside}
		if λ_all[int(i_inside)] == λcurrent {
			F = divide_list_float(F, i_inside)
		} else {
			F = divide_list_float(F, i_outside)
		}
		// F changed here too

		_, row := I_F_MAT.Dims()
		col, _ := COV_INV.Dims()

		// A := 1  / 1 ′ F Σ −1 F 1F
		A1_matrix := mat.NewDense(row, col, Multiply_mat_dense(I_F_MAT.T(), COV_INV)) // To be verified
		// fmt.Println("DEBUG A1 ", A1_matrix)
		A1 := mat.NewDense(row, col, Multiply_dense_dense(*A1_matrix, I_F_MAT)) // To be verified - scalar value required
		fmt.Println(A1)
		// A = 1 / A1

		// B := ( (1 ′ F Σ −1 F μF) / (1 ′ F Σ −1 F 1F ) λcurrent)
		B1_matrix := mat.NewDense(row, col, Multiply_mat_dense(I_F_MAT.T(), COV_INV)) // To be verified
		B1 := mat.NewDense(row, col, Multiply_dense_list(*B1_matrix, expected_returns))
		fmt.Println(B1)

		// B = B1 / A1

		// γ := 1  / 1 ′ F Σ −1 F 1F  - ( (1 ′ F Σ −1 F μF) / (1 ′ F Σ −1 F 1F ) λcurrent) => A - B

		// γ = A - B

		// W_f_t := λcurrent * (Σ −1 F μF) + γ * (Σ −1 F 1F)

		turning_point_weights = []float64{}
		for ind_i := 0; ind_i < 1; ind_i++ {
			turning_point_weights = append(turning_point_weights, 1)
		}

	}

	// ------------------------------------------------------------------------------------

	turning_points_weights = append(turning_points_weights, turning_point_weights)

	fmt.Println("------------------------------------------------------------------------")
	return turning_points_weights

}

// Main Code
func main() {
	fmt.Println("------------------------------------------------------------------------")
	returns := read_csv("../returns.csv")
	returns_mean := df_mean(returns)

	cov_list := cov(returns)

	turning_points_weights := CLA(cov_list, returns_mean)
	fmt.Println(turning_points_weights)

	// a := mat.NewDense(3, 2, []float64{
	// 	1, 2,
	// 	3, 4,
	// 	5, 6,
	// })
	// b := mat.NewDense(2, 2, []float64{
	// 	1, 2,
	// 	3, 4,
	// })

	// // weights := []float64{}
	// weights := Multiply_mat_mat(a, b)
	// fmt.Println(weights)

}
