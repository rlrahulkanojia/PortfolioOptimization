package main

// Imports
import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"

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

func contains(s []int, e int) bool {
	for _, a := range s {
		if a == e {
			return true
		}
	}
	return false
}

func getB(length int, f []int) []int {
	ret := []int{}
	ranged_array := []int{}
	for ind_i := 0; ind_i < length; ind_i++ {
		ranged_array = append(ranged_array, ind_i)
	}

	for ind_i := 0; ind_i < length; ind_i++ {
		if contains(f, ranged_array[ind_i]) == false {
			ret = append(ret, ranged_array[ind_i])
		}
	}
	return ret
}
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
	for ind_i := 0; ind_i < F*F; ind_i++ {
		I_F = append(I_F, 1)
	}
	I_F_MAT := mat.NewDense(F, F, I_F)
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

func Multiply_list_list(a []float64, b []float64) []float64 {
	aa := []float64{}
	ans := 0.0
	for row := 0; row < len(a); row++ {
		ans += a[row] * b[row]
	}
	return append(aa, ans)
}

func Multiply_matrix_list(a mat.Matrix, b []float64) []float64 {
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

type initAlgo_type struct {
	index int
	mu    float64
}

func divide_list_float(l []float64, den float64) []float64 {
	weights := []float64{}
	for ind_i := 0; ind_i < len(l); ind_i++ {
		weights = append(weights, l[ind_i]/den)
	}
	return weights
}

func sum_f64(array []float64) float64 {
	result := 0.0
	for _, v := range array {
		result += v
	}
	return result
}

func initAlgo(mean []float64, lb []float64, uB []float64) ([]int, []float64) {

	temp := []initAlgo_type{}
	for ind_i := 0; ind_i < len(mean); ind_i++ {
		temp = append(temp, initAlgo_type{ind_i, mean[ind_i]})
	}
	sort.Slice(temp, func(i, j int) bool { return temp[i].mu < temp[j].mu })

	i := len(mean)
	w := lb
	for {
		if (i == 0) || sum_f64(w) <= 0 {
			break
		}
		i = i - 1
		w[temp[i].index] = uB[temp[i].index]
	}
	f := []int{}
	w[temp[i].index] += 1 - sum_f64(w)
	return append(f, temp[i].index), w
}

func getMatrices(covar []float64, f []int, num_assets int, mean []float64, turning_points_weights [][]float64) ([]float64, []float64, []float64, []float64) {

	covarF := reduceMatrix(covar, f, f, num_assets)
	// meanF := reduceMatrix(mean, f, temp, num_assets)
	meanF := reduceMean(mean, f)
	b := getB(len(mean), f)
	covarFB := reduceMatrix(covar, f, b, num_assets)
	wB := reduceMean(turning_points_weights[len(turning_points_weights)-1], b)
	return covarF, covarFB, meanF, wB

}

func reduceMean(mean []float64, f []int) []float64 {
	out := []float64{}
	for ind_i := 0; ind_i < len(f); ind_i++ {
		out = append(out, mean[f[ind_i]])
	}
	return out
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

// l,bi=self.computeLambda(covarF_inv,covarFB,meanF,wB,meanF.shape[0]-1, self.w[-1][i])
func computeLambda(covarF_inv mat.Dense,
	covarFB []float64,
	meanF []float64,
	wB []float64,
	meanF_len int,
	turning_points_weights [][]float64,
	b_index int) (float64, float64) {

	fmt.Println(" Initializing lambda computation ")
	// 1) C
	// onesF=np.ones(meanF.shape)
	onesF := Identity_F(len(meanF))
	// c1 = np.dot(np.dot(onesF.T, covarF_inv), onesF)

	_, row := onesF.Dims()
	col, _ := covarF_inv.Dims()
	// np.dot(onesF.T, covarF_inv)
	c1_matrix := mat.NewDense(row, col, Multiply_mat_dense(onesF.T(), covarF_inv)) // To be verified
	fmt.Println("C1 Debug :")
	// c1 := mat.NewDense(row, col, Multiply_dense_dense(*c1_matrix, onesF))
	c1 := Multiply_dense_dense(*c1_matrix, onesF) // To be verified - scalar value required
	// fmt.Println("C1 : ", c1)

	// c2=np.dot(covarF_inv,meanF)
	// c2 := mat.NewDense(row, col, Multiply_dense_list(covarF_inv, meanF)
	c2 := Multiply_dense_list(covarF_inv, meanF)
	// fmt.Println("C2 :", c2)

	// c3=np.dot(np.dot(onesF.T,covarF_inv),meanF)
	// c3 := mat.NewDense(row, col, Multiply_dense_list(*c1_matrix, meanF))
	c3 := Multiply_dense_list(*c1_matrix, meanF)
	// fmt.Println("C3 :", c3)

	// c4=np.dot(covarF_inv,onesF)
	// c4 := mat.NewDense(row, col, Multiply_dense_dense(covarF_inv, onesF))
	c4 := Multiply_dense_dense(covarF_inv, onesF)
	fmt.Println(" --- C4 --- :", c4)
	return 1.0, 1.0

	// c=-c1*c2[i]+c3*c4[i]
	c := -1*c1[0]*c2[meanF_len] + c3[0]*c4[meanF_len]
	if c == 0 {
		return math.NaN(), math.NaN()
	}
	// temp := turning_points_weights[len(turning_points_weights)-1]
	// if reflect.ValueOf(temp).Kind() == "[]float64" {
	// if type(bi)==list:
	// bi=self.computeBi(c,bi)

	// FIX BI logic
	bi := 1.0
	if len(wB) == 0 {
		return (c4[meanF_len] - c1[0]*bi) / c, bi // to be verified
	} else {
		// onesB=np.ones(wB.shape)
		onesB := Identity_F(len(wB))
		// l1=np.dot(onesB.T,wB)
		l1 := Multiply_matrix_list(onesB.T(), wB)
		// l2=np.dot(covarF_inv,covarFB)
		l2 := Multiply_dense_list(covarF_inv, covarFB)
		// l3=np.dot(l2,wB)
		l3 := Multiply_list_list(l2, wB)
		// l2=np.dot(onesF.T,l3)
		l2 = Multiply_matrix_list(onesF.T(), l3)
		bi := []float64{1.0}
		fmt.Println(bi, l1)
		// return ((1-l1[0]+l2[0])*c4[meanF_len] - c1*(bi+l3[meanF_len])) / c, bi

	}

	// }

	return 1.0, 1.0

}

func Zeros(a int) []float64 {
	out := []float64{}
	for ind_i := 0; ind_i < a; ind_i++ {
		out = append(out, 0)
	}
	return out
}

func remove_from_list(slice []float64, s float64) []float64 {
	out := []float64{}
	for ind_i := 0; ind_i < len(slice); ind_i++ {
		if slice[ind_i] != s {
			out = append(out, 0)
		}
	}
	return out
}

func CLA(covar []float64, mean []float64, lower_bound []float64, upper_bound []float64, NUM_ASSETS int) [][]float64 {

	// Global Variables

	/*
		Notes :
		1) Fix covar_inverse incorrect references ( covarF_inv_list_pre > covarF_inv_list)
	*/

	// // Initialization
	Y := []float64{} // gammas
	λ := []float64{} // Lambdas
	// free_weights := []float64{}
	turning_point_weights := []float64{0, 0, 0, 0} // temporary varaible for weights
	turning_points_weights := [][]float64{}        // solution

	F := [][]int{}
	// flag := true

	// λcurrent := math.Inf(1)
	// i_inside := math.NaN()
	// i_outside := math.NaN()
	// iteration := 0

	f, turning_point_weights := initAlgo(mean, lower_bound, upper_bound)
	turning_points_weights = append(turning_points_weights, turning_point_weights)
	λ = append(λ, math.NaN()) // l.append(None)
	Y = append(Y, math.NaN()) // g.append(None)
	F = append(F, f)          // self.f.append(f[:])

	for {

		if Y[len(Y)-1] == 0 {
			break
		}
		l_in := math.NaN()
		// 1) case a): Bound one free weight
		if len(f) > 1 {
			fmt.Println(" -- DEBUG --")
			covarF, covarFB, meanF, wB := getMatrices(covar, f, NUM_ASSETS, mean, turning_points_weights)
			fmt.Println("COVAR F ", covarF)
			fmt.Println("COVAR FB ", covarFB)
			fmt.Println("Mean F ", meanF)
			fmt.Println("wB", wB)
			break
		}
		// 2) case b): Free one bounded weight
		l_out := math.NaN()
		i_out := 0
		if len(f) < len(mean) {
			b := getB(len(mean), f)
			fmt.Println("B : ", b)
			for ind_i := 0; ind_i < len(b); ind_i++ {
				// covarF,covarFB,meanF,wB=self.getMatrices(f+[i])
				covarF, covarFB, meanF, wB := getMatrices(covar, append(f, b[ind_i]), NUM_ASSETS, mean, turning_points_weights)
				// covarF_inv=np.linalg.inv(covarF)

				fmt.Println("covarF  : ", covarF)
				fmt.Println("covarFB : ", covarFB)
				fmt.Println("meanF   : ", meanF)
				fmt.Println("wB      : ", wB)

				var covarF_inv mat.Dense
				// Row and cols are tricky here
				rows := len(covarF) / 2
				cols := len(covarF) / 2
				covarF_inv_list := mat.NewDense(rows, cols, covarF)
				_ = covarF_inv.Inverse(covarF_inv_list)
				fmt.Println(covarF_inv)
				fmt.Println(" ----- Computing Lambda ------")
				// fix references to index and variable in computelambda functions.
				l, _ := computeLambda(covarF_inv, covarFB, meanF, wB, len(meanF), turning_points_weights, b[ind_i])
				if (λ[len(λ)-1] == math.NaN() || l < λ[len(λ)-1]) && l > l_out {
					l_out = l
					i_out = b[ind_i]
				}
			}

			// 3) compute minimum variance solution
			if (l_in == math.NaN() || l_in < 0) && (l_out == math.NaN() || l_out < 0) {
				λ = append(λ, 0)
				covarF, _, meanF, _ := getMatrices(covar, f, NUM_ASSETS, mean, turning_points_weights)
				// covarF_inv=np.linalg.inv(covarF)
				var covarF_inv mat.Dense
				// Row and cols are tricky here
				rows := len(covarF) / 2
				cols := len(covarF) / 2
				covarF_inv_list := mat.NewDense(rows, cols, covarF)
				_ = covarF_inv.Inverse(covarF_inv_list)
				meanF = Zeros(len(meanF))
			} else {
				// #) decide lambda
				if l_in > l_out {
					λ = append(λ, l_in)
					// f = remove_from_list(f, i_in) // to be added after fixing case 1
					// turning_point_weights[i_in] = bi_in // to be added after fixing case 1
				} else {
					λ = append(λ, l_out)
					f = append(f, i_out)
					covarF, _, _, _ := getMatrices(covar, f, NUM_ASSETS, mean, turning_points_weights)
					// covarF_inv=np.linalg.inv(covarF)
					var covarF_inv mat.Dense
					// Row and cols are tricky here
					rows := len(covarF) / 2
					cols := len(covarF) / 2
					covarF_inv_list_pre := mat.NewDense(rows, cols, covarF)
					_ = covarF_inv.Inverse(covarF_inv_list_pre)

				}
			}

			// 5) compute solution vector
			// wF,g=self.computeW(covarF_inv,covarFB,meanF,wB)

		}

	}

	for ind_i := 0; ind_i < 1; ind_i++ {
		turning_point_weights = append(turning_point_weights, 1)
	}

	// ------------------------------------------------------------------------------------

	turning_points_weights = append(turning_points_weights, turning_point_weights)
	return turning_points_weights

}

// Main Code
func main() {
	fmt.Println("------------------------------------------------------------------------")
	returns := read_csv("../returns.csv")
	returns_mean := df_mean(returns)
	lb := []float64{0.1, 0.1, 0.1, 0.1} // n size list containing lower bound of n assets
	ub := []float64{0.4, 0.2, 0.2, 0.2} // n size list containing upper bound of n assets
	NUM_ASSETS := 4
	covar := cov(returns)

	fmt.Println(initAlgo(returns_mean, lb, ub))

	turning_points_weights := CLA(covar, returns_mean, lb, ub, NUM_ASSETS)
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
