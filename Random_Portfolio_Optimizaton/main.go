package main

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

type weights_record_struct struct {
	weights []float64
}

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

// add edge case of multplile min values
func findMinElement(arr []float64) (float64, int) {
	min_num := arr[0]
	min_ind := 0

	for i := 0; i < len(arr); i++ {
		if arr[i] < min_num {
			min_num = arr[i]
			min_ind = i
		}
	}
	return min_num, min_ind
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

type Results struct {
	Deviation    []float64
	Returns      []float64
	Sharpe_Ratio []float64
}

// **************************************************************************************************************

func portfolio_annualised_performance(weights []float64, mean_returns []float64, cov_matrix [][]float64) (float64, float64) {

	//fmt.Println("    In Portfolio Annualized ")
	sum := 0.0

	//fmt.Println("    Calculating Sum and intermidiate 2d array")

	for ind_i := 0; ind_i < 4; ind_i++ {
		for ind_j := 0; ind_j < 4; ind_j++ {
			sum += mean_returns[ind_j] * cov_matrix[ind_j][ind_j]
			// fmt.Println("Debug 2", ind_i, ind_j)
			// TEMP1[ind_j][ind_j] = mean_returns[ind_j] * cov_matrix[ind_j][ind_j]
		}
	}

	//fmt.Println("    Length of weights :", len(weights))
	//fmt.Println("    Calculating weghts matrix")

	weights_vector := mat.NewVecDense(4, weights)

	//fmt.Println("    Calculating COV matrix")

	//fmt.Println("    Calculating DOT Product")

	DOT_SUM := []float64{}

	for ind_i := 0; ind_i < 4; ind_i++ {
		// fmt.Println(" Debug 1 ", len(cov_matrix[ind_i]))
		cov_vector := mat.NewVecDense(4, cov_matrix[ind_i])
		dotProduct := mat.Dot(weights_vector, cov_vector)
		DOT_SUM = append(DOT_SUM, dotProduct)
		// fmt.Println(" Sum dot ------> ", dotProduct)
	}

	//fmt.Println("    Calculations Done")
	// fmt.Println("    ***************")
	// fmt.Println("    Testing ", weights_vector)
	// fmt.Println("    Testing ", cov_vector)
	// fmt.Println("    Testing ", DOT_SUM)

	DOT_SUM_vector := mat.NewVecDense(4, DOT_SUM)
	// weights_matrix = mat.NewDense(4,1, weights)

	dotProduct := mat.Dot(weights_vector, DOT_SUM_vector)
	// fmt.Println(" Sum dot 2 ------> ", dotProduct)

	std := math.Sqrt(dotProduct) * math.Sqrt(252)
	//fmt.Println("    Standard Deviation  ------> ", std)

	mean_returns_vector := mat.NewVecDense(4, mean_returns)
	returns := mat.Dot(weights_vector, mean_returns_vector) * 252
	//fmt.Println("    Returns             ------> ", returns)

	return std, returns
}

// Partial
func random_portfolios(num_portfolios int, mean_returns []float64, cov_matrix [][]float64, risk_free_rate float64) (Results, [][]float64) { //, ) {
	// num_portfolios, mean_returns, cov_matrix, risk_free_rate
	//fmt.Println("    In Random Portfolio ")
	// var we ights_record []weights_record_struct
	// weights := make([][]float64, num_portfolios, num_portfolios)
	// var weights [][]float64
	// fmt.Println("    ", len(weights), weights[0][0])
	weights := [][]float64{}
	var result Results

	//fmt.Println("    Starting For loop for Num Portfolios ")
	for ind_i := 0; ind_i < num_portfolios; ind_i++ {
		weight := []float64{}

		// raw_weights := randFloats(0, 1, 5)
		// sum_weights := sum_float64(raw_weights)
		// weights := make([]float64, len(raw_weights))
		// for weights_ind := 0; weights_ind < len(raw_weights); weights_ind++ {
		// 	weights[weights_ind] = raw_weights[weights_ind] / sum_weights
		// }

		//fmt.Println("    Creating Weights")
		random_weights := randFloats(0, 1, 4)
		sum_weights := sum_float64(random_weights)
		for ind_j := 0; ind_j < 4; ind_j++ {
			//fmt.Println(" DEBUG ", ind_j)
			// weights = append(weight, random_weights[ind_j]/sum_weights)
			weight = append(weight, random_weights[ind_j]/sum_weights)
		}
		weights = append(weights, weight)

		//fmt.Println("\n***** Weights Record ***********")
		//fmt.Println("    ", weight[2])

		var portfolio_std_dev, portfolio_return = portfolio_annualised_performance(weight, mean_returns, cov_matrix)
		//fmt.Println("\n    Loop iterated", portfolio_std_dev, portfolio_return)

		result.Deviation = append(result.Deviation, portfolio_std_dev)
		result.Returns = append(result.Returns, portfolio_return)
		result.Sharpe_Ratio = append(result.Sharpe_Ratio, (portfolio_return-risk_free_rate)/portfolio_std_dev)

	}
	return result, weights
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

func main() {
	// fmt.Println("Reading CSV .......")
	returns := read_csv("../returns.csv")

	NUM_PORTFOLIO := 25000

	// fmt.Println(returns.Dims())
	// fmt.Println("Calculating Returns Mean ")
	returns_mean := df_mean(returns)
	// fmt.Println(returns_mean.DataFrame.dataframe())
	// fmt.Println(returns)

	//fmt.Println("Calculating Covariance Matrix ")
	cov_matrix := cov_matrix(returns)
	// fmt.Println("****** Covariance Matrix **********\n", cov_matrix)
	// res := reshape(cov_matrix, 4, 4)
	// fmt.Println(res)
	// fmt.Println("****************")

	// returns_matrix := matrix{returns}
	// fmt.Print(returns_matrix.T)

	// stat.CovarianceMatrix(returns_matrix)
	risk_free_rate := 0.0017

	//fmt.Println("Calculating Random Portfolio ")
	var results Results
	var weights [][]float64
	results, weights = random_portfolios(NUM_PORTFOLIO, returns_mean, cov_matrix, risk_free_rate)
	// fmt.Println(len(results.Sharpe_Ratio), len(weights))

	// max_sharpe_allocation
	max_sharpe_ratio, max_sharpe_idx := findMaxElement(results.Sharpe_Ratio)
	// fmt.Println("Maximum Sharpe Ratio ", max_sharpe_ratio)
	sdp, rp := results.Deviation[max_sharpe_idx], results.Returns[max_sharpe_idx]
	max_sharpe_allocation := weights[max_sharpe_idx]
	for ind_i := 0; ind_i < 4; ind_i++ {
		max_sharpe_allocation[ind_i] = max_sharpe_allocation[ind_i] * 100
	}
	// fmt.Println(sdp, rp, max_sharpe_allocation)

	min_deviation, min_vol_idx := findMinElement(results.Deviation)
	// fmt.Println("Minimum Deviation ", min_deviation)
	sdp_min, rp_min := results.Deviation[min_vol_idx], results.Returns[min_vol_idx]
	min_vol_allocation := weights[min_vol_idx]
	for ind_i := 0; ind_i < 4; ind_i++ {
		min_vol_allocation[ind_i] = min_vol_allocation[ind_i] * 100
	}
	// fmt.Println(sdp_min, rp_min, min_vol_allocation)

	fmt.Println("---------------------------------------------")
	fmt.Println("Maximum Sharpe Ratio Portfolio Allocation : ", max_sharpe_ratio)
	fmt.Println("Annualised Return:", rp)
	fmt.Println("Annualised Volatility:", sdp)
	fmt.Println(max_sharpe_allocation)
	fmt.Println("---------------------------------------------")
	fmt.Println("Minimum Volatility Portfolio Allocation : ", min_deviation)
	fmt.Println("Annualised Return:", rp_min)
	fmt.Println("Annualised Volatility:", sdp_min)
	fmt.Println(min_vol_allocation)
	// AAPL  AMZN    FB  GOOGL

	// num_portfolios := 25000
	// risk_free_rate := 0.0178

}
