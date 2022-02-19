package main

import (
	"fmt"
	"log"
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

// **************************************************************************************************************

func portfolio_annualised_performance(weights []float64, mean_returns []float64, cov_matrix [][]float64) (float64, float64) {

	fmt.Println("    In Portfolio Annualized ")
	sum := 0.0
	// TEMP1 := make([][]float64, 4, 4)

	fmt.Println("    Calculating Sum and intermidiate 2d array")

	for ind_i := 0; ind_i < 4; ind_i++ {
		for ind_j := 0; ind_j < 4; ind_j++ {
			sum += mean_returns[ind_j] * cov_matrix[ind_j][ind_j]
			// fmt.Println("Debug 2", ind_i, ind_j)
			// TEMP1[ind_j][ind_j] = mean_returns[ind_j] * cov_matrix[ind_j][ind_j]
		}
	}

	fmt.Print(len(weights))
	fmt.Println("    Calculating mean matrix")

	a := mat.NewDense(4, 4, weights)
	fmt.Println("***************")
	fmt.Println("TESTTING ", a)
	// fmt.Println("***************")
	// fmt.Println("***** Portfolio Sum ***********")
	// fmt.Println("SUM :", sum)
	// sum = sum * 252
	// fmt.Println("SUM :", sum)
	// fmt.Println("***************")

	// fmt.Println(len(weights))

	// TEMP2 := 0
	// for ind_i := 0; ind_i < 4; ind_i++ {
	// 	for ind_j := 0; ind_j < 4; ind_j++ {
	// 		sum += mean_returns[ind_j] * cov_matrix[ind_j][ind_j]
	// 		TEMP1[ind_j][ind_j] = mean_returns[ind_j] * cov_matrix[ind_j][ind_j]
	// 	}
	// }

	// std = np.sqrt(np.dot(weights.T, np.dot(cov_matrix, weights))) * np.sqrt(252)
	// return std, returns
	return 1.0, 1.0
}

// }

// Partial
func random_portfolios(num_portfolios int, mean_returns []float64, cov_matrix [][]float64, risk_free_rate float64) { //, ) {
	// num_portfolios, mean_returns, cov_matrix, risk_free_rate
	fmt.Println("    In Random Portfolio ")
	// var weights_record []weights_record_struct
	// weights := make([][]float64, num_portfolios, num_portfolios)
	// var weights [][]float64
	// fmt.Println("    ", len(weights), weights[0][0])
	weights := []float64{}
	// results := mat.NewDense(3, num_portfolios, nil)

	fmt.Println("    Starting For loop for Num Portfolios ")
	for i := 0; i < 1; i++ {

		// raw_weights := randFloats(0, 1, 5)
		// sum_weights := sum_float64(raw_weights)
		// weights := make([]float64, len(raw_weights))
		// for weights_ind := 0; weights_ind < len(raw_weights); weights_ind++ {
		// 	weights[weights_ind] = raw_weights[weights_ind] / sum_weights
		// }

		fmt.Println("    Creating Weights")
		for ind_i := 0; ind_i < num_portfolios; ind_i++ {
			// random_weights_normalised := make([]float64, num_portfolios*num_portfolios)
			random_weights := randFloats(0, 1, 5)
			sum_weights := sum_float64(random_weights)
			for ind_j := 0; ind_j < num_portfolios; ind_j++ {
				weights = append(weights, random_weights[ind_j]/sum_weights)
			}
			// weights = append(weights, random_weights_normalised)

		}

		// record := weights_record_struct{weights: weights}
		// weights_record = append(weights_record, record)
		fmt.Println("***** Weights Record ***********")
		fmt.Println(weights[2])

		var portfolio_std_dev, portfolio_return = portfolio_annualised_performance(weights, mean_returns, cov_matrix)
		fmt.Println(portfolio_std_dev, portfolio_return)

	}

	//     results[0,i] = portfolio_std_dev # volatility
	//     results[1,i] = portfolio_return
	//     results[2,i] = (portfolio_return - risk_free_rate) / portfolio_std_dev # sharpe Ratio
	// return results, weights_record

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
	// fmt.Println(COV_DATA)
	return COV_DATA
}

func main() {
	fmt.Println("Reading CSV .......")
	returns := read_csv("../returns.csv")

	// fmt.Println(returns.Dims())
	fmt.Println("Calculating Returns Mean ")
	returns_mean := df_mean(returns)
	// fmt.Println(returns_mean.DataFrame.dataframe())
	// fmt.Println(returns)

	fmt.Println("Calculating Covariance Matrix ")
	cov_matrix := cov_matrix(returns)
	// fmt.Println("****** Covariance Matrix **********\n", cov_matrix)
	// res := reshape(cov_matrix, 4, 4)
	// fmt.Println(res)
	// fmt.Println("****************")

	// returns_matrix := matrix{returns}
	// fmt.Print(returns_matrix.T)

	// stat.CovarianceMatrix(returns_matrix)
	risk_free_rate := 0.0017

	fmt.Println("Calculating Random Portfolio ")
	random_portfolios(4, returns_mean, cov_matrix, risk_free_rate)

	// num_portfolios := 25000
	// risk_free_rate := 0.0178

}
