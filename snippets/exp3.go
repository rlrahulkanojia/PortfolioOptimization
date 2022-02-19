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

func read_csv(PATH string) dataframe.DataFrame {
	file, err := os.Open(PATH)
	if err != nil {
		log.Fatal(err)
	}
	df := dataframe.ReadCSV(file)
	// fmt.Println(df)
	return df
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

func main() {
	fmt.Println("hello world")
	PATH := "table_filtered.csv"
	DF := read_csv(PATH)
	// fmt.Println(DF.Dims())

	// fmt.Println(DF.)
	// fmt.Println(COV_DATA)
	// fmt.Println()

	COV_DATA := cov(DF)
	// fmt.Println(COV_DATA)
	COV_MATRIX := mat.NewDense(4, 4, COV_DATA)
	fmt.Println(COV_MATRIX))

}
