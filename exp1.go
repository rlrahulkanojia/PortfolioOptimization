package main

import (
	"fmt"
	"log"
	"math/rand"
	"os"

	"github.com/go-gota/gota/dataframe"
	"gonum.org/v1/gonum/mat"
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

func main() {
	fmt.Println("hello world")
	PATH := "table.csv"
	DF := read_csv(PATH)
	// fmt.Println(DF.Dims())
	// fmt.Println(DF.Nrow())
	// fmt.Println(DF.Ncol())
	// fmt.Println(DF.Names())
	// fmt.Println(DF.Types())
	// fmt.Println(DF.Describe())
	// fmt.Println(DF.Select("GOOGL"))
	// fmt.Println(DF.Select(0))
	// fmt.Println(DF.Select([]string{"GOOGL", "AAPL"}))
	// fmt.Println(DF.Subset(0))

	// Apply conditions
	df_series := DF.Col("GOOGL")
	fmt.Printf("%T \n", df_series)

	// df_series_mean   := df_series.Mean()
	// df_series_median := df_series.Median()
	// fmt.Println(df_series_mean, df_series_median)

	// Check for missing values
	// fmt.Println(df_series.IsNaN())

	// gmean := stat.Mean(df_series.Float(), nil)
	// fmt.Println(gmean)

	// type F struct {
	// 	Colname    string
	// 	Comparator series.Comparator
	// 	Comparando interface{}
	// }

	istest := DF.Filter(dataframe.F{
		Colidx:     0,
		Colname:    "GOOGL",
		Comparator: "!=",
		Comparando: "759.440000",
	})
	fmt.Println(istest)

	// mat.Matrix((1))

	// Allocate a zeroed real matrix of size 3×5
	zero := mat.NewDense(3, 5, nil)
	fmt.Println(zero.Dims())

	// Generate a 6×6 matrix of random values.
	data := make([]float64, 36)
	for i := range data {
		data[i] = rand.NormFloat64()
	}
	a := mat.NewDense(6, 6, data)
	fmt.Println(a.Dims())

	// Operations involving matrix data are implemented as functions when the values of the matrix remain unchanged
	tr := mat.Trace(a)
	fmt.Println(tr)

	// and are implemented as methods when the operation modifies the receiver.
	zero.Copy(a)

}
