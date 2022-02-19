// Case a) Free assets move to its bound
if len(F) > 1 {
	for ind_i := 0; ind_i < len(F); ind_i++ {
		// do Ci ← −(1 ′ F Σ −1 F 1F)(Σ −1 F μF)i + (1 ′ F Σ −1 F μF)(Σ −1 F 1F)i
		// λi ← −(Σ −1 F 1F)i/Ci

	}
	// i inside ← argmaxi∈F{λi|λi < λcurrent}
}

// Case b) Asset on its bound becomes free
if len(F) < len(expected_returns) {
	// for i not in F
	// Do
	// Fi ← F ∪ {i}
	// Ci ← −(1 ′ Fi Σ −1 Fi 1Fi)(Σ −1 Fi μFi)i + (1 ′ Fi Σ −1 Fi μFi)(Σ −1 Fi 1Fi)i
	// λi ← (Σ −1 Fi 1Fi)i/Ci

	// i outside ← argmaxi/ ∈F {λi|λi < λcurre
}

// Find Turning Points by comparing cases
// if i inside 6= nil or i outside 6= nil
// 		then t ← t + 1
// 			λcurrent ← max{λi inside,λi outside}
// 			if λi inside = max{λi inside,λi outside}
// 				then F ← F\{i inside}
// 				else F ← F ∪ {i outside}
// 			γ ← 1 1 ′ F Σ −1 F 1F − 1 ′ F Σ −1 F μF 1 ′ F Σ −1 F 1F λcurrent
// 			w (t) F ← λcurrentΣ −1 F μF + γΣ −1 F 1F
//

}

	// 	Cov Inverse method

	// cov_matrix := mat.NewDense(4, 4, cov_list)
	// var cov_matrix_inverse mat.Dense
	// _ = cov_matrix_inverse.Inverse(cov_matrix)
	// fmt.Println("cov_matrix_inverse ", cov_matrix_inverse)
