#pragma once
#include "library/matrix/matrix.hpp"
// {rank, det(非正方行列の場合は未定義)} を返す
// 掃き出しまでしたい場合は、GaussElimination(a,-1,true)
template <typename T>
pair<int, T> GaussElimination(Matrix<T> &a, int pivot_end = -1,
                                   bool diagonalize = false) {
    if (a.H() == 0) return {0, 1};
    int H = a.size(), W = a[0].size(), rank = 0;
    if (pivot_end == -1) pivot_end = W;
    T det = 1;
    matrix_internal::EliminationTolerance<T> tolerance(a, H, W);
    for (int j = 0; j < pivot_end; j++) {
		int idx = Matrix<T>::choose_pivot(a, rank, j);
		if (idx == -1) {
			det = 0;
			continue;
		}
		if (rank != idx) {
			det = -det;
			swap(a[rank], a[idx]);
			tolerance.swap_rows(rank, idx);
		}
		det *= a[rank][j];
		if (diagonalize && a[rank][j] != T(1)) {
			T coeff = T(1) / a[rank][j];
			T coeff_scale = tolerance.quotient_scale(coeff, a[rank][j], T(0), rank, j);
			for (int k = j; k < W; k++) {
				tolerance.multiply(a[rank][k], rank, k, coeff, coeff_scale);
			}
		}
		int is = diagonalize ? 0 : rank + 1;
		for (int i = is; i < H; i++) {
			if (i == rank) continue;
			if (a[i][j] != T(0)) {
				T coeff = a[i][j] / a[rank][j];
				T coeff_scale = tolerance.quotient_scale(coeff, a[rank][j], tolerance.entry_scale(i, j), rank, j);
				for (int k = j; k < W; k++) tolerance.subtract_product(a[i][k], a[rank][k], coeff, coeff_scale, i, k, rank);
			}
		}
		rank++;
    }
    return make_pair(rank, det);
}
