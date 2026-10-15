#pragma once
#include "library/matrix/matrix-array.hpp"
// {rank, det}
// det は基本的に H == W かつ pivot_end == W の正方行列の場合のみ意味を持つ
template <typename T, int H, int W>
pair<int, T> GaussEliminationArray(MatrixArray<T, H, W> &a,
                              int pivot_end = -1,
                              bool diagonalize = false) {
    if constexpr (H == 0) {
        return {0, T(1)};
    }

    if (pivot_end == -1) pivot_end = W;
    assert(0 <= pivot_end && pivot_end <= W);

    int rank = 0;
    T det = T(1);
    matrix_internal::EliminationTolerance<T> tolerance(a, H, W);

    for (int j = 0; j < pivot_end; j++) {
        int idx = MatrixArray<T, H, W>::choose_pivot(a, rank, j);

        if (idx == -1) {
            det = T(0);
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
                for (int k = j; k < W; k++) {
                    tolerance.subtract_product(a[i][k], a[rank][k], coeff, coeff_scale, i, k, rank);
                }
            }
        }

        rank++;
    }

    return {rank, det};
}
