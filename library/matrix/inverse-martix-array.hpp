#pragma once
#include "library/matrix/gauss-elimination-array.hpp"
//逆行列が村座しない場合は、nullopt
template <typename mint, int N>
optional<MatrixArray<mint, N, N>> inverse_matrix_array(const MatrixArray<mint, N, N> &a) {
    static_assert(N > 0);

    MatrixArray<mint, N, 2 * N> m;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            m[i][j] = a[i][j];
        }
        m[i][N + i] = mint(1);
    }

    auto [rank, det] = GaussEliminationArray(m, N, true);

    if (rank != N) return nullopt;

    MatrixArray<mint, N, N> b;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b[i][j] = m[i][j + N];
        }
    }

    return b;
}