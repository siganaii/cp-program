#pragma once
#include "library/matrix/gauss-elimination.hpp"
template <typename mint>
Matrix<mint> inverse_matrix(const Matrix<mint> & a) {
    int N = a.size();
    assert(N > 0);
    assert(N == (int)a[0].size());

    Matrix<mint> m(N,2 * N);
    for (int i = 0; i < N; i++) {
		copy(begin(a[i]), end(a[i]), begin(m[i]));
		m[i][N + i] = 1;
    }

    auto [rank, det] = GaussElimination(m, N, true);
    if (rank != N) return {};

    Matrix<mint> b(N);
    for (int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			b[i][j] = m[i][j+N];
		}
    }
    return b;
}