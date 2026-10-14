#pragma once
#include "library/matrix/f2-matrix.hpp"

template <size_t H_MAX, size_t W_MAX>
vector<bitset<W_MAX>> LinearEquationMod2(
    const F2_Matrix<H_MAX, W_MAX>& A,
    const bitset<H_MAX>& b
) {
    int H = A.H;
    int W = A.W;
    F2_Matrix<H_MAX, W_MAX + 1> C(H, W + 1);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            C[i][j] = A[i][j];
        }
        C[i][W] = b[i];
    }
    vector<int> pivot(W, -1);
    int rank = 0;
    for (int col = 0; col < W; col++) {
        int idx = -1;
        for (int row = rank; row < H; row++) {
            if (C[row][col]) {
                idx = row;
                break;
            }
        }
        if (idx == -1) continue;
        swap(C[idx], C[rank]);
        for (int row = 0; row < H; row++) {
            if (row != rank && C[row][col]) {
                C[row] ^= C[rank];
            }
        }
        pivot[col] = rank;
        rank++;
    }
    for (int row = rank; row < H; row++) {
        if (C[row][W]) {
            return {};
        }
    }
    vector<bitset<W_MAX>> res;
    bitset<W_MAX> v;
    for (int j = 0; j < W; j++) {
        if (pivot[j] != -1) {
            v[j] = C[pivot[j]][W];
        }
    }
    res.push_back(v);
    for (int free_col = 0; free_col < W; free_col++) {
        if (pivot[free_col] != -1) continue;
        bitset<W_MAX> x;
        x[free_col] = 1;
        for (int j = 0; j < W; j++) {
            if (pivot[j] != -1) {
                x[j] = C[pivot[j]][free_col];
            }
        }
        res.push_back(x);
    }
    return res;
}