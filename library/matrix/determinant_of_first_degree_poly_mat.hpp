#pragma once
#include "library/matrix/CharateristicPolynominal.hpp"
//det(M_0 + M_1 x)
//O(n^3)
template <class mint>
vector<mint> determinant_of_first_degree_poly_mat(Matrix<mint> M0,Matrix<mint> M1) {
    const int N = M0.size();
    assert(M0.H() == M0.W() && M1.H() == M1.W() && M0.H() == M1.H());
    int multiply_by_x = 0;
    mint detAdetBinv = 1;
    for (int p = 0;p < N;++p) {
        int pivot = -1;
        for (int row = p; row < N; ++row) {
            if (M1[row][p] != mint()) {
                pivot = row;
                break;
            }
        }
        if(pivot < 0) {
            ++multiply_by_x;
            if (multiply_by_x > N) return vector<mint>(N + 1);
            for (int row = 0; row < p; ++row) {
                mint v = M1[row][p];
                M1[row][p] = 0;
                for (int i = 0; i < N; ++i) M0[i][p] -= v * M0[i][row];
            }
            for (int i = 0; i < N; ++i) swap(M0[i][p], M1[i][p]);
            --p;
            continue;
        }
        if(pivot != p) {
            M1[pivot].swap(M1[p]);
            M0[pivot].swap(M0[p]);
            detAdetBinv *= -1;
        }
        mint v = M1[p][p], vinv = v.inverse();
        detAdetBinv *= v;
        for (int col = 0; col < N; ++col) {
            M0[p][col] *= vinv;
            M1[p][col] *= vinv;
        }
        for (int row = 0; row < N; ++row) {
            if (row == p) continue;
            mint v = M1[row][p];
            for (int col = 0; col < N; ++col) {
                M0[row][col] -= M0[p][col] * v;
                M1[row][col] -= M1[p][col] * v;
            }
        }
    }
    for (int i = 0;i < N;++i) {
        for (auto &x : M0[i]) x = -x;
    }
    auto poly = CharacteristicPolynomial(M0);
    if(N % 2) detAdetBinv *= -1;
    for (auto &x : poly) x *= detAdetBinv;
    poly.erase(poly.begin(), poly.begin() + multiply_by_x);
    poly.resize(N + 1);
    return poly;
}