#pragma once
#include <cmath>
#include <type_traits>
#include "library/matrix/elimination-utils.hpp"
template <class T>
struct Matrix {
    vector<vector<T>> A;
    // 行が0本の場合の列数を保持する。
    int empty_width = 0;
    Matrix() = default;
    Matrix(int n, int m) : A(n, vector<T>(m, T())), empty_width(m) {}
    Matrix(int n) : Matrix(n, n) {};
    int H() const { return A.size(); }
    int W() const { return A.empty() ? empty_width : A[0].size(); }
    int size() const { return A.size(); }
    inline const vector<T> &operator[](int k) const { return A[k]; }
    inline vector<T> &operator[](int k) { return A[k]; }
    static Matrix I(int n) {
        Matrix mat(n);
        for (int i = 0; i < n; i++) mat[i][i] = 1;
        return (mat);
    }
    Matrix &operator+=(const Matrix &B) {
        int n = H(), m = W();
        assert(n == B.H() && m == B.W());
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) (*this)[i][j] += B[i][j];
        return (*this);
    }
    Matrix &operator-=(const Matrix &B) {
        int n = H(), m = W();
        assert(n == B.H() && m == B.W());
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) (*this)[i][j] -= B[i][j];
        return (*this);
    }
    Matrix &operator*=(const Matrix &B) {
        int n = H(), m = B.W(), p = W();
        assert(p == B.H());
        vector<vector<T>> C(n, vector<T>(m, T{}));
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < p; k++) {
                if((*this)[i][k] == 0) continue;
                for (int j = 0; j < m; j++) C[i][j] += (*this)[i][k] * B[k][j];
            }
        }
        A.swap(C);
        empty_width = m;
        return (*this);
    }
    Matrix &operator^=(long long k) {
        assert(H() == W());
        Matrix B = Matrix::I(H());
        while (k > 0) {
            if (k & 1) B *= *this;
            k >>= 1LL;
            if(k) *this *= *this;
        }
        empty_width = B.W();
        A.swap(B.A);
        return (*this);
    }
    Matrix operator+(const Matrix &B) const { return (Matrix(*this) += B); }
    Matrix operator-(const Matrix &B) const { return (Matrix(*this) -= B); }
    Matrix operator*(const Matrix &B) const { return (Matrix(*this) *= B); }
    Matrix operator^(const long long k) const { return (Matrix(*this) ^= k); }
    bool operator==(const Matrix &B) const {
        if (H() != B.H() || W() != B.W()) return false;
        for (int i = 0; i < H(); i++) for (int j = 0; j < W(); j++) {
            if (A[i][j] != B[i][j]) return false;
        }
        return true;
    }
    bool operator!=(const Matrix &B) const {
        return !(*this == B);
    }
    friend ostream &operator<<(ostream &os, const Matrix &p) {
        int n = p.H(), m = p.W();
        for (int i = 0; i < n; i++) {
            os << (i ? "   " : "") << "[";
            for (int j = 0; j < m; j++) {
                os << p[i][j] << (j + 1 == m ? "]\n" : ",");
            }
        }
        return (os);
    }
    template <typename T2, typename enable_if<is_floating_point<T2>::value>::type * = nullptr>
    static int choose_pivot(const Matrix<T2> &mtr, int h, int c) noexcept {
        int piv = -1;
        for (int j = h; j < mtr.H(); j++) {
            if (mtr[j][c] != T2(0) && (piv < 0 || std::abs(mtr[j][c]) > std::abs(mtr[piv][c])))
                piv = j;
        }
        return piv;
    }
    template <typename T2, typename enable_if<!is_floating_point<T2>::value>::type * = nullptr>
    static int choose_pivot(const Matrix<T2> &mtr, int h, int c) noexcept {
        for (int j = h; j < mtr.H(); j++) {
            if (mtr[j][c] != T2()) return j;
        }
        return -1;
    }
    T determinant() const {
        assert(H() == W());
        Matrix B(*this);
        int n = H();
        T ret = T(1);
        matrix_internal::EliminationTolerance<T> tolerance(B, n, n);
        for (int i = 0; i < n; i++) {
            int idx = choose_pivot(B, i, i);
            if (idx == -1) return T(0);
            if (i != idx) {
                ret *= T(-1);
                swap(B[i], B[idx]);
                tolerance.swap_rows(i, idx);
            }
            ret *= B[i][i];
            T inv = T(1) / B[i][i];
            T inv_scale = tolerance.quotient_scale(inv, B[i][i], T(0), i, i);
            for (int j = i; j < n; j++) {
                tolerance.multiply(B[i][j], i, j, inv, inv_scale);
            }
            for (int j = i + 1; j < n; j++) {
                T a = B[j][i];
                T a_scale = tolerance.entry_scale(j, i);
                if (a == T(0)) continue;
                for (int k = i; k < n; k++) {
                    tolerance.subtract_product(B[j][k], B[i][k], a, a_scale, j, k, i);
                }
            }
        }
        return ret;
    }
};
