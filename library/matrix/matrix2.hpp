#pragma once
#include <cmath>
#include <type_traits>
template <class T,T(*zero)(),T(*one)(),T(*Add)(T,T),T(*Mul)(T,T)>
struct Matrix {
    vector<vector<T>> A;
    // 行が0本の場合の列数を保持する。
    int empty_width = 0;

    Matrix() = default;
    Matrix(int n, int m) : A(n, vector<T>(m, zero())), empty_width(m) {}
    Matrix(int n) : Matrix(n, n) {};

    int H() const { return A.size(); }

    int W() const { return A.empty() ? empty_width : A[0].size(); }

    int size() const { return A.size(); }

    inline const vector<T> &operator[](int k) const { return A[k]; }

    inline vector<T> &operator[](int k) { return A[k]; }

    static Matrix I(int n) {
        Matrix mat(n);
        for (int i = 0; i < n; i++) {
            for(int j = 0; j < n;j++) {
                if(i == j) mat[i][j] = one();
                else mat[i][j] = zero();
            }
        }
        return (mat);
    }

    Matrix &operator+=(const Matrix &B) {
        int n = H(), m = W();
        assert(n == B.H() && m == B.W());
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) (*this)[i][j] = Add((*this)[i][j],B[i][j]);
        return (*this);
    }

    /*Matrix &operator-=(const Matrix &B) {
        int n = H(), m = W();
        assert(n == B.H() && m == B.W());
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) (*this)[i][j] -= B[i][j];
        return (*this);
    }*/

    Matrix &operator*=(const Matrix &B) {
        int n = H(), m = B.W(), p = W();
        assert(p == B.H());
        vector<vector<T>> C(n, vector<T>(m, zero()));
        for (int i = 0; i < n; i++)
            for (int k = 0; k < p; k++)
                for (int j = 0; j < m; j++) C[i][j] = Add(C[i][j],Mul((*this)[i][k],B[k][j]));
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

    //Matrix operator-(const Matrix &B) const { return (Matrix(*this) -= B); }

    Matrix operator*(const Matrix &B) const { return (Matrix(*this) *= B); }

    Matrix operator^(const long long k) const { return (Matrix(*this) ^= k); }

    bool operator==(const Matrix &B) const {
        if(H() != B.H() || W() != B.W()) return false;
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

    T determinant() const {
        Matrix B(*this);
        assert(H() == W());
        T ret = 1;
        for (int i = 0; i < H(); i++) {
            int idx = -1;
            for (int j = i; j < W(); j++) {
                if (B[j][i] != 0) {
                    if constexpr (is_floating_point<T>::value) {
                        if(idx == -1 || std::abs(B[j][i]) > std::abs(B[idx][i])) idx = j;
                    } else {
                        idx = j;
                        break;
                    }
                }
            }
            if (idx == -1) return 0;
            if (i != idx) {
                ret *= T(-1);
                swap(B[i], B[idx]);
            }
            ret *= B[i][i];
            T inv = T(1) / B[i][i];
            for (int j = 0; j < W(); j++) {
                B[i][j] *= inv;
            }
            for (int j = i + 1; j < H(); j++) {
                T a = B[j][i];
                if (a == 0) continue;
                for (int k = i; k < W(); k++) {
                    B[j][k] -= B[i][k] * a;
                }
            }
        }
    return ret;
    }
};

