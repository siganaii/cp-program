#pragma once
#include <cmath>
#include <type_traits>
template<typename T,int H,int W>
struct MatrixArray {
    using Array = array<array<T,W>,H>;
    Array A;
    MatrixArray():A() {
        for(int i = 0;i < H;i++) {
            for(int j = 0;j < W;j++) {
                (*this)[i][j] = T();
            }
        }
    }
    int height() const {return H;}
    int width() const {return W;}
    inline const array<T,W> &operator[](int k) const {return A[k];}
    inline array<T,W> &operator[](int k){return A[k];}
    static MatrixArray I() {
        assert(H == W);
        MatrixArray mat;
        for(int i = 0;i < H;i++) mat[i][i] = 1;
        return (mat);
    }
    MatrixArray &operator+=(const MatrixArray &B) {
        for(int i = 0;i < H;i++) {
            for(int j = 0;j < W;j++) {
                A[i][j] += B[i][j];
            }
        }
        return (*this);
    }
    MatrixArray &operator-=(const MatrixArray &B) {
        for(int i = 0;i < H;i++) {
            for(int j = 0;j < W;j++) {
                A[i][j] -= B[i][j];
            }
        }
        return (*this);
    }
    MatrixArray &operator*=(const MatrixArray &B) {
        assert(H == W);
        MatrixArray C;
        for(int i = 0;i < H;i++) {
            for(int k = 0;k < H;k++) {
                for(int j = 0;j < H;j++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        A.swap(C.A);
        return (*this);
    }
    MatrixArray &operator^=(long long k) {
        MatrixArray B = MatrixArray::I();
        while(k > 0) {
            if(k & 1) B *= *this;
            k >>= 1LL;
            if(k) *this *= *this;
        }
        A.swap(B.A);
        return (*this);
    }
    MatrixArray operator+(const MatrixArray &B) const{return (MatrixArray(*this) += B);}
    MatrixArray operator-(const MatrixArray &B) const{return (MatrixArray(*this) -= B);}
    MatrixArray operator*(const MatrixArray &B) const{return (MatrixArray(*this) *= B);}
    MatrixArray operator^(const long long k) const{return (MatrixArray(*this) ^= k);}
    bool operator==(const MatrixArray &B) const {
        for(int i = 0;i < H;i++) {
            for(int j = 0;j < W;j++) {
                if(A[i][j] != B[i][j]) return false;
            }
        }
        return true;
    }
    bool operator!=(const MatrixArray &B) const {
        for(int i = 0;i < H;i++) {
            for(int j = 0;j < W;j++) {
                if(A[i][j] != B[i][j]) return true;
            }
        }
        return false;
    }
    friend ostream &operator<<(ostream &os,const MatrixArray &p) {
        for(int i = 0;i < H;i++) {
            os << "[";
            for(int j = 0;j < W;j++) {
                os << p[i][j] << (j+1 == W ?"]\n":",");
            }
        }
        return(os);
    }
    static int choose_pivot(const MatrixArray &mat, int h, int c, int row_end = H) {
        int pivot = -1;
        for(int i = h; i < row_end; i++) {
            if(mat[i][c] == T(0)) continue;
            if constexpr (is_floating_point<T>::value) {
                if(pivot == -1 || std::abs(mat[i][c]) > std::abs(mat[pivot][c])) pivot = i;
            } else {
                return i;
            }
        }
        return pivot;
    }
    T determinant(int n = -1) {
        if(n == -1) n = H;
        MatrixArray B(*this);
        T ret = 1;
        for(int i = 0;i < n;i++) {
            int idx = choose_pivot(B, i, i, n);
            if(idx == -1) return 0;
            if(i != idx) {
                ret *= T(-1);
                swap(B[i],B[idx]);
            }
            ret *= B[i][i];
            T inv = T(1) / B[i][i];
            for(int j = 0;j < n;j++) {
                B[i][j] *= inv;
            }
            for(int j = i+1;j < n;j++) {
                T a = B[j][i];
                if(a == 0) continue;
                for(int k = i;k < n;k++) {
                    B[j][k] -= B[i][k] * a;
                }
            }
        }
        return(ret);
    }
};
