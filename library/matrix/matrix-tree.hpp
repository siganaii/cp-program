#pragma once
#include "library/matrix/matrix.hpp"
//行列木定理で全域木の個数を数える
//有向木の場合はis_directed = true
template <typename T,bool is_directed = false>
struct MatrixTree {
    int n;
    vector<vector<T>> M1,M2;
    MatrixTree() = default;
    MatrixTree(int _n) : n(_n), M1(n,vector<T>(n)),M2(n,vector<T>(n)) { assert(n > 0); }
    void add(int i, int j, const T& x) {
        if(i == j) return;
        M1[i][i] += x;
        M1[i][j] -= x;
        M2[j][j] += x;
        M2[j][i] -= x;
    }
    //有向木の場合、
    //根から全ての点に行ける場合は、is_start = true
    //全ての点から根に行ける場合は、is_start = false
    T calc(int root = -1,bool is_start = true) {
        if(root == -1) root = 0;
        Matrix<T> M(n-1);
        for(int i = 0;i < n;++i) {
            if(i == root) continue;
            for(int j = 0;j < n;++j) {
                if(j == root) continue;
                if(!is_directed) M[i - (i > root)][j - (j > root)] = M1[i][j] + M2[i][j];
                else if(is_start) M[i - (i > root)][j - (j > root)] = M2[i][j];
                else M[i - (i > root)][j - (j > root)] = M1[i][j];
            }
        }
        return M.determinant();
    }
};
