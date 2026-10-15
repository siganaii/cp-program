#pragma once
#include "library/modulo/binomial.hpp"
#include "library/matrix/matrix-tree.hpp"
template <typename mint>
struct Count_Euler_Circuit {
    int n;
    vector<int> outdeg,indeg;
    vector<vector<int>> G;
    Count_Euler_Circuit() = default;
    Count_Euler_Circuit(int _n) : n(_n), G(n,vector<int>(n)), outdeg(n), indeg(n) {assert(n > 0);}
    void add(int i, int j, const int& x) {
        G[i][j] += x;
        outdeg[i] += x;
        indeg[j] += x;
    }
    mint count(Binomial<mint> &C) {
        for(int i = 0;i < n;++i) {
            if(indeg[i] != outdeg[i]) return 0;
        }
        vector<int> seen(n,-1);
        for(int i = 0;i < n;++i) {
            for(int j = 0;j < n;++j) {
                if(G[i][j] != 0) {
                    seen[i] = 1;
                    seen[j] = 1;
                }
            }
        }
        int idx = 0;
        for(int i = 0;i < n;++i) if(seen[i] == 1) seen[i] = idx++;
        // 辺がない場合は、空の巡回路を1通りと数える。
        if(idx == 0) return 1;
        MatrixTree<mint,true> MT(idx);
        for (int i = 0; i < n; i++) if (seen[i] != -1) {
            for (int j = 0; j < n; j++) if (seen[j] != -1) {
                MT.add(seen[i],seen[j],G[i][j]);
            }
        }
        mint ret = 1;
        for(int i = 0;i < n;++i) {
            if(seen[i] != -1) ret *= C.fac(outdeg[i] - 1);
        }
        return ret * MT.calc(-1,false);
    }
};
