#pragma once
#include "library/matrix/count_Euler_circuit.hpp"
//オイラー路の数え上げ
//O(N^3 + M)
template<class mint>
struct Count_Eulerian_Trail {
    int n;
    Count_Euler_Circuit<mint> CE;
    Count_Eulerian_Trail(int _n) : n(_n), CE(n) {assert(n > 0);}
    void add(int i,int j,int x) {
        CE.add(i,j,x);
    }
    mint count(Binomial<mint> &C) {
        bool ok = 1;
        if(n == 0) return 1;
        int tot = 0;
        for(int i = 0;i < n;++i) {
            tot += CE.outdeg[i];
        }
        int st = -1, ed = -1;
        for(int i = 0;i < n;++i) if(CE.outdeg[i] != CE.indeg[i]) {
            if(abs(CE.outdeg[i] - CE.indeg[i]) > 1) {
                ok = 0;
                break;
            }
            if(CE.outdeg[i] > CE.indeg[i]) {
                if(ed != -1) ok = 0;
                ed = i;
            }
            else {
                if(st != -1) ok = 0;
                st = i;
            }
        }
        if(!ok) return 0;
        if((st == -1) ^ (ed == -1)) return 0;
        if(st == -1) {
            mint tmp = CE.count(C);
            return tmp * tot;
        }
        CE.add(st,ed,1);
        mint ret = CE.count(C);
        CE.add(st,ed,-1);
        return ret;
    }
};