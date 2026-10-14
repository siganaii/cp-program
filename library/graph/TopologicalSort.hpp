#pragma once
#include "library/graph/graph-template.hpp"
template<typename T>
vector<int>TopologicalSort(T &g) {
    int N = g.size();
    vector<int> marked(N),tmp(N),v;
    auto visit = [&](auto &f,int i) -> bool {
        if(tmp[i] == 1) return false;
        if(marked[i] == 0) {
            tmp[i] = 1;
            for(auto &e:g[i]) {
                if(f(f,e) == false) return false;
            }
            marked[i] = 1;
            v.push_back(i);
            tmp[i] = 0;
        }
        return true;
    };
    for(int i = 0;i < N;i++) {
        if(marked[i] == 0) {
            if(visit(visit,i) == false) return vector<int>();
        }
    }
    reverse(v.begin(),v.end());
    return v;
}