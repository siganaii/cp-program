#pragma once
#include "library/graph/graph-template.hpp"
//自己ループは判定しない
template <typename G>
vector<pair<int,int>> cycleditection(const G &g,bool directed = true) {
    vector<int> pidx(g.size(),-1),vis(g.size());
    vector<pair<int,int>> cycle;
    int fin = 0;
    auto dfs = [&](auto &rec,int cur,int pval,int par) -> int {
        pidx[cur] = pval;
        vis[cur] = 1;
        for(auto &dst:g[cur]) {
            if(fin) return -1;
            if(!directed && dst == par) continue;
            if(pidx[dst] == pidx[cur]) {
                cycle.emplace_back(cur,dst);
                return dst;
            }
            if(vis[dst]) continue;
            int nx = rec(rec,dst,pval,cur);
            if(nx != -1) {
                cycle.emplace_back(cur,dst);
                if(cur == nx) {
                    fin = 1;
                    return -1;
                }
                return nx;
            }
        } 
        pidx[cur] = -1;
        return -1;
    };
    for(int i = 0;i < (int)g.size();i++) {
        if(vis[i])continue;
        dfs(dfs,i,i,-1);
    }
    if(fin) {
        reverse(begin(cycle),end(cycle));
        return cycle;
    }
    return vector<pair<int,int>> {};
}
