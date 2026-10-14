#pragma once
#include "library/graph/graph-template.hpp"
template<typename T>
struct TreeDiameter_Weighted {
    private:
    Wgraph<T> &g;
    vector<int> to;
    pair<T,int> dfs(int cur,int par) {
        pair<T,int> ret(0,cur);
        for(auto &ch : g[cur]) {
            if(ch == par) continue;
            pair<T,int> cost = dfs(ch,cur);
            cost.first += ch.cost;
            if(ret < cost) {
                ret = cost;
                to[cur] = ch;
            }
        }
        return ret;
    }
    public:
    vector<int> path;
    T distance;
    TreeDiameter_Weighted(Wgraph<T> &_g):g(_g){build();}
    void build() {
        to.assign(g.size(),-1);
        pair<T,int> p = dfs(0,-1);
        pair<T,int> q = dfs(p.second,-1);
        distance = q.first;
        int cur = p.second;
        path.emplace_back(cur);
        while(cur != q.second) {
            path.emplace_back(to[cur]);
            cur = to[cur];
        }
    }
};
