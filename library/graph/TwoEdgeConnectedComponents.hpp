#pragma once
#include "library/graph/graph-template.hpp"
#include "library/graph/LowLink.hpp"
//橋がないグラフに分解
//groups それぞれのグループの頂点が入る
//comp その点が属するグループが入る
//tree 橋となる辺だけを残した木になる
template<typename G>
struct TwoEdgeConnectedComponents {
    const G &g;
    LowLink<G> low;
    vector<int> comp;
    int k;
    vector<vector<int>> groups,tree;
    TwoEdgeConnectedComponents(const G &_g):g(_g),low(_g),comp(_g.size(),-1),k(0) {
        for(int i = 0;i < (int)g.size();i++) {
            if(comp[i] == -1) dfs(i,-1);
        }
        groups.resize(k);
        tree.resize(k);
        for(int i = 0;i < (int)g.size();i++) {
            groups[comp[i]].push_back(i);
        }
        for(auto &e:low.bridge) {
            int u = comp[e.first],v = comp[e.second];
            tree[u].emplace_back(v);
            tree[v].emplace_back(u);
        }
    };
    int operator[](const int &k) const{return comp[k];}
    void dfs(int i,int p) {
        if(p >= 0 && low.ord[p] >= low.low[i]) comp[i] = comp[p];
        else comp[i] = k++;
        for(auto &d:g[i]) if(comp[d] == -1) dfs(d,i);
    }
};