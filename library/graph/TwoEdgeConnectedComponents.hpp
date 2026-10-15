#pragma once
#include "library/graph/graph-template.hpp"
#include "library/graph/LowLink.hpp"
//橋がないグラフに分解
//groups それぞれのグループの頂点が入る
//comp その点が属するグループが入る
//tree 橋となる辺だけを残した森になる
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
        vector<pair<int,int>> st{{i, p}};
        while(!st.empty()) {
            auto [v, par] = st.back();
            st.pop_back();
            if(par >= 0 && low.ord[par] >= low.low[v]) comp[v] = comp[par];
            else comp[v] = k++;
            // LowLink と同じ DFS 木をたどり、橋を越えるときに成分を分ける。
            for(auto it = low.dfstree[v].rbegin(); it != low.dfstree[v].rend(); ++it)
                st.emplace_back(*it, v);
        }
    }
};
