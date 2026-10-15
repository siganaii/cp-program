#pragma once
#include "library/graph/graph-template.hpp"
// 自己ループは無視する。無向グラフの多重辺による長さ2の閉路は検出する。
template <typename G>
vector<pair<int,int>> cycleditection(const G &g,bool directed = true) {
    // 0: 未訪問、1: 探索中、2: 探索終了
    vector<int> state(g.size()), parent(g.size(), -1);
    vector<size_t> next(g.size());
    vector<bool> skipped_parent(g.size());
    vector<int> st;
    for(int i = 0;i < (int)g.size();i++) {
        if(state[i]) continue;
        state[i] = 1;
        st.push_back(i);
        while(!st.empty()) {
            int cur = st.back();
            if(next[cur] == g[cur].size()) {
                state[cur] = 2;
                st.pop_back();
                continue;
            }
            int to = g[cur][next[cur]++];
            if(to == cur) continue;
            if(!directed && to == parent[cur] && !skipped_parent[cur]) {
                skipped_parent[cur] = true;
                continue;
            }
            if(state[to] == 1) {
                vector<pair<int,int>> cycle{{cur, to}};
                for(int v = cur; v != to; v = parent[v])
                    cycle.emplace_back(parent[v], v);
                reverse(cycle.begin(), cycle.end());
                return cycle;
            }
            if(state[to] == 0) {
                parent[to] = cur;
                state[to] = 1;
                st.push_back(to);
            }
        }
    }
    return {};
}
