#pragma once
#include "library/graph/graph-template.hpp"
template<typename T>
vector<int>TopologicalSort(T &g) {
    int N = g.size();
    // 0: 未訪問、1: 探索中、2: 探索終了
    vector<int> state(N),v;
    vector<pair<int, size_t>> st;
    for(int i = 0;i < N;i++) {
        if(state[i]) continue;
        state[i] = 1;
        st.emplace_back(i, 0);
        while(!st.empty()) {
            int cur = st.back().first;
            auto &next = st.back().second;
            if(next == g[cur].size()) {
                state[cur] = 2;
                v.push_back(cur);
                st.pop_back();
                continue;
            }
            int to = g[cur][next++];
            if(state[to] == 1) return {};
            if(state[to] == 0) {
                state[to] = 1;
                st.emplace_back(to, 0);
            }
        }
    }
    reverse(v.begin(),v.end());
    return v;
}
