#pragma once
#include "library/graph/graph-template.hpp"
// bridge:橋,articulation:関節点
template<typename G>
struct LowLink {
    int N;
    const G &g;
    vector<int> ord, low, articularion, parent;
    vector<pair<int,int>> bridge;
    vector<vector<int>> dfstree;
    LowLink(const G&_g):g(_g) {
        N = g.size();
        ord.resize(N,-1);
        low.resize(N,-1);
        parent.resize(N,-1);
        dfstree.resize(N);
        int k = 0;
        for(int i = 0;i < N;i++) {
            if(ord[i] == -1) k = dfs(i,k,-1);
        }
    }
    int dfs(int idx,int k,int par) {
        struct Frame {
            int v, par;
            size_t next = 0;
            int children = 0;
            bool skipped_parent = false, is_arti = false;
        };
        vector<Frame> st{{idx, par}};
        parent[idx] = par;
        low[idx] = (ord[idx] = k++);
        while(!st.empty()) {
            auto &frame = st.back();
            int v = frame.v;
            if(frame.next == g[v].size()) {
                bool is_arti = frame.par == -1 ? frame.children > 1 : frame.is_arti;
                if(is_arti) articularion.push_back(v);
                st.pop_back();
                if(!st.empty()) {
                    auto &p = st.back();
                    low[p.v] = min(low[p.v], low[v]);
                    p.is_arti |= p.par != -1 && low[v] >= ord[p.v];
                    if(ord[p.v] < low[v]) bridge.emplace_back(minmax(p.v, v));
                }
                continue;
            }
            int to = g[v][frame.next++];
            // 無向多重辺では、親へ戻る最初の1辺だけを除外する。
            if(to == frame.par && !exchange(frame.skipped_parent, true)) continue;
            if(ord[to] == -1) {
                frame.children++;
                dfstree[v].push_back(to);
                parent[to] = v;
                low[to] = (ord[to] = k++);
                st.push_back({to, v});
            } else low[v] = min(low[v], ord[to]);
        }
        return k;
    }
};
