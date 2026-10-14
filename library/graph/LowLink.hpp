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
        parent[idx] = par;
        low[idx] = (ord[idx] = k++);
        int cnt = 0;
        bool is_arti = false,flg = false;
        for(auto &to:g[idx]) {
            if(ord[to] == -1) {
                cnt++;
                dfstree[idx].emplace_back(to);
                k = dfs(to,k,idx);
                low[idx] = min(low[idx],low[to]);
                is_arti |= (par != -1) && (low[to] >= ord[idx]);
                if(ord[idx] < low[to]) bridge.emplace_back(minmax(idx,(int)to));
            }
            else if(to != par || exchange(flg,true)) low[idx] = min(low[idx],ord[to]);
        }
        is_arti |= par == -1 && cnt > 1;
        if(is_arti) articularion.emplace_back(idx);
        return k;
    }
};