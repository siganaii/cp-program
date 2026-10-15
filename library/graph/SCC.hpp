#pragma once
#include "library/graph/graph-template.hpp"
// scc[k] kが属するグループのindex
// belong(i) i番目のグループのリスト
// dag グループを縮約した時のグラフ
template <typename G>
struct SCC {
    private:
        const G &g;
        vector<vector<int>> rg;
        vector<int> comp,order;
        vector<char> used;
        vector<vector<int>> blng;
    public:
        vector<vector<int>> dag;
        SCC(G &_g):g(_g),used(g.size(),0){ build(); }
        int operator[](int k){return comp[k]; }
        vector<int> &belong(int i){return blng[i]; }
    private:
        void dfs(int idx) {
            if(used[idx]) return;
            used[idx] = true;
            vector<pair<int, size_t>> st{{idx, 0}};
            while(!st.empty()) {
                int v = st.back().first;
                auto &next = st.back().second;
                if(next == g[v].size()) {
                    order.push_back(v);
                    st.pop_back();
                    continue;
                }
                int to = g[v][next++];
                if(!used[to]) {
                    used[to] = true;
                    st.emplace_back(to, 0);
                }
            }
        }
        void rdfs(int idx,int cnt) {
            if(comp[idx] != -1) return;
            comp[idx] = cnt;
            vector<int> st{idx};
            while(!st.empty()) {
                int v = st.back();
                st.pop_back();
                for(int to:rg[v]) if(comp[to] == -1) {
                    comp[to] = cnt;
                    st.push_back(to);
                }
            }
        }
        void build() {
            for(int i = 0;i < (int)g.size(); i++) dfs(i);
            reverse(begin(order),end(order));
            used.clear();
            used.shrink_to_fit();
            comp.resize(g.size(),-1);
            rg.resize(g.size());
            for(int i = 0;i < (int)g.size(); i++) {
                for(auto e:g[i]) {
                    rg[(int)e].emplace_back(i);
                }
            }
            int ptr = 0;
            for(int i:order) if(comp[i] == -1) rdfs(i,ptr),ptr++;
            rg.clear();
            rg.shrink_to_fit();
            order.clear();
            order.shrink_to_fit();
            dag.resize(ptr);
            blng.resize(ptr);
            for(int i = 0;i < (int)g.size();i++) {
                blng[comp[i]].push_back(i);
                for(auto &to:g[i]) {
                    int x = comp[i],y = comp[to];
                    if(x == y)continue;
                    dag[x].push_back(y);
                }
            }
        }
};
