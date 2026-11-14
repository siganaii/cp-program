#pragma once
#include "library/graph/graph-template.hpp"
// scc[k] kが属するグループのindex（トポロジカル順）
// belong(i) i番目のグループのリスト（頂点番号の昇順）
// dag グループを縮約した時のグラフ（多重辺は残す）
// 非再帰のTarjan法。計算量 O(V + E)。逆向きのグラフは構築しない。
// SCC(g, false) はdagの構築を省略する。size() は強連結成分の個数。
template <typename G>
struct SCC {
    private:
        const G &g;
        vector<int> comp;
        vector<vector<int>> blng;
    public:
        vector<vector<int>> dag;
        SCC(G &_g, bool build_dag = true):g(_g){ build(build_dag); }
        int size() const { return blng.size(); }
        int operator[](int k){return comp[k]; }
        vector<int> &belong(int i){return blng[i]; }
    private:
        void build(bool build_dag) {
            const int n = g.size();
            comp.resize(n);
            vector<size_t> start(n + 1);
            for (int v = 0; v < n; ++v) start[v + 1] = start[v] + g[v].size();
            vector<int> edges(start.back());
            for (int v = 0; v < n; ++v) {
                size_t next = start[v];
                for (const auto &to : g[v]) edges[next++] = to;
            }
            int groups = 0;
            {
                vector<int> ord(n, -1), active;
                struct Frame { int v, low; size_t next; };
                vector<Frame> st;
                active.reserve(n);
                st.reserve(n);
                int timer = 0;
                for (int root = 0; root < n; ++root) {
                    if (ord[root] != -1) continue;
                    ord[root] = timer++;
                    active.push_back(root);
                    st.push_back({root, ord[root], start[root]});
                    while (!st.empty()) {
                        int v = st.back().v;
                        if (st.back().next < start[v + 1]) {
                            int to = edges[st.back().next++];
                            if (ord[to] == -1) {
                                ord[to] = timer++;
                                active.push_back(to);
                                st.push_back({to, ord[to], start[to]});
                            } else {
                                st.back().low = min(st.back().low, ord[to]);
                            }
                            continue;
                        }
                        int low = st.back().low;
                        if (low == ord[v]) {
                            while (true) {
                                int u = active.back();
                                active.pop_back();
                                comp[u] = groups;
                                ord[u] = n; // 確定した成分への辺はlowに影響させない。
                                if (u == v) break;
                            }
                            ++groups;
                        }
                        st.pop_back();
                        if (!st.empty()) {
                            st.back().low = min(st.back().low, low);
                        }
                    }
                }
            }
            if (build_dag) dag.resize(groups);
            blng.resize(groups);
            if (groups == 0) return;
            if (groups == 1) {
                blng[0].resize(n);
                for (int v = 0; v < n; ++v) blng[0][v] = v;
                return;
            }
            vector<int> sizes(groups);
            for (int v = 0; v < n; ++v) {
                comp[v] = groups - 1 - comp[v];
                ++sizes[comp[v]];
            }
            for (int c = 0; c < groups; ++c) blng[c].reserve(sizes[c]);
            for (int v = 0; v < n; ++v) blng[comp[v]].push_back(v);
            if (!build_dag) return;
            vector<size_t> degrees(groups);
            for (int v = 0; v < n; ++v) {
                int c = comp[v];
                for (size_t i = start[v]; i < start[v + 1]; ++i) {
                    if (c != comp[edges[i]]) ++degrees[c];
                }
            }
            for (int c = 0; c < groups; ++c) dag[c].reserve(degrees[c]);
            for (int v = 0; v < n; ++v) {
                int c = comp[v];
                for (size_t i = start[v]; i < start[v + 1]; ++i) {
                    int to = comp[edges[i]];
                    if (c != to) dag[c].push_back(to);
                }
            }
        }
};
