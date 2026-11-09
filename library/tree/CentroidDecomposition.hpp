#pragma once
// root：根 tree：重心分解の探索する際に根となる順番 sub：探索する時の頂点数
template <typename G>
struct CentroidDecomposition {
    const G &g;
    vector<int> sub;
    vector<bool> v;
    vector<vector<int>> tree;
    int root = -1;
    CentroidDecomposition(const G &g_, int isbuild = true) : g(g_) {
        sub.resize(g.size(), 0);
        v.resize(g.size(), false);
        if (isbuild) build();
    }
    void build() {
        tree.assign(g.size(), {});
        sub.assign(g.size(), 0);
        v.assign(g.size(), false);
        root = g.empty() ? -1 : build_dfs(0);
    }
    int get_size(int cur, int par) {
        vector<pair<int, int>> order{{cur, par}};
        for (int i = 0; i < (int)order.size(); ++i) {
            auto [u, p] = order[i];
            sub[u] = 1;
            for (int dst : g[u]) {
                if (dst != p && !v[dst]) order.emplace_back(dst, u);
            }
        }
        for (int i = (int)order.size() - 1; i > 0; --i) {
            auto [u, p] = order[i];
            sub[p] += sub[u];
        }
        return sub[cur];
    }
    int get_centroid(int cur, int par, int mid) {
        while (true) {
            int next = -1;
            for (int dst : g[cur]) {
                if (dst != par && !v[dst] && sub[dst] > mid) {
                    next = dst;
                    break;
                }
            }
            if (next == -1) return cur;
            par = cur;
            cur = next;
        }
    }
    int build_dfs(int cur) {
        // 重心で成分が半分以下になるため、この再帰の深さは O(log N)。
        int centroid = get_centroid(cur, -1, get_size(cur, -1) / 2);
        v[centroid] = true;
        for (auto &dst : g[centroid]) {
            if (!v[dst]) {
                int nxt = build_dfs(dst);
                if (centroid != nxt) tree[centroid].emplace_back(nxt);
            }
        }
        v[centroid] = false;
        return centroid;
    }
};
