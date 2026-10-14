#pragma once
#include "library/graph/graph-template.hpp"
template <typename G = vector<vector<int>>>
struct HLD {
private:
    void dfs_sz(int root) {
        vector<int> vs;
        vs.reserve(g.size());
        par[root] = -1;
        depth[root] = 0;
        vector<int> st;
        st.push_back(root);
        while (!st.empty()) {
            int cur = st.back();
            st.pop_back();
            vs.push_back(cur);
            for (int dst : g[cur]) {
                if (dst == par[cur]) continue;
                par[dst] = cur;
                depth[dst] = depth[cur] + 1;
                st.push_back(dst);
            }
        }
        for (int i = int(vs.size()) - 1; i >= 0; i--) {
            int cur = vs[i];
            size[cur] = 1;
            int heavy_pos = -1;
            int heavy_size = -1;
            for (int j = 0; j < int(g[cur].size()); j++) {
                int dst = g[cur][j];
                if (dst == par[cur]) continue;
                size[cur] += size[dst];
                if (size[dst] > heavy_size) {
                    heavy_size = size[dst];
                    heavy_pos = j;
                }
            }
            if (heavy_pos != -1) {
                swap(g[cur][0], g[cur][heavy_pos]);
            }
        }
    }
    void dfs_hld(int root) {
        vector<pair<int, int>> st;
        st.emplace_back(root, 0);
        nxt[root] = root;
        while (!st.empty()) {
            auto [cur, state] = st.back();
            st.pop_back();
            if (state == 0) {
                ord[id] = cur;
                down[cur] = id++;
                st.emplace_back(cur, 1);
                for (int i = int(g[cur].size()) - 1; i >= 0; i--) {
                    int dst = g[cur][i];
                    if (dst == par[cur]) continue;
                    if (dst == g[cur][0]) {
                        nxt[dst] = nxt[cur];
                    } else {
                        nxt[dst] = dst;
                    }
                    st.emplace_back(dst, 0);
                }
            } else {
                up[cur] = id;
            }
        }
    }
public:
    // [u, v)
    vector<pair<int, int>> ascend(int u, int v) const {
        vector<pair<int, int>> res;
        while (nxt[u] != nxt[v]) {
            res.emplace_back(down[u], down[nxt[u]]);
            u = par[nxt[u]];
        }
        if (u != v) res.emplace_back(down[u], down[v] + 1);
        return res;
    }
    // (u, v]
    vector<pair<int, int>> descend(int u, int v) const {
        if (u == v) return {};
        if (nxt[u] == nxt[v]) return {{down[u] + 1, down[v]}};
        auto res = descend(u, par[nxt[v]]);
        res.emplace_back(down[nxt[v]], down[v]);
        return res;
    }
    G g;
    int root, id;
    vector<int> size, depth, down, up, ord, nxt, par;
    HLD() = default;
    HLD(G& _g, int _root = 0)
        : g(_g),
          root(_root),
          id(0),
          size(g.size(), 0),
          depth(g.size(), 0),
          down(g.size(), -1),
          up(g.size(), -1),
          ord(g.size(), 0),
          nxt(g.size(), root),
          par(g.size(), -1) {
        dfs_sz(root);
        dfs_hld(root);
    }
    void build(int _root) {
        root = _root;
        id = 0;
        fill(size.begin(), size.end(), 0);
        fill(depth.begin(), depth.end(), 0);
        fill(down.begin(), down.end(), -1);
        fill(up.begin(), up.end(), -1);
        fill(ord.begin(), ord.end(), 0);
        fill(nxt.begin(), nxt.end(), root);
        fill(par.begin(), par.end(), -1);
        dfs_sz(root);
        dfs_hld(root);
    }
    pair<int, int> idx(int i) const {
        return make_pair(down[i], up[i]);
    }
    template <typename F>
    void path_query(int u, int v, bool vertex, const F& f) {
        int l = lca(u, v);
        for (auto&& [a, b] : ascend(u, l)) {
            int s = a + 1, t = b;
            if (s > t) f(t, s);
            else f(s, t);
        }
        if (vertex) f(down[l], down[l] + 1);
        for (auto&& [a, b] : descend(l, v)) {
            int s = a, t = b + 1;
            if (s > t) f(t, s);
            else f(s, t);
        }
    }
    template <typename F>
    void path_noncommutative_query(int u, int v, bool vertex, const F& f) {
        int l = lca(u, v);
        for (auto&& [a, b] : ascend(u, l)) {
            f(a + 1, b);
        }
        if (vertex) f(down[l], down[l] + 1);
        for (auto&& [a, b] : descend(l, v)) {
            f(a, b + 1);
        }
    }
    template <typename F>
    void subtree_query(int u, bool vertex, const F& f) {
        f(down[u] + int(!vertex), up[u]);
    }
    int lca(int a, int b) const {
        while (nxt[a] != nxt[b]) {
            if (down[a] < down[b]) swap(a, b);
            a = par[nxt[a]];
        }
        return depth[a] < depth[b] ? a : b;
    }
    int dist(int a, int b) const {
        return depth[a] + depth[b] - depth[lca(a, b)] * 2;
    }
    int kth_ancestor(int u, int k) const {
        if (k < 0) return -1;
        while (u >= 0) {
            int h = nxt[u];
            if (down[u] - k >= down[h]) {
                return ord[down[u] - k];
            }
            k -= down[u] - down[h] + 1;
            u = par[h];
        }
        return -1;
    }
    int next(int s, int t) const {
        assert(s != t);
        assert(0 <= s && s < int(g.size()));
        assert(0 <= t && t < int(g.size()));
        if (depth[s] >= depth[t]) return par[s];
        int u = kth_ancestor(t, depth[t] - depth[s] - 1);
        return par[u] == s ? u : par[s];
    }
    // s - t 間のパス上の頂点のうち s から距離 d の頂点
    // dist(s, t) < d のとき -1
    int jump(int s, int t, int d) const {
        int lc = lca(s, t);
        int d1 = depth[s] - depth[lc];
        if (d <= d1) {
            return kth_ancestor(s, d);
        }
        int d2 = d1 + depth[t] - depth[lc];
        if (d <= d2) {
            return kth_ancestor(t, d2 - d);
        }
        return -1;
    }
    vector<int> path(int s, int t) const {
        vector<int> pre, suf;
        while (depth[s] > depth[t]) {
            pre.emplace_back(s);
            s = par[s];
        }
        while (depth[s] < depth[t]) {
            suf.emplace_back(t);
            t = par[t];
        }
        while (s != t) {
            pre.emplace_back(s);
            suf.emplace_back(t);
            s = par[s];
            t = par[t];
        }
        pre.push_back(s);
        reverse(begin(suf), end(suf));
        copy(begin(suf), end(suf), back_inserter(pre));
        return pre;
    }
    // u is in subtree of v
    bool in_subtree(int u, int v) const {
        return down[v] <= down[u] && down[u] < up[v];
    }
};