#pragma once
#include "library/graph/graph-template.hpp"
// 辺の重みは非負で、経路の重みの和が T に収まることを前提とする。
template<typename T>
struct TreeDiameter_Weighted {
    private:
    Wgraph<T> &g;
    vector<int> to;
    pair<T,int> dfs(int cur,int par) {
        // 親を記録しながら反復探索し、最遠頂点を求める。
        pair<T,int> ret(0,cur);
        vector<pair<int,T>> stack{{cur, 0}};
        to[cur] = par;
        while (!stack.empty()) {
            auto [u, dist] = stack.back();
            stack.pop_back();
            ret = max(ret, make_pair(dist, u));
            for (auto &ch : g[u]) {
                if (ch == to[u]) continue;
                to[ch] = u;
                stack.emplace_back(ch, dist + ch.cost);
            }
        }
        return ret;
    }
    public:
    vector<int> path;
    T distance;
    TreeDiameter_Weighted(Wgraph<T> &_g):g(_g){build();}
    void build() {
        path.clear();
        to.assign(g.size(),-1);
        distance = 0;
        if (g.empty()) return;
        pair<T,int> p = dfs(0,-1);
        pair<T,int> q = dfs(p.second,-1);
        distance = q.first;
        for (int cur = q.second; cur != -1; cur = to[cur]) path.push_back(cur);
        reverse(path.begin(), path.end());
    }
};
