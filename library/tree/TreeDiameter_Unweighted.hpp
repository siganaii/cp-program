#pragma once
template <typename G>
struct TreeDiameter_Unweighted {
    private:
    G &g;
    vector<int> to;
    pair<int,int> dfs(int cur,int par) {
        // 親を記録しながら反復探索し、最遠頂点を求める。
        pair<int,int> ret(0,cur);
        vector<pair<int,int>> stack{{cur, 0}};
        to[cur] = par;
        while (!stack.empty()) {
            auto [u, dist] = stack.back();
            stack.pop_back();
            ret = max(ret, make_pair(dist, u));
            for (auto &ch : g[u]) {
                if (ch == to[u]) continue;
                to[ch] = u;
                stack.emplace_back(ch, dist + 1);
            }
        }
        return ret;
    }
    public:
    vector<int> path;
    int distance;
    TreeDiameter_Unweighted(G &_g):g(_g){build();}
    void build() {
        path.clear();
        to.assign(g.size(),-1);
        distance = 0;
        if (g.empty()) return;
        pair<int,int> p = dfs(0,-1);
        pair<int,int> q = dfs(p.second,-1);
        distance = q.first;
        for (int cur = q.second; cur != -1; cur = to[cur]) path.push_back(cur);
        reverse(path.begin(), path.end());
    }
};
