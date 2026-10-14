#pragma once
template <typename G>
struct TreeDiameter_Unweighted {
    private:
    G &g;
    vector<int> to;
    pair<int,int> dfs(int cur,int par) {
        pair<int,int> ret(0,cur);
        for(auto &ch : g[cur]) {
            if(ch == par) continue;
            pair<int,int> cost = dfs(ch,cur);
            cost.first++;
            if(ret < cost) {
                ret = cost;
                to[cur] = ch;
            }
        }
        return ret;
    }
    public:
    vector<int> path;
    int distance;
    TreeDiameter_Unweighted(G &_g):g(_g){build();}
    void build() {
        to.assign(g.size(),-1);
        pair<int,int> p = dfs(0,-1);
        pair<int,int> q = dfs(p.second,-1);
        distance = q.first;
        int cur = p.second;
        path.emplace_back(cur);
        while(cur != q.second) {
            path.emplace_back(to[cur]);
            cur = to[cur];
        }
    }
};