#pragma once
template<bool directed>
struct EulerianTrail {
    vector<vector<pair<int,int>>> g;
    vector<pair<int,int>> es;
    int M;
    vector<int> used_edge,used_vertex,deg;
    explicit EulerianTrail(int N) : g(N),M(0),used_vertex(N),deg(N) {}
    void add_edge(int u,int v) {
        es.emplace_back(u,v);
        g[u].emplace_back(v,M);
        if(directed) {
            deg[u]++;
            deg[v]--;
        }
        else {
            g[v].emplace_back(u,M);
            deg[u]++;
            deg[v]++;
        }
        M++;
    }
    pair<int,int> get_edge(int idx) {
        assert(0 <= idx && idx < M);
        return es[idx];
    }
    // pair<頂点の順番，辺の順番>
    pair<vector<int>,vector<int>> enumerate_eulerian_trail(int start = -1) {
        if(start == -1) {
            if(directed) {
                start = max_element(deg.begin(),deg.end()) - deg.begin();
                if(deg[start] == 0) {
                    start = (M == 0 ? 0 : es[0].first);
                }
            }
            else {
                for(int i = 0;i < g.size();i++) {
                    if(deg[i] & 1) {
                        start = i;
                        break;
                    }
                }
                if(start == -1) {
                    start = (M == 0 ? 0 : es[0].first);
                }
            }
        }
        if(M == 0) {
            return {{start},{}};
        }
        //辞書順最小
        // for(int i = 0;i < g.size();i++) {
        //     sort(g[i].begin(),g[i].end(),greater());
        // }
        used_edge.assign(M,0);
        vector<pair<int,int>> st;
        vector<int> ord_edge,ord_vertex,D(g.size());
        st.emplace_back(start,-1);
        D[start]++;
        while(!st.empty()) {
            auto [idx,edge_idx] = st.back();
            used_vertex[idx] = 1;
            if(g[idx].empty()) {
                ord_edge.emplace_back(edge_idx);
                ord_vertex.emplace_back(idx);
                st.pop_back();
            }
            else {
                auto e = g[idx].back();
                g[idx].pop_back();
                if(used_edge[e.second]) continue;
                used_edge[e.second] = 1;
                st.emplace_back(e);
                D[idx]--;
                D[e.first]++;
            }
        }
        ord_edge.pop_back();
        reverse(ord_edge.begin(),ord_edge.end());
        reverse(ord_vertex.begin(),ord_vertex.end());
        for(auto &x : D) if(x < 0) {return {{},{}};}
        if(ord_edge.size() != M) {
            return {{},{}};
        }
        return {ord_vertex,ord_edge};
    }   
};