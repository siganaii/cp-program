#pragma once
#include "library/graph/graph-template.hpp"
template <typename T>
vector<T> dijkstra(Wgraph<T> &g,int start = 0) {
    int N = (int)g.size();
    using P = pair<T,int>;
    vector<T> d(N,T(-1));
    priority_queue<P,vector<P>,greater<P>> pq;
    d[start] = 0;
    pq.emplace(0,start);
    while(!pq.empty()) {
        P p = pq.top();
        pq.pop();
        int now = p.second;
        if (d[now] < p.first) continue;
        for (auto nex: g[now]) {
            if (d[nex] == T(-1) || d[nex] > d[now] + nex.cost) {
                d[nex] = d[now] + nex.cost;
                pq.emplace(d[nex],nex); 
            }
        }
    }
    return d;
}
template <typename T>
vector<int> restore_shortest_path(Wgraph<T> &g,vector<T> &dist,int start,int goal) {
    vector<int> res;
    if (start == goal) {
        res.push_back(start);
        return res;
    }
    int N = g.size();
    Wgraph<T> rg(N);
    for (int i = 0;i < N;i++) {
        for (auto &edge:g[i]) rg[edge.to].emplace_back(edge.to,i,edge.cost);
    }

    res.push_back(goal);
    int p = goal;
    T dd = dist[goal];
    vector<int> vis(N);
    vis[goal] = 1;
    const T INF = numeric_limits<T>::max() / 2;
    while(p != start) {
        int nex = -1;
        T nval = INF;
        for(auto &edge:rg[p]) {
            if (vis[edge] || dist[edge] + edge.cost != dd) continue;
            if (dist[edge] != -1 && dist[edge] < nval) {
                nval = dist[edge];
                nex = edge;
                break;
            }
        }
        res.push_back((vis[nex] = 1,dd = nval, p = nex));
    }
    reverse(begin(res),end(res));
    return res;
};