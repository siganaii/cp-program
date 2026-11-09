#pragma once
#include <limits>
#include "library/graph/graph-template.hpp"
// 辺の重みは非負。整数型では最短距離が T に収まることを前提とする。
// 未到達は T(-1) で表すため、符号なし整数型では最大値を距離に使えない。
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
            // 型の範囲を超える遠回りの候補は加算しない。
            if constexpr (numeric_limits<T>::is_integer) {
                if (nex.cost > numeric_limits<T>::max() - d[now]) continue;
            }
            T nd = d[now] + nex.cost;
            if (nd == T(-1)) continue;
            if (d[nex] == T(-1) || d[nex] > nd) {
                d[nex] = nd;
                pq.emplace(d[nex],nex); 
            }
        }
    }
    return d;
}
template <typename T>
vector<int> restore_shortest_path(Wgraph<T> &g,vector<T> &dist,int start,int goal) {
    vector<int> res;
    if (dist[goal] == T(-1)) return res;
    int N = g.size();
    // 最短距離を満たす辺を幅優先探索し、重み 0 の閉路でも復元できるようにする。
    vector<int> prev(N, -1);
    queue<int> q;
    prev[start] = start;
    q.push(start);
    while (!q.empty() && prev[goal] == -1) {
        int now = q.front();
        q.pop();
        for (auto &edge : g[now]) {
            if (prev[edge.to] != -1 || dist[edge.to] == T(-1)) continue;
            if constexpr (numeric_limits<T>::is_integer) {
                if (edge.cost > numeric_limits<T>::max() - dist[now]) continue;
            }
            if (dist[now] + edge.cost != dist[edge.to]) continue;
            prev[edge.to] = now;
            q.push(edge.to);
        }
    }
    if (prev[goal] == -1) return res;
    for (int p = goal; p != start; p = prev[p]) res.push_back(p);
    res.push_back(start);
    reverse(begin(res),end(res));
    return res;
};
