#pragma once
template <typename T> 
struct Edge {
	int from, to;
	T cost;
	Edge() = default;
	Edge(int _to, T _cost) : from(-1), to(_to), cost(_cost) {}
	Edge(int _from, int _to, T _cost) : from(_from), to(_to), cost(_cost) {}
	bool operator < (const Edge &a) const { return cost < a.cost; }
	bool operator > (const Edge &a) const { return cost > a.cost; }
    Edge &operator = (const int &x) {
        to = x;
        return *this;
    }
    operator int() const { return to; }
    friend ostream& operator<<(ostream &os, const Edge &edge) {
        return os << edge.to;
    }
};
 
template <typename T>
using Edges = vector<Edge<T>>;
template <typename T>
using Wgraph = vector<Edges<T>>;
using Ugraph = vector<vector<int>>;
// M を省略した場合は木の辺数。空グラフでは入力を消費しない。
Ugraph uinput(int N, int M = -1, bool is_directed = false, int origin = 1) {
    assert(N >= 0 && M >= -1 && (N > 0 || M <= 0));
    Ugraph g(N);
    if (M == -1) M = max(0, N - 1);
    while(M--) {
        int a,b;
        cin >> a >> b;
        a -= origin, b -= origin;
        g[a].push_back(b);
        if(!is_directed) g[b].push_back(a);
    }
    return g;
}
template <typename T>
Wgraph<T> winput(int N, int M = -1, bool is_directed = false,int origin = 1) {
    assert(N >= 0 && M >= -1 && (N > 0 || M <= 0));
    Wgraph<T> g(N);
    if (M == -1) M = max(0, N - 1);
    while(M--) {
        int a,b;
        T c;
        cin >> a >> b >> c;
        a -= origin, b -= origin;
        g[a].emplace_back(b,c);
        if(!is_directed) g[b].emplace_back(a,c);
    }
    return g;
}
