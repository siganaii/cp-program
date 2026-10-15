#pragma once
#include "library/data-structure/unionfind.hpp"
#include "library/tree/HLD.hpp"
struct FunctionalGraph {
    int N,M;
    vector<int> TO;
    vector<int> root;
    vector<vector<int>> G;
    HLD<vector<vector<int>>> hld;
    FunctionalGraph() {}
    FunctionalGraph(int N) : N(N), M(0), TO(N, -1), root(N, -1) {}
    void add_edge(int from,int to) {
        assert(0 <= from && from < N);
        assert(TO[from] == -1);
        ++M;
        TO[from] = to;
    }
    void build() {
        assert(N == M);
        unionfind uf(N);
        for(int i = 0;i < N;i++) {
            if(uf.same(i,TO[i])) root[i] = i;
            else uf.merge(i,TO[i]);
        }
        for(int i = 0;i < N;i++) {
            if(root[i] == i) {
                root[uf.find(i)] = i;
            }
        }
        for(int i = 0;i < N;i++) {
            root[i] = root[uf.find(i)];
        }
        G = vector<vector<int>>(N + 1);
        for(int i = 0;i < N;i++) {
            if(root[i] == i) {
                G[N].emplace_back(i);
            }
            else {
                G[TO[i]].emplace_back(i);
            }
        }
        hld = HLD<vector<vector<int>>>(G,N);
    }
    // from->toの距離、不可能なら-1
    int dist(int from,int to) {
        if(hld.down[to] <= hld.down[from] && hld.down[from] < hld.up[to]) {
            return hld.depth[from] - hld.depth[to];
        } 
        int r = root[from];
        int btm = TO[r];
        if(hld.down[to] <= hld.down[btm] && hld.down[btm] < hld.up[to]) {
            int x = hld.depth[from] - hld.depth[r];
            x++;
            x += hld.depth[btm] - hld.depth[to];
            return x;
        }
        return -1;
    }
    // functional graph に向かって進む
    int jump(int v, long long step) {
        int d = hld.depth[v];
        if (step <= d - 1) return hld.jump(v, N, step);
        v = root[v];
        step -= d - 1;
        int bottom = TO[v];
        int c = hld.depth[bottom];
        step %= c;
        if (step == 0) return v;
        return hld.jump(bottom, N, step - 1);
    }
    vector<int> jump_all(long long step) {
        vector<int> res(N, -1);
        // v の k 個先を res[w] に入れる
        vector<vector<pair<int,int>>> query(N);
        for(int v = 0;v < N;v++) {
            int d = hld.depth[v];
            int r = root[v];
            if (d - 1 > step) { query[v].emplace_back(v, step); }
            if (d - 1 <= step) {
                long long k = step - (d - 1);
                int bottom = TO[r];
                int c = hld.depth[bottom];
                k %= c;
                if (k == 0) {
                    res[v] = r;
                    continue;
                }
                query[bottom].emplace_back(v, k - 1);
            }
        }
        vector<int> path;
        // HLD の DFS 順では、深さに合わせて縮めると祖先の列が残る。
        for(int v: hld.ord) {
            if(v == N) continue;
            path.resize(hld.depth[v] - 1);
            path.push_back(v);
            for(auto [w, k]: query[v]) res[w] = path[(int)path.size() - 1 - k];
        }
        return res;
    }
};
