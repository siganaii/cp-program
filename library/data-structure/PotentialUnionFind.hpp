#pragma once
template<class T>
struct PotentialUnionFind {
    vector<int> dat;
    vector<T> pot;
    int _N;
    PotentialUnionFind(int N): dat(N,-1),pot(N,T()),_N(N) {}
    int find(int x) {
        if(dat[x] < 0) return x;
        int r = find(dat[x]);
        pot[x] += pot[dat[x]];
        return dat[x] = r;
    }
    // return P(x) - P(find(x))
    T potential(int x) {
        find(x);
        return pot[x];
    }
    bool same(int x,int y) {return find(x) == find(y);}
    // return P(x) - P(y)
    T diff(int x,int y) {
        return potential(x) - potential(y);
    }
    // P(x) = P(y) + p
    // return Satisfiability
    bool merge(int x,int y,T p) {
        p += potential(y) - potential(x);
        x = find(x),y = find(y);
        if(x == y) return p == T();
        if(dat[x] < dat[y]) swap(x,y),p = -p;
        dat[y] += dat[x];
        dat[x] = y;
        pot[x] = p;
        return true;
    }
    int size(int x) {return -dat[find(x)];}
    vector<vector<int>> groups() {
        vector<int> leader_buf(_N), group_size(_N);
        for (int i = 0; i < _N; i++) {
            leader_buf[i] = find(i);
            group_size[leader_buf[i]]++;
        }
        vector<vector<int>> result(_N);
        for (int i = 0; i < _N; i++) {
            result[i].reserve(group_size[i]);
        }
        for (int i = 0; i < _N; i++) {
            result[leader_buf[i]].push_back(i);
        }
        result.erase(
            remove_if(result.begin(), result.end(),
                           [&](const vector<int>& v) { return v.empty(); }),
            result.end());
        return result;
    }
};
template<class T>
struct PotentialUnionFindXOR {
    vector<int> dat;
    vector<T> pot;
    int _N;
    PotentialUnionFindXOR(int N): dat(N,-1),pot(N,T()),_N(N) {}
    int find(int x) {
        if(dat[x] < 0) return x;
        int r = find(dat[x]);
        pot[x] ^= pot[dat[x]];
        return dat[x] = r;
    }
    // return P(x) ^ P(find(x))
    T potential(int x) {
        find(x);
        return pot[x];
    }
    bool same(int x,int y) {return find(x) == find(y);}
    // return P(x) ^ P(y)
    T diff(int x,int y) {
        return potential(x) ^ potential(y);
    }
    // P(x) = P(y) ^ p
    // return Satisfiability
    bool merge(int x,int y,T p) {
        p ^= potential(y) ^ potential(x);
        x = find(x),y = find(y);
        if(x == y) return p == T();
        if(dat[x] < dat[y]) {
			swap(x,y);
		}
        dat[y] += dat[x];
        dat[x] = y;
        pot[x] = p;
        return true;
    }
    int size(int x) {return -dat[find(x)];}
    vector<vector<int>> groups() {
        vector<int> leader_buf(_N), group_size(_N);
        for (int i = 0; i < _N; i++) {
            leader_buf[i] = find(i);
            group_size[leader_buf[i]]++;
        }
        vector<vector<int>> result(_N);
        for (int i = 0; i < _N; i++) {
            result[i].reserve(group_size[i]);
        }
        for (int i = 0; i < _N; i++) {
            result[leader_buf[i]].push_back(i);
        }
        result.erase(
            remove_if(result.begin(), result.end(),
                           [&](const vector<int>& v) { return v.empty(); }),
            result.end());
        return result;
    }
};