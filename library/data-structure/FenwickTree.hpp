#pragma once
template <typename T>
struct FenwickTree{
    int N;
    T total = 0;
    vector<T> data;
    FenwickTree() = default;
    FenwickTree(int size) {init(size);}
    FenwickTree(vector<T> &v) {
        N = v.size() + 2;
        total = 0;
        data.reserve(N + 1);

        data.emplace_back(0);
        for (const auto &e : v) {
            total += e;
            data.emplace_back(e);
        }
        data.emplace_back(0);
        data.emplace_back(0);

        for (int i = 1; i < N; ++i) {
            int j = i + (i & -i);
            if (j < N) data[j] += data[i];
        }
    }
    void init(int size) {
        N = size + 2;
        total = T{};
        data.assign(N + 1,{});
    }
    // get sum of [0,k]
    T prod(int k) const {
        if (k < 0) return T{};
        T ret{};
        for (++k;k > 0;k -= k & -k) ret += data[k];
        return ret;
    }
    // get sum of [l,r)
    inline T prod(int l,int r) const {return prod(r - 1) - prod(l - 1);}
    // get value of k
    inline T get(int k) const {return prod(k) - prod(k - 1); }
    T all_prod() const {return total;}
    void add(int k, T x) { 
        total += x;
        for(++k;k < N;k += k & -k) data[k] += x;
    }
    // minimize i s.t. sum(i) >= w
    int lower_bound(T w) {
        if (w <= 0) return 0;
        int x = 0;
        for(int k = 1 <<__lg(N);k;k >>= 1) {
            if (x + k <= N - 1 && data[x + k] < w) {
                w -= data[x + k];
                x += k;
            }
        }
        if(x > N - 2) return N - 2;
        return x;
    }
    // minimize i s.t. sum(i) > w
    int upper_bound(T w) {
        if (w < 0) return 0;
        int x = 0;
        for(int k = 1 <<__lg(N);k;k >>= 1) {
            if (x + k <= N - 1 && data[x + k] <= w) {
                w -= data[x + k];
                x += k;
            }
        }
        if(x > N - 2) return N - 2;
        return x;
    }
};
