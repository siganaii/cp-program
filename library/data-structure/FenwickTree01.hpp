#pragma once
#include "library/data-structure/FenwickTree.hpp"
struct FenwickTree_01 {
    int N, n;
    using u64 = unsigned long long;
    vector<u64> dat;
    FenwickTree<int> bit;
    FenwickTree_01() {}
    FenwickTree_01(int n) { build(n); }
    void build(int m) {
        N = m;
        n = N / 64 + 1;
        dat.assign(n, u64(0));
        bit = FenwickTree<int>(n);
    }
    void set_all_one() {
        dat.assign(n, ~u64(0));

        vector<int> cnt(n, 64);

        int rem = N % 64;
        if (rem == 0) {
            dat[n - 1] = 0;
            cnt[n - 1] = 0;
        } else {
            dat[n - 1] = (u64(1) << rem) - 1;
            cnt[n - 1] = rem;
        }

        bit = FenwickTree<int>(cnt);
    }

    void set_all_zero() {
        dat.assign(n, u64(0));
        bit = FenwickTree<int>(n);
    }
    int all_prod() const { return bit.all_prod();}
    // get sum of [0,k]
    int prod(int k) const {
        if (k < 0) return 0;
        if (k >= N) k = N - 1;
        if (N == 0) return 0;
        k++;
        int ans = bit.prod(k / 64 - 1);
        ans += __builtin_popcountll(dat[k / 64] & ((u64(1) << (k % 64)) - 1));
        return ans;
    }
    // get sum of [l,r)
    int prod(int l,int r) const {
    
        if(l == 0) return prod(r-1);
        int ans = 0;
        ans -= __builtin_popcountll(dat[l / 64] & ((u64(1) << (l % 64)) - 1));
        ans += __builtin_popcountll(dat[r / 64] & ((u64(1) << (r % 64)) - 1));
        ans += bit.prod(l / 64, r / 64);
        return ans;
    }
    // get value of k
    int get(int k) const {
        return dat[k / 64] >> (k % 64) & 1;
    }
    // x == 1 or -1
    void add(int k,int x) {
        if(x == 1) add(k);
        if(x == -1) remove(k);
    }
    void add(int k) {
        if(!get(k)) {
            dat[k / 64] |= u64(1) << (k % 64);
            bit.add(k / 64, 1);
        }
    }
    void remove(int k) {
        if(get(k)) {
            dat[k / 64] &= ~(u64(1) << (k % 64));
            bit.add(k / 64, -1);
        }
    }
    // L以上でk(0-indexed)番目に小さい数
    // 存在しない時はN
    int kth(int k, int L = 0) {
        k += prod(L-1);
        if (k >= all_prod()) return N;
        int idx = bit.lower_bound(k+1);
        if (idx >= n) return N;
        k -= bit.prod(idx-1);
        u64 x = dat[idx];
        int p = __builtin_popcountll(x);
        if (p <= k) return N;
        int ok = 0;
        int ng = 64;
        while(ng - ok > 1) {
            int mid = (ng + ok) / 2;
            if(p - __builtin_popcountll(x >> mid) <= k) ok = mid;
            else ng = mid; 
        }
        return 64 * idx + ok;
    }
    // k以上で最小の数
    // 存在しない時はN
    int next(int k) {
        if (k < 0) k = 0;
        if (k >= N) return N;
        int idx = k / 64;
        k %= 64;
        u64 x = dat[idx] & ~((u64(1) << k) - 1);
        if (x) return 64 * idx + __builtin_ctzll(x);
        idx = bit.lower_bound(1 + bit.prod(idx));
        if (idx >= n || !dat[idx]) return N;
        return 64 * idx + __builtin_ctzll(dat[idx]);
    }
    // k以下で最大の数
    // 存在しない時は-1
    int prev(int k) {
        if (k < 0) return -1;
        if (k >= N) k = N - 1;
        if (N == 0) return -1;
        int idx = k / 64;
        k %= 64;
        u64 x = dat[idx];
        if (k < 63) x &= (u64(1) << (k + 1)) - 1;
        if (x) return 64 * idx + 63 - __builtin_clzll(x);
        int val = bit.prod(idx - 1);
        if(val == 0) return -1;
        idx = bit.lower_bound(val);
        return 64 * idx + 63 - __builtin_clzll(dat[idx]);
    }
};