#pragma once
struct Fast01 {
    using u64 = unsigned long long;

    int N, B, SB;
    vector<u64> dat; // 実データ
    vector<u64> sup; // dat の各 block が非空か
    vector<vector<u64>> upper; // sup より上の階層。最上位は 1 word。

private:
    vector<u64>& level(int h) {
        return h == 0 ? dat : (h == 1 ? sup : upper[h - 2]);
    }
    const vector<u64>& level(int h) const {
        return h == 0 ? dat : (h == 1 ? sup : upper[h - 2]);
    }
    int height() const { return 2 + (int)upper.size(); }

public:

    Fast01() : N(0), B(0), SB(0) {}
    Fast01(int n) { build(n); }

    void build(int n) {
        N = n;
        B = (N >> 6) + ((N & 63) != 0);
        SB = (B + 63) >> 6;
        dat.assign(B, 0ULL);
        sup.assign(SB, 0ULL);
        upper.clear();
        for (int sz = SB; sz > 1;) {
            sz = (sz + 63) >> 6;
            upper.emplace_back(sz, 0ULL);
        }
    }

    void rebuild_sup() {
        for (int h = 1; h < height(); h++) {
            const auto& lower = level(h - 1);
            auto& higher = level(h);
            fill(higher.begin(), higher.end(), 0ULL);
            for (int i = 0; i < (int)lower.size(); i++) {
                if (lower[i]) higher[i >> 6] |= 1ULL << (i & 63);
            }
        }
    }

    void set_all_zero() {
        for (int h = 0; h < height(); h++) {
            auto& v = level(h);
            fill(v.begin(), v.end(), 0ULL);
        }
    }

    void set_all_one() {
        fill(dat.begin(), dat.end(), ~0ULL);

        // 末尾の余計な bit を消す
        if (N == 0) return;

        int rem = N & 63;
        if (rem) {
            dat[B - 1] = (1ULL << rem) - 1;
        }

        rebuild_sup();
    }

    inline int get(int k) const {
        return (dat[k >> 6] >> (k & 63)) & 1;
    }

    inline void set1(int k) {
        for (int h = 0; h < height(); h++) {
            auto& v = level(h);
            int b = k >> 6;
            u64 old = v[b];
            v[b] |= 1ULL << (k & 63);
            if (old) break;
            k = b;
        }
    }

    inline void set0(int k) {
        for (int h = 0; h < height(); h++) {
            auto& v = level(h);
            int b = k >> 6;
            u64 old = v[b];
            v[b] &= ~(1ULL << (k & 63));
            if (!old || v[b]) break;
            k = b;
        }
    }

    inline void add(int k, int x) {
        if (x == 1) set1(k);
        else if (x == -1) set0(k);
    }

    // k 以上で最小の 1 の位置
    // 存在しないなら N
    // 更新・next・prev は O(log_64(max(N, 2)))。
    int next(int k) const {
        if (k < 0) k = 0;
        if (k >= N) return N;

        for (int h = 0; h < height(); h++) {
            const auto& v = level(h);
            int b = k >> 6;
            if (b >= (int)v.size()) return N;
            u64 x = v[b] >> (k & 63);
            if (x) {
                k += __builtin_ctzll(x);
                for (int j = h - 1; j >= 0; j--) {
                    k = (k << 6) + __builtin_ctzll(level(j)[k]);
                }
                return k;
            }
            k = b + 1;
        }

        return N;
    }

    // k 以下で最大の 1 の位置
    // 存在しないなら -1
    int prev(int k) const {
        if (k < 0) return -1;
        if (N == 0) return -1;
        if (k >= N) k = N - 1;

        for (int h = 0; h < height(); h++) {
            if (k < 0) return -1;
            const auto& v = level(h);
            int b = k >> 6;
            u64 x = v[b] & (~0ULL >> (63 - (k & 63)));
            if (x) {
                k = (b << 6) + 63 - __builtin_clzll(x);
                for (int j = h - 1; j >= 0; j--) {
                    k = (k << 6) + 63 - __builtin_clzll(level(j)[k]);
                }
                return k;
            }
            k = b - 1;
        }

        return -1;
    }
};
