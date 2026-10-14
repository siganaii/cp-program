#pragma once
struct Fast01 {
    using u64 = unsigned long long;

    int N, B, SB;
    vector<u64> dat; // 実データ
    vector<u64> sup; // dat の各 block が非空か

    Fast01() : N(0), B(0), SB(0) {}
    Fast01(int n) { build(n); }

    void build(int n) {
        N = n;
        B = (N + 63) >> 6;
        SB = (B + 63) >> 6;
        dat.assign(B, 0ULL);
        sup.assign(SB, 0ULL);
    }

    void rebuild_sup() {
        fill(sup.begin(), sup.end(), 0ULL);
        for (int i = 0; i < B; i++) {
            if (dat[i]) sup[i >> 6] |= 1ULL << (i & 63);
        }
    }

    void set_all_zero() {
        fill(dat.begin(), dat.end(), 0ULL);
        fill(sup.begin(), sup.end(), 0ULL);
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
        int b = k >> 6;
        int r = k & 63;

        u64 old = dat[b];
        dat[b] |= 1ULL << r;

        if (!old && dat[b]) {
            sup[b >> 6] |= 1ULL << (b & 63);
        }
    }

    inline void set0(int k) {
        int b = k >> 6;
        int r = k & 63;

        dat[b] &= ~(1ULL << r);

        if (!dat[b]) {
            sup[b >> 6] &= ~(1ULL << (b & 63));
        }
    }

    inline void add(int k, int x) {
        if (x == 1) set1(k);
        else if (x == -1) set0(k);
    }

    // k 以上で最小の 1 の位置
    // 存在しないなら N
    int next(int k) const {
        if (k < 0) k = 0;
        if (k >= N) return N;

        int b = k >> 6;
        int r = k & 63;

        u64 x = dat[b] & (~0ULL << r);
        if (x) {
            int res = (b << 6) + __builtin_ctzll(x);
            return res < N ? res : N;
        }

        int sb = b >> 6;
        int br = (b & 63) + 1;

        // br == 64 のときに << 64 してはいけない
        if (br < 64) {
            u64 y = sup[sb] & (~0ULL << br);
            if (y) {
                int nb = (sb << 6) + __builtin_ctzll(y);
                int res = (nb << 6) + __builtin_ctzll(dat[nb]);
                return res < N ? res : N;
            }
        }

        for (int s = sb + 1; s < SB; s++) {
            if (sup[s]) {
                int nb = (s << 6) + __builtin_ctzll(sup[s]);
                int res = (nb << 6) + __builtin_ctzll(dat[nb]);
                return res < N ? res : N;
            }
        }

        return N;
    }

    // k 以下で最大の 1 の位置
    // 存在しないなら -1
    int prev(int k) const {
        if (k < 0) return -1;
        if (N == 0) return -1;
        if (k >= N) k = N - 1;

        int b = k >> 6;
        int r = k & 63;

        u64 x = dat[b];
        if (r != 63) {
            x &= (1ULL << (r + 1)) - 1;
        }

        if (x) {
            return (b << 6) + 63 - __builtin_clzll(x);
        }

        int sb = b >> 6;
        int br = b & 63;

        u64 y;
        if (br == 0) {
            y = 0ULL;
        } else {
            y = sup[sb] & ((1ULL << br) - 1);
        }

        if (y) {
            int nb = (sb << 6) + 63 - __builtin_clzll(y);
            return (nb << 6) + 63 - __builtin_clzll(dat[nb]);
        }

        for (int s = sb - 1; s >= 0; s--) {
            if (sup[s]) {
                int nb = (s << 6) + 63 - __builtin_clzll(sup[s]);
                return (nb << 6) + 63 - __builtin_clzll(dat[nb]);
            }
        }

        return -1;
    }
};