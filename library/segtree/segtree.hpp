#pragma once
template <typename T, typename OP>
struct segtree {
    int N;
    int size;
    vector<T> seg;
    const OP op;
    const T I;
    segtree(OP _op, const T I_) : segtree(0, _op, I_) {}
    explicit segtree(int _N, const OP _op, const T I_) : op(_op), I(I_) { init(_N); }
    explicit segtree(const vector<T> &v, const OP _op, const T I_) : op(_op), I(I_) {
        init(v.size());
        for (int i = 0; i < (int)v.size(); i++) {
            seg[i + size] = v[i];
        }
        build();
    }
    void init(int _N) {
        N = _N;
        size = 1;
        while (size < N) size <<= 1;
        seg.assign(2 * size, I);
    }
    void build() {
        for (int k = size - 1; k > 0; k--) {
            seg[k] = op(seg[2 * k], seg[2 * k + 1]);
        }
    }
    void set(int p, T x) {
        assert(0 <= p && p < N);
        p += size;
        seg[p] = x;
        while (p >>= 1) {
            seg[p] = op(seg[2 * p], seg[2 * p + 1]);
        }
    }
    void add(int p, T x) {
        assert(0 <= p && p < N);
        p += size;
        seg[p] += x;
        while (p >>= 1) {
            seg[p] = op(seg[2 * p], seg[2 * p + 1]);
        }
    }
    T get(int p) const {
        assert(0 <= p && p < N);
        return seg[p + size];
    }
    // query to [l, r)
    T prod(int l, int r) {
        assert(0 <= l && l <= r && r <= N);
        T L = I, R = I;
        for (l += size, r += size; l < r; l >>= 1, r >>= 1) {
            if (l & 1) L = op(L, seg[l++]);
            if (r & 1) R = op(seg[--r], R);
        }
        return op(L, R);
    }
    T all_prod() {return seg[1];}
    // check(a[l] * ...  * a[r-1]) が true となる最大の r
    // (右端まですべて true なら N を返す)
    template <class C>
    int max_right(int l, C check) {
        assert(0 <= l && l <= N);
        assert(check(I) == true);
        if (l == N) return N;
        l += size;
        T sm = I;
        do {
            while (l % 2 == 0) l >>= 1;
            if (!check(op(sm, seg[l]))) {
                while (l < size) {
                    l = (2 * l);
                    if (check(op(sm, seg[l]))) {
                        sm = op(sm, seg[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, seg[l]);
            l++;
        } while ((l & -l) != l);
        return N;
    }
    // check(a[l] * ... * a[r-1]) が true となる最小の l
    // (左端まで true なら 0 を返す)
    template <typename C>
    int min_left(int r, C check) {
        assert(0 <= r && r <= N);
        assert(check(I) == true);
        if (r == 0) return 0;
        r += size;
        T sm = I;
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!check(op(seg[r], sm))) {
                while (r < size) {
                    r = (2 * r + 1);
                    if (check(op(seg[r], sm))) {
                        sm = op(seg[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(seg[r], sm);
        } while ((r & -r) != r);
        return 0;
    }
};
