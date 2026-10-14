#pragma once
//https://atcoder.jp/contests/abc467/submissions/77698223
#include <immintrin.h>
struct wm_bit_vector {
    using u32 = uint32_t;
    using u64 = uint64_t;
    static constexpr u32 w = 64;
    vector<u64> block;
    vector<u32> count;
    u32 n = 0, zeros = 0;
    wm_bit_vector() = default;
    wm_bit_vector(int _n) { init(_n); }
    __attribute__((optimize("O3,unroll-loops")))
    void init(int _n) {
        n = zeros = _n;
        block.assign(n / w + 1, 0);
        count.assign(block.size(), 0);
    }
    inline u32 get(u32 i) const {
        return u32(block[i / w] >> (i % w)) & 1u;
    }
    inline void set(u32 i) {
        block[i / w] |= 1ULL << (i % w);
    }
    __attribute__((target("popcnt")))
    void build() {
        for (u32 i = 1; i < block.size(); ++i) {
            count[i] = count[i - 1] + _mm_popcnt_u64(block[i - 1]);
        }
        zeros = rank0(n);
    }
    inline u32 rank0(u32 i) const {
        return i - rank1(i);
    }
    __attribute__((target("bmi2,popcnt")))
    inline u32 rank1(u32 i) const {
        return count[i / w] + _mm_popcnt_u64(_bzhi_u64(block[i / w], i % w));
    }
};

// X: x座標型
// Y: y座標・値型
// S: Fenwick に載せる値型
//
// op, e, inv は ACL 風。
// ただし Fenwick なので inv が必要。
// op は可換であることを想定。
template <class X, class Y, class S, S (*op)(S, S), S (*e)(), S (*inv)(S)>
struct WaveletMatrixFenwick {
private:
    using u32 = uint32_t;
    using P = pair<X, Y>;
    struct BIT {
        int N = 0;
        vector<S> data;
        BIT() = default;
        BIT(int n) { init(n); }
        void init(int n) {
            N = n;
            data.assign(N + 1, e());
        }
        void add(int k, S x) {
            for (++k; k <= N; k += k & -k) {
                data[k] = op(data[k], x);
            }
        }
        S prod(int r) const {
            S ret = e();
            for (; r > 0; r -= r & -r) {
                ret = op(ret, data[r]);
            }
            return ret;
        }
        // [l, r)
        S prod(int l, int r) const {
            S ret = e();
            while (l != r) {
                if (l < r) {
                    ret = op(ret, data[r]);
                    r -= r & -r;
                } else {
                    ret = op(ret, inv(data[l]));
                    l -= l & -l;
                }
            }
            return ret;
        }
    };
    int n = 0;
    int lg = 0;
    vector<P> ps;
    vector<Y> ys;
    vector<wm_bit_vector> bv;
    vector<BIT> bit;
    BIT all_bit;
    // 各 compressed point の現在値
    vector<S> cur;
    int point_id(X x, Y y) const {
        P p{x, y};
        auto it = lower_bound(ps.begin(), ps.end(), p);
        assert(it != ps.end() && *it == p);
        return int(it - ps.begin());
    }
    void apply_by_point_id(int i, S delta) {
        cur[i] = op(cur[i], delta);
        all_bit.add(i, delta);

        for (int h = lg - 1; h >= 0; h--) {
            int i0 = bv[h].rank0(i);

            if (bv[h].get(i)) {
                i += bv[h].zeros - i0;
            } else {
                i = i0;
            }

            bit[h].add(i, delta);
        }
    }

    // 高速版: y-id [0, ?)
    template <class F>
    pair<int, S> max_right_prefix_by_xid(int xl, int xr, F f) const {
        S acc = e();
        assert(f(acc));

        S total = all_bit.prod(xl, xr);
        if (f(total)) {
            return {(int)ys.size(), total};
        }

        int l = xl;
        int r = xr;
        int y = 0;

        for (int h = lg - 1; h >= 0; h--) {
            int l0 = bv[h].rank0(l);
            int r0 = bv[h].rank0(r);

            int l1 = l + bv[h].zeros - l0;
            int r1 = r + bv[h].zeros - r0;

            S val0 = bit[h].prod(l0, r0);
            S nxt = op(acc, val0);

            if (f(nxt)) {
                acc = nxt;
                y |= 1 << h;
                l = l1;
                r = r1;
            } else {
                l = l0;
                r = r0;
            }
        }

        return {min<int>(y, ys.size()), acc};
    }

    // 高速版: y-id [?, size_y)
    template <class F>
    pair<int, S> min_left_suffix_by_xid(int xl, int xr, F f) const {
        S acc = e();
        assert(f(acc));

        S total = all_bit.prod(xl, xr);
        if (f(total)) {
            return {0, total};
        }

        int l = xl;
        int r = xr;
        int y = 1 << lg;

        for (int h = lg - 1; h >= 0; h--) {
            int l0 = bv[h].rank0(l);
            int r0 = bv[h].rank0(r);

            int l1 = l + bv[h].zeros - l0;
            int r1 = r + bv[h].zeros - r0;

            S val1 = bit[h].prod(l1, r1);
            S nxt = op(acc, val1);

            if (f(nxt)) {
                acc = nxt;
                y -= 1 << h;
                l = l0;
                r = r0;
            } else {
                l = l1;
                r = r1;
            }
        }

        return {min<int>(y, ys.size()), acc};
    }
public:
    WaveletMatrixFenwick() = default;
    void reserve(int m) {
        ps.reserve(m);
        ys.reserve(m);
    }
    // build 前に、将来 add / set する可能性のある点を全て登録する
    void add_point(X x, Y y) {
        ps.emplace_back(x, y);
        ys.emplace_back(y);
    }
    int size_y() const {
        return (int)ys.size();
    }

    Y get_y(int i) const {
        assert(0 <= i && i < (int)ys.size());
        return ys[i];
    }

    int xid(X x) const {
        return lower_bound(
            ps.begin(), ps.end(), x,
            [](const P& a, const X& b) {
                return a.first < b;
            }
        ) - ps.begin();
    }

    int yid(Y y) const {
        return lower_bound(ys.begin(), ys.end(), y) - ys.begin();
    }

    void build() {
        sort(ps.begin(), ps.end());
        ps.erase(unique(ps.begin(), ps.end()), ps.end());

        sort(ys.begin(), ys.end());
        ys.erase(unique(ys.begin(), ys.end()), ys.end());

        n = ps.size();

        int m = max<int>(1, ys.size());
        lg = 0;
        while ((1u << lg) < (u32)m) lg++;
        if (lg == 0) lg = 1;

        bv.assign(lg, wm_bit_vector());
        bit.assign(lg, BIT());

        for (int h = 0; h < lg; h++) {
            bv[h].init(n);
            bit[h].init(n);
        }

        all_bit.init(n);
        cur.assign(n, e());

        vector<u32> a(n), b(n);

        for (int i = 0; i < n; i++) {
            a[i] = yid(ps[i].second);
        }

        for (int h = lg - 1; h >= 0; h--) {
            for (int i = 0; i < n; i++) {
                if ((a[i] >> h) & 1u) {
                    bv[h].set(i);
                }
            }

            bv[h].build();

            auto it0 = b.begin();
            auto it1 = b.begin() + bv[h].zeros;

            for (int i = 0; i < n; i++) {
                if (bv[h].get(i)) {
                    *it1++ = a[i];
                } else {
                    *it0++ = a[i];
                }
            }

            swap(a, b);
        }
    }
    // 点 (x, y) に delta を加える
    void add(X x, Y y, S delta) {
        apply_by_point_id(point_id(x, y), delta);
    }

    // 点 (x, y) の値を v にする
    void set(X x, Y y, S v) {
        int i = point_id(x, y);
        S delta = op(inv(cur[i]), v);
        apply_by_point_id(i, delta);
    }

    S get(X x, Y y) const {
        return cur[point_id(x, y)];
    }

    // x-id in [l, r)
    S prod_xid(int l, int r) const {
        return all_bit.prod(l, r);
    }

    // x in [L, R)
    S prod_x(X L, X R) const {
        return prod_xid(xid(L), xid(R));
    }

    S all_prod() const {
        return all_bit.prod(0, n);
    }

    // x-id in [l, r), y-id in [0, upper)
    S prod_yid_lt_by_xid(int l, int r, int upper) const {
        if (upper <= 0) return e();
        if (upper >= (int)ys.size()) return all_bit.prod(l, r);

        S res = e();

        for (int h = lg - 1; h >= 0; h--) {
            int l0 = bv[h].rank0(l);
            int r0 = bv[h].rank0(r);

            if ((upper >> h) & 1u) {
                res = op(res, bit[h].prod(l0, r0));

                l += bv[h].zeros - l0;
                r += bv[h].zeros - r0;
            } else {
                l = l0;
                r = r0;
            }
        }

        return res;
    }

    // x-id in [l, r), y-id in [lower, upper)
    S prod_yid_by_xid(int l, int r, int lower, int upper) const {
        lower = max(lower, 0);
        upper = min<int>(upper, ys.size());
        if (lower >= upper) return e();

        S res = prod_yid_lt_by_xid(l, r, upper);
        res = op(res, inv(prod_yid_lt_by_xid(l, r, lower)));
        return res;
    }

    // x in [XL, XR), y in [YL, YR)
    S prod(X XL, X XR, Y YL, Y YR) const {
        int l = xid(XL);
        int r = xid(XR);
        return prod_yid_by_xid(l, r, yid(YL), yid(YR));
    }

    // ACL segtree::max_right 風。
    //
    // 仮想配列 B:
    //   B[i] = prod(x-id in [xl,xr), y-id == i)
    //
    // 戻り値 ret:
    //   f(B[yl] + ... + B[ret-1]) == true
    // を満たす最大の ret。
    //
    // pair の second は B[yl] + ... + B[ret-1]。
    template <class F>
    pair<int, S> max_right_info_yid_by_xid(
        int xl, int xr, int yl, F f, int yr = -1
    ) const {
        if (yr == -1) yr = ys.size();

        assert(0 <= xl && xl <= xr && xr <= n);
        assert(0 <= yl && yl <= yr && yr <= (int)ys.size());

        if (yl == 0 && yr == (int)ys.size()) {
            return max_right_prefix_by_xid(xl, xr, f);
        }

        S acc = e();
        assert(f(acc));

        int dom = 1 << lg;

        struct Node {int h;int l, r;int nl, nr;S val;};

        vector<Node> st;
        st.reserve(4 * lg + 10);
        st.push_back(Node{lg - 1, xl, xr, 0, dom, all_bit.prod(xl, xr)});

        while (!st.empty()) {
            Node cur_node = st.back();
            st.pop_back();

            int h = cur_node.h;
            int l = cur_node.l;
            int r = cur_node.r;
            int nl = cur_node.nl;
            int nr = cur_node.nr;
            S val = cur_node.val;

            if (r <= l) continue;
            if (nr <= yl || yr <= nl) continue;

            if (yl <= nl && nr <= yr) {
                S nxt = op(acc, val);

                if (f(nxt)) {
                    acc = nxt;
                    continue;
                }

                if (h < 0) {
                    return {nl, acc};
                }
            }

            if (h < 0) continue;

            int mid = (nl + nr) >> 1;

            int l0 = bv[h].rank0(l);
            int r0 = bv[h].rank0(r);

            int l1 = l + bv[h].zeros - l0;
            int r1 = r + bv[h].zeros - r0;

            S val0 = bit[h].prod(l0, r0);
            S val1 = bit[h].prod(l1, r1);

            // 左から見たいので、右を先に積む
            st.push_back(Node{h - 1, l1, r1, mid, nr, val1});
            st.push_back(Node{h - 1, l0, r0, nl, mid, val0});
        }

        return {yr, acc};
    }

    template <class F>
    int max_right_yid_by_xid(int xl, int xr, int yl, F f, int yr = -1) const {
        return max_right_info_yid_by_xid(xl, xr, yl, f, yr).first;
    }

    template <class F>
    pair<int, S> max_right_info_yid(X L, X R, int yl, F f, int yr = -1) const {
        return max_right_info_yid_by_xid(xid(L), xid(R), yl, f, yr);
    }

    template <class F>
    int max_right_yid(X L, X R, int yl, F f, int yr = -1) const {
        return max_right_yid_by_xid(xid(L), xid(R), yl, f, yr);
    }

    // ACL segtree::min_left 風。
    //
    // 戻り値 ret:
    //   f(B[ret] + ... + B[yr-1]) == true
    // を満たす最小の ret。
    //
    // pair の second は B[ret] + ... + B[yr-1]。
    template <class F>
    pair<int, S> min_left_info_yid_by_xid(int xl, int xr, int yr, F f, int yl = 0) const {
        assert(0 <= xl && xl <= xr && xr <= n);
        assert(0 <= yl && yl <= yr && yr <= (int)ys.size());

        if (yl == 0 && yr == (int)ys.size()) {
            return min_left_suffix_by_xid(xl, xr, f);
        }

        S acc = e();
        assert(f(acc));

        int dom = 1 << lg;

        struct Node {
            int h;
            int l, r;
            int nl, nr;
            S val;
        };

        vector<Node> st;
        st.reserve(4 * lg + 10);
        st.push_back(Node{lg - 1, xl, xr, 0, dom, all_bit.prod(xl, xr)});

        while (!st.empty()) {
            Node cur_node = st.back();
            st.pop_back();

            int h = cur_node.h;
            int l = cur_node.l;
            int r = cur_node.r;
            int nl = cur_node.nl;
            int nr = cur_node.nr;
            S val = cur_node.val;

            if (r <= l) continue;
            if (nr <= yl || yr <= nl) continue;

            if (yl <= nl && nr <= yr) {
                S nxt = op(acc, val);

                if (f(nxt)) {
                    acc = nxt;
                    continue;
                }

                if (h < 0) {
                    return {nr, acc};
                }
            }

            if (h < 0) continue;

            int mid = (nl + nr) >> 1;

            int l0 = bv[h].rank0(l);
            int r0 = bv[h].rank0(r);

            int l1 = l + bv[h].zeros - l0;
            int r1 = r + bv[h].zeros - r0;

            S val0 = bit[h].prod(l0, r0);
            S val1 = bit[h].prod(l1, r1);

            // 右から見たいので、左を先に積む
            st.push_back(Node{h - 1, l0, r0, nl, mid, val0});
            st.push_back(Node{h - 1, l1, r1, mid, nr, val1});
        }

        return {yl, acc};
    }

    template <class F>
    int min_left_yid_by_xid(int xl, int xr, int yr, F f, int yl = 0) const {
        return min_left_info_yid_by_xid(xl, xr, yr, f, yl).first;
    }

    template <class F>
    pair<int, S> min_left_info_yid(X L, X R, int yr, F f, int yl = 0) const {
        return min_left_info_yid_by_xid(xid(L), xid(R), yr, f, yl);
    }

    template <class F>
    int min_left_yid(X L, X R, int yr, F f, int yl = 0) const {
        return min_left_yid_by_xid(xid(L), xid(R), yr, f, yl);
    }
};

/*
| 名前の一部       | 意味                      |
| ----------- | ----------------------- |
| `prod`      | 指定範囲の積を返す               |
| `max_right` | 左から条件を満たす最大境界を探す        |
| `min_left`  | 右から条件を満たす最小境界を探す        |
| `_info`     | 境界だけでなく、採用部分の積 `S` も返す  |
| `_yid`      | `y` 範囲を圧縮 ID で指定する      |
| `_by_xid`   | `x` 範囲も点列上の ID で指定する    |
| `_lt`       | 指定境界未満 `[0, upper)` を表す |
*/
