#pragma once
template <class S, class F, auto mapping, auto composition, auto e, auto id>
struct dual_segtree {
private:
    int n, sz, lg2;
    vector<F> lazy;
    vector<unsigned char> has_lazy;
    vector<S> val;
    void all_apply(int i, const F& f) {
        if (i < sz) {
            if (has_lazy[i]) {
                lazy[i] = composition(f, lazy[i]);
            } else {
                lazy[i] = f;
                has_lazy[i] = 1;
            }
        } else {
            val[i - sz] = mapping(f, val[i - sz]);
        }
    }
    void push(int i) {
        if (!has_lazy[i]) return;
        all_apply(2 * i, lazy[i]);
        all_apply(2 * i + 1, lazy[i]);
        has_lazy[i] = 0;
    }
public:
    dual_segtree() : dual_segtree(0) {}
    dual_segtree(int _n)
        : dual_segtree(vector<S>(_n, e())) {}
    dual_segtree(const vector<S>& a) {
        n = (int)a.size();
        sz = bit_ceil((unsigned int)max(1, n));
        lg2 = countr_zero((unsigned int)sz);

        lazy.assign(sz, id());
        has_lazy.assign(sz, 0);

        val = a;
        val.resize(sz, e());
    }
    int size() const {
        return n;
    }
    void apply(int l, int r, const F& f) {
        assert(0 <= l && l <= r && r <= n);
        if (l == r) return;
        l += sz;
        r += sz;
        for (int i = lg2; i >= 1; --i) {
            int li = -1;
            if (((l >> i) << i) != l) {
                li = l >> i;
                push(li);
            }
            if (((r >> i) << i) != r) {
                int ri = (r - 1) >> i;
                if (ri != li) push(ri);
            }
        }
        while (l < r) {
            if (l & 1) all_apply(l++, f);
            if (r & 1) all_apply(--r, f);
            l >>= 1;
            r >>= 1;
        }
    }
    void apply_all(const F& f) {
        if (n == 0) return;
        // 木全体が実要素ならO(1)、端にダミー要素がある場合はO(log n)で更新する。
        if (n == sz) all_apply(1, f);
        else apply(0, n, f);
    }
    S get(int p) {
        assert(0 <= p && p < n);
        p += sz;
        for (int i = lg2; i >= 1; --i) {
            push(p >> i);
        }
        return val[p - sz];
    }
    void set(int p, const S& x) {
        assert(0 <= p && p < n);
        p += sz;
        for (int i = lg2; i >= 1; --i) {
            push(p >> i);
        }
        val[p - sz] = x;
    }
    vector<S> materialize() {
        for (int i = 1; i < sz; ++i) push(i);
        vector<S> res(n);
        for (int i = 0; i < n; ++i) res[i] = val[i];
        return res;
    }
};
