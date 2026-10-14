#pragma once
template <typename T>
struct SequenceRBST {
    struct Node {
        int l, r;
        T key;
        int cnt;

        Node() : l(0), r(0), key(), cnt(0) {}
        Node(const T& key) : l(0), r(0), key(key), cnt(1) {}
    };
    vector<Node> tr;
    vector<int> stock;
    int root = 0;
    SequenceRBST() {
        tr.emplace_back();
    }
    SequenceRBST(const vector<T>& v) {
        tr.reserve(v.size() + 1);
        tr.emplace_back();
        root = build(0, (int)v.size(), v);
    }
    int size() const {
        return count(root);
    }
    bool empty() const {
        return root == 0;
    }
    // k 番目に x を挿入
    // 0 <= k <= size()
    void insert(int k, const T& x) {
        assert(0 <= k && k <= size());
        root = insert(root, k, new_node(x));
    }
    // k 番目を削除
    // 0 <= k < size()
    void erase(int k) {
        assert(0 <= k && k < size());
        root = erase(root, k);
    }
    // k 番目を取得
    // 0 <= k < size()
    T& get(int k) {
        assert(0 <= k && k < size());
        return tr[kth_element(k)].key;
    }
    const T& get(int k) const {
        assert(0 <= k && k < size());
        return tr[kth_element(k)].key;
    }
    T& operator[](int k) {
        return get(k);
    }
    const T& operator[](int k) const {
        return get(k);
    }
private:
    static uint32_t rng() {
        static uint64_t x = 88172645463325252ULL;
        x ^= x << 7;
        x ^= x >> 9;
        return x & 0xFFFFFFFFu;
    }
    int count(int t) const {
        return t ? tr[t].cnt : 0;
    }
    void update(int t) {
        tr[t].cnt = 1 + count(tr[t].l) + count(tr[t].r);
    }
    int new_node(const T& x) {
        if (!stock.empty()) {
            int t = stock.back();
            stock.pop_back();
            tr[t] = Node(x);
            return t;
        }
        tr.emplace_back(x);
        return (int)tr.size() - 1;
    }
    void del_node(int t) {
        stock.push_back(t);
    }
    pair<int, int> split(int t, int k) {
        if (!t) return {0, 0};
        if (k <= count(tr[t].l)) {
            auto [a, b] = split(tr[t].l, k);
            tr[t].l = b;
            update(t);
            return {a, t};
        } else {
            auto [a, b] =
                split(tr[t].r, k - count(tr[t].l) - 1);
            tr[t].r = a;
            update(t);
            return {t, b};
        }
    }
    int merge(int l, int r) {
        if (!l || !r) return l ? l : r;
        uint64_t x = (uint64_t)rng() * (count(l) + count(r));
        if ((int)(x >> 32) < count(l)) {
            tr[l].r = merge(tr[l].r, r);
            update(l);
            return l;
        } else {
            tr[r].l = merge(l, tr[r].l);
            update(r);
            return r;
        }
    }
    int insert(int t, int k, int x) {
        if (!t) return x;
        int n = count(t);
        // 確率 1/(n+1) で新ノードをこの部分木の根にする
        if ((uint64_t)rng() * (n + 1) >> 32 == 0) {
            auto [l, r] = split(t, k);
            tr[x].l = l;
            tr[x].r = r;
            update(x);
            return x;
        }
        int ls = count(tr[t].l);
        if (k <= ls) {
            tr[t].l = insert(tr[t].l, k, x);
        } else {
            tr[t].r = insert(tr[t].r, k - ls - 1, x);
        }
        update(t);
        return t;
    }

    int erase(int t, int k) {
        int ls = count(tr[t].l);
        if (k < ls) {
            tr[t].l = erase(tr[t].l, k);
            update(t);
            return t;
        }
        if (k > ls) {
            tr[t].r = erase(tr[t].r, k - ls - 1);
            update(t);
            return t;
        }
        int res = merge(tr[t].l, tr[t].r);
        del_node(t);
        return res;
    }

    int kth_element(int k) const {
        int t = root;
        while (true) {
            int ls = count(tr[t].l);
            if (k < ls) {
                t = tr[t].l;
            } else if (k == ls) {
                return t;
            } else {
                k -= ls + 1;
                t = tr[t].r;
            }
        }
    }
    int build(int l, int r, const vector<T>& v) {
        if (l == r) return 0;
        int m = (l + r) >> 1;
        int t = new_node(v[m]);
        tr[t].l = build(l, m, v);
        tr[t].r = build(m + 1, r, v);
        update(t);
        return t;
    }
};