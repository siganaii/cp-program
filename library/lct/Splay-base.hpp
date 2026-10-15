#pragma once
template <typename Node>
struct SplayTreeBase {
    using Ptr = Node *;
    template <typename... Args>
    Ptr my_new(const Args &...args) {
        return new Node(args...);
    }
    void my_del(Ptr p) { delete p; }
    bool is_root(Ptr t) { return !(t->p) || (t->p->l != t && t->p->r != t); }
    int size(Ptr t) const { return count(t); }
    virtual void splay(Ptr t) {
        if (!t) return;
        push(t);
        while (!is_root(t)) {
            Ptr q = t->p;
            if (is_root(q)) {
                push(q), push(t);
                rot(t);
            } 
            else {
                Ptr r = q->p;
                push(r), push(q), push(t);
                if (pos(q) == pos(t))
                rot(q), rot(t);
                else
                rot(t), rot(t);
            }
        }
    }
    Ptr get_left(Ptr t) {
        while (t->l) push(t), t = t->l;
        return t;
    }
    Ptr get_right(Ptr t) {
        while (t->r) push(t), t = t->r;
        return t;
    }
    // 根tの木を先頭k要素と残りに分割する。0 <= k <= size(t) を要求する。
    pair<Ptr, Ptr> split(Ptr t, int k) {
        if (!t) return {nullptr, nullptr};
        if (k == 0) return {nullptr, t};
        if (k == count(t)) return {t, nullptr};
        // 右側の先頭を根に上げて切り離す。深い木でも再帰しない。
        while (true) {
            push(t);
            int left = count(t->l);
            if (k < left) t = t->l;
            else if (k > left) {
                k -= left + 1;
                t = t->r;
            } else break;
        }
        splay(t);
        Ptr l = t->l;
        t->l = nullptr;
        l->p = nullptr;
        return {l, update(t)};
    }
    Ptr merge(Ptr l, Ptr r) {
        if (!l && !r) return nullptr;
        if (!l) return splay(r), r;
        if (!r) return splay(l), l;
        splay(l), splay(r);
        l = get_right(l);
        splay(l);
        l->r = r;
        r->p = l;
        update(l);
        return l;
    }
    using Key = decltype(Node::key);
    Ptr build(const vector<Key> &v) { return build(0, v.size(), v); }
    Ptr build(int l, int r, const vector<Key> &v) {
        if (l == r) return nullptr;
        if (l + 1 == r) return my_new(v[l]);
        return merge(build(l, (l + r) >> 1, v), build((l + r) >> 1, r, v));
    }
    template <typename... Args>
    void insert(Ptr &t, int k, const Args &...args) {
        splay(t);
        auto x = split(t, k);
        t = merge(merge(x.first, my_new(args...)), x.second);
    }
    void erase(Ptr &t, int k) {
        splay(t);
        auto x = split(t, k);
        auto y = split(x.second, 1);
        my_del(y.first);
        t = merge(x.first, y.second);
    }
    virtual Ptr update(Ptr t) = 0;
    protected:
    inline int count(Ptr t) const { return t ? t->cnt : 0; }
    virtual void push(Ptr t) = 0;
    Ptr build(const vector<Ptr> &v) { return build(0, v.size(), v); }
    Ptr build(int l, int r, const vector<Ptr> &v) {
        if (l + 1 >= r) return v[l];
        return merge(build(l, (l + r) >> 1, v), build((l + r) >> 1, r, v));
    }
    inline int pos(Ptr t) {
        if (t->p) {
            if (t->p->l == t) return -1;
            if (t->p->r == t) return 1;
        }
        return 0;
    }
    virtual void rot(Ptr t) {
        Ptr x = t->p, y = x->p;
        if (pos(t) == -1) {
            if ((x->l = t->r)) t->r->p = x;
            t->r = x, x->p = t;
        } 
        else {
            if ((x->r = t->l)) t->l->p = x;
            t->l = x, x->p = t;
        }
        update(x), update(t);
        if ((t->p = y)) {
            if (y->l == x) y->l = t;
            if (y->r == x) y->r = t;
        }
    }
};
