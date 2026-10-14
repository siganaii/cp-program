#pragma once
template <typename Tree, typename Node, typename T, T (*f)(T, T), T (*ts)(T)>
struct ReversibleBBST : Tree {
    using Tree::merge;
    using Tree::split;
    using typename Tree::Ptr;
    ReversibleBBST() = default;
    virtual void toggle(Ptr t) {
        if(!t) return;
        swap(t->l, t->r);
        t->SUM = ts(t->SUM);
        t->rev ^= true;
    }
    T fold(Ptr &t, int a, int b) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        auto ret = SUM(y.first);
        t = merge(x.first, merge(y.first, y.second));
        return ret;
    }
    void reverse(Ptr &t, int a, int b) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        toggle(y.first);
        t = merge(x.first, merge(y.first, y.second));
    }
    Ptr update(Ptr t) override {
        if (!t) return t;
        t->cnt = 1;
        t->SUM = t->key;
        if (t->l) t->cnt += t->l->cnt, t->SUM = f(t->l->SUM, t->SUM);
        if (t->r) t->cnt += t->r->cnt, t->SUM = f(t->SUM, t->r->SUM);
        return t;
    }
    protected:
    inline T SUM(const Ptr t) { return t ? t->SUM : T(); }
    void push(Ptr t) override {
        if (!t) return;
        if (t->rev) {
        if (t->l) toggle(t->l);
        if (t->r) toggle(t->r);
        t->rev = false;
        }
    }
};
