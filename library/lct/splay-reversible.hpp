#pragma once
#include "library/lct/reversible-bbst-base.hpp"
#include "library/lct/Splay-base.hpp"
template <typename T>
struct ReversibleSplayTreeNode {
    using Ptr = ReversibleSplayTreeNode *;
    Ptr l, r, p;
    T key, SUM;
    int cnt;
    bool rev;
    ReversibleSplayTreeNode(const T &t = T())
        : l(), r(), p(), key(t), SUM(t), cnt(1), rev(false) {}
};

template <typename T, T (*f)(T, T), T (*ts)(T)>
struct ReversibleSplayTree
    : ReversibleBBST<SplayTreeBase<ReversibleSplayTreeNode<T>>,
                     ReversibleSplayTreeNode<T>, T, f, ts> {
  using Node = ReversibleSplayTreeNode<T>;
};

/**
 * @brief 反転可能Splay Tree
 */