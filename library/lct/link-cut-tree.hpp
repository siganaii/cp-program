#pragma once
#include "library/lct/splay-reversible.hpp"
#include "library/lct/link-cut-base.hpp"
template <typename T, T (*f)(T, T), T (*ts)(T)>
struct LinkCutTree : LinkCutBase<ReversibleSplayTree<T, f, ts>> {};

/**
 * @brief Link/Cut Tree
 * @docs docs/lct/link-cut-tree.md
 */