#pragma once
#include "library/tree/RootedTreeHash.hpp"
#include "library/tree/centroid.hpp"
template <typename G>
vector<typename RootedTreeHash<G>::Hash> tree_hash(const G& g) {
    using Hash = typename RootedTreeHash<G>::Hash;
    auto cs = centroid(g);
    vector<Hash> hs;
    for (auto& c : cs) hs.push_back(RootedTreeHash{g, c}.hash[c]);
    // 重心が一つの場合も返り値は二要素とし、計算済みの値を複製する。
    if (hs.size() == 1) hs.push_back(hs[0]);
    sort(begin(hs), end(hs));
    return hs;
}
