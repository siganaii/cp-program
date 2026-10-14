#pragma once
#include "library/tree/RootedTreeHash.hpp"
#include "library/tree/centroid.hpp"
template <typename G>
vector<typename RootedTreeHash<G>::Hash> tree_hash(const G& g) {
    using Hash = typename RootedTreeHash<G>::Hash;
    auto cs = centroid(g);
    if ((int)cs.size() == 1) cs.push_back(cs[0]);
    vector<Hash> hs;
    for (auto& c : cs) hs.push_back(RootedTreeHash{g, c}.hash[c]);
    sort(begin(hs), end(hs));
    return hs;
}