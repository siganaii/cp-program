#pragma once
#include "library/inner/inner-hash.hpp"
template <typename G>
struct RootedTreeHash {
	using Hash = inner::Hash<1>;
	const G& g;
	int n;
	vector<Hash> hash;
	vector<int> depth;
	static vector<Hash>& xs() {
		static vector<Hash> _xs;
		return _xs;
	}
	RootedTreeHash(const G& _g, int root = 0) : g(_g), n(g.size()) {
		hash.resize(n);
		depth.resize(n, 0);
		if (n != 0) dfs(root, -1);
	}
	private:
	int dfs(int c, int p) {
		vector<pair<int, int>> order{{c, p}};
		for (int i = 0; i < (int)order.size(); ++i) {
			auto [u, par] = order[i];
			for (int d : g[u]) if (d != par) order.emplace_back(d, u);
		}
		// 子から親の順に高さとハッシュを計算する。
		for (int i = (int)order.size() - 1; i >= 0; --i) {
			auto [u, par] = order[i];
			int dep = 0;
			for (int d : g[u]) if (d != par) dep = max(dep, depth[d] + 1);
			// 必要になった高さの乱数だけを生成し、別の木でも共有する。
			while ((int)xs().size() <= dep) xs().push_back(Hash::get_basis());
			Hash x = xs()[dep], h = Hash::set(1);
			for (int d : g[u]) if (d != par) h = h * (x + hash[d]);
			hash[u] = h;
			depth[u] = dep;
		}
		return depth[c];
	}
};

/**
 * @brief 根付き木のハッシュ
 */
