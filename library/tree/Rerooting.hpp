#pragma once
// f1 は結合則を満たす演算。可換性は不要で、g[i] の隣接順に併合する。
// f2(memo[i], chd, par) は子から親へ渡す値を計算する。
// leaf は併合する隣接頂点がない場合の値で、f1 の単位元でなくてもよい。
// memo[i] は根を 0 とした部分木、memo2[i] は親側からの寄与、dp[i] は i を根とした結果。
template <typename T, typename G, typename F1, typename F2>
struct Rerooting {
  const G &g;
  const F1 f1;
  const F2 f2;
  vector<T> memo,memo2,dp;
  T leaf;

  Rerooting(const G &_g, const F1 _f1, const F2 _f2, const T &_leaf)
      : g(_g), f1(_f1), f2(_f2), memo(g.size()), memo2(g.size()), dp(g.size()), leaf(_leaf) {
    if (g.empty()) return;
    dfs(0, -1);
    dfs2(0, -1, T{});
  }

  const T &operator[](int i) const { return dp[i]; }

  void dfs(int cur, int par) {
    // 親より子を後に並べ、逆順に部分木の値を計算する。
    vector<pair<int, int>> order{{cur, par}};
    for (int i = 0; i < (int)order.size(); ++i) {
      auto [u, p] = order[i];
      for (int v : g[u]) if (v != p) order.emplace_back(v, u);
    }
    for (int i = (int)order.size() - 1; i >= 0; --i) {
      auto [u, p] = order[i];
      bool first = true;
      memo[u] = leaf;
      for (auto &dst : g[u]) {
        if (dst == p) continue;
        T val = f2(memo[dst], dst, u);
        memo[u] = first ? val : f1(memo[u], val);
        first = false;
      }
    }
  }

  void dfs2(int cur, int par, const T &pval) {
    memo2[cur] = pval;
    vector<pair<int, int>> stack{{cur, par}};
    vector<T> buf, head, tail;
    while (!stack.empty()) {
      auto [u, p] = stack.back();
      stack.pop_back();
      buf.clear();
      // 親側の寄与も隣接リスト上の位置に入れ、順序を変えない。
      for (auto &dst : g[u]) {
        buf.push_back(dst == p ? memo2[u] : f2(memo[dst], dst, u));
      }
      int n = buf.size();
      head.resize(n);
      tail.resize(n);
      if (n != 0) {
        head[0] = buf[0];
        for (int i = 1; i < n; ++i) head[i] = f1(head[i - 1], buf[i]);
        tail[n - 1] = buf[n - 1];
        for (int i = n - 2; i >= 0; --i) tail[i] = f1(buf[i], tail[i + 1]);
      }
      dp[u] = n == 0 ? leaf : head.back();
      for (int i = 0; i < n; ++i) {
        int dst = g[u][i];
        if (dst == p) continue;
        T val = n == 1 ? leaf
              : i == 0 ? tail[1]
              : i == n - 1 ? head[i - 1]
              : f1(head[i - 1], tail[i + 1]);
        memo2[dst] = f2(val, u, dst);
        stack.emplace_back(dst, u);
      }
    }
  }
};

/**
 * @brief Rerooting(全方位木DP)
 * @docs docs/tree/rerooting.md
 */
