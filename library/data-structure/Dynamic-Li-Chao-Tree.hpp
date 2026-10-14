#pragma once

// 直線追加かつ非永続なら空間 O(Q)。
// 関数は ll -> T。[L, R) 上 f が overflow しないように注意。
// evaluate を書き変えると、totally monotone な関数群にも使える。
//
// 例:
// struct F {
//   using value_type = ll;
//   int a;
//   ll b;
//   ll operator()(ll x) { return a * x + b; }
// };
//
// Node は int 3つなので通常 12 byte。
// pool は vector で必要に応じて自動拡張される。
template <typename FUNC, bool PERSISTENT, bool MINIMIZE, long long INFINTY = 1LL << 60>
struct Dynamic_LiChao_Tree {
  using T = typename FUNC::value_type;
  using np = int;

  vector<FUNC> funcs;
  static inline T evaluate(FUNC &f, long long x) { return f(x); }

  struct Node {
    int fid;
    int l, r;
  };

  vector<Node> pool;
  long long L, R;

  Dynamic_LiChao_Tree(long long L, long long R) : L(L), R(R) {}
  Dynamic_LiChao_Tree() {}

  void reserve_nodes(size_t n) { pool.reserve(n); }
  void reserve_funcs(size_t n) { funcs.reserve(n); }
  int node_count() const { return (int)pool.size(); }
  size_t node_capacity() const { return pool.capacity(); }
  int func_count() const { return (int)funcs.size(); }

  void reset() {
    funcs.clear();
    pool.clear();
  }

  np new_root() const { return -1; }

  np chmin_line(np root, FUNC f) {
    static_assert(MINIMIZE);
    int fid = funcs.size();
    funcs.emplace_back(std::move(f));
    if (root == -1) root = new_node();
    return add_line_rec(root, fid, L, R);
  }

  np chmax_line(np root, FUNC f) {
    static_assert(!MINIMIZE);
    int fid = funcs.size();
    funcs.emplace_back(std::move(f));
    if (root == -1) root = new_node();
    return add_line_rec(root, fid, L, R);
  }

  // [xl, xr)
  np chmin_segment(np root, long long xl, long long xr, FUNC f) {
    static_assert(MINIMIZE);
    int fid = funcs.size();
    funcs.emplace_back(std::move(f));
    if (root == -1) root = new_node();
    return add_segment_rec(root, xl, xr, fid, L, R);
  }

  // [xl, xr)
  np chmax_segment(np root, long long xl, long long xr, FUNC f) {
    static_assert(!MINIMIZE);
    int fid = funcs.size();
    funcs.emplace_back(std::move(f));
    if (root == -1) root = new_node();
    return add_segment_rec(root, xl, xr, fid, L, R);
  }

  // (値・関数番号)
  pair<T, int> query(np root, long long x) {
    assert(L <= x && x < R);
    if (root == -1) {
      if constexpr (MINIMIZE) return {INFINTY, -1};
      else return {-INFINTY, -1};
    }
    return query_rec(root, x, L, R);
  }

private:
  np new_node(int fid = -1, int l = -1, int r = -1) {
    int id = pool.size();
    pool.push_back({fid, l, r});
    return id;
  }

  np copy_node(np c) {
    if constexpr (!PERSISTENT) return c;
    else {
      Node t = pool[c];
      int id = pool.size();
      pool.push_back(t);
      return id;
    }
  }

  inline T evaluate_inner(int fid, long long x) {
    if (fid == -1) {
      if constexpr (MINIMIZE) return INFINTY;
      else return -INFINTY;
    }
    return evaluate(funcs[fid], x);
  }

  np add_segment_rec(np c, long long xl, long long xr, int fid, long long node_l, long long node_r) {
    chmax(xl, node_l), chmin(xr, node_r);
    if (xl >= xr) return c;
    if (node_l < xl || xr < node_r) {
      c = copy_node(c);
      long long node_m = node_l + (node_r - node_l) / 2;
      int l = pool[c].l, r = pool[c].r;
      if (l == -1) l = new_node();
      if (r == -1) r = new_node();
      l = add_segment_rec(l, xl, xr, fid, node_l, node_m);
      r = add_segment_rec(r, xl, xr, fid, node_m, node_r);
      pool[c].l = l, pool[c].r = r;
      return c;
    }
    return add_line_rec(c, fid, node_l, node_r);
  }

  np add_line_rec(np c, int fid, long long node_l, long long node_r) {
    int gid = pool[c].fid;
    T fl = evaluate_inner(fid, node_l), fr = evaluate_inner(fid, node_r - 1);
    T gl = evaluate_inner(gid, node_l), gr = evaluate_inner(gid, node_r - 1);
    bool bl = (MINIMIZE ? fl < gl : fl > gl);
    bool br = (MINIMIZE ? fr < gr : fr > gr);

    if (bl && br) {
      c = copy_node(c);
      pool[c].fid = fid;
      return c;
    }
    if (!bl && !br) return c;

    c = copy_node(c);
    long long node_m = node_l + (node_r - node_l) / 2;
    T fm = evaluate_inner(fid, node_m), gm = evaluate_inner(gid, node_m);
    bool bm = (MINIMIZE ? fm < gm : fm > gm);

    if (bm) {
      pool[c].fid = fid;
      if (bl) {
        int r = pool[c].r;
        if (r == -1) r = new_node();
        pool[c].r = add_line_rec(r, gid, node_m, node_r);
      } else {
        int l = pool[c].l;
        if (l == -1) l = new_node();
        pool[c].l = add_line_rec(l, gid, node_l, node_m);
      }
    } else {
      if (!bl) {
        int r = pool[c].r;
        if (r == -1) r = new_node();
        pool[c].r = add_line_rec(r, fid, node_m, node_r);
      } else {
        int l = pool[c].l;
        if (l == -1) l = new_node();
        pool[c].l = add_line_rec(l, fid, node_l, node_m);
      }
    }
    return c;
  }

  pair<T, int> query_rec(np c, long long x, long long node_l, long long node_r) {
    int fid = pool[c].fid;
    pair<T, int> res = {evaluate_inner(fid, x), fid};
    long long node_m = node_l + (node_r - node_l) / 2;

    if (x < node_m) {
      int l = pool[c].l;
      if (l != -1) {
        auto res1 = query_rec(l, x, node_l, node_m);
        if constexpr (MINIMIZE) res = min(res, res1);
        else res = max(res, res1);
      }
    } else {
      int r = pool[c].r;
      if (r != -1) {
        auto res1 = query_rec(r, x, node_m, node_r);
        if constexpr (MINIMIZE) res = min(res, res1);
        else res = max(res, res1);
      }
    }
    return res;
  }
};