#pragma once
// https://atcoder.jp/contests/abc448/submissions/78649355
// 実数座標 Dynamic Li Chao Tree
//
// FUNC:
// struct F {
//   using value_type = double;
//   double a, b;
//   double operator()(double x) { return a * x + b; }
// };
//
// PERSISTENT : 永続化するか
// MINIMIZE   : true -> min, false -> max
// DEPTH      : 最大深さ
// X          : x 座標の型 (double / long double)

//
// DEPTH で打ち切るため、実数区間は最終的に
// およそ (R-L)/2^DEPTH の幅で離散化される。
// 特に segment add の境界はこの精度で近似される。
template <
    typename FUNC,
    bool PERSISTENT,
    bool MINIMIZE,
    int DEPTH = 60,
    typename X = double
>
struct Dynamic_LiChao_Tree_Real {
  using T = typename FUNC::value_type;
  using np = int;
  static_assert(std::is_floating_point_v<X>);
  struct Node {
    int fid;
    int l, r;
  };
  vector<FUNC> funcs;
  vector<Node> pool;
  X L, R;
  T INF;
  static inline T evaluate(FUNC &f, X x) {
    return f(x);
  }
  static T default_infinity() {
    if constexpr (numeric_limits<T>::has_infinity) {
      return numeric_limits<T>::infinity();
    } else {
      return numeric_limits<T>::max();
    }
  }
  Dynamic_LiChao_Tree_Real(
      X L,
      X R,
      T INF = default_infinity()
  ) : L(L), R(R), INF(INF) {
    assert(L < R);
  }
  Dynamic_LiChao_Tree_Real() {}
  // 必須ではない。再確保を減らしたいときだけ使う。
  void reserve_nodes(size_t n) {
    pool.reserve(n);
  }
  void reserve_funcs(size_t n) {
    funcs.reserve(n);
  }
  int node_count() const {
    return (int)pool.size();
  }
  size_t node_capacity() const {
    return pool.capacity();
  }
  int func_count() const {
    return (int)funcs.size();
  }
  void reset() {
    funcs.clear();
    pool.clear();
  }
  np new_root() const {
    return -1;
  }
  np chmin_line(np root, FUNC f) {
    static_assert(MINIMIZE);
    int fid = (int)funcs.size();
    funcs.emplace_back(std::move(f));
    return add_line_rec(root, fid, L, R, 0);
  }
  np chmax_line(np root, FUNC f) {
    static_assert(!MINIMIZE);
    int fid = (int)funcs.size();
    funcs.emplace_back(std::move(f));
    return add_line_rec(root, fid, L, R, 0);
  }
  // [xl, xr)
  np chmin_segment(np root, X xl, X xr, FUNC f) {
    static_assert(MINIMIZE);
    int fid = (int)funcs.size();
    funcs.emplace_back(std::move(f));
    return add_segment_rec(root, xl, xr, fid, L, R, 0);
  }
  // [xl, xr)
  np chmax_segment(np root, X xl, X xr, FUNC f) {
    static_assert(!MINIMIZE);
    int fid = (int)funcs.size();
    funcs.emplace_back(std::move(f));
    return add_segment_rec(root, xl, xr, fid, L, R, 0);
  }
  // (値, 関数番号)
  pair<T, int> query(np root, X x) {
    assert(L <= x && x < R);
    if (root == -1) {
      if constexpr (MINIMIZE) {
        return {INF, -1};
      } else {
        return {-INF, -1};
      }
    }
    return query_rec(root, x, L, R, 0);
  }
private:
  static inline X midpoint(X l, X r) {
    return l + (r - l) / X(2);
  }
  static inline bool better(const T &a, const T &b) {
    if constexpr (MINIMIZE) {
      return a < b;
    } else {
      return a > b;
    }
  }
  np new_node(
      int fid = -1,
      int l = -1,
      int r = -1
  ) {
    int id = (int)pool.size();
    pool.push_back({fid, l, r});
    return id;
  }
  np copy_node(np c) {
    if constexpr (!PERSISTENT) {
      return c;
    } else {
      Node t = pool[c];
      int id = (int)pool.size();
      pool.push_back(t);
      return id;
    }
  }
  inline T evaluate_inner(int fid, X x) {
    if (fid == -1) {
      if constexpr (MINIMIZE) {
        return INF;
      } else {
        return -INF;
      }
    }
    return evaluate(funcs[fid], x);
  }
  np add_segment_rec(
      np c,
      X xl,
      X xr,
      int fid,
      X node_l,
      X node_r,
      int depth
  ) {
    // 交差なし
    if (xr <= node_l || node_r <= xl) {
      return c;
    }
    // 完全被覆
    if (xl <= node_l && node_r <= xr) {
      return add_line_rec(
          c,
          fid,
          node_l,
          node_r,
          depth
      );
    }
    X node_m = midpoint(node_l, node_r);
    // これ以上分割できない。
    //
    // leaf を midpoint で代表させ、
    // midpoint が [xl,xr) 内ならこの leaf 全体に追加する。
    if (depth == DEPTH ||
        node_m == node_l ||
        node_m == node_r) {
      if (xl <= node_m && node_m < xr) {
        return add_line_rec(
            c,
            fid,
            node_l,
            node_r,
            depth
        );
      }
      return c;
    }
    int old_l = (c == -1 ? -1 : pool[c].l);
    int old_r = (c == -1 ? -1 : pool[c].r);
    int new_l = old_l;
    int new_r = old_r;
    if (xl < node_m) {
      new_l = add_segment_rec(
          old_l,
          xl,
          xr,
          fid,
          node_l,
          node_m,
          depth + 1
      );
    }
    if (node_m < xr) {
      new_r = add_segment_rec(
          old_r,
          xl,
          xr,
          fid,
          node_m,
          node_r,
          depth + 1
      );
    }
    if (new_l == old_l && new_r == old_r) {
      return c;
    }
    if (c == -1) {
      c = new_node();
    } else {
      c = copy_node(c);
    }
    pool[c].l = new_l;
    pool[c].r = new_r;
    return c;
  }
  np add_line_rec(
      np c,
      int fid,
      X node_l,
      X node_r,
      int depth
  ) {
    if (c == -1) {
      return new_node(fid);
    }
    int gid = pool[c].fid;
    if (gid == -1) {
      c = copy_node(c);
      pool[c].fid = fid;
      return c;
    }
    T fl = evaluate_inner(fid, node_l);
    T gl = evaluate_inner(gid, node_l);
    T fr = evaluate_inner(fid, node_r);
    T gr = evaluate_inner(gid, node_r);
    bool bl = better(fl, gl);
    bool br = better(fr, gr);
    // fid が区間全体で良い
    if (bl && br) {
      c = copy_node(c);
      pool[c].fid = fid;
      return c;
    }
    // gid が区間全体で良い
    if (!bl && !br) {
      return c;
    }
    X node_m = midpoint(node_l, node_r);
    // これ以上分割できない場合は midpoint の良い方を採用
    if (depth == DEPTH ||
        node_m == node_l ||
        node_m == node_r) {
      T fm = evaluate_inner(fid, node_m);
      T gm = evaluate_inner(gid, node_m);
      if (better(fm, gm)) {
        c = copy_node(c);
        pool[c].fid = fid;
      }
      return c;
    }
    T fm = evaluate_inner(fid, node_m);
    T gm = evaluate_inner(gid, node_m);
    bool bm = better(fm, gm);
    if (bm) {
      // midpoint では fid の方が良いので
      // fid をこの node に置き、gid を片側へ送る。
      if (bl) {
        int old_r = pool[c].r;
        int new_r = add_line_rec(
            old_r,
            gid,
            node_m,
            node_r,
            depth + 1
        );
        c = copy_node(c);
        pool[c].fid = fid;
        pool[c].r = new_r;
      } else {
        int old_l = pool[c].l;
        int new_l = add_line_rec(
            old_l,
            gid,
            node_l,
            node_m,
            depth + 1
        );
        c = copy_node(c);
        pool[c].fid = fid;
        pool[c].l = new_l;
      }
    } else {
      // midpoint では gid の方が良いので
      // gid をこの node に残し、fid を片側へ送る。
      if (!bl) {
        int old_r = pool[c].r;
        int new_r = add_line_rec(
            old_r,
            fid,
            node_m,
            node_r,
            depth + 1
        );
        c = copy_node(c);
        pool[c].r = new_r;
      } else {
        int old_l = pool[c].l;
        int new_l = add_line_rec(
            old_l,
            fid,
            node_l,
            node_m,
            depth + 1
        );
        c = copy_node(c);
        pool[c].l = new_l;
      }
    }
    return c;
  }
  pair<T, int> query_rec(
      np c,
      X x,
      X node_l,
      X node_r,
      int depth
  ) {
    int fid = pool[c].fid;
    pair<T, int> res = {
        evaluate_inner(fid, x),
        fid
    };
    if (depth == DEPTH) {
      return res;
    }
    X node_m = midpoint(node_l, node_r);
    if (node_m == node_l ||
        node_m == node_r) {
      return res;
    }
    int child;
    X next_l, next_r;
    if (x < node_m) {
      child = pool[c].l;
      next_l = node_l;
      next_r = node_m;
    } else {
      child = pool[c].r;
      next_l = node_m;
      next_r = node_r;
    }
    if (child != -1) {
      auto res1 = query_rec(
          child,
          x,
          next_l,
          next_r,
          depth + 1
      );
      if (better(res1.first, res.first)) {
        res = res1;
      }
    }
    return res;
  }
};