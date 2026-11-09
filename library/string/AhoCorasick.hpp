#pragma once
#include "library/string/trie.hpp"
template <size_t X = 26, char margin = 'a'>
struct AhoCorasick : Trie<X + 1, margin> {
  using TRIE = Trie<X + 1, margin>;
  using TRIE::next;
  using TRIE::st;
  using TRIE::TRIE;
  vector<int> cnt;
  vector<int> last;
  vector<int> BFS;
  vector<int> suffix_link;
  bool built = false;
  // 構築時に補完した遷移を取り除き、元の Trie の辺に戻す。O(頂点数 * X)。
  void restore_trie() {
    if (!built) return;
    for (auto &node : st) node.nxt.fill(-1);
    for (int i = 1; i < (int)st.size(); ++i) {
      next(st[i].par, st[i].key - margin) = i;
    }
    built = false;
  }
  // 構築後に追加した場合は、検索前にもう一度 build() を呼ぶ。
  int add(const string &s, int id) {
    restore_trie();
    return TRIE::add(s, id);
  }
  void add(const string &s) {
    int pos = add(s,(int)last.size());
    last.emplace_back(pos);
  }
  void build() {
    restore_trie();
    int n = st.size();
    cnt.resize(n);
    suffix_link.assign(n, 0);
    for (int i = 0; i < n; i++) {
      cnt[i] = (int)st[i].idxs.size();
    }
    queue<int> que;
    BFS.clear();
    next(0, X) = 0;
    for (int i = 0; i < (int)X; i++) {
      if (~next(0, i)) {
        next(next(0, i), X) = 0;
        que.emplace(next(0, i));
        BFS.emplace_back(next(0, i));
      } else {
        next(0, i) = 0;
      }
    }
    while (!que.empty()) {
      auto &x = st[que.front()];
      int fail = x.nxt[X];
      cnt[que.front()] += cnt[fail];
      que.pop();
      for (int i = 0; i < (int)X; i++) {
        int &nx = x.nxt[i];
        if (nx < 0) {
          nx = next(fail, i);
          continue;
        }
        BFS.emplace_back(nx);
        suffix_link[nx] = next(fail,i);
        que.emplace(nx);
        next(nx, X) = next(fail, i);
      }
    }
    reverse(BFS.begin(),BFS.end());
    built = true;
  }

  // build() 後に呼ぶ。空パターンは文字列の各境界に一致し、出現数は長さ + 1。
  vector<int> match(string &s,int heavy = true) {
    vector<int> res;
    if(heavy) {
        vector<int> visit_cnt(TRIE::size());
        visit_cnt[0] = 1;
        int pos = 0;
        for(auto &c:s) {
            pos = next(pos, c - margin);
            visit_cnt[pos]++;
        }
        res.assign(last.size(),0);
        for(auto &x:BFS) {
            visit_cnt[suffix_link[x]] += visit_cnt[x];
        }
        for(int i = 0,size = last.size();i < size;++i) {
            res[i] = visit_cnt[last[i]];
        }
    }
    else {
        res.assign(1,cnt[0]);
        int pos = 0;
        for(auto &c:s) {
            pos = next(pos, c - margin);
            res[0] += cnt[pos];
        }
    }
    return res;
  }

  int count(int pos) { return cnt[pos]; }
};
