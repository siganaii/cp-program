#pragma once
#include "library/string/trie.hpp"
template <size_t X = 26, char margin = 'a'>
struct AhoCorasick : Trie<X + 1, margin> {
  using TRIE = Trie<X + 1, margin>;
  using TRIE::next;
  using TRIE::st;
  using TRIE::add;
  using TRIE::TRIE;
  vector<int> cnt;
  vector<int> last;
  vector<int> BFS;
  vector<int> suffix_link;
  void add(const string &s) {
    int pos = add(s,(int)last.size());
    last.emplace_back(pos);
  }
  void build() {
    int n = st.size();
    cnt.resize(n);
    suffix_link.resize(n);
    for (int i = 0; i < n; i++) {
      cnt[i] = (int)st[i].idxs.size();
    }
    queue<int> que;
    for (int i = 0; i < (int)X; i++) {
      if (~next(0, i)) {
        next(next(0, i), X) = 0;
        que.emplace(next(0, i));
      } else {
        next(0, i) = 0;
      }
    }
    BFS.clear();
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
  }

  vector<int> match(string &s,int heavy = true) {
    vector<int> res;
    if(heavy) {
        vector<int> visit_cnt(TRIE::size());
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
        res.assign(1,0);
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