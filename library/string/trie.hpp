#pragma once
template <size_t X = 26, char margin = 'a'>
struct Trie {
    struct Node {
        array<int, X> nxt;
        vector<int> idxs;
        int idx;
        int par;
        char key;
        Node(char c = '$', int p = -1) : idx(-1), par(p), key(c) {
            fill(nxt.begin(), nxt.end(), -1);
        }
    };

    vector<Node> st;
    Trie(char c = '$') { st.emplace_back(c, -1); }

    inline int &next(int i, int j) { return st[i].nxt[j]; }

    int add(const string &s, int x) {
        int pos = 0;
        for (int i = 0; i < (int)s.size(); i++) {
            int k = s[i] - margin;
            if (~next(pos, k)) {
                pos = next(pos, k);
                continue;
            }
            int npos = (int)st.size();
            next(pos, k) = npos;
            st.emplace_back(s[i], pos);
            pos = npos;
        }
        st[pos].idx = x;
        st[pos].idxs.emplace_back(x);
        return pos;
    }

    int find(const string &s) const {
        int pos = 0;
        for (int i = 0; i < (int)s.size(); i++) {
            int k = s[i] - margin;
            if (st[pos].nxt[k] < 0) return -1;
            pos = st[pos].nxt[k];
        }
        return pos;
    }

    int move(int pos, char c) const {
        if (pos < 0) return -1;
        return st[pos].nxt[c - margin];
    }

    int size() const { return (int)st.size(); }
    int idx(int pos) const { return pos < 0 ? -1 : st[pos].idx; }
    vector<int> idxs(int pos) const { return pos < 0 ? vector<int>() : st[pos].idxs; }
    int par(int pos) const { return pos <= 0 ? -1 : st[pos].par; }
};