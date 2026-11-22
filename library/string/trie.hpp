#pragma once

template <size_t X = 26, char margin = 'a'>
struct Trie {
    struct Node {
        array<int, X> nxt;
        int idx;
        int par;
        char key;
        vector<int> idxs;
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

// 実在する辺と登録済みの ID を格納する Trie。文字は [margin, margin + X) を使う。
template <size_t X = 26, char margin = 'a'>
class SparseTrie {
    static_assert(1 <= X && X <= 32);

    struct Node {
        uint32_t symbols = 0;
        int children = 0;  // 子が1個なら状態番号、2個以上なら edges 内の先頭位置。
        int parent = -1;
        int last_pattern = -1;

        explicit Node(int p = -1) : parent(p) {}
    };

    struct Pattern {
        int id;
        int previous;  // 同じ終端状態に登録した1つ前のパターン。
    };

    vector<Node> nodes{Node{}};
    // 各状態の子を文字順に格納する。ブロックが満杯になるたびに新しい領域へ移す。
    vector<int> edges;
    vector<Pattern> patterns;

protected:
    static uint32_t symbol_bit(char c) {
        return uint32_t{1} << (c - margin);
    }

    int child(int state, uint32_t bit) const {
        const auto& node = nodes[state];
        if ((node.symbols & bit) == 0) return 0;
        if (node.symbols == bit) return node.children;
        int index = __builtin_popcount(node.symbols & (bit - 1));
        return edges[node.children + index];
    }

    template <typename Visit>
    void for_each_child(int state, Visit visit) const {
        const auto& node = nodes[state];
        auto symbols = node.symbols;
        if (symbols == 0) return;
        const int* edge = (symbols & (symbols - 1)) == 0 ? &node.children : edges.data() + node.children;
        for (; symbols != 0; symbols &= symbols - 1) {
            // 最下位のビットが、現在の辺の文字を表す。
            uint32_t bit = symbols & -symbols;
            visit(bit, *edge++);
        }
    }

    int terminal_count(int state) const {
        int result = 0;
        for (int pattern = nodes[state].last_pattern; pattern != -1; pattern = patterns[pattern].previous) {
            ++result;
        }
        return result;
    }

public:
    int add(string_view pattern, int id) {
        int state = 0;
        for (char c : pattern) {
            auto bit = symbol_bit(c);
            uint32_t symbols = nodes[state].symbols;
            int first = nodes[state].children;
            if (symbols == bit) {
                state = first;
                continue;
            }
            int index = __builtin_popcount(symbols & (bit - 1));
            if ((symbols & bit) != 0) {
                state = edges[first + index];
                continue;
            }

            int degree = __builtin_popcount(symbols);
            if (degree != 0 && (degree & (degree - 1)) == 0) {
                int next_first = static_cast<int>(edges.size());
                edges.resize(edges.size() + degree * 2);
                if (degree == 1) edges[next_first] = first;
                else copy_n(edges.begin() + first, degree, edges.begin() + next_first);
                first = next_first;
            }

            int next_state = size();
            nodes.emplace_back(state);
            if (degree == 0) {
                first = next_state;
            } else {
                auto row = edges.begin() + first;
                move_backward(row + index, row + degree, row + degree + 1);
                edges[first + index] = next_state;
            }
            nodes[state].children = first;
            nodes[state].symbols = symbols | bit;
            state = next_state;
        }

        int pattern_index = static_cast<int>(patterns.size());
        patterns.push_back({id, nodes[state].last_pattern});
        nodes[state].last_pattern = pattern_index;
        return state;
    }

    int move(int state, char c) const {
        if (state < 0) return -1;
        int target = child(state, symbol_bit(c));
        return target == 0 ? -1 : target;
    }

    int find(string_view pattern) const {
        int state = 0;
        for (char c : pattern) {
            state = move(state, c);
            if (state < 0) break;
        }
        return state;
    }

    int size() const { return static_cast<int>(nodes.size()); }
    int par(int state) const { return state <= 0 ? -1 : nodes[state].parent; }

    int idx(int state) const {
        if (state < 0 || nodes[state].last_pattern == -1) return -1;
        return patterns[nodes[state].last_pattern].id;
    }

    vector<int> idxs(int state) const {
        vector<int> ids;
        if (state < 0) return ids;
        for (int pattern = nodes[state].last_pattern; pattern != -1; pattern = patterns[pattern].previous) {
            ids.push_back(patterns[pattern].id);
        }
        reverse(ids.begin(), ids.end());
        return ids;
    }
};
