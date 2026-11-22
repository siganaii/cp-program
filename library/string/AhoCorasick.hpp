#pragma once
#include "library/string/trie.hpp"

template <size_t X = 26, char margin = 'a'>
struct AhoCorasick : Trie<X + 1, margin> {
    using TRIE = Trie<X + 1, margin>;
    using TRIE::next;
    using TRIE::st;
    using TRIE::TRIE;
    vector<int> cnt;  // 各状態で末尾に一致するパターン数。
    vector<int> last;  // 自動採番で追加した各パターンの終端状態。
    vector<int> bfs_order;  // 根を除いた幅優先順。
    vector<int> suffix_link;
    bool built = false;
    // 構築時に補完した遷移を取り除き、元の Trie の辺に戻す。O(頂点数 * X)。
    void restore_trie() {
        if (!built) return;
        for (auto &node : st) node.nxt.fill(-1);
        for (int state = 1; state < static_cast<int>(st.size()); ++state) {
            const auto &node = st[state];
            next(node.par, node.key - margin) = state;
        }
        built = false;
    }
    // 構築後に追加した場合は、検索前にもう一度 build() を呼ぶ。
    int add(const string &pattern, int id) {
        restore_trie();
        return TRIE::add(pattern, id);
    }
    void add(const string &pattern) {
        int state = add(pattern, static_cast<int>(last.size()));
        last.push_back(state);
    }
    // add() または restore_trie() が呼ばれるまで、構築済みの遷移を再利用する。
    void build() {
        if (built) return;

        constexpr bool copy_sparse_rows = X >= 4;
        if constexpr (copy_sparse_rows) {
            // 幅優先で一致数を書き込むまで、cnt に唯一の子を記録する。0 は子なし、-1 は複数。
            cnt.assign(st.size(), 0);
            for (int child = 1; child < static_cast<int>(st.size()); ++child) {
                int &entry = cnt[st[child].par];
                entry = entry == 0 ? child : -1;
            }
        } else {
            cnt.resize(st.size());
        }

        suffix_link.assign(st.size(), 0);
        bfs_order.clear();
        bfs_order.reserve(st.size() - 1);
        cnt[0] = static_cast<int>(st[0].idxs.size());
        next(0, X) = 0;
        for (int symbol = 0; symbol < static_cast<int>(X); ++symbol) {
            int &child = next(0, symbol);
            if (child < 0) {
                child = 0;
            } else {
                bfs_order.push_back(child);
            }
        }
        for (size_t i = 0; i < bfs_order.size(); ++i) {
            int state = bfs_order[i];
            int failure = suffix_link[state];
            auto &transitions = st[state].nxt;
            const auto &fallbacks = st[failure].nxt;
            int only_child = copy_sparse_rows ? cnt[state] : -1;
            cnt[state] = static_cast<int>(st[state].idxs.size()) + cnt[failure];

            if (only_child >= 0) {
                transitions = fallbacks;
                if (only_child != 0) {
                    int symbol = st[only_child].key - margin;
                    suffix_link[only_child] = transitions[symbol];
                    transitions[symbol] = only_child;
                    bfs_order.push_back(only_child);
                }
            } else {
                for (int symbol = 0; symbol < static_cast<int>(X); ++symbol) {
                    int &child = transitions[symbol];
                    int fallback = fallbacks[symbol];
                    if (child < 0) {
                        child = fallback;
                    } else {
                        suffix_link[child] = fallback;
                        bfs_order.push_back(child);
                    }
                }
            }
            transitions[X] = failure;  // Trie の末尾スロットにも失敗リンクを格納する。
        }
        built = true;
    }
    // build() 後に呼ぶ。heavy が真なら自動採番したパターンごとの出現数、偽なら全パターンの合計を返す。
    // 空パターンは文字列の各境界に一致し、出現数は長さ + 1。
    vector<int> match(string &text, int heavy = true) {
        int state = 0;
        if (!heavy) {
            int total = cnt[0];
            for (char c : text) {
                state = next(state, c - margin);
                total += cnt[state];
            }
            return {total};
        }
        if (last.empty()) return {};
        vector<int> visits(st.size());
        visits[0] = 1;
        for (char c : text) {
            state = next(state, c - margin);
            ++visits[state];
        }
        // 深い状態から失敗リンクをたどり、短い接尾辞にも出現数を加える。
        for (auto it = bfs_order.rbegin(); it != bfs_order.rend(); ++it) {
            int node = *it;
            int occurrences = visits[node];
            if (occurrences != 0) visits[suffix_link[node]] += occurrences;
        }
        vector<int> result;
        result.reserve(last.size());
        for (int node : last) result.push_back(visits[node]);
        return result;
    }
    int count(int state) { return cnt[state]; }
};
