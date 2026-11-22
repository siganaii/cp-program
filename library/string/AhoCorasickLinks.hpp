#pragma once
#include "library/string/trie.hpp"

// 実在する辺を格納する Aho-Corasick。文字は [margin, margin + X) を使う。
template <size_t X = 26, char margin = 'a'>
class AhoCorasickLinks : private SparseTrie<X, margin> {
    using TRIE = SparseTrie<X, margin>;
    using TRIE::child;
    using TRIE::for_each_child;
    using TRIE::symbol_bit;
    using TRIE::terminal_count;

    int transition(int state, uint32_t bit) const {
        int target = child(state, bit);
        while (state != 0 && target == 0) {
            state = suffix_link[state];
            target = child(state, bit);
        }
        return target;
    }

public:
    vector<int> cnt;
    vector<int> last;
    vector<int> bfs_order;
    vector<int> suffix_link;
    bool built = false;

    using TRIE::idx;
    using TRIE::idxs;
    using TRIE::par;
    using TRIE::size;

    int count(int state) const { return cnt[state]; }

    int move(int state, char c) const {
        if (!built) return TRIE::move(state, c);
        if (state < 0) return -1;
        return transition(state, symbol_bit(c));
    }

    int find(string_view pattern) const {
        if (!built) return TRIE::find(pattern);
        int state = 0;
        for (char c : pattern) state = transition(state, symbol_bit(c));
        return state;
    }

    void restore_trie() { built = false; }

    // ID を指定したパターンは終端状態に記録する。last には自動採番の分を登録する。
    int add(string_view pattern, int id) {
        restore_trie();
        return TRIE::add(pattern, id);
    }

    int add(string_view pattern) {
        int id = static_cast<int>(last.size());
        last.push_back(0);
        try {
            last.back() = add(pattern, id);
        } catch (...) {
            last.pop_back();
            throw;
        }
        return last.back();
    }

    // 検索と失敗リンク・一致数の参照は build() 後に行う。
    void build() {
        if (built) return;
        int n = size();
        cnt.resize(n);
        suffix_link.assign(n, 0);
        bfs_order.clear();
        bfs_order.reserve(n - 1);
        cnt[0] = terminal_count(0);
        for_each_child(0, [&](uint32_t, int child) {
            bfs_order.push_back(child);
        });

        for (size_t i = 0; i < bfs_order.size(); ++i) {
            int state = bfs_order[i];
            int failure = suffix_link[state];
            cnt[state] = terminal_count(state) + cnt[failure];
            for_each_child(state, [&](uint32_t bit, int child) {
                suffix_link[child] = transition(failure, bit);
                bfs_order.push_back(child);
            });
        }
        built = true;
    }

    // build() 後に呼ぶ。heavy が真なら自動採番したパターンごとの出現数、偽なら全パターンの合計を返す。
    // 空パターンは文字列の各境界に一致し、出現数は長さ + 1。
    vector<int> match(string_view text, int heavy = true) const {
        int state = 0;
        if (!heavy) {
            int total = cnt[0];
            for (char c : text) {
                state = transition(state, symbol_bit(c));
                total += cnt[state];
            }
            return {total};
        }
        if (last.empty()) return {};
        vector<int> visits(size());
        visits[0] = 1;
        for (char c : text) {
            state = transition(state, symbol_bit(c));
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
};
