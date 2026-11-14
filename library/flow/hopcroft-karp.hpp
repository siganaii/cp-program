#pragma once

// 容量1・重みなしの二部最大マッチング。左右の頂点番号は別々に0-indexed。
// V = L + R, E = 辺数として O((V + E) sqrt(V)), 空間 O(V + E)。
// 辺の追加後も再計算可能。max_matching() は追加分ではなく現在の総マッチング数。
// matching() は最大マッチングを計算して (左頂点, 右頂点) の一覧を返す。
// 多重辺・孤立頂点・空グラフに対応。探索は非再帰。
struct HopcroftKarp {
    explicit HopcroftKarp(int L, int R) : left_size(L), right_size(R) {
        assert(L >= 0 && R >= 0);
        match_left.assign(L, -1);
        match_right.assign(R, -1);
        dist.assign(L, -1);
        iter.resize(L);
        que.reserve(L);
        path.reserve(L);
    }

    void add_edge(int l, int r) {
        assert(0 <= l && l < left_size);
        assert(0 <= r && r < right_size);
        input_edges.emplace_back(l, r);
        solved = false;
    }

    int max_matching() {
        if (solved) return matching_size;
        build_graph();
        greedy_matching();
        initial_matching();
        roots.clear();
        for (int l = 0; l < left_size; ++l) {
            if (match_left[l] == -1 && start[l] != start[l + 1]) roots.push_back(l);
        }
        while (bfs()) {
            size_t remaining = 0;
            for (int l : roots) {
                if (augment(l)) ++matching_size;
                else roots[remaining++] = l;
            }
            roots.resize(remaining);
        }
        solved = true;
        return matching_size;
    }

    vector<pair<int, int>> matching() {
        max_matching();
        vector<pair<int, int>> ret;
        ret.reserve(matching_size);
        for (int l = 0; l < left_size; ++l) {
            if (match_left[l] != -1) ret.emplace_back(l, match_left[l]);
        }
        return ret;
    }

private:
    int left_size, right_size;
    int matching_size = 0, shortest = -1;
    bool solved = false;
    vector<pair<int, int>> input_edges;
    vector<int> match_left, match_right, dist, que, path;
    vector<size_t> iter, start, right_start;
    vector<int> edges, reverse_edges, roots;
    void build_graph() {
        start.assign(left_size + 1, 0);
        right_start.assign(right_size + 1, 0);
        for (auto [l, r] : input_edges) { ++start[l + 1]; ++right_start[r + 1]; }
        for (int l = 0; l < left_size; ++l) start[l + 1] += start[l];
        for (int r = 0; r < right_size; ++r) right_start[r + 1] += right_start[r];
        edges.resize(input_edges.size());
        reverse_edges.resize(input_edges.size());
        auto left_cursor = start, right_cursor = right_start;
        for (auto [l, r] : input_edges) {
            edges[left_cursor[l]++] = r;
            reverse_edges[right_cursor[r]++] = l;
        }
    }
    void greedy_matching() {
        vector<int> degree_left(left_size), degree_right(right_size);
        for (auto [l, r] : input_edges) {
            if (match_left[l] == -1 && match_right[r] == -1) {
                ++degree_left[l];
                ++degree_right[r];
            }
        }
        vector<pair<int, bool>> leaves;
        for (int l = 0; l < left_size; ++l) if (degree_left[l] == 1) leaves.emplace_back(l, false);
        for (int r = 0; r < right_size; ++r) if (degree_right[r] == 1) leaves.emplace_back(r, true);
        int next_left = 0;
        while (true) {
            int l = -1, r = -1;
            if (!leaves.empty()) {
                auto [v, right] = leaves.back();
                leaves.pop_back();
                if (right) {
                    if (match_right[v] != -1 || degree_right[v] == 0) continue;
                    r = v;
                    for (size_t i = right_start[r]; i < right_start[r + 1]; ++i) {
                        if (match_left[reverse_edges[i]] == -1) { l = reverse_edges[i]; break; }
                    }
                } else {
                    if (match_left[v] != -1 || degree_left[v] == 0) continue;
                    l = v;
                }
            } else {
                while (next_left < left_size &&
                       (match_left[next_left] != -1 || degree_left[next_left] == 0)) ++next_left;
                if (next_left == left_size) break;
                l = next_left++;
            }
            if (r == -1) {
                for (size_t i = start[l]; i < start[l + 1]; ++i) {
                    if (match_right[edges[i]] == -1) { r = edges[i]; break; }
                }
            }
            match_left[l] = r;
            match_right[r] = l;
            ++matching_size;
            for (size_t i = start[l]; i < start[l + 1]; ++i) {
                int v = edges[i];
                if (match_right[v] == -1 && --degree_right[v] == 1) leaves.emplace_back(v, true);
            }
            for (size_t i = right_start[r]; i < right_start[r + 1]; ++i) {
                int u = reverse_edges[i];
                if (match_left[u] == -1 && --degree_left[u] == 1) leaves.emplace_back(u, false);
            }
        }
    }
    void initial_matching() {
        const int infinity = left_size;
        vector<int> label(left_size), bfs_queue, pending;
        bfs_queue.reserve(left_size);
        pending.reserve(right_size);
        for (int r = 0; r < right_size; ++r) {
            if (match_right[r] == -1 && right_start[r] != right_start[r + 1]) pending.push_back(r);
        }
        size_t budget = 8 * (size_t(left_size) + right_size + reverse_edges.size());
        size_t head = 0, since_relabel = 0;
        while (head < pending.size()) {
            if (since_relabel == 0) {
                fill(label.begin(), label.end(), infinity);
                bfs_queue.clear();
                for (int l = 0; l < left_size; ++l) if (match_left[l] == -1) {
                    label[l] = 0;
                    bfs_queue.push_back(l);
                }
                for (size_t i = 0; i < bfs_queue.size(); ++i) {
                    int l = bfs_queue[i];
                    for (size_t j = start[l]; j < start[l + 1]; ++j) {
                        int r = edges[j];
                        if (budget == 0) return;
                        --budget;
                        int next = match_right[r];
                        if (next != -1 && label[next] == infinity) {
                            label[next] = label[l] + 1;
                            bfs_queue.push_back(next);
                        }
                    }
                }
            }
            int r = pending[head++], best = -1, distance = infinity;
            for (size_t i = right_start[r]; i < right_start[r + 1]; ++i) {
                if (budget == 0) return;
                --budget;
                int l = reverse_edges[i];
                if (label[l] < distance) { distance = label[l]; best = l; }
            }
            if (best != -1) {
                int old = match_left[best];
                if (old == -1) ++matching_size;
                else {
                    match_right[old] = -1;
                    pending.push_back(old);
                }
                match_right[r] = best;
                match_left[best] = r;
                label[best] = distance + 1;
            }
            if (++since_relabel == size_t(left_size) + right_size) since_relabel = 0;
        }
    }
    bool bfs() {
        for (int l : que) dist[l] = -1;
        que.clear();
        shortest = -1;
        for (int l : roots) {
            dist[l] = 0;
            iter[l] = start[l];
            que.push_back(l);
        }
        for (size_t i = 0; i < que.size(); ++i) {
            int l = que[i];
            if (shortest != -1 && dist[l] >= shortest) continue;
            for (size_t j = start[l]; j < start[l + 1]; ++j) {
                int r = edges[j];
                int next = match_right[r];
                if (next == -1) {
                    shortest = dist[l] + 1;
                } else if (dist[next] == -1) {
                    dist[next] = dist[l] + 1;
                    iter[next] = start[next];
                    que.push_back(next);
                }
            }
        }
        return shortest != -1;
    }

    bool augment(int root) {
        path.clear();
        path.push_back(root);
        while (!path.empty()) {
            int l = path.back();
            bool advanced = false;
            while (iter[l] < start[l + 1]) {
                int r = edges[iter[l]++];
                int next = match_right[r];
                if (next == -1) {
                    if (dist[l] + 1 != shortest) continue;
                    // 末端から交互路を反転する。以前の相手が一つ前の辺に対応する。
                    while (!path.empty()) {
                        int u = path.back();
                        path.pop_back();
                        int old = match_left[u];
                        match_left[u] = r;
                        match_right[r] = u;
                        r = old;
                    }
                    return true;
                }
                if (dist[next] == dist[l] + 1 && dist[next] < shortest) {
                    path.push_back(next);
                    advanced = true;
                    break;
                }
            }
            if (!advanced) {
                dist[l] = -1;
                path.pop_back();
            }
        }
        return false;
    }
};
