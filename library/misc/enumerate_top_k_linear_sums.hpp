#pragma once
// c * x + y の大きい順に K 個選んだときの
// {sum_x, sum_y} の候補を列挙する。
//
// 返り値は c = +inf から c = -inf へ動かした順。
// 計算量: O(N^2 log N)
// 空間計算量: O(N^2)
//
// 座標差同士の積は __int128_t に収まることを仮定する。
// sum_x, sum_y は long long に収まることを仮定する。
vector<pair<long long, long long>>
enumerate_top_k_linear_sums(
    vector<pair<long long, long long>> xy,
    int k
) {
    using ll = long long;
    using i128 = __int128_t;
    using P = pair<ll, ll>;
    const int n = (int)xy.size();
    assert(0 <= k && k <= n);
    if (k == 0) return {{0, 0}};
    // c = +inf では x の降順。
    // x が等しい場合は、常に y の大きい方が上になる。
    sort(xy.begin(), xy.end(), [](const P& a, const P& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second > b.second;
    });
    if (k == n) {
        // 最終的な合計が収まっていても、途中の部分和は超えることがある。
        i128 sx = 0, sy = 0;
        for (auto [x, y] : xy) {
            sx += x;
            sy += y;
        }
        return {{(ll)sx, (ll)sy}};
    }
    struct Event {
        int i, j;
    };
    vector<Event> events;
    events.reserve((size_t)n * (n - 1) / 2);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // x が等しい2点は、cを動かしても順序が変化しない。
            if (xy[i].first == xy[j].first) continue;
            events.push_back({i, j});
        }
    }
    // i,j の順位が入れ替わる c は
    //
    // c * x_i + y_i = c * x_j + y_j
    // c = -(y_i-y_j)/(x_i-x_j)
    //
    // c を降順に見ることは、
    // (y_i-y_j)/(x_i-x_j) を昇順に見ることに相当する。
    auto slope_less = [&](const Event& a, const Event& b) {
        i128 adx = (i128)xy[a.i].first  - xy[a.j].first;
        i128 ady = (i128)xy[a.i].second - xy[a.j].second;
        i128 bdx = (i128)xy[b.i].first  - xy[b.j].first;
        i128 bdy = (i128)xy[b.i].second - xy[b.j].second;
        return ady * bdx < bdy * adx;
    };
    // 同じ傾きのイベントは、生成時の辞書順を維持して処理する。
    // 3点以上が同時に交差する場合に必要。
    stable_sort(events.begin(), events.end(), slope_less);
    vector<int> pos(n);
    iota(pos.begin(), pos.end(), 0);
    i128 sx = 0;
    i128 sy = 0;
    for (int i = 0; i < k; ++i) {
        sx += xy[i].first;
        sy += xy[i].second;
    }
    vector<P> result;
    result.emplace_back((ll)sx, (ll)sy);
    for (const auto& event : events) {
        int i = event.i;
        int j = event.j;
        // 同じ傾きのイベントを順に処理すると、
        // 既にこの2点の順序が反転していることがある。
        if (pos[i] > pos[j]) continue;
        bool in_i = pos[i] < k;
        bool in_j = pos[j] < k;
        // 片方だけが上位K個に含まれている場合、
        // 選択集合が i -> j に変化する。
        if (in_i != in_j) {
            // 座標差も long long の範囲を超えうるので、引く前に拡張する。
            sx += (i128)xy[j].first - xy[i].first;
            sy += (i128)xy[j].second - xy[i].second;
            P sum{(ll)sx, (ll)sy};
            if (result.back() != sum) {
                result.push_back(sum);
            }
        }
        swap(pos[i], pos[j]);
    }
    return result;
}
