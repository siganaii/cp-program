#pragma once
// MST を含む候補辺を重み順に返す。木を得るには Kruskal 等を適用する。
// 整数座標は64bitまで。中間計算を128bitで行う。
// 浮動小数点は有限の2進数、最近接丸め、fast-math 無効が前提。
// 中間計算は有限で、候補辺の距離が T に収まることを要求する。
// 候補辺の選択は丸め誤差を補償する。返す重みには T の丸め誤差がある。
template <typename T>
vector<tuple<T, int, int>> manhattan_mst(vector<T> xs, vector<T> ys) {
    assert(xs.size() == ys.size());
    using Wide = conditional_t<numeric_limits<T>::is_integer, __int128_t, T>;
    auto sum_key = [](Wide a, Wide b) {
        if constexpr (numeric_limits<T>::is_integer) return a + b;
        else {
            Wide hi = a + b;
            // TwoSum: hi + lo は a + b の厳密な値。辞書順で比較できる。
            Wide z = hi - a;
            Wide lo = (a - (hi - z)) + (b - z);
            return pair<Wide, Wide>{hi, lo};
        }
    };
    vector<Wide> x(xs.begin(), xs.end()), y(ys.begin(), ys.end());
    const int n = xs.size();
    using Key = decltype(sum_key(Wide{}, Wide{}));
    vector<Key> sums(n), differences(n);
    vector<int> idx(n);
    iota(idx.begin(),idx.end(),0);
    vector<tuple<T,int,int>> ret;
    for (int s = 0; s < 2; s++) {
        for (int t = 0; t < 2; t++) {
            for(int i = 0; i < n; ++i) {
                sums[i] = sum_key(x[i], y[i]);
                differences[i] = sum_key(x[i], -y[i]);
            }
            auto cmp = [&](int i, int j) {return sums[i] < sums[j];};
            sort(idx.begin(),idx.end(),cmp);
            map<Wide,int> sweep;
            for (int i:idx) {
                for (auto it = sweep.lower_bound(-y[i]); it != sweep.end(); it = sweep.erase(it)) {
                    int j = it->second;
                    if (differences[i] < differences[j]) break;
                    Wide dx = x[i] - x[j], dy = y[i] - y[j];
                    Wide distance = (dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy);
                    ret.emplace_back(static_cast<T>(distance),i,j);
                }
                sweep[-y[i]] = i;
            }
            swap(x,y);
        }
        for (auto &value:x) value = -value;
    }
    sort(ret.begin(),ret.end());
    return ret;
}
