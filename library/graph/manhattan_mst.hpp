#pragma once
// MST を含む候補辺を重み順に返す。木を得るには Kruskal 等を適用する。
// 整数座標は64bitまで。中間計算を128bitで行う。
// 候補辺の距離が返り値の型 T に収まらない場合は overflow_error。
template <typename T>
vector<tuple<T, int, int>> manhattan_mst(vector<T> xs, vector<T> ys) {
    static_assert(!numeric_limits<T>::is_integer || sizeof(T) <= 8,
                  "manhattan_mst supports integer coordinates up to 64 bits");
    assert(xs.size() == ys.size());
    using Wide = conditional_t<numeric_limits<T>::is_integer, __int128_t, T>;
    vector<Wide> x(xs.begin(), xs.end()), y(ys.begin(), ys.end());
    const int n = xs.size();
    vector<int> idx(n);
    iota(idx.begin(),idx.end(),0);
    vector<tuple<T,int,int>> ret;
    for (int s = 0; s < 2; s++) {
        for (int t = 0; t < 2; t++) {
            auto cmp = [&](int i, int j) {return x[i] + y[i] < x[j] + y[j];};
            sort(idx.begin(),idx.end(),cmp);
            map<Wide,int> sweep;
            for (int i:idx) {
                for (auto it = sweep.lower_bound(-y[i]); it != sweep.end(); it = sweep.erase(it)) {
                    int j = it->second;
                    if (x[i] - x[j] < y[i] - y[j]) break;
                    Wide dx = x[i] - x[j], dy = y[i] - y[j];
                    Wide distance = (dx < 0 ? -dx : dx) + (dy < 0 ? -dy : dy);
                    if constexpr (numeric_limits<T>::is_integer) {
                        if(distance > numeric_limits<T>::max())
                            throw overflow_error("manhattan_mst: edge weight does not fit T");
                    }
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
