#pragma once
// 偶数頂点の無向完全グラフをハミルトンパスに分解
vector<vector<int>> to_paths(int N) {
    assert(N >= 0 && N % 2 == 0);
    int K = N / 2;
    vector<vector<int>> paths(K,vector<int>(N));
    for(int k = 0;k < K;k++) {
        for(int i = 0;i < N;i++) {
            auto &x = paths[k][i];
            x = k + (i & 1 ? (i + 1) / 2 : -(i / 2));
            if(x < 0) x += N;
            if(x >= N) x -= N;
        }
    }
    return paths;
}

// 奇数頂点の無向完全グラフをハミルトンサイクルに分解
vector<vector<int>> to_cycles(int N) {
    assert(N >= 0 && N % 2 == 1);
    auto paths = to_paths(N - 1);
    for (auto &&P: paths) P.emplace_back(N - 1);
    return paths;
}
// 偶数頂点の無向完全グラフをマッチングに分解
// 辺彩色しているといってもよい
// 異なる2個の和は頂点素な閉路の集合になるが、単一のハミルトンサイクルとは限らない。
vector<vector<pair<int,int>>> to_matchings(int N) {
    assert(N > 0 && N % 2 == 0);
    vector<vector<pair<int,int>>> res(N - 1);
    const int mod = N - 1;
    for(int i = 0;i < mod;i++) {
        res[i].reserve(N / 2);
        res[i].emplace_back(N - 1,i);
        int x = i, y = i;
        for(int j = 0;j < N / 2 - 1;j++) {
            --x, ++y;
            if (x < 0) x += N - 1;
            if (y >= N - 1) y -= N - 1;
            res[i].emplace_back(x,y);
        }
    }
    return res;
}
// 偶数頂点の無向完全グラフを、ハミルトンサイクルと 1 つのマッチングに分解
pair<vector<vector<int>>, vector<int>> to_cycles_and_matching(ll N) {
    assert(N > 0 && N % 2 == 0);
    vector<vector<int>> cycles(N/2-1,vector<int>(N));
    int mod = N - 1;
    for(int i = 0;i < mod - 1;i++) if (i % 2 == 0) {
        auto &C = cycles[i/2];
        // 和が 2a, 2a+2 のどちらかになるときに結ぶ
        C[0] = i;
        for(int j = 0;j < N - 2;j++) {
            int nxt = (j % 2 == 0 ? 2 * i + 2 - C[j] : 2 * i - C[j]);
            if (nxt < 0) nxt += mod;
            if (nxt >= mod) nxt -= mod;
            C[j + 1] = nxt;
        }
        C[N - 1] = mod;
    }
    vector<int> match(N);
    for(int i = 0;i < mod / 2;i++) {
        int j = mod - 2 - i;
        match[2 * i] = i, match[2 * i + 1] = j;
    }
    match[N - 2] = N - 2, match[N - 1] = N - 1;
    return {cycles, match};
}
