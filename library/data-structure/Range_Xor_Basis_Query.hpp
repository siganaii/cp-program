#pragma once
// 列もクエリも static
template <typename INT_TYPE, int LOG>
struct Range_Xor_Basis_Query {
    int n;
    int Q;
    vector<INT_TYPE>& A;
    vector<vector<pair<int, int>>> QUERY;

    Range_Xor_Basis_Query(vector<INT_TYPE>& A) : n(A.size()), Q(0), A(A) {
        QUERY.resize(n + 1);
    }

    // クエリ追加
    void query(int L, int R) {
        assert(0 <= L && L <= R && R <= n);
        QUERY[R].emplace_back(Q++, L);
    }

    // 関数 F(qid, basis) を渡す
    template <typename CALC>
    void calc(CALC F) {
        using P = pair<INT_TYPE, int>;
        array<P, LOG> d;
        for(int k = 0;k < LOG;++k) d[k] = {INT_TYPE(1) << k, -1};
        vector<INT_TYPE> basis;
        for(int i = 0;i < n + 1;++i) {
            for (auto&& [qid, L]: QUERY[i]) {
                basis.clear();
                for(int k = LOG - 1;k >= 0;--k) if (d[k].second >= L) basis.emplace_back(d[k].first);
                F(qid, basis);
            }
            if (i == n) break;
            P p = {A[i], i};
            for(int k = LOG - 1;k >= 0;--k) {
                if (!(p.first >> k & 1)) continue;
                if (p.second > d[k].second) swap(p, d[k]);
                p.first ^= d[k].first;
            }
        }
    }
};