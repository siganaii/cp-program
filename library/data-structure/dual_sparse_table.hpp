#pragma once
// composition は結合則・可換性・冪等性を満たし、e() はその単位元であること。
// 例: min / max / bitwise OR。加算は更新が重複するため使えない。
// apply をすべて実行してから build し、get で結果を取得する。
template<class S,S (*e)(),S (*composition)(S, S)>
struct dual_sparse_table {
    std::vector<std::vector<S>> table;
    int n2;
    dual_sparse_table (): dual_sparse_table(0) {}
    dual_sparse_table (const int _n): dual_sparse_table(vector<S>(_n, e())) {}
    dual_sparse_table(const vector<S>& v) {
        n2 = 1; // 空配列でも table[0] を確保する。
        int n = v.size();
        while ((1LL<<n2) <= n) n2++;
        table.resize(n2);
        table[0] = v;
        for (int i = 0; i < n2-1; i++){
            int nsz = table[i].size() - (1 << i);
            table[i + 1].resize(nsz,e());
        }
    }
    void apply(int l, int r,S x) {
        assert(0 <= l && l < r && r <= (int)(table[0].size()));
        int lgs = 31 - __builtin_clz((unsigned int)(r-l));
        table[lgs][l] = composition(x,table[lgs][l]);
        table[lgs][r - (1 << lgs)] = composition(x,table[lgs][r - (1 << lgs)]);
    }
    void build() {
        for (int i = n2 - 2;i >= 0;i--) {
            for(int j = 0,si = table[i+1].size();j < si;j++) {
                table[i][j] = composition(table[i][j],table[i + 1][j]);
                table[i][j + (1 << i)] = composition(table[i][j + (1 << i)],table[i + 1][j]);
            }
        }
    }
    S get(int x) {
        return table[0][x];
    }
};
