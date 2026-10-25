#pragma once
#include "library/modulo/binomial.hpp"
#include "library/misc/Mo.hpp"

// 各 (n, m) に対して sum_{i=0}^m C(n, i) を返す。0 <= m <= n とする。
// mod は素数とし、mod != 2 の場合は n < mod とする。
template <typename mint>
vector<mint> multipoint_binomial_sum(const vector<pair<int, int>>& qs) {
    if (qs.empty()) return {};
    int N = 0;
    for (auto& p : qs) {
        assert(0 <= p.second && p.second <= p.first);
        N = max(N, p.first);
    }
    if (mint::get_mod() == 2) {
        // n > 0 なら部分和の偶奇は C(n-1, m) と一致する。
        // C(a, b) が奇数である条件は、b の立つビットがすべて a に含まれること。
        vector<mint> ans;
        ans.reserve(qs.size());
        for (auto [n, m] : qs) ans.emplace_back(n == 0 || (m & (n - 1)) == m);
        return ans;
    }
    Binomial<mint> b(N);
    int Q = qs.size();
    Mo mo(N, Q);
    for (auto& p : qs) {
        assert(p.second <= p.first);
        assert(p.first <= N);
        mo.insert(p.second, p.first);
    }
    vector<mint> ans(Q);
    mint cur = 1;
    int n = 0, m = 0;
    auto al = [&](int) { cur -= b.C(n, m--); };
    auto ar = [&](int) { cur += cur - b.C(n++, m); };
    auto el = [&](int) { cur += b.C(n, ++m); };
    auto er = [&](int) { cur = (cur + b.C(--n, m)) * b.inv(2); };
    auto q = [&](int i) { ans[i] = cur; };
    mo.run(al, ar, el, er, q);
    return ans;
}
