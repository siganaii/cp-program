#pragma once
//https://259-momone.hatenablog.com/entry/2025/04/28/232456
//https://atcoder.jp/contests/abc429/submissions/77489599
template <class S, S (*op)(S, S), S (*e)()>
struct floor_prod {
    using ull = unsigned long long;
    using u128 = __uint128_t;
    S x, y;
    floor_prod(S x_, S y_) : x(x_), y(y_) {}
    S pow_monoid(S a, ull k) const {
        S res = e();
        while (k) {
            if (k & 1) res = op(res, a);
            k >>= 1;
            if (k) a = op(a, a);
        }
        return res;
    }
    S solve(ull n, ull m, ull a, u128 b, S x, S y) const {
        assert(m > 0);
        x = op(x, pow_monoid(y, a / m));
        a %= m;
        S prefix = pow_monoid(y, b / m);
        b %= m;
        u128 y_max = u128(a) * n + b;
        if (y_max < m) {
            return op(prefix, pow_monoid(x, n));
        }
        ull q = ull(y_max / m);
        ull r = ull(y_max % m);
        return op(op(op(prefix,solve(q - 1,a,m,u128(m) + a - b - 1,y,x)),y),pow_monoid(x, r / a));
    }

    S operator()(ull n, ull m, ull a, ull b) const {
        return solve(n, m, a, b, x, y);
    }
};
