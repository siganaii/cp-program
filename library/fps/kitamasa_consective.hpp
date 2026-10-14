// return A[i] if exists, otherwise 0
template <typename mint>
mint fps_get(const FormalPowerSeries<mint>& A, int i) {
    if (i < 0 || i >= (int)A.size()) return mint(0);
    return A[i];
}

// A(x) -> A_even(x) = A_0 + A_2 x + A_4 x^2 + ...
template <typename mint>
FormalPowerSeries<mint> fps_even(const FormalPowerSeries<mint>& A) {
    FormalPowerSeries<mint> B;
    for (int i = 0; i < (int)A.size(); i += 2) {
        B.push_back(A[i]);
    }
    B.shrink();
    return B;
}

// A(x) -> A_odd(x) = A_1 + A_3 x + A_5 x^2 + ...
template <typename mint>
FormalPowerSeries<mint> fps_odd(const FormalPowerSeries<mint>& A) {
    FormalPowerSeries<mint> B;
    for (int i = 1; i < (int)A.size(); i += 2) {
        B.push_back(A[i]);
    }
    B.shrink();
    return B;
}

// A(x) -> A(-x)
template <typename mint>
FormalPowerSeries<mint> fps_neg_x(FormalPowerSeries<mint> A) {
    for (int i = 1; i < (int)A.size(); i += 2) {
        A[i] = -A[i];
    }
    return A;
}

// small r 用の愚直展開
// return [x^l], ..., [x^{r-1}] P(x) / Q(x)
// assumes 0 <= l <= r
template <typename mint>
FormalPowerSeries<mint> NaiveRangeRational(
    FormalPowerSeries<mint> P,
    FormalPowerSeries<mint> Q,
    long long l,
    long long r
) {
    vector<mint> coef(r);
    mint inv_q0 = Q[0].inverse();

    for (int i = 0; i < r; i++) {
        mint s = fps_get(P, i);

        for (int j = 1; j < (int)Q.size() && j <= i; j++) {
            s -= Q[j] * coef[i - j];
        }

        coef[i] = s * inv_q0;
    }

    FormalPowerSeries<mint> res;
    for (long long i = l; i < r; i++) {
        res.push_back(coef[i]);
    }
    return res;
}

// deg(P) < deg(Q) を仮定した区間版 Bostan-Mori
// return [x^l], ..., [x^{r-1}] P(x) / Q(x)
template <typename mint>
FormalPowerSeries<mint> RangeRationalProper(
    FormalPowerSeries<mint> P,
    FormalPowerSeries<mint> Q,
    long long l,
    long long r
) {
    if (l >= r) return {};

    // 負次数の係数は 0
    if (r <= 0) {
        return FormalPowerSeries<mint>(r - l, mint(0));
    }

    if (l < 0) {
        FormalPowerSeries<mint> res(-l, mint(0));
        auto right = RangeRationalProper(P, Q, 0, r);
        res.insert(res.end(), right.begin(), right.end());
        return res;
    }

    P.shrink();
    Q.shrink();

    if (P.empty()) {
        return FormalPowerSeries<mint>(r - l, mint(0));
    }

    assert(!Q.empty());
    assert(Q[0] != mint(0));

    static constexpr long long THRESHOLD = 64;

    if (r <= THRESHOLD) {
        return NaiveRangeRational(P, Q, l, r);
    }

    FormalPowerSeries<mint>::set_fft();

    // Q2(x) = Q(-x)
    auto Q2 = fps_neg_x(Q);

    // S(x) = P(x) Q(-x)
    auto S = P * Q2;

    // T(x) = Q(x) Q(-x) = V(x^2)
    auto T = Q * Q2;

    // V(x)
    auto V = fps_even(T);

    // P_even(x), P_odd(x)
    auto Pe = fps_even(S);
    auto Po = fps_odd(S);

    // even index: i = 2m
    // 2m in [l, r) <=> m in [ceil(l/2), ceil(r/2))
    long long el = (l + 1) / 2;
    long long er = (r + 1) / 2;

    // odd index: i = 2m + 1
    // 2m+1 in [l, r) <=> m in [floor(l/2), floor(r/2))
    long long ol = l / 2;
    long long orr = r / 2;

    auto even_coef = RangeRationalProper(Pe, V, el, er);
    auto odd_coef  = RangeRationalProper(Po, V, ol, orr);

    FormalPowerSeries<mint> res(r - l);

    int ei = 0;
    int oi = 0;

    for (long long i = l; i < r; i++) {
        if ((i & 1) == 0) {
            res[i - l] = even_coef[ei++];
        } else {
            res[i - l] = odd_coef[oi++];
        }
    }

    return res;
}

// return [x^l], ..., [x^{r-1}] P(x) / Q(x)
template <typename mint>
FormalPowerSeries<mint> LinearRecurrenceRange(
    long long l,
    long long r,
    FormalPowerSeries<mint> Q,
    FormalPowerSeries<mint> P
) {
    if (l >= r) return {};

    Q.shrink();
    P.shrink();

    assert(!Q.empty());
    assert(Q[0] != mint(0));

    FormalPowerSeries<mint> ret(r - l, mint(0));

    // improper fraction の多項式部分を先に処理
    if ((int)P.size() >= (int)Q.size()) {
        auto R = P / Q;
        P -= R * Q;
        P.shrink();

        for (long long i = l; i < r; i++) {
            if (0 <= i && i < (long long)R.size()) {
                ret[i - l] += R[i];
            }
        }
    }

    if (P.empty()) {
        return ret;
    }

    auto frac = RangeRationalProper(P, Q, l, r);

    for (int i = 0; i < (int)frac.size(); i++) {
        ret[i] += frac[i];
    }

    return ret;
}


