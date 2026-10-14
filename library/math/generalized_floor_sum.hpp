#pragma once
using i128 = __int128_t;
struct FloorSum2 {
    i128 sum_floor;    // f_{0,1} = Σ floor((a*i+b)/m)
    i128 sum_i_floor;  // f_{1,1} = Σ i * floor((a*i+b)/m)
    i128 sum_floor2;   // f_{0,2} = Σ floor((a*i+b)/m)^2
};

i128 sum1_i128(i128 n) {
    return n * (n - 1) / 2;
}

i128 sum2_i128(i128 n) {
    return n * (n - 1) * (2 * n - 1) / 6;
}
// 0 <= n, 1 <= m, 0 <= a, 0 <= b を想定
FloorSum2 floor_sum2(long long n, long long m, long long a, long long b) {
    auto dfs = [&](auto&& dfs, long long n, long long m, long long a, long long b) -> FloorSum2 {
        if (n == 0) return {0, 0, 0};
        i128 N = n;
        i128 S1 = sum1_i128(N); 
        i128 S2 = sum2_i128(N); 
        long long qa = a / m;
        long long qb = b / m;
        long long ra = a % m;
        long long rb = b % m;
        if (qa != 0 || qb != 0) {
            auto g = dfs(dfs, n, m, ra, rb);
            i128 A = qa;
            i128 B = qb;
            FloorSum2 res;
            res.sum_floor = A * S1 + B * N + g.sum_floor;
            res.sum_i_floor = A * S2 + B * S1 + g.sum_i_floor;
            res.sum_floor2 =
                A * A * S2
                + 2 * A * B * S1
                + B * B * N
                + 2 * A * g.sum_i_floor
                + 2 * B * g.sum_floor
                + g.sum_floor2;
            return res;
        }
        if (a == 0) {
            return {0, 0, 0};
        }
        long long k = (long long)(((i128)a * (n - 1) + b) / m);
        if (k == 0) {
            return {0, 0, 0};
        }
        long long nb = m - b + a - 1;
        auto g = dfs(dfs, k, a, m, nb);
        i128 K = k;
        FloorSum2 res;
        res.sum_floor = N * K - g.sum_floor;
        // f_{1,1}(N,M,A,B)
        // = S1(N)*K + (f_{0,1}' - f_{0,2}') / 2
        res.sum_i_floor = S1 * K + (g.sum_floor - g.sum_floor2) / 2;
        // f_{0,2}(N,M,A,B)
        // = N*K^2 - f_{0,1}' - 2*f_{1,1}'
        res.sum_floor2 = N * K * K - g.sum_floor - 2 * g.sum_i_floor;
        return res;
    };

    return dfs(dfs, n, m, a, b);
}


// using ll = long long;
// using i128 = __int128_t;

// // T は atcoder::modint998244353 などを想定
// template <class T>
// struct GeneralizedFloorSum {
//     int D;
//     vector<vector<T>> C;
//     vector<vector<T>> s_coef;

//     GeneralizedFloorSum(int max_degree) : D(max_degree) {
//         build_binom();
//         build_faulhaber();
//     }

//     void build_binom() {
//         int N = D + 2;
//         C.assign(N + 1, vector<T>(N + 1, T(0)));

//         for (int i = 0; i <= N; i++) {
//             C[i][0] = C[i][i] = T(1);
//             for (int j = 1; j < i; j++) {
//                 C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
//             }
//         }
//     }

//     // S_p(n) = sum_{i=0}^{n-1} i^p = sum_j s_coef[p][j] n^j
//     void build_faulhaber() {
//         s_coef.assign(D + 1, vector<T>(D + 2, T(0)));

//         for (int p = 0; p <= D; p++) {
//             auto &c = s_coef[p];

//             // coefficient of n^p in S_p(n+1)-S_p(n) is 1
//             c[p + 1] = T(1) / T(p + 1);

//             for (int d = p - 1; d >= 0; d--) {
//                 T sum = T(0);

//                 for (int j = d + 2; j <= p + 1; j++) {
//                     sum += c[j] * C[j][d];
//                 }

//                 c[d + 1] = -sum / T(d + 1);
//             }

//             c[0] = T(0);
//         }
//     }

//     T pow_int(T a, int e) const {
//         T r = T(1);
//         while (e > 0) {
//             if (e & 1) r *= a;
//             a *= a;
//             e >>= 1;
//         }
//         return r;
//     }

//     vector<T> powers(T x, int n) const {
//         vector<T> p(n + 1);
//         p[0] = T(1);
//         for (int i = 0; i < n; i++) p[i + 1] = p[i] * x;
//         return p;
//     }

//     T sum_power(int p, ll n) const {
//         vector<T> pw(p + 2);
//         pw[0] = T(1);

//         T N = T(n);
//         for (int i = 0; i <= p; i++) {
//             pw[i + 1] = pw[i] * N;
//         }

//         T res = T(0);
//         for (int j = 0; j <= p + 1; j++) {
//             res += s_coef[p][j] * pw[j];
//         }
//         return res;
//     }

//     vector<vector<T>> zero_table() const {
//         return vector<vector<T>>(D + 1, vector<T>(D + 1, T(0)));
//     }

//     // ret[p][q] = sum_{i=0}^{n-1} i^p floor((a*i+b)/m)^q
//     // ただし p+q <= D の範囲を計算する
//     vector<vector<T>> calc(ll n, ll m, ll a, ll b) const {
//         assert(n >= 0);
//         assert(m > 0);
//         assert(a >= 0);
//         assert(b >= 0);

//         auto res = zero_table();

//         vector<T> S(D + 1);
//         for (int p = 0; p <= D; p++) {
//             S[p] = sum_power(p, n);
//             res[p][0] = S[p];
//         }

//         if (n == 0) return res;

//         // floor((a*i+b)/m)
//         // = qa*i + qb + floor((ra*i+rb)/m)
//         if (a >= m || b >= m) {
//             ll qa = a / m;
//             ll qb = b / m;
//             ll ra = a % m;
//             ll rb = b % m;

//             auto g = calc(n, m, ra, rb);

//             auto pow_qa = powers(T(qa), D);
//             auto pow_qb = powers(T(qb), D);

//             for (int p = 0; p <= D; p++) {
//                 for (int q = 1; p + q <= D; q++) {
//                     T val = T(0);

//                     for (int r = 0; r <= q; r++) {
//                         // y^r
//                         for (int s = 0; s <= q - r; s++) {
//                             // (qa*i)^s * qb^(q-r-s)
//                             val += C[q][r]
//                                  * C[q - r][s]
//                                  * pow_qa[s]
//                                  * pow_qb[q - r - s]
//                                  * g[p + s][r];
//                         }
//                     }

//                     res[p][q] = val;
//                 }
//             }

//             return res;
//         }

//         // ここから 0 <= a,b < m
//         if (a == 0) {
//             // floor(b/m) = 0
//             return res;
//         }

//         ll k = (ll)(((i128)a * (n - 1) + b) / m);

//         if (k == 0) {
//             return res;
//         }

//         // d(j+1) = ceil((m*(j+1)-b)/a)
//         //        = floor((m*j + m + a - b - 1) / a)
//         ll nb = m + a - b - 1;

//         auto g = calc(k, a, m, nb);

//         auto pow_k = powers(T(k), D);

//         for (int p = 0; p <= D; p++) {
//             for (int q = 1; p + q <= D; q++) {
//                 T val = pow_k[q] * S[p];

//                 // subtract:
//                 // sum_{r=0}^{p+1} s_p(r)
//                 //   sum_{t=0}^{q-1} C(q,t) f_{t,r}(k,a,m,nb)
//                 for (int r = 0; r <= p + 1; r++) {
//                     T coef = s_coef[p][r];

//                     for (int t = 0; t < q; t++) {
//                         val -= coef * C[q][t] * g[t][r];
//                     }
//                 }

//                 res[p][q] = val;
//             }
//         }

//         return res;
//     }
// };

// // ret[p][q] を 0 <= p <= max_p, 0 <= q <= max_q だけ返す
// template <class T>
// vector<vector<T>> generalized_floor_sum(ll n, ll m, ll a, ll b, int max_p, int max_q) {
//     int D = max_p + max_q;

//     GeneralizedFloorSum<T> solver(D);
//     auto all = solver.calc(n, m, a, b);

//     vector<vector<T>> ret(max_p + 1, vector<T>(max_q + 1));

//     for (int p = 0; p <= max_p; p++) {
//         for (int q = 0; q <= max_q; q++) {
//             ret[p][q] = all[p][q];
//         }
//     }

//     return ret;
// }