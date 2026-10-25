#pragma once
#include <type_traits>
// 入力長は2の冪とする。空入力は何もしない。
// modintの逆変換では入力長の逆元が必要。長さ2以上の偶数modには対応しない。
template <typename T>
void walsh_hadamard_transform(vector<T>& f, bool inverse = false) {
    if (f.empty()) return;
    int n = f.size();
    assert((n & (n - 1)) == 0);
    for (int len = 1; len < n; len <<= 1) {
        for (int l = 0; l < n; l += len << 1) {
            for (int i = 0; i < len; i++) {
                T x = f[l + i];
                T y = f[l + i + len];
                f[l + i] = x + y;
                f[l + i + len] = x - y;
            }
        }
    }
    if (inverse) {
        // GCCの標準モードでも拡張整数型を整数として扱う。
        if constexpr (is_integral<T>::value || is_same<T, __int128_t>::value ||
                      is_same<T, __uint128_t>::value) {
            for (auto& x : f) x /= n;
        } 
        else {
            T invn = T(1) / T(f.size());
            for (auto& x : f) x *= invn;
        }
    }
}
