#pragma once
template <typename T>
void walsh_hadamard_transform(vector<T>& f, bool inverse = false) {
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
        if constexpr (is_integral<T>::value) {
            for (auto& x : f) x /= n;
        } 
        else {
            T invn = T(1) / T(f.size());
            for (auto& x : f) x *= invn;
        }
    }
}
