#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

namespace matrix_internal {
// 浮動小数点の消去で生じる丸め誤差の目安を、要素ごとに引き継ぐ。
// 誤差と区別できない残差はゼロとする。厳密なランク判定を保証するものではない。
// modint などの型では誤差用の領域を確保せず、通常の四則演算だけを行う。
template <class T>
struct EliminationTolerance {
    std::vector<std::vector<T>> scale;
    T relative = T(0);

    template <class Mat>
    EliminationTolerance(const Mat &a, int h, int w) {
        if constexpr (std::is_floating_point<T>::value) {
            relative = (T(8) * std::numeric_limits<T>::epsilon()) * T(h);
            scale.assign(h, std::vector<T>(w));
            for(int i = 0; i < h; ++i) for(int j = 0; j < w; ++j)
                scale[i][j] = std::abs(a[i][j]);
        }
    }
    void swap_rows(int i, int j) {
        if constexpr (std::is_floating_point<T>::value) std::swap(scale[i], scale[j]);
    }
    T entry_scale(int row, int col) const {
        if constexpr (std::is_floating_point<T>::value) return scale[row][col];
        else return T(0);
    }
    T quotient_scale(const T &quotient, const T &denominator,
                     const T &numerator_scale, int row, int col) const {
        if constexpr (std::is_floating_point<T>::value) {
            // 割る値の誤差も、消去に使う係数へ伝わる。
            return (numerator_scale + std::abs(quotient) * scale[row][col]) / std::abs(denominator);
        } else return T(0);
    }
    void multiply(T &value, int row, int col, const T &factor, const T &factor_scale) {
        if constexpr (std::is_floating_point<T>::value)
            scale[row][col] = scale[row][col] * std::abs(factor) + std::abs(value) * factor_scale;
        value *= factor;
    }
    void subtract_product(T &value, const T &a, const T &b, const T &factor_scale,
                          int row, int col, int pivot) {
        value -= a * b;
        if constexpr (std::is_floating_point<T>::value) {
            // 過去の打ち消し合いも含め、演算に使った値の大きさを引き継ぐ。
            // 要素ごとの基準なので、未変更の小さな非ゼロ値は消さない。
            scale[row][col] += scale[pivot][col] * std::abs(b) + std::abs(a) * factor_scale;
            if(std::abs(value) <= relative * scale[row][col]) value = T(0);
        }
    }
};
}
