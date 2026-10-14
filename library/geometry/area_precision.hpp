#pragma once
#include "library/geometry/geometry-base.hpp"

namespace geometry_detail {
// 面積の再計算用。小さい成分から順に2個または4個のlong doubleで値を保持する。
// 積・和の丸め残差を残し、大きな円弧同士の相殺にも余裕を持たせる。
template <int Count>
struct AreaNumber {
    static_assert(Count == 2 || Count == 4);
    array<Real, Count> parts{};
    AreaNumber(Real value = 0) { parts[Count - 1] = value; }

    template <size_t N>
    static void append(array<Real, N>& terms, int& count, Real carry) {
        if (carry == 0) return;
        int size = 0;
        for (int i = 0; i < count; i++) {
            Real sum = carry + terms[i], v = sum - carry;
            Real error = (carry - (sum - v)) + (terms[i] - v);
            if (error != 0) terms[size++] = error;
            carry = sum;
        }
        if (carry != 0) terms[size++] = carry;
        count = size;
    }
    template <size_t N>
    static AreaNumber pack(array<Real, N> terms, int count) {
        AreaNumber result;
        if (count == 0) return result;
        // 残差を詰め直してから、必要な数の上位成分を取り出す。
        Real q = terms[count - 1];
        int bottom = count - 1;
        for (int i = count - 2; i >= 0; i--) {
            Real sum = q + terms[i], error = terms[i] - (sum - q);
            if (error != 0) { terms[bottom--] = sum; q = error; }
            else q = sum;
        }
        terms[bottom] = q;
        int size = 0;
        for (int i = bottom + 1; i < count; i++) {
            Real sum = terms[i] + q, error = q - (sum - terms[i]);
            if (error != 0) terms[size++] = error;
            q = sum;
        }
        if (q != 0) terms[size++] = q;
        int take = min(size, Count);
        for (int i = 0; i < take; i++) result.parts[Count - take + i] = terms[size - take + i];
        return result;
    }
    explicit operator Real() const {
        Real result = 0;
        for (Real part : parts) result += part;
        return result;
    }
    AreaNumber operator-() const {
        AreaNumber result;
        for (int i = 0; i < Count; i++) result.parts[i] = -parts[i];
        return result;
    }
    friend AreaNumber operator+(const AreaNumber& a, const AreaNumber& b) {
        array<Real, 2 * Count> terms{};
        int count = 0;
        for (Real part : a.parts) append(terms, count, part);
        for (Real part : b.parts) append(terms, count, part);
        return pack(terms, count);
    }
    friend AreaNumber operator-(const AreaNumber& a, const AreaNumber& b) { return a + (-b); }
    friend AreaNumber operator*(const AreaNumber& a, const AreaNumber& b) {
        array<Real, 2 * Count * Count> terms{};
        int count = 0;
        for (Real x : a.parts) for (Real y : b.parts) {
            if (x == 0 || y == 0) continue;
            Real hi = x * y;
            append(terms, count, fmal(x, y, -hi));
            append(terms, count, hi);
        }
        return pack(terms, count);
    }
    friend AreaNumber operator/(AreaNumber a, const AreaNumber& b) {
        assert(b.parts[Count - 1] != 0);
        AreaNumber result;
        for (int i = 0; i < Count && a.parts[Count - 1] != 0; i++) {
            AreaNumber q(a.parts[Count - 1] / b.parts[Count - 1]);
            result = result + q;
            a = a - b * q;
        }
        return result;
    }
    friend bool operator<(const AreaNumber& a, const AreaNumber& b) { return (a - b).parts[Count - 1] < 0; }
    friend bool operator>(const AreaNumber& a, const AreaNumber& b) { return b < a; }
    friend bool operator<=(const AreaNumber& a, const AreaNumber& b) { return !(b < a); }
    friend bool operator>=(const AreaNumber& a, const AreaNumber& b) { return !(a < b); }
    friend bool operator==(const AreaNumber& a, const AreaNumber& b) { return (a - b).parts[Count - 1] == 0; }
    friend bool operator!=(const AreaNumber& a, const AreaNumber& b) { return !(a == b); }
};

Real area_sqrt(Real value) { return sqrtl(value); }
template <int Count>
AreaNumber<Count> area_sqrt(const AreaNumber<Count>& value) {
    using Q = AreaNumber<Count>;
    if (value == 0) return 0;
    Q result(sqrtl(static_cast<Real>(value)));
    for (int i = 0; i < (Count == 2 ? 2 : 3); i++) result = (result + value / result) * Q(0.5L);
    return result;
}
// 引数の絶対値を0.2以下にしてから使う。
template <int Count>
AreaNumber<Count> area_atan_series(const AreaNumber<Count>& value) {
    using Q = AreaNumber<Count>;
    Q square = value * value, power = value, result = value;
    Real tolerance = ldexpl(1.L, -Count * numeric_limits<Real>::digits - 8);
    for (int i = 1; i <= Count * numeric_limits<Real>::digits; i++) {
        power = -power * square;
        Q term = power / Q(2 * i + 1);
        result = result + term;
        if (abs(static_cast<Real>(term)) <= tolerance * abs(static_cast<Real>(result))) break;
    }
    return result;
}
Real area_atan2(Real y, Real x) { return atan2l(y, x); }
template <int Count>
AreaNumber<Count> area_atan2(AreaNumber<Count> y, AreaNumber<Count> x) {
    using Q = AreaNumber<Count>;
    if (x == 0 && y == 0) return 0;
    // long doubleの円周率を拡張せず、同じ精度で一度だけ計算する。
    static const Q pi = area_atan_series(Q(1) / Q(5)) * Q(16)
                              - area_atan_series(Q(1) / Q(239)) * Q(4);
    bool negative_x = x < 0, negative_y = y < 0;
    if (negative_x) x = -x;
    if (negative_y) y = -y;
    bool swapped = x < y;
    if (swapped) swap(x, y);
    Q z = y / x;
    // 小さい角度はそのまま使い、大きい場合だけ半角公式で収束を速める。
    int halves = 0;
    while (static_cast<Real>(z) > 0.2L) {
        z = z / (Q(1) + area_sqrt(Q(1) + z * z));
        ++halves;
    }
    Q result = area_atan_series(z) * Q(ldexpl(1.L, halves));
    if (swapped) result = pi * Q(0.5L) - result;
    if (negative_x) result = pi - result;
    return negative_y ? -result : result;
}

// 三角形・扇形の符号付き面積の2倍を足す。Real版は誤差の目安も返す。
template <class T>
pair<T, Real> circle_polygon_integral(const Polygon& polygon, const Circle& circle) {
    using P = PointBase<T>;
    T radius(circle.r), area_sum(0), correction(0);
    Real scale_sum = 0;
    auto add = [&](T value) {
        if constexpr (is_same_v<T, Real>) {
            T next = area_sum + value;
            correction += abs(area_sum) >= abs(value) ? (area_sum - next) + value : (value - next) + area_sum;
            area_sum = next;
        } else area_sum = area_sum + value;
    };
    auto sector = [&](const P& a, const P& b) {
        add(radius * radius * area_atan2(cross(a, b), dot(a, b)));
    };
    for (int i = 0; i < (int)polygon.size(); i++) {
        const Point& from = polygon[i];
        const Point& to = polygon[(i + 1) % polygon.size()];
        if (from == to) continue;
        // 座標の引き算より前に精度を広げ、入力の小さい成分を失わない。
        P a(T(from.x) - T(circle.p.x), T(from.y) - T(circle.p.y));
        P b(T(to.x) - T(circle.p.x), T(to.y) - T(circle.p.y)), d = b - a;
        if constexpr (is_same_v<T, Real>) {
            Real scale = max({abs(a.x), abs(a.y), abs(b.x), abs(b.y), circle.r});
            scale_sum += scale * scale;
        }
        T H = cross(a, b), rr = radius * radius;
        if (dot(a, a) <= rr && dot(b, b) <= rr) { add(H); continue; }
        T len = area_sqrt(dot(d, d)), h = H / len;
        if (h >= radius || h <= -radius) { sector(a, b); continue; }
        P e = d / len, normal(-e.y, e.x);
        T xa = dot(e, a), xb = dot(e, b);
        T extent = area_sqrt((radius - h) * (radius + h));
        T left = max(xa, -extent), right = min(xb, extent);
        if (right <= left) { sector(a, b); continue; }
        // 交点は端点からの比率ではなく、円の中心付近で構成する。
        // 辺が長くても、円内に残る短い区間の幅を保てる。
        add(h * (right - left));
        if (xa < -extent) sector(a, e * (-extent) - normal * h);
        if (xb > extent) sector(e * extent - normal * h, b);
    }
    return {(area_sum + correction) / T(2), scale_sum * (32 * numeric_limits<Real>::epsilon())};
}
}
