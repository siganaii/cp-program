#pragma once
template <typename R>
struct PointBase {
    using P = PointBase;
    using C = conditional_t<is_integral_v<R>, __int128_t, R>;
    R x, y;
    PointBase() : x(0), y(0) {}
    PointBase(R _x, R _y) : x(_x), y(_y) {}
    template <typename T, typename U>
    PointBase(const pair<T, U>& p) : x(p.first), y(p.second) {}
    // 通常の座標演算は結果がRの範囲に収まることを要求する。
    // 範囲を超えうる差分を幾何判定に使う場合はdifferenceを使う。
    P operator+(const P& r) const { return {x + r.x, y + r.y}; }
    P operator-(const P& r) const { return {x - r.x, y - r.y}; }
    P operator*(R r) const { return {x * r, y * r}; }
    P operator/(R r) const { return {x / r, y / r}; }
    P operator-() const { return {-x, -y}; }
    P& operator+=(const P& r) { return (*this) = (*this) + r; }
    P& operator-=(const P& r) { return (*this) = (*this) - r; }
    P& operator*=(R r) { return (*this) = (*this) * r; }
    P& operator/=(R r) { return (*this) = (*this) / r; }
    bool operator<(const P& r) const { return x != r.x ? x < r.x : y < r.y; }
    bool operator>(const P& r) const { return x != r.x ? x > r.x : y > r.y; }
    bool operator==(const P& r) const { return x == r.x and y == r.y; }
    bool operator!=(const P& r) const { return !((*this) == r); }
    P rotate(long double rad) const {
        static_assert(!is_integral<R>::value);
        return {x * cos(rad) - y * sin(rad), x * sin(rad) + y * cos(rad)};
    }
    P conj() const {return {x,-y};}
    R real() const { return x; }
    R imag() const { return y; }
    friend R real(const P& p) { return p.x; }
    friend R imag(const P& p) { return p.y; }
    friend C dot(const P& l, const P& r) {return C(l.x) * r.x + C(l.y) * r.y;}
    friend C cross(const P& l, const P& r) {return C(l.x) * r.y - C(l.y) * r.x;}
    // 幾何判定用の差分。整数では引き算の前に型を広げる。
    friend PointBase<C> difference(const P& l, const P& r) {
        return {C(l.x) - C(r.x), C(l.y) - C(r.y)};
    }
    friend long double abs(const P& p) { return hypotl((long double)p.x, (long double)p.y); }
    friend C norm(const P& p) { return dot(p, p); }
    friend long double arg(const P& p) { return atan2(p.y, p.x); }
    friend istream& operator>>(istream& is, P& p) {
        R a, b;
        is >> a >> b;
        p = P{a, b};
        return is;
    }
    friend ostream& operator<<(ostream& os, const P& p) {
        return os << p.x << " " << p.y;
    }
};

namespace geometry_detail {
// 積同士が近いときの桁落ちをFMAで補正する。整数では元の拡張型を使う。
template <typename T>
auto cross_value(const PointBase<T>& a, const PointBase<T>& b) {
    if constexpr (is_floating_point_v<T>) {
        T p = a.y * b.x;
        return fma(a.x, b.y, -p) + fma(-a.y, b.x, p);
    } else {
        return cross(a, b);
    }
}
// 方向の符号にはEPSを使わず、異なる方向を同一視しない。
template <typename T>
int cross_sign(const PointBase<T>& a, const PointBase<T>& b) {
    auto c = cross_value(a, b);
    return c > 0 ? 1 : c < 0 ? -1 : 0;
}
// 2本の差分ベクトルの外積。引き算で失われた下位桁も加える。
template <typename T>
auto difference_cross(const PointBase<T>& a, const PointBase<T>& b,
                      const PointBase<T>& c, const PointBase<T>& d) {
    if constexpr (is_floating_point_v<T>) {
        using P = PointBase<long double>;
        auto error = [](const P& x, const P& y, const P& delta) {
            P v = x - delta;
            return (x - (delta + v)) + (v - y);
        };
        P pa(a.x, a.y), pb(b.x, b.y), pc(c.x, c.y), pd(d.x, d.y);
        P u = pb - pa, v = pd - pc;
        P ue = error(pb, pa, u), ve = error(pd, pc, v);
        return cross_value(u, v) + cross_value(u, ve)
             + cross_value(ue, v) + cross_value(ue, ve);
    } else {
        return cross(difference(b, a), difference(d, c));
    }
}
template <typename T>
auto line_cross(const PointBase<T>& a, const PointBase<T>& b,
                const PointBase<T>& p) {
    return difference_cross(a, b, a, p);
}
// 頂点の向きは、差分を作るときの丸め誤差も補正して判定する。
template <typename T>
int orientation_sign(const PointBase<T>& a, const PointBase<T>& b,
                     const PointBase<T>& c) {
    auto v = line_cross(a, b, c);
    return v > 0 ? 1 : v < 0 ? -1 : 0;
}
}
