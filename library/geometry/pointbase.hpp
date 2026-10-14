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
    friend long double abs(const P& p) { return sqrtl(p.x * p.x + p.y * p.y); }
    friend R norm(const P& p) { return p.x * p.x + p.y * p.y; }
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