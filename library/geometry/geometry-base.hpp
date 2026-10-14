#pragma once
#include "library/geometry/pointbase.hpp"
using Real = long double;
using Point = PointBase<Real>;
long double EPS = 1e-11;
inline int sign(const Real &r) { return r <= -EPS ? -1 : r >= EPS ? 1 : 0; }
template <class T, class U>
bool equals(T a, U b) {
    if constexpr (is_same_v<T, long double> || is_same_v<U, long double>) {
        return fabsl((long double)a - (long double)b) < EPS;
    } else {
        return a == b;
    }
}
template <class T, class U>
bool equals(const PointBase<T>& a, const PointBase<U>& b) {
    return equals(a.x, b.x) && equals(a.y, b.y);
}
namespace geometry_detail {
// 外積・内積の符号には、距離のEPSではなく演算の丸め誤差を使う。
template <typename T>
int cross_sign(const PointBase<T>& a, const PointBase<T>& b) {
    auto c = cross(a, b);
    if constexpr (is_floating_point_v<T>) {
        auto tol = 16 * numeric_limits<T>::epsilon()
            * (abs(a.x * b.y) + abs(a.y * b.x));
        return c > tol ? 1 : c < -tol ? -1 : 0;
    } else {
        return c > 0 ? 1 : c < 0 ? -1 : 0;
    }
}
Real dot_error(const Point& a, const Point& b) {
    return 16 * numeric_limits<Real>::epsilon()
        * (abs(a.x * b.x) + abs(a.y * b.y));
}
// 積同士が近いときの桁落ちをFMAで補正する。
Real cross_value(const Point& a, const Point& b) {
    Real p = a.y * b.x;
    return fmal(a.x, b.y, -p) + fmal(-a.y, b.x, p);
}
// 2本の差分ベクトルの外積。引き算で失われた下位桁も加える。
Real difference_cross(const Point& a, const Point& b,
                      const Point& c, const Point& d) {
    auto error = [](const Point& x, const Point& y, const Point& d) {
        Point v = x - d;
        return (x - (d + v)) + (v - y);
    };
    Point u = b - a, v = d - c;
    Point ue = error(b, a, u), ve = error(d, c, v);
    return cross_value(u, v) + cross_value(u, ve)
         + cross_value(ue, v) + cross_value(ue, ve);
}
Real line_cross(const Point& a, const Point& b, const Point& p) {
    return difference_cross(a, b, a, p);
}
Real radial_tolerance(Real a, Real b) {
    return EPS + 16 * numeric_limits<Real>::epsilon() * max(abs(a), abs(b));
}
// 点と直線の位置関係。外積を方向ベクトルの長さで割り、距離で比較する。
int line_side(const Point& a, const Point& b, const Point& p) {
    Point d = b - a, v = p - a;
    Real len = abs(d);
    if (len == 0) return 0;
    Real h = cross_value(d, v) / len;
    Real tol = EPS + 16 * numeric_limits<Real>::epsilon()
        * (abs(d.x * v.y) + abs(d.y * v.x)) / len;
    return h > tol ? 1 : h < -tol ? -1 : 0;
}
}
Real radian_to_degree(Real r) { return (r * 180.0 / PI); }
Real degree_to_radian(Real d) { return (d * PI / 180.0); }
int ccw(const Point &a,Point b,Point c) {
    Point d = b - a;
    Real len = abs(d);
    if (len == 0) return abs(c - a) <= EPS ? 0 : -2;
    int side = geometry_detail::line_side(a, b, c);
    if (side > 0) return +1; // 反時計回り
    if (side < 0) return -1; // 時計回り
    d /= len;
    Point from_a = c - a, from_b = c - b;
    if (dot(d, from_a) < -EPS - geometry_detail::dot_error(d, from_a)) return +2; //b-a-c
    if (dot(d, from_b) > EPS + geometry_detail::dot_error(d, from_b)) return -2; //a-b-c
    return 0;                        //a-c-b
}
// a-bベクトルとb-cベクトルのなす角
long double get_angle(const Point &a,const Point &b,const Point &c) {
    const Point v(b - a), w(c - b);
    long double alpha = arg(v), beta = arg(w);
    if(alpha > beta) swap(alpha,beta);
    long double theta = beta - alpha;
    return min(theta,2 * acosl(-1) - theta);
}
struct Circle {
    Point p;
    Real r;
    Circle() = default;
    Circle(Point _p,Real _r):p(_p),r(_r){}
};
struct Line {
    Point a, b;
    Line() = default;
    Line(Point a, Point b) : a(a), b(b) {}
    Line(Real A, Real B, Real C) // Ax + By = C
    {
        // 係数の一律な拡大・縮小で直線を変えないため、ゼロは厳密に判定する。
        assert(A != 0 || B != 0);
        if(A == 0) a = Point(0, C / B), b = Point(1, C / B);
        else if(B == 0) a = Point(C / A, 0), b = Point(C / A, 1);
        else if(C == 0) {
            Real scale = max(abs(A), abs(B));
            a = Point(0, 0), b = Point(B / scale, -A / scale);
        }
        else a = Point(0, C / B), b = Point(C / A, 0);
    }
    friend ostream &operator<<(ostream &os, Line &p) {
        return os << p.a << " to " << p.b;
    }
    friend istream &operator>>(istream &is, Line &a) {
        return is >> a.a >> a.b;
    }
};
struct Segment : Line {
    Segment() = default;
    Segment(Point a, Point b) : Line(a, b) {}
};
//aを始点としてbに向けて伸びる半直線
struct Ray : Line {
    Ray() = default;
    Ray(Point a, Point b) : Line(a, b) {}
};
bool parallel(const Line &a, const Line &b) {
    return geometry_detail::cross_sign(a.b - a.a, b.b - b.a) == 0;
}
bool orthogonal(const Line &a, const Line &b) {
    Point d = a.b - a.a, e = b.b - b.a;
    return abs(dot(d, e)) <= geometry_detail::dot_error(d, e);
}
Point projection(const Line &l, const Point &p) {
    Point d = l.b - l.a;
    Real len = abs(d);
    assert(len > 0);
    Point n(-d.y / len, d.x / len);
    Real h = geometry_detail::line_cross(l.a, l.b, p) / len;
    // 遠い端点に戻すのではなく、対象の点から垂直方向に移動する。
    return {fmal(-n.x, h, p.x), fmal(-n.y, h, p.y)};
}
Point reflection(const Line &l, const Point &p) {
    return p + (projection(l, p) - p) * 2.0;
}
using Points = vector<Point>;
using Polygon = vector<Point>;
using Lines = vector<Line>;
using Segments = vector<Segment>;
using Rays = vector<Ray>;
using Circles = vector<Circle>;
