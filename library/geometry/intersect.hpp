#pragma once
#include "library/geometry/geometry-base.hpp"
bool intersect(const Line &l, const Point &p) {
    return abs(ccw(l.a, l.b, p)) != 1;
}
bool intersect(const Line &l, const Line &m) {
    return !parallel(l, m) || geometry_detail::line_side(l.a, l.b, m.a) == 0;
}
bool intersect(const Segment &s, const Point &p) {
    return ccw(s.a, s.b, p) == 0;
}
bool intersect(const Line &l, const Segment &s) {
    return geometry_detail::line_side(l.a, l.b, s.a)
         * geometry_detail::line_side(l.a, l.b, s.b) <= 0;
}
bool intersect(const Segment &s,const Line &l) {
    return intersect(l,s);   
}
bool intersect(const Circle &c, const Point &p) {
    Real d = abs(p - c.p);
    return abs(d - c.r) <= geometry_detail::radial_tolerance(d, c.r);
}
namespace geometry_detail {
// 交点は円の中心から計算し、長い線分の始点との加減算による桁落ちを避ける。
pair<int, array<Point, 2>> circle_line_intersections(const Circle& c, const Line& l) {
    array<Point, 2> ps{};
    Point d = l.b - l.a;
    Real len = abs(d);
    if (len == 0) {
        ps[0] = l.a;
        return {intersect(c, l.a) ? 1 : 0, ps};
    }
    Point e = d / len, n(-e.y, e.x);
    Real signed_h = line_cross(l.a, l.b, c.p) / len;
    Real h = abs(signed_h), tol = radial_tolerance(h, c.r);
    if (h - c.r > tol) return {0, ps};
    Point foot(fmal(-n.x, signed_h, c.p.x), fmal(-n.y, signed_h, c.p.y));
    ps[0] = foot;
    if (abs(h - c.r) <= tol) return {1, ps};
    Real offset = sqrtl(max(Real(0), (c.r - h) * (c.r + h)));
    ps[0] = foot - e * offset;
    ps[1] = foot + e * offset;
    return {2, ps};
}
}
bool intersect(const Circle &c, const Line &l) {
    return geometry_detail::circle_line_intersections(c, l).first != 0;
}
bool intersect(const Segment &s, const Segment &t) {
    // 長さ0の線分は点として判定する。
    if (s.a == s.b) return intersect(t, s.a);
    if (t.a == t.b) return intersect(s, t.a);
    return ccw(s.a, s.b, t.a) * ccw(s.a, s.b, t.b) <= 0 && ccw(t.a, t.b, s.a) * ccw(t.a, t.b, s.b) <= 0;
}
bool intersect(const Ray &r, const Point &p) {
    Point d = r.b - r.a;
    Real len = abs(d);
    if (len == 0) return equals(r.a, p);
    d /= len;
    if (geometry_detail::line_side(r.a, r.b, p) != 0) return false;
    return dot(d, p - r.a) >= -EPS - geometry_detail::dot_error(d, p - r.a);
}

bool intersect(const Point &p, const Ray &r) {
    return intersect(r, p);
}
bool intersect(const Ray &r, const Line &l) {
    Point d = r.b - r.a;
    Point e = l.b - l.a;
    Real len_d = abs(d), len_e = abs(e);
    if (len_e == 0) return intersect(r, l.a);
    if (len_d == 0) return intersect(l, r.a);
    if (parallel(r, l)) {
        return geometry_detail::line_side(r.a, r.b, l.a) == 0;
    }
    // 方向は正規化せず、パラメータ側のEPSを長さで換算する。
    Real det = geometry_detail::difference_cross(r.a, r.b, l.a, l.b);
    Real t = geometry_detail::line_cross(l.a, l.b, r.a) / det;
    return t >= -EPS / len_d;
}
bool intersect(const Line &l, const Ray &r) {
    return intersect(r, l);
}
bool intersect(const Ray &r, const Segment &s) {
    Point d = r.b - r.a;
    Point e = s.b - s.a;
    Real len_d = abs(d), len_e = abs(e);
    if (len_d == 0) return intersect(s, r.a);
    if (len_e == 0) return intersect(r, s.a);
    if (parallel(r, s)) {
        if (geometry_detail::line_side(r.a, r.b, s.a) != 0) return false;
        return intersect(r, s.a) || intersect(r, s.b);
    }
    Real det = geometry_detail::difference_cross(r.a, r.b, s.a, s.b);
    Real t = geometry_detail::line_cross(s.a, s.b, r.a) / det;
    Real u = -geometry_detail::line_cross(r.a, r.b, s.a) / det;
    // u - 1 の引き算を避け、終点側のずれも元の座標から計算する。
    Real v = -geometry_detail::line_cross(r.a, r.b, s.b) / det;
    return t >= -EPS / len_d && u >= -EPS / len_e && v <= EPS / len_e;
}

bool intersect(const Segment &s, const Ray &r) {
    return intersect(r, s);
}
bool intersect(const Ray &r, const Ray &s) {
    Point d = r.b - r.a;
    Point e = s.b - s.a;
    Real len_d = abs(d), len_e = abs(e);
    if (len_d == 0) return intersect(s, r.a);
    if (len_e == 0) return intersect(r, s.a);
    if (parallel(r, s)) {
        if (geometry_detail::line_side(r.a, r.b, s.a) != 0) return false;
        return intersect(r, s.a) || intersect(s, r.a);
    }
    Real det = geometry_detail::difference_cross(r.a, r.b, s.a, s.b);
    Real t = geometry_detail::line_cross(s.a, s.b, r.a) / det;
    Real u = -geometry_detail::line_cross(r.a, r.b, s.a) / det;
    return t >= -EPS / len_d && u >= -EPS / len_e;
}
namespace geometry_detail {
// 始点から終点の順に交点を返す。端点を含め、動的確保は行わない。
pair<int, array<Point, 2>> circle_segment_intersections(const Circle& c,
                                                       const Segment& s) {
    auto [n, ps] = circle_line_intersections(c, s);
    if (s.a == s.b) return {n, ps};
    Point e = (s.b - s.a) / abs(s.b - s.a);
    int count = 0;
    for (int i = 0; i < n; i++) {
        Point from_a = ps[i] - s.a, from_b = ps[i] - s.b;
        Real a = dot(e, from_a), b = dot(e, from_b);
        if (a < -EPS - dot_error(e, from_a) || b > EPS + dot_error(e, from_b)) continue;
        Point p = a < 0 ? s.a : b > 0 ? s.b : ps[i];
        if (count == 0 || p != ps[count - 1]) ps[count++] = p;
    }
    return {count, ps};
}
}
// 円周との異なる交点数。接点や端点はそれぞれ1点として数える。
int intersect(const Circle &c, const Segment &l) {
    return geometry_detail::circle_segment_intersections(c, l).first;
}
int intersect(const Circle &c, const Ray &r) {
    Point d = r.b - r.a;
    Real len = abs(d);
    if (len == 0) return intersect(c, r.a) ? 1 : 0;
    d /= len;
    auto [n, ps] = geometry_detail::circle_line_intersections(c, r);
    int res = 0;
    for (int i = 0; i < n; i++) {
        Point v = ps[i] - r.a;
        if (dot(d, v) >= -EPS - geometry_detail::dot_error(d, v)) ++res;
    }
    return res;
}
int intersect(const Ray &r, const Circle &c) {
    return intersect(c, r);
}
int intersect(Circle c1,Circle c2) {
    if(c1.r < c2.r) swap(c1,c2);
    Real d = abs(c1.p - c2.p);
    Real outer = c1.r + c2.r, inner = c1.r - c2.r;
    Real outer_tol = geometry_detail::radial_tolerance(d, outer);
    Real inner_tol = geometry_detail::radial_tolerance(d, inner);
    if (abs(d - outer) <= outer_tol) return 3; //外接
    if (d > outer) return 4;                  //円が離れている
    if (abs(d - inner) <= inner_tol) return 1; //内接
    if (d > inner) return 2;                  //2点で交わる
    return 0;                       //一方の円の内部にもう一方の円がある
}
