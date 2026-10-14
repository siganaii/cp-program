#pragma once
#include "library/geometry/distance.hpp"
#include "library/geometry/intersect.hpp"
Point crosspoint(const Line &l, const Line &m) {
    Real A = cross(l.b - l.a, m.b - m.a);
    Real B = cross(l.b - l.a, l.b - m.a);
    if(equals(abs(A), 0.0) && equals(abs(B), 0.0)) return m.a;
    return (m.a * (A - B) + m.b * B) / A;
}
Point crosspoint(const Segment &l, const Segment &m) {
    return crosspoint(Line(l), Line(m));
}
pair<Point,Point> crosspoint(const Circle &c, const Line &l) {
    Point pr = projection(l, c.p);
    Point e = (l.b - l.a) / abs(l.b - l.a);
    if(equals(distance(l, c.p), c.r)) return {pr, pr};
    Real base = sqrt(c.r * c.r - norm(pr - c.p));
    return {pr - e * base, pr + e * base};
}
pair<Point, Point> crosspoint(const Circle &c, const Segment &l) {
    auto [count, ts] = geometry_detail::circle_segment_parameters(c, l);
    // 交点が存在することを要求する。1点の場合は同じ点を2回返す。
    assert(count > 0);
    Point d = l.b - l.a;
    auto at = [&](Real t) {
        if (t == 0) return l.a;
        if (t == 1) return l.b;
        return l.a + d * t;
    };
    return {at(ts[0]), at(ts[count == 2 ? 1 : 0])};
}
pair<Point,Point> crosspoint(const Circle &c1,const Circle &c2) {
    Real d = abs(c1.p - c2.p);
    Real x = (c1.r * c1.r + d * d - c2.r * c2.r) / (2 * c1.r * d);
    if (abs(x) > 1) x = (x > 0 ? 1.0 : -1.0);
    Real a = acos(x);
    Real t = atan2(c2.p.imag() - c1.p.imag(), c2.p.real() - c1.p.real());
    Point p1 = c1.p + Point(cos(t + a) * c1.r, sin(t + a) * c1.r);
    Point p2 = c1.p + Point(cos(t - a) * c1.r, sin(t - a) * c1.r);
    return {p1, p2};
}
