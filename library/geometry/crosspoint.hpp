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
    Line aa = Line(l.a, l.b);
    if (intersect(c, l) == 2) return crosspoint(c, aa);
    auto ret = crosspoint(c, aa);
    if (dot(l.a - ret.first, l.b - ret.first) < 0)
        ret.second = ret.first;
    else
        ret.first = ret.second;
    return ret;
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