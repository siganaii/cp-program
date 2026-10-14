#pragma once
#include "library/geometry/distance.hpp"
#include "library/geometry/intersect.hpp"
Point crosspoint(const Line &l, const Line &m) {
    Point d = l.b - l.a, e = m.b - m.a;
    if (parallel(l, m)) {
        assert(geometry_detail::line_side(l.a, l.b, m.a) == 0);
        return m.a;
    }
    Real t = geometry_detail::cross_value(m.a - l.a, e)
           / geometry_detail::cross_value(d, e);
    return {fmal(d.x, t, l.a.x), fmal(d.y, t, l.a.y)};
}
Point crosspoint(const Segment &l, const Segment &m) {
    return crosspoint(Line(l), Line(m));
}
pair<Point,Point> crosspoint(const Circle &c, const Line &l) {
    auto [count, ps] = geometry_detail::circle_line_intersections(c, l);
    assert(count > 0);
    return {ps[0], ps[count == 2 ? 1 : 0]};
}
pair<Point, Point> crosspoint(const Circle &c, const Segment &l) {
    auto [count, ps] = geometry_detail::circle_segment_intersections(c, l);
    // 交点が存在することを要求する。1点の場合は同じ点を2回返す。
    assert(count > 0);
    return {ps[0], ps[count == 2 ? 1 : 0]};
}
pair<Point,Point> crosspoint(const Circle &c1,const Circle &c2) {
    // 交点が存在することを要求する。半径0の円は1点として扱う。
    if (c1.r == 0) {
        assert(intersect(c2, c1.p));
        return {c1.p, c1.p};
    }
    if (c2.r == 0) {
        assert(intersect(c1, c2.p));
        return {c2.p, c2.p};
    }
    Real d = abs(c1.p - c2.p);
    assert(d > 0);
    // 半径の二乗同士を引かず、中心線上の距離を直接求める。
    Real a = d * 0.5L + ((c1.r - c2.r) / d) * (c1.r * 0.5L + c2.r * 0.5L);
    a = clamp(a, -c1.r, c1.r);
    Real h = sqrtl(max(Real(0), (c1.r - a) * (c1.r + a)));
    Point e = (c2.p - c1.p) / d, n(-e.y, e.x);
    Point foot = c1.p + e * a;
    int relation = intersect(c1, c2);
    if (relation == 1 || relation == 3) return {foot, foot};
    return {foot + n * h, foot - n * h};
}
