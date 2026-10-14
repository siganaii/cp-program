#pragma once
#include "library/geometry/distance.hpp"
#include "library/geometry/intersect.hpp"
Point crosspoint(const Line &l, const Line &m) {
    Point d = l.b - l.a, e = m.b - m.a;
    if (parallel(l, m)) {
        assert(geometry_detail::line_side(l.a, l.b, m.a) == 0);
        return m.a;
    }
    Real det = geometry_detail::difference_cross(l.a, l.b, m.a, m.b);
    Real t = geometry_detail::line_cross(m.a, m.b, l.a) / det;
    Point p(fmal(d.x, t, l.a.x), fmal(d.y, t, l.a.y));
    // 遠い端点から求めた座標を、両直線からのずれに基づいて局所的に補正する。
    Real along_d = geometry_detail::line_cross(m.a, m.b, p) / det;
    Real along_e = -geometry_detail::line_cross(l.a, l.b, p) / det;
    return {fmal(d.x, along_d, fmal(e.x, along_e, p.x)),
            fmal(d.y, along_d, fmal(e.y, along_e, p.y))};
}
Point crosspoint(const Segment &l, const Segment &m) {
    // 交点の存在を要求する。重なる場合は両方の線分に属する端点を返す。
    if (parallel(l, m)) {
        if (intersect(m, l.a)) return l.a;
        if (intersect(m, l.b)) return l.b;
        if (intersect(l, m.a)) return m.a;
        if (intersect(l, m.b)) return m.b;
        assert(false);
    }
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
    // 小さい円を基準にし、大きい半径と中心線上の距離の差が潰れるのを避ける。
    bool swapped = c1.r > c2.r;
    const Circle& small = swapped ? c2 : c1;
    const Circle& large = swapped ? c1 : c2;
    Real a;
    if (small.r < large.r * 0.5L) {
        // 半径差が大きい場合は d^2 - large.r^2 を先に因数分解する。
        a = (d - large.r) * (0.5L + 0.5L * (large.r / d))
          + small.r * (0.5L * (small.r / d));
    } else {
        // 半径が近い場合は半径同士の差を使い、近接した中心にも対応する。
        a = d * 0.5L + ((small.r - large.r) / d) * (small.r * 0.5L + large.r * 0.5L);
    }
    a = clamp(a, -small.r, small.r);
    Real h = sqrtl(max(Real(0), (small.r - a) * (small.r + a)));
    Point e = (large.p - small.p) / d;
    // 返す2点の順序は、元の c1 -> c2 の向きにそろえる。
    Point n = swapped ? Point(e.y, -e.x) : Point(-e.y, e.x);
    Point foot(fmal(e.x, a, small.p.x), fmal(e.y, a, small.p.y));
    int relation = intersect(c1, c2);
    if (relation == 1 || relation == 3) return {foot, foot};
    return {foot + n * h, foot - n * h};
}
