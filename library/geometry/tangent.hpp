#pragma once
#include "library/geometry/geometry-base.hpp"
#include "library/geometry/crosspoint.hpp"
pair<Point, Point> tangent(const Circle &c1, const Point &p2) {
    if (c1.r == 0) return {c1.p, c1.p};
    Point v = p2 - c1.p;
    Real d = abs(v);
    // 円の外側または円周上の点を要求する。
    assert(d > 0 && d >= c1.r - geometry_detail::radial_tolerance(d, c1.r));
    Point e = v / d, n(-e.y, e.x);
    Real ratio = min(Real(1), c1.r / d);
    // 補助円との交点に変換せず、接点の中心線方向・垂直方向の成分を求める。
    Point foot = c1.p + e * (c1.r * ratio);
    Real h = c1.r * sqrtl(max(Real(0), ((d - c1.r) / d) * (1 + ratio)));
    return {foot + n * h, foot - n * h};
}
Lines tangent(Circle c1, Circle c2) {
    Lines ret;
    if (c1.r < c2.r) swap(c1, c2);
    Real d = abs(c1.p - c2.p);
    if (d == 0) return ret;
    if (c1.r == 0) return {Line(c1.p, c2.p)};
    Point u = (c2.p - c1.p) / d;
    Point v(-u.y, u.x);
    for (int s : {-1, 1}) {
        // 半径0の円では内公接線と外公接線が一致する。
        if (s == 1 && c2.r == 0) continue;
        Real r = c1.r + s * c2.r;
        Real tol = geometry_detail::radial_tolerance(r, d);
        if (r - d > tol) continue;
        Real h = clamp(r / d, Real(0), Real(1));
        // 二乗値ではなく、距離と半径の差を比較する。
        if (abs(d - r) <= tol) {
            Point p = c1.p + u * c1.r;
            // 大きい座標でも方向ベクトルが丸めで消えない長さを使う。
            Real step = max({Real(1), abs(p.x), abs(p.y)});
            Point q(fmal(v.x, step, p.x), fmal(v.y, step, p.y));
            ret.emplace_back(p, q);
        } 
        else {
            Point uu = u * h, vv = v * sqrtl(max(Real(0), (1 - h) * (1 + h)));
            ret.emplace_back(c1.p + (uu + vv) * c1.r, c2.p - (uu + vv) * c2.r * s);
            ret.emplace_back(c1.p + (uu - vv) * c1.r, c2.p - (uu - vv) * c2.r * s);
        }
    }
    return ret;
}
