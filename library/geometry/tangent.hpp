#pragma once
#include "library/geometry/geometry-base.hpp"
#include "library/geometry/crosspoint.hpp"
pair<Point, Point> tangent(const Circle &c1, const Point &p2) {
    return crosspoint(c1, Circle(p2, sqrt(norm(c1.p - p2) - c1.r * c1.r)));
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
        if (r - d >= EPS) continue;
        Real h = clamp(r / d, Real(0), Real(1));
        // 二乗値ではなく、距離と半径の差を比較する。
        if (abs(d - r) < EPS) {
            Point p = c1.p + u * c1.r;
            ret.emplace_back(p, p + v);
        } 
        else {
            Point uu = u * h, vv = v * sqrtl(max(Real(0), (1 - h) * (1 + h)));
            ret.emplace_back(c1.p + (uu + vv) * c1.r, c2.p - (uu + vv) * c2.r * s);
            ret.emplace_back(c1.p + (uu - vv) * c1.r, c2.p - (uu - vv) * c2.r * s);
        }
    }
    return ret;
}
