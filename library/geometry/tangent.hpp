#pragma once
#include "library/geometry/geometry-base.hpp"
#include "library/geometry/crosspoint.hpp"
pair<Point, Point> tangent(const Circle &c1, const Point &p2) {
    return crosspoint(c1, Circle(p2, sqrt(norm(c1.p - p2) - c1.r * c1.r)));
}
Lines tangent(Circle c1, Circle c2) {
    Lines ret;
    if (c1.r < c2.r) swap(c1, c2);
    Real g = norm(c1.p - c2.p);
    if (equals(g, 0)) return ret;
    Point u = (c2.p - c1.p) / sqrt(g);
    Point v = u.rotate(PI * 0.5L);
    for (int s : {-1, 1}) {
        Real h = (c1.r + s * c2.r) / sqrt(g);
        if (equals(1 - h * h, 0)) {
            ret.emplace_back(c1.p + u * c1.r, c1.p + (u + v) * c1.r);
        } 
        else if (1 - h * h > 0) {
            Point uu = u * h, vv = v * sqrt(1 - h * h);
            ret.emplace_back(c1.p + (uu + vv) * c1.r, c2.p - (uu + vv) * c2.r * s);
            ret.emplace_back(c1.p + (uu - vv) * c1.r, c2.p - (uu - vv) * c2.r * s);
        }
    }
    return ret;
}
