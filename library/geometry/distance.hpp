#pragma once
#include "library/geometry/geometry-base.hpp"
#include "library/geometry/intersect.hpp"
Real distance(const Point &a, const Point &b) { return abs(a - b); }
Real distance(const Line &l, const Point &p) {
    return abs(p - projection(l, p));
}
Real distance(const Line &l, const Line &m) {
    return intersect(l, m) ? 0 : distance(l, m.a);
}
Real distance(const Segment &s, const Point &p) {
    Point d = s.b - s.a;
    Real len2 = norm(d);
    if (len2 == 0) return abs(s.a - p);
    // 垂線の足の位置は内積で判定し、外積の丸め誤差を避ける。
    Real t = clamp(dot(p - s.a, d) / len2, Real(0), Real(1));
    return abs((p - s.a) - d * t);
}
Real distance(const Segment &a, const Segment &b) {
    if (intersect(a, b)) return 0;
    return min({distance(a, b.a), distance(a, b.b), distance(b, a.a), distance(b, a.b)});
}
Real distance(const Line &l, const Segment &s) {
    if (intersect(l, s)) return 0;
    return min(distance(l, s.a), distance(l, s.b));
}
