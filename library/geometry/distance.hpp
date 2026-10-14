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
    Point r = projection(s, p);
    if (intersect(s, r)) return abs(r - p);
    return min(abs(s.a - p), abs(s.b - p));
}
Real distance(const Segment &a, const Segment &b) {
    if (intersect(a, b)) return 0;
    return min({distance(a, b.a), distance(a, b.b), distance(b, a.a), distance(b, a.b)});
}
Real distance(const Line &l, const Segment &s) {
    if (intersect(l, s)) return 0;
    return min(distance(l, s.a), distance(l, s.b));
}