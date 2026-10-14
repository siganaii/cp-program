#pragma once
#include "library/geometry/geometry-base.hpp"
bool intersect(const Line &l, const Point &p) {
    return abs(ccw(l.a, l.b, p)) != 1;
}
bool intersect(const Line &l, const Line &m) {
    return abs(cross(l.b - l.a, m.b - m.a)) > EPS || abs(cross(l.b - l.a, m.b - l.a)) < EPS;
}
bool intersect(const Segment &s, const Point &p) {
    return ccw(s.a, s.b, p) == 0;
}
bool intersect(const Line &l, const Segment &s) {
    return cross(l.b - l.a, s.a - l.a) * cross(l.b - l.a, s.b - l.a) < EPS;
}
bool intersect(const Segment &s,const Line &l) {
    return intersect(l,s);   
}
bool intersect(const Circle &c, const Line &l) {
    return abs(c.p - projection(l, c.p)) <= c.r + EPS;
}
bool intersect(const Circle &c, const Point &p) {
    return abs(abs(p - c.p) - c.r) < EPS;
}
bool intersect(const Segment &s, const Segment &t) {
    return ccw(s.a, s.b, t.a) * ccw(s.a, s.b, t.b) <= 0 && ccw(t.a, t.b, s.a) * ccw(t.a, t.b, s.b) <= 0;
}
bool intersect(const Ray &r, const Point &p) {
    if (abs(cross(r.b - r.a, p - r.a)) > EPS) return false;
    return dot(r.b - r.a, p - r.a) >= -EPS;
}

bool intersect(const Point &p, const Ray &r) {
    return intersect(r, p);
}
bool intersect(const Ray &r, const Line &l) {
    Point d = r.b - r.a;
    Point e = l.b - l.a;
    Real det = cross(d, e);
    if (abs(det) < EPS) {
        return abs(cross(d, l.a - r.a)) < EPS;
    }
    Real t = cross(l.a - r.a, e) / det;
    return t >= -EPS;
}
bool intersect(const Line &l, const Ray &r) {
    return intersect(r, l);
}
bool intersect(const Ray &r, const Segment &s) {
    Point d = r.b - r.a;
    Point e = s.b - s.a;
    Real det = cross(d, e);
    if (abs(det) < EPS) {
        if (abs(cross(d, s.a - r.a)) > EPS) return false;
        return dot(d, s.a - r.a) >= -EPS
            || dot(d, s.b - r.a) >= -EPS;
    }
    Real t = cross(s.a - r.a, e) / det;
    Real u = cross(s.a - r.a, d) / det;
    return t >= -EPS && u >= -EPS && u <= 1 + EPS;
}

bool intersect(const Segment &s, const Ray &r) {
    return intersect(r, s);
}
bool intersect(const Ray &r, const Ray &s) {
    Point d = r.b - r.a;
    Point e = s.b - s.a;
    Real det = cross(d, e);
    if (abs(det) < EPS) {
        if (abs(cross(d, s.a - r.a)) > EPS) return false;
        return intersect(r, s.a) || intersect(s, r.a);
    }
    Real t = cross(s.a - r.a, e) / det;
    Real u = cross(s.a - r.a, d) / det;
    return t >= -EPS && u >= -EPS;
}
int intersect(const Circle &c, const Segment &l) {
    if (norm(projection(l, c.p) - c.p) - c.r * c.r > EPS) return 0;
    auto d1 = abs(c.p - l.a), d2 = abs(c.p - l.b);
    if (d1 < c.r + EPS && d2 < c.r + EPS) return 0;
    if (d1 < c.r - EPS && d2 > c.r + EPS || d1 > c.r + EPS && d2 < c.r - EPS)
        return 1;
    const Point h = projection(l, c.p);
    if (dot(l.a - h, l.b - h) < 0) return 2;
    return 0;
}
int intersect(const Circle &c, const Ray &r) {
    Point d = r.b - r.a;
    Real t = dot(c.p - r.a, d) / norm(d);
    Point h = r.a + d * t;
    Real h2 = c.r * c.r - norm(h - c.p);
    if (h2 < -EPS) return 0;
    if (abs(h2) < EPS) {
        return t >= -EPS ? 1 : 0;
    }
    Real dt = sqrtl(max((Real)0, h2) / norm(d));
    Real t1 = t - dt;
    Real t2 = t + dt;
    int res = 0;
    if (t1 >= -EPS) res++;
    if (t2 >= -EPS) res++;
    return res;
}
int intersect(const Ray &r, const Circle &c) {
    return intersect(c, r);
}
int intersect(Circle c1,Circle c2) {
    if(c1.r < c2.r) swap(c1,c2);
    Real d = abs(c1.p - c2.p);
    if(c1.r + c2.r < d) return 4;   //円が離れている
    if(equals(c1.r + c2.r,d)) return 3; //円の外部と外部が接している
    if(c1.r - c2.r < d) return 2;   //円が交わっている
    if(equals(c1.r - c2.r,d)) return 1; //一方の円の内部ともう一方の円の外部と接している
    return 0;                       //一方の円の内部にもう一方の円がある
}