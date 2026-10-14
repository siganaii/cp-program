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
    return sign(cross(l.b - l.a, s.a - l.a))
         * sign(cross(l.b - l.a, s.b - l.a)) <= 0;
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
namespace geometry_detail {
// 方向の平行判定には長さのEPSを使わず、外積の丸め誤差を見込む。
bool parallel_directions(const Point& d, const Point& e) {
    Real tol = 16 * numeric_limits<Real>::epsilon()
        * (abs(d.x * e.y) + abs(d.y * e.x));
    return abs(cross(d, e)) <= tol;
}
}
bool intersect(const Ray &r, const Point &p) {
    Point d = r.b - r.a;
    Real len = abs(d);
    if (len == 0) return equals(r.a, p);
    d /= len;
    if (abs(cross(d, p - r.a)) > EPS) return false;
    return dot(d, p - r.a) >= -EPS;
}

bool intersect(const Point &p, const Ray &r) {
    return intersect(r, p);
}
bool intersect(const Ray &r, const Line &l) {
    Point d = r.b - r.a;
    Point e = l.b - l.a;
    Real len_d = abs(d), len_e = abs(e);
    if (len_e == 0) return intersect(r, l.a);
    if (len_d == 0) return intersect(l, r.a);
    // tを始点からの距離にし、方向ベクトルの長さによらずEPSを適用する。
    d /= len_d;
    e /= len_e;
    Real det = cross(d, e);
    if (geometry_detail::parallel_directions(d, e)) {
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
    Real len_d = abs(d), len_e = abs(e);
    if (len_d == 0) return intersect(s, r.a);
    if (len_e == 0) return intersect(r, s.a);
    d /= len_d;
    e /= len_e;
    Real det = cross(d, e);
    if (geometry_detail::parallel_directions(d, e)) {
        if (abs(cross(d, s.a - r.a)) > EPS) return false;
        return dot(d, s.a - r.a) >= -EPS
            || dot(d, s.b - r.a) >= -EPS;
    }
    Real t = cross(s.a - r.a, e) / det;
    Real u = cross(s.a - r.a, d) / det;
    return t >= -EPS && u >= -EPS && u <= len_e + EPS;
}

bool intersect(const Segment &s, const Ray &r) {
    return intersect(r, s);
}
bool intersect(const Ray &r, const Ray &s) {
    Point d = r.b - r.a;
    Point e = s.b - s.a;
    Real len_d = abs(d), len_e = abs(e);
    if (len_d == 0) return intersect(s, r.a);
    if (len_e == 0) return intersect(r, s.a);
    d /= len_d;
    e /= len_e;
    Real det = cross(d, e);
    if (geometry_detail::parallel_directions(d, e)) {
        if (abs(cross(d, s.a - r.a)) > EPS) return false;
        return intersect(r, s.a) || intersect(s, r.a);
    }
    Real t = cross(s.a - r.a, e) / det;
    Real u = cross(s.a - r.a, d) / det;
    return t >= -EPS && u >= -EPS;
}
namespace geometry_detail {
// 円周と線分の交点を a + (b-a)*t の t で返す。端点を含め、昇順に並ぶ。
// 固定長配列を使い、交点数の判定でも動的確保を行わない。
pair<int, array<Real, 2>> circle_segment_parameters(const Circle& c,
                                                   const Segment& s) {
    array<Real, 2> ts{};
    int count = 0;
    Point d = s.b - s.a;
    Real len = abs(d);
    if (len == 0) return {intersect(c, s.a) ? 1 : 0, ts};
    Point e = d / len, v = c.p - s.a;
    Real h = abs(cross(e, v)), along = dot(e, v);
    const Real roundoff = 16 * numeric_limits<Real>::epsilon();
    Real tol = EPS + roundoff * max(c.r, abs(v));
    if (h > c.r + tol) return {0, ts};
    Real dt = abs(h - c.r) <= tol ? 0
        : sqrtl(max(Real(0), (c.r - h) * (c.r + h)));
    Real t_tol = EPS / len + roundoff;
    auto add = [&](Real t) {
        if (t < -t_tol || t > 1 + t_tol) return;
        t = clamp(t, Real(0), Real(1));
        if (count == 0 || t != ts[count - 1]) ts[count++] = t;
    };
    add((along - dt) / len);
    if (dt != 0) add((along + dt) / len);
    return {count, ts};
}
}
// 円周との異なる交点数。接点や端点はそれぞれ1点として数える。
int intersect(const Circle &c, const Segment &l) {
    return geometry_detail::circle_segment_parameters(c, l).first;
}
int intersect(const Circle &c, const Ray &r) {
    Point d = r.b - r.a;
    Real len = abs(d);
    if (len == 0) return intersect(c, r.a) ? 1 : 0;
    d /= len;
    Point v = c.p - r.a;
    Real t = dot(v, d), h = abs(cross(d, v));
    Real tol = EPS + 16 * numeric_limits<Real>::epsilon() * max(c.r, abs(v));
    if (h > c.r + tol) return 0;
    if (abs(h - c.r) <= tol) {
        return t >= -EPS ? 1 : 0;
    }
    Real dt = sqrtl(max(Real(0), (c.r - h) * (c.r + h)));
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
