#pragma once
#include "library/geometry/pointbase.hpp"
using Real = long double;
using Point = PointBase<Real>;
long double EPS = 1e-11;
inline int sign(const Real &r) { return r <= -EPS ? -1 : r >= EPS ? 1 : 0; }
template <class T, class U>
bool equals(T a, U b) {
    if constexpr (is_same_v<T, long double> || is_same_v<U, long double>) {
        return fabsl((long double)a - (long double)b) < EPS;
    } else {
        return a == b;
    }
}
template <class T, class U>
bool equals(const PointBase<T>& a, const PointBase<U>& b) {
    return equals(a.x, b.x) && equals(a.y, b.y);
}
Real radian_to_degree(Real r) { return (r * 180.0 / PI); }
Real degree_to_radian(Real d) { return (d * PI / 180.0); }
int ccw(const Point &a,Point b,Point c) {
    b = b - a, c = c - a;
    if(cross(b,c) > EPS) return +1;  //反時計回り
    if(cross(b,c) < -EPS) return -1; //時計回り
    if(dot(b,c) < -EPS) return +2;      //b-a-c
    if(norm(b) < norm(c)) return -2; //a-b-c
    return 0;                        //a-c-b
}
// a-bベクトルとb-cベクトルのなす角
long double get_angle(const Point &a,const Point &b,const Point &c) {
    const Point v(b - a), w(c - b);
    long double alpha = arg(v), beta = arg(w);
    if(alpha > beta) swap(alpha,beta);
    long double theta = beta - alpha;
    return min(theta,2 * acosl(-1) - theta);
}
struct Circle {
    Point p;
    Real r;
    Circle() = default;
    Circle(Point _p,Real _r):p(_p),r(_r){}
};
struct Line {
    Point a, b;
    Line() = default;
    Line(Point a, Point b) : a(a), b(b) {}
    Line(Real A, Real B, Real C) // Ax + By = C
    {
        if(equals(A,0)) a = Point(0, C / B), b = Point(1, C / B);
        else if(equals(B,0)) a = Point(C / A, 0), b = Point(C / A, 1);
        else if(equals(C,0)) a = Point(0, C / B), b = Point(1, (C - A) / B);
        else a = Point(0, C / B), b = Point(C / A, 0);
    }
    friend ostream &operator<<(ostream &os, Line &p) {
        return os << p.a << " to " << p.b;
    }
    friend istream &operator>>(istream &is, Line &a) {
        return is >> a.a >> a.b;
    }
};
struct Segment : Line {
    Segment() = default;
    Segment(Point a, Point b) : Line(a, b) {}
};
//aを始点としてbに向けて伸びる半直線
struct Ray : Line {
    Ray() = default;
    Ray(Point a, Point b) : Line(a, b) {}
};
bool parallel(const Line &a, const Line &b) {
    return equals(cross(a.b - a.a, b.b - b.a),0);
}
bool orthogonal(const Line &a, const Line &b) {
    return equals(dot(a.a - a.b, b.a - b.b),0);
}
Point projection(const Line &l, const Point &p) {
    Real t = dot(p - l.a, l.a - l.b) / norm(l.a - l.b);
    return l.a + (l.a - l.b) * t;
}
Point reflection(const Line &l, const Point &p) {
    return p + (projection(l, p) - p) * 2.0;
}
using Points = vector<Point>;
using Polygon = vector<Point>;
using Lines = vector<Line>;
using Segments = vector<Segment>;
using Rays = vector<Ray>;
using Circles = vector<Circle>;
