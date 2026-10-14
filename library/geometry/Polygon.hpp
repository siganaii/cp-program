#pragma once
#include "library/geometry/geometry-base.hpp"
#include "library/geometry/crosspoint.hpp"
#include "library/geometry/intersect.hpp"

// 単純多角形を前提とする。時計回り・反時計回りの両方に対応。
template <typename T>
bool is_convex(const vector<PointBase<T>>& p) {
    int n = p.size(), dir = 0;
    for (int i = 0; i < n; i++) {
        auto c = cross(difference(p[(i + 1) % n], p[i]),
                       difference(p[(i + 2) % n], p[(i + 1) % n]));
        if (equals(c, 0)) continue;
        int d = c > 0 ? 1 : -1;
        if (dir != 0 && dir != d) return false;
        dir = d;
    }
    return true;
}

// 入力点の添字を返す。同一座標は1点にまとめる。
// full: 辞書順最小の点から反時計回り。
// lower/upper: 辞書順の左から右。
// 全点共線なら、allow_180=true は全点、false は両端点を返す。
// sorted=true は XY が辞書順にソート済みであることを要求する。
template <typename T, bool allow_180 = false>
vector<int> convex_hull(const vector<PointBase<T>>& XY,
                        string mode = "full", bool sorted = false) {
    assert(mode == "full" || mode == "lower" || mode == "upper");
    vector<int> I(XY.size());
    iota(I.begin(), I.end(), 0);
    if (!sorted) {
        sort(I.begin(), I.end(), [&](int i, int j) {
            if (XY[i] != XY[j]) return XY[i] < XY[j];
            return i < j;
        });
    }
    I.erase(unique(I.begin(), I.end(), [&](int i, int j) {
        return XY[i] == XY[j];
    }), I.end());

    int n = I.size();
    if (n <= 2) return I;
    bool collinear = true;
    for (int k = 1; k + 1 < n; k++) {
        if (cross(difference(XY[I.back()], XY[I[0]]),
                  difference(XY[I[k]], XY[I[0]])) != 0) {
            collinear = false;
            break;
        }
    }
    if (collinear) {
        if constexpr (allow_180) return I;
        return {I.front(), I.back()};
    }

    auto calc = [&](bool upper) {
        vector<int> H;
        for (int t = 0; t < n; t++) {
            int k = I[upper ? n - 1 - t : t];
            while (H.size() >= 2) {
                auto d = cross(difference(XY[H.back()], XY[H[H.size() - 2]]),
                               difference(XY[k], XY[H.back()]));
                if constexpr (allow_180) {
                    if (d >= 0) break;
                } else {
                    if (d > 0) break;
                }
                H.pop_back();
            }
            H.push_back(k);
        }
        return H;
    };
    if (mode == "lower") return calc(false);
    if (mode == "upper") {
        auto H = calc(true);
        reverse(H.begin(), H.end());
        return H;
    }
    auto L = calc(false), U = calc(true);
    L.pop_back();
    U.pop_back();
    L.insert(L.end(), U.begin(), U.end());
    return L;
}

template <typename T, bool allow_180 = false>
vector<int> convex_hull(const vector<pair<T, T>>& XY,
                        string mode = "full", bool sorted = false) {
    vector<PointBase<T>> points;
    points.reserve(XY.size());
    for (auto [x, y] : XY) points.emplace_back(x, y);
    return convex_hull<T, allow_180>(points, mode, sorted);
}

enum { OUT, ON, IN };

// 単純多角形の内外判定。O(N)。
int contains(const Polygon& Q, const Point& p) {
    bool in = false;
    for (int i = 0; i < (int)Q.size(); i++) {
        Point a = Q[i] - p, b = Q[(i + 1) % Q.size()] - p;
        if (a.imag() > b.imag()) swap(a, b);
        if (equals(cross(a, b), 0) && dot(a, b) <= 0) return ON;
        // 境界判定以外では高さの基準を0に統一し、頂点を二重に数えない。
        if (a.imag() <= 0 && 0 < b.imag() && cross(a, b) < 0)
            in = !in;
    }
    return in ? IN : OUT;
}

// 凸多角形の内外判定。時計回り・反時計回りの両方に対応。
// 余分な同一直線上の頂点がなければ O(log N)。
// 退化ケースや一部の境界ケースでは contains による O(N) 判定。
int convex_contains(const Polygon& Q, const Point& p) {
    int n = Q.size();
    if (n == 0) return OUT;
    if (n == 1) return equals(Q[0], p) ? ON : OUT;
    auto on_segment = [&](const Point& a, const Point& b) {
        return equals(cross(b - a, p - a), 0)
            && dot(p - a, p - b) <= 0;
    };
    if (n == 2) return on_segment(Q[0], Q[1]) ? ON : OUT;
    auto sgn = [](Real x) {
        return equals(x, 0) ? 0 : (x > 0 ? 1 : -1);
    };
    int dir = sgn(cross(Q[1] - Q[0], Q.back() - Q[0]));
    if (dir == 0) return contains(Q, p);

    auto side = [&](const Point& a, const Point& b) {
        return dir * sgn(cross(a - Q[0], b - Q[0]));
    };
    int a = side(Q[1], p), b = side(Q.back(), p);
    if (a < 0 || b > 0) return OUT;
    if (a == 0) {
        if (equals(cross(Q[1] - Q[0], Q[2] - Q[0]), 0))
            return contains(Q, p);
        return on_segment(Q[0], Q[1]) ? ON : OUT;
    }
    if (b == 0) {
        if (equals(cross(Q.back() - Q[0], Q[n - 2] - Q[0]), 0))
            return contains(Q, p);
        return on_segment(Q[0], Q.back()) ? ON : OUT;
    }

    int l = 1, r = n - 1;
    while (r - l > 1) {
        int m = (l + r) / 2;
        if (side(Q[m], p) >= 0) l = m;
        else r = m;
    }
    int s = dir * sgn(cross(Q[r] - Q[l], p - Q[l]));
    return s < 0 ? OUT : s == 0 ? ON : IN;
}

// 同一直線上で重なる・接する線分を併合する。
// O(N^2 + N log N)。長さ0の線分にも対応。
void merge_segments(vector<Segment>& segs) {
    vector<Segment> S, result;
    vector<Point> points;
    for (auto s : segs) {
        if (s.b < s.a) swap(s.a, s.b);
        if (s.a == s.b) points.push_back(s.a);
        else S.push_back(s);
    }

    vector<bool> used(S.size());
    for (int i = 0; i < (int)S.size(); i++) {
        if (used[i]) continue;
        Point d = S[i].b - S[i].a;
        vector<Segment> group;
        for (int j = i; j < (int)S.size(); j++) {
            if (used[j]) continue;
            if (!equals(cross(d, S[j].b - S[j].a), 0)) continue;
            if (!equals(cross(d, S[j].a - S[i].a), 0)) continue;
            used[j] = true;
            group.push_back(S[j]);
        }
        sort(group.begin(), group.end(), [](const Segment& a,
                                            const Segment& b) {
            return a.a != b.a ? a.a < b.a : a.b < b.b;
        });
        Segment cur = group[0];
        for (int j = 1; j < (int)group.size(); j++) {
            if (!(cur.b < group[j].a) || equals(cur.b, group[j].a)) {
                cur.b = max(cur.b, group[j].b);
            } else {
                result.push_back(cur);
                cur = group[j];
            }
        }
        result.push_back(cur);
    }

    sort(points.begin(), points.end());
    points.erase(unique(points.begin(), points.end()), points.end());
    int k = result.size();
    for (const auto& p : points) {
        bool covered = false;
        for (int i = 0; i < k; i++) {
            if (intersect(result[i], p)) {
                covered = true;
                break;
            }
        }
        if (!covered) result.emplace_back(p, p);
    }
    segs = move(result);
}

// 線分の端点・交点を頂点とするグラフを構築する。
// ps に既存の点がある場合、それらも頂点候補に含める。
// 浮動小数点の同一交点の統合は厳密比較。
// 重なり合う線分を渡すと多重辺が生じる場合がある。
vector<vector<int>> segment_arrangement(vector<Segment>& segs,
                                         vector<Point>& ps) {
    int N = segs.size();
    vector<vector<Point>> on(N);
    for (int i = 0; i < N; i++) {
        ps.emplace_back(segs[i].a);
        ps.emplace_back(segs[i].b);
        on[i].push_back(segs[i].a);
        on[i].push_back(segs[i].b);
        for (int j = i + 1; j < N; j++) {
            Point p1 = segs[i].b - segs[i].a;
            Point p2 = segs[j].b - segs[j].a;
            if (cross(p1, p2) == 0) continue;
            if (intersect(segs[i], segs[j])) {
                Point p = crosspoint(segs[i], segs[j]);
                ps.push_back(p);
                // 算出元の線分への所属は、丸め誤差のある座標で再判定しない。
                on[i].push_back(p);
                on[j].push_back(p);
            }
        }
    }
    sort(ps.begin(), ps.end());
    ps.erase(unique(ps.begin(), ps.end()), ps.end());
    int M = ps.size();
    vector<vector<int>> g(M);
    for (int i = 0; i < N; i++) {
        vector<bool> known(M);
        for (const auto& p : on[i])
            known[lower_bound(ps.begin(), ps.end(), p) - ps.begin()] = true;
        vector<int> vec;
        for (int j = 0; j < M; j++) {
            if (known[j] || intersect(segs[i], ps[j])) vec.emplace_back(j);
        }
        for (int j = 1; j < (int)vec.size(); j++) {
            g[vec[j - 1]].push_back(vec[j]);
            g[vec[j]].push_back(vec[j - 1]);
        }
    }
    return g;
}

// 凸多角形を直線 l で切断し、l.a -> l.b の左側を返す。
// l.a != l.b を要求する。
Polygon convex_cut(const Polygon& U, Line l) {
    Polygon ret;
    for (int i = 0; i < (int)U.size(); i++) {
        const Point& now = U[i];
        const Point& nxt = U[(i + 1) % U.size()];
        auto cf = cross(l.a - now, l.b - now);
        auto cs = cross(l.a - nxt, l.b - nxt);
        if (sign(cf) >= 0) ret.emplace_back(now);
        if (sign(cf) * sign(cs) < 0)
            ret.emplace_back(crosspoint(Line(now, nxt), l));
    }
    return ret;
}

// 符号付き面積。反時計回りで正、時計回りで負。
long double area(const Polygon& p) {
    Real A = 0;
    for (int i = 0; i < (int)p.size(); i++)
        A += cross(p[i], p[(i + 1) % p.size()]);
    return A * 0.5L;
}

// 多角形と円の共通部分の符号付き面積。
// 多角形が反時計回りで正、時計回りで負。
// 各辺を交点で高々3区間に分ける。再帰を使わず O(N)。
long double area(const Polygon& p, const Circle& c) {
    if (p.size() < 3 || c.r == 0) return 0;
    Real A = 0;
    Circle origin(Point(0, 0), c.r);
    for (int i = 0; i < (int)p.size(); i++) {
        Point a = p[i] - c.p, b = p[(i + 1) % p.size()] - c.p;
        Point d = b - a;
        if (a == b) continue;
        auto [count, ts] = geometry_detail::circle_segment_parameters(origin, Segment(a, b));
        array<Real, 4> cuts{};
        for (int j = 0; j < count; j++) cuts[j + 1] = ts[j];
        cuts[count + 1] = 1;
        for (int j = 0; j <= count; j++) {
            if (cuts[j] == cuts[j + 1]) continue;
            Point u = a + d * cuts[j], v = a + d * cuts[j + 1];
            Point mid = a + d * ((cuts[j] + cuts[j + 1]) * 0.5L);
            if (norm(mid) <= c.r * c.r) A += cross(u, v);
            else A += c.r * c.r * atan2l(cross(u, v), dot(u, v));
        }
    }
    return A * 0.5L;
}

// 最遠点対の添字。異なる2点を選ぶ。2点未満なら {-1,-1}。
// 凸包構築 O(N log N)、回転キャリパー O(N)。
template <typename T>
pair<int, int> furthest_pair(const vector<PointBase<T>>& p) {
    if (p.size() < 2) return {-1, -1};
    using D = decltype(dot(p[0], p[0]));
    D best = -1;
    pair<int, int> ans = {0, 1};
    auto upd = [&](int i, int j) {
        auto delta = difference(p[i], p[j]);
        auto d = dot(delta, delta);
        if (d > best) {
            best = d;
            ans = {i, j};
        }
    };

    auto I = convex_hull(p);
    int n = I.size();
    if (n == 1) return {0, 1};
    if (n == 2) return {I[0], I[1]};
    int j = 1;
    for (int i = 0; i < n; i++) {
        int ni = (i + 1) % n;
        auto edge = difference(p[I[ni]], p[I[i]]);
        while (cross(edge, difference(p[I[(j + 1) % n]], p[I[j]])) > 0)
            j = (j + 1) % n;
        upd(I[i], I[j]);
        upd(I[ni], I[j]);
        if (cross(edge, difference(p[I[(j + 1) % n]], p[I[j]])) == 0) {
            upd(I[i], I[(j + 1) % n]);
            upd(I[ni], I[(j + 1) % n]);
        }
    }
    return ans;
}

// 最近点対の距離。2点以上を要求する。
long double closest_pair(Points ps) {
    if (ps.size() <= 1) throw(0);
    sort(ps.begin(), ps.end());
    auto compare_y = [](const Point& a, const Point& b) {
        return imag(a) < imag(b);
    };
    vector<Point> buf(ps.size());
    const Real INF = numeric_limits<Real>::infinity();
    auto rec = [&](auto& self, int left, int right) -> Real {
        if (right - left <= 1) return INF;
        int mid = (left + right) / 2;
        Real x = real(ps[mid]);
        Real ret = min(self(self, left, mid), self(self, mid, right));
        inplace_merge(ps.begin() + left, ps.begin() + mid,
                      ps.begin() + right, compare_y);
        int ptr = 0;
        for (int i = left; i < right; i++) {
            if (abs(real(ps[i]) - x) >= ret) continue;
            for (int j = ptr - 1; j >= 0; j--) {
                Point d = ps[i] - buf[j];
                if (imag(d) >= ret) break;
                ret = min(ret, abs(d));
            }
            buf[ptr++] = ps[i];
        }
        return ret;
    };
    return rec(rec, 0, ps.size());
}
