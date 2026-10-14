#pragma once
#include "library/geometry/geometry-base.hpp"
#include "library/geometry/crosspoint.hpp"
#include "library/geometry/intersect.hpp"

// 単純多角形を前提とする。時計回り・反時計回りの両方に対応。
template <typename T>
bool is_convex(const vector<PointBase<T>>& p) {
    int n = p.size(), dir = 0;
    for (int i = 0; i < n; i++) {
        int d = geometry_detail::orientation_sign(p[i], p[(i + 1) % n], p[(i + 2) % n]);
        if (d == 0) continue;
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
        if (geometry_detail::orientation_sign(XY[I[0]], XY[I.back()], XY[I[k]]) != 0) {
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
                int d = geometry_detail::orientation_sign(XY[H[H.size() - 2]], XY[H.back()], XY[k]);
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
        Point a = Q[i], b = Q[(i + 1) % Q.size()];
        if (a.y > b.y) swap(a, b);
        if (ccw(Q[i], Q[(i + 1) % Q.size()], p) == 0) return ON;
        // 高さは半開区間で比較し、向きは座標差分の丸め誤差も補正する。
        if (a.y <= p.y && p.y < b.y && geometry_detail::line_cross(a, b, p) < 0)
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
        return ccw(a, b, p) == 0;
    };
    if (n == 2) return on_segment(Q[0], Q[1]) ? ON : OUT;
    int dir = geometry_detail::orientation_sign(Q[0], Q[1], Q.back());
    if (dir == 0) return contains(Q, p);

    auto side = [&](const Point& a, const Point& b) {
        return dir * geometry_detail::line_side(Q[0], a, b);
    };
    int a = side(Q[1], p), b = side(Q.back(), p);
    if (on_segment(Q[0], Q[1]) || on_segment(Q.back(), Q[0])) return ON;
    // 扇形の境界付近では、端点周辺も含めて通常の内外判定にそろえる。
    if (a == 0 || b == 0) return contains(Q, p);
    if (a < 0 || b > 0) return OUT;

    int l = 1, r = n - 1;
    while (r - l > 1) {
        int m = (l + r) / 2;
        // 内部の対角線ではEPSを使わず、実際の向きから調べる辺を選ぶ。
        if (dir * geometry_detail::line_cross(Q[0], Q[m], p) >= 0) l = m;
        else r = m;
    }
    // 頂点の近くでは、隣の辺の許容誤差内に入る場合もある。
    if (on_segment(Q[l - 1], Q[l]) || on_segment(Q[l], Q[r])
        || on_segment(Q[r], Q[(r + 1) % n])) return ON;
    int s = dir * geometry_detail::line_side(Q[l], Q[r], p);
    return s <= 0 ? OUT : IN;
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
        vector<Segment> group;
        for (int j = i; j < (int)S.size(); j++) {
            if (used[j]) continue;
            if (!parallel(S[i], S[j])) continue;
            if (geometry_detail::line_side(S[i].a, S[i].b, S[j].a) != 0) continue;
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
            if (parallel(segs[i], segs[j])) continue;
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
        int cf = geometry_detail::line_side(l.a, l.b, now);
        int cs = geometry_detail::line_side(l.a, l.b, nxt);
        if (cf >= 0) ret.emplace_back(now);
        if (cf * cs < 0)
            ret.emplace_back(crosspoint(Line(now, nxt), l));
    }
    return ret;
}

// 符号付き面積。反時計回りで正、時計回りで負。
long double area(const Polygon& p) {
    if (p.size() < 3) return 0;
    Real A = 0, correction = 0;
    // 1頂点を基準にし、座標差分・外積・面積の加算それぞれの丸め誤差を補正する。
    for (int i = 1; i + 1 < (int)p.size(); i++) {
        Real value = geometry_detail::line_cross(p[0], p[i], p[i + 1]);
        Real next = A + value;
        correction += abs(A) >= abs(value) ? (A - next) + value : (value - next) + A;
        A = next;
    }
    return (A + correction) * 0.5L;
}

// 多角形と円の共通部分の符号付き面積。
// 多角形が反時計回りで正、時計回りで負。
// 各辺を交点で高々3区間に分ける。再帰を使わず O(N)。
long double area(const Polygon& p, const Circle& c) {
    if (p.size() < 3 || c.r == 0) return 0;
    bool inside = true;
    for (const auto& v : p) {
        if (abs(v - c.p) > c.r) {
            inside = false;
            break;
        }
    }
    // 円は凸なので、全頂点が内部なら多角形全体も内部にある。
    if (inside) return area(p);
    Real A = 0, correction = 0;
    auto add = [&](Real value) {
        Real next = A + value;
        correction += abs(A) >= abs(value)
            ? (A - next) + value : (value - next) + A;
        A = next;
    };
    auto on_circle = [&](Point v) {
        Real len = abs(v);
        return len == 0 ? v : v * (c.r / len);
    };
    Point reference = p[0] - c.p;
    if (abs(reference) > c.r) reference = on_circle(reference);
    auto arc_remainder = [](Real theta) {
        if (abs(theta) >= 0.1L) return theta - sinl(theta);
        // 小さい角度でも theta - sin(theta) の有効桁を失わないようにする。
        Real t2 = theta * theta, term = theta * t2 / 6, sum = term;
        for (int k = 5; k <= 13; k += 2) {
            term *= -t2 / ((k - 1) * k);
            sum += term;
        }
        return sum;
    };
    Circle origin(Point(0, 0), c.r);
    for (int i = 0; i < (int)p.size(); i++) {
        Point a = p[i] - c.p, b = p[(i + 1) % p.size()] - c.p;
        if (a == b) continue;
        auto [count, ps] = geometry_detail::circle_segment_intersections(origin, Segment(a, b));
        array<Point, 4> cuts{};
        cuts[0] = a;
        for (int j = 0; j < count; j++) cuts[j + 1] = ps[j];
        cuts[count + 1] = b;
        for (int j = 0; j <= count; j++) {
            if (cuts[j] == cuts[j + 1]) continue;
            Point u = cuts[j], v = cuts[j + 1];
            Point mid = (u + v) * 0.5L;
            if (abs(mid) <= c.r) {
                add(geometry_detail::cross_value(u - reference, v - reference));
            } else {
                Real theta = atan2l(geometry_detail::cross_value(u, v), dot(u, v));
                Point a = on_circle(u), b = on_circle(v);
                // 円外の部分は円弧に写す。近い基準点からの三角形と弓形に分ける。
                add(geometry_detail::cross_value(a - reference, b - reference));
                add(c.r * c.r * arc_remainder(theta));
            }
        }
    }
    return (A + correction) * 0.5L;
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
        // 凸包と同様に、座標差分と外積の丸め誤差を補正して探索方向を判定する。
        auto edge_cross = [&](int k) {
            return geometry_detail::difference_cross(
                p[I[i]], p[I[ni]], p[I[k]], p[I[(k + 1) % n]]);
        };
        auto turn = edge_cross(j);
        while (turn > 0) {
            j = (j + 1) % n;
            turn = edge_cross(j);
        }
        upd(I[i], I[j]);
        upd(I[ni], I[j]);
        if (turn == 0) {
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
