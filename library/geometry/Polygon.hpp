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
// 各辺を円の内外で高々3区間に分ける。再帰を使わず O(N)。
long double area(const Polygon& p, const Circle& c) {
    if (p.size() < 3 || c.r == 0) return 0;
    bool inside = true, has_reference = false;
    Point reference_a, reference_b;
    Real reference_t = 0;
    for (const auto& v : p) {
        if (abs(v - c.p) > c.r) {
            inside = false;
        } else if (!has_reference || v < reference_a) {
            reference_a = reference_b = v;
            has_reference = true;
        }
    }
    // 円は凸なので、全頂点が内部なら多角形全体も内部にある。
    if (inside) return area(p);
    if (!has_reference) {
        // 円内の頂点がなければ、円と交わる辺上の垂線の足を基準にする。
        // 辺の選択と向きをそろえ、開始頂点や列挙方向に依存させない。
        for (int i = 0; i < (int)p.size(); i++) {
            Point a = p[i], b = p[(i + 1) % p.size()];
            if (b < a) swap(a, b);
            if (a == b) continue;
            if (has_reference && !(a < reference_a || (a == reference_a && b < reference_b))) continue;
            Point d = b - a;
            Real len = abs(d), along = dot(d / len, c.p - a);
            if (along < 0 || along > len) continue;
            if (abs(geometry_detail::line_cross(a, b, c.p)) / len > c.r) continue;
            reference_a = a;
            reference_b = b;
            reference_t = clamp(along / len, Real(0), Real(1));
            has_reference = true;
        }
    }
    if (!has_reference) {
        // 境界が円と交わらなければ、共通部分は円全体か空集合。
        // 小さい円も区別できるよう、中心の内外判定にはEPSを使わない。
        bool contains_center = false;
        for (int i = 0; i < (int)p.size(); i++) {
            Point a = p[i], b = p[(i + 1) % p.size()];
            if (a.y > b.y) swap(a, b);
            if (a.y <= c.p.y && c.p.y < b.y && geometry_detail::line_cross(a, b, c.p) < 0)
                contains_center = !contains_center;
        }
        return contains_center ? copysignl(PI * c.r * c.r, area(p)) : 0;
    }
    // 基準点の座標を作ると、細い多角形の外へ丸められる場合がある。
    // 元の端点との外積を補間し、基準点を辺上の点として扱う。
    auto reference_cross = [&](const Point& a, const Point& b, const Point& origin) {
        Real x = geometry_detail::difference_cross(a, b, origin, reference_a);
        if (reference_t == 0) return x;
        Real y = geometry_detail::difference_cross(a, b, origin, reference_b);
        Real left = (1 - reference_t) * x, right = reference_t * y;
        Real value = fmal(reference_t, y, left);
        if (abs(value) >= (abs(left) + abs(right)) * 0.5L) return value;

        // 外積同士が相殺するときだけ、積と和の丸め残差を保持して補間する。
        // 長い辺と小さい円でも、各外積を先に丸めて小さい成分を失わない。
        // 16個の積について、積と重み付けの残差を含めても高々64項。
        array<Real, 64> terms{};
        int count = 0;
        auto append = [&](Real carry) {
            if (carry == 0) return;
            int size = 0;
            for (int i = 0; i < count; i++) {
                Real sum = carry + terms[i], v = sum - carry;
                Real error = (carry - (sum - v)) + (terms[i] - v);
                if (error != 0) terms[size++] = error;
                carry = sum;
            }
            if (carry != 0) terms[size++] = carry;
            count = size;
        };
        auto product = [&](Real u, Real v, Real weight) {
            Real hi = u * v, lo = fmal(u, v, -hi);
            for (Real part : {hi, lo}) {
                Real scaled = part * weight;
                append(scaled);
                append(fmal(part, weight, -scaled));
            }
        };
        auto weighted_cross = [&](const Point& from, const Point& to, Real weight) {
            product(b.x, to.y, weight);   product(b.x, from.y, -weight);
            product(a.x, to.y, -weight);  product(a.x, from.y, weight);
            product(b.y, to.x, -weight);  product(b.y, from.x, weight);
            product(a.y, to.x, weight);   product(a.y, from.x, -weight);
        };
        // cross(b-a, reference_a-origin) + t * cross(b-a, reference_b-reference_a)。
        weighted_cross(origin, reference_a, 1);
        weighted_cross(reference_a, reference_b, reference_t);
        Real result = 0;
        for (int i = 0; i < count; i++) result += terms[i];
        return result;
    };
    Real A = 0, correction = 0;
    auto add = [&](Real value) {
        Real next = A + value;
        correction += abs(A) >= abs(value)
            ? (A - next) + value : (value - next) + A;
        A = next;
    };
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
    for (int i = 0; i < (int)p.size(); i++) {
        const Point& a = p[i];
        const Point& b = p[(i + 1) % p.size()];
        if (a == b) continue;
        Point u = a - c.p, v = b - c.p, d = b - a;
        Real la = abs(u), lb = abs(v), len = abs(d);
        // K は基準点から見た三角形の符号付き面積の2倍。
        Real K = reference_cross(a, b, a);
        if (la <= c.r && lb <= c.r) {
            add(K);
            continue;
        }
        // R は円の中心から基準点へのベクトル。H = cross(u, v), X = cross(u, R), Y = cross(R, v)。
        Real H = geometry_detail::line_cross(a, b, c.p);
        Real X = reference_cross(c.p, a, c.p);
        Real Y = -reference_cross(c.p, b, c.p);
        auto arc = [&](Real theta) { add(c.r * c.r * arc_remainder(theta)); };
        auto outside = [&]() {
            Real sa = c.r / la, sb = c.r / lb;
            // 拡大率から1を引かず、半径との差から補正項を求める。
            Real da = (c.r - la) / la, db = (c.r - lb) / lb;
            // 円周に写した端点と基準点の外積を、元の座標の外積で展開する。
            add(sa * sb * K);
            add(sa * db * X);
            add(da * sb * Y);
            arc(atan2l(H, dot(u, v)));
        };
        Real h = H / len;
        if (abs(h) >= c.r) {
            outside();
            continue;
        }
        Point e = d / len;
        // 円の中心から辺への垂線の足を原点とする、辺方向の座標。
        Real xa = dot(e, u), xb = dot(e, v);
        Real extent = sqrtl(max(Real(0), (c.r - abs(h)) * (c.r + abs(h))));
        if (xa >= extent || xb <= -extent) {
            outside();
            continue;
        }
        // 交点座標を作らず、区間の長さと元の辺の外積から面積を求める。
        // 長い辺を小さい円で切る場合も、近いパラメータ同士を引かない。
        Real fraction = clamp((min(xb, extent) - max(xa, -extent)) / len, Real(0), Real(1));
        add(fraction * K);
        Real S = -reference_cross(a, b, c.p);
        // 垂線の足までの比率と、そこから交点までの短い距離を分ける。
        // 長い辺では交点までの比率を1値に丸めると円内の幅が消えてしまう。
        auto clipped_cross = [&](Real along, Real endpoint_cross, Real t) {
            // 交点が端点に近い場合は、短い区間の比率をそのまま使う。
            if (extent > along * 0.5L) return fmal(t, S, -endpoint_cross);
            Real center = along / len;
            Real rest = (fmal(-center, len, along) - extent) / len;
            return fmal(center, S, -endpoint_cross) + rest * S;
        };
        if (xa < -extent) {
            Real t = clamp((-extent - xa) / len, Real(0), Real(1));
            Real scale = c.r / la, delta = (c.r - la) / la;
            add(scale * t * K);
            add(delta * clipped_cross(-xa, X, t));
            arc(atan2l(t * H, fmal(xa, -extent, h * h)));
        }
        if (xb > extent) {
            Real t = clamp((xb - extent) / len, Real(0), Real(1));
            Real scale = c.r / lb, delta = (c.r - lb) / lb;
            add(scale * t * K);
            add(delta * clipped_cross(xb, Y, t));
            arc(atan2l(t * H, fmal(extent, xb, h * h)));
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
