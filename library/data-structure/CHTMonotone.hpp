#pragma once
template<typename T,bool is_min>
struct CHTMonotone {
    deque<pair<T,T>> deq;
    CHTMonotone() = default;
    bool empty() const {return deq.empty();}
    void clear() {deq.clear();}
    inline int sgn(T x) {return x == 0 ? 0 : (x < 0 ? -1 : 1);}
    inline bool check(const pair<T,T> &a,const pair<T,T> &b,const pair<T,T> &c) {
        if constexpr(is_integral<T>::value) {
            // 傾きは a.first > b.first > c.first。64-bit 整数の差も拡張してから取る。
            __int128 x = (__int128)b.second - a.second;
            __int128 y = (__int128)c.second - b.second;
            if ((x < 0) != (y < 0)) return x >= 0;
            // 絶対値の積は unsigned 128-bit に収め、符号付き積の overflow を避ける。
            __uint128_t lhs = (__uint128_t)(x < 0 ? -x : x) *
                              (__uint128_t)((__int128)b.first - c.first);
            __uint128_t rhs = (__uint128_t)(y < 0 ? -y : y) *
                              (__uint128_t)((__int128)a.first - b.first);
            return x < 0 ? lhs <= rhs : lhs >= rhs;
        }
        else {
            if(b.second == a.second || c.second == b.second) return sgn(b.first - a.first) * sgn(c.second - b.second) >= sgn(b.second - a.second) * sgn(c.first - b.first);
            return (b.first - a.first) * sgn(c.second - b.second) / abs(b.second - a.second) >= (c.first - b.first) * sgn(b.second - a.second) / abs(c.second - b.second);
        }
    }
    void add(T a,T b) {
        if(!is_min) a *= -1,b *= -1;
        pair<T,T> line(a,b);
        if(empty()) {
            deq.emplace_front(line);
            return;
        }
        if(deq.front().first <= a) {
            if(deq.front().first == a) {
                if(deq.front().second <= b) return;
                deq.pop_front();
            }
            while(deq.size() >= 2 && check(line,deq.front(),deq[1])) deq.pop_front();
            deq.emplace_front(line);
        }
        else {
            assert(a <= deq.back().first);
            if(deq.back().first == a) {
                if(deq.back().second <= b) return;
                deq.pop_back();
            }
            while(deq.size() >= 2 && check(deq[deq.size() - 2],deq.back(),line)) deq.pop_back();
            deq.emplace_back(line);
        }
    }
    inline T get_y(const pair<T,T> &a,const T &x) {
        return a.first * x + a.second;
    }
    T query(T x) {
        assert(!empty());
        int l = -1,r = deq.size() - 1;
        while(l + 1 < r) {
            int m = (l + r) >> 1;
            if(get_y(deq[m],x) >= get_y(deq[m+1],x)) l = m;
            else r = m;
        }
        if(is_min) return get_y(deq[r],x);
        else return -get_y(deq[r],x);
    }
    T query_monotone_inc(T x) {
        assert(!empty());
        while(deq.size() >= 2 && get_y(deq.front(),x) >= get_y(deq[1],x)) deq.pop_front();
        if(is_min) return get_y(deq.front(),x);
        else return -get_y(deq.front(),x);
    }
    T query_monotone_dec(T x) {
        assert(!empty());
        while(deq.size() >= 2 && get_y(deq.back(),x) >= get_y(deq[deq.size() - 2],x)) deq.pop_back();
        if(is_min) return get_y(deq.back(),x);
        else return -get_y(deq.back(),x);
    }
};
