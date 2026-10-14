#pragma once
#include "library/modint/barrett-reduction.hpp"
struct ArbitraryModint {
    int x;
    ArbitraryModint():x(0) {}
    ArbitraryModint(int64_t y) {
        int z = y % get_mod();
        if(z < 0) z += get_mod();
        x = z;
    }
    ArbitraryModint &operator+=(const ArbitraryModint &p) {
        if((x += p.x) >= get_mod()) x -= get_mod();
        return *this;
    }
    ArbitraryModint &operator-=(const ArbitraryModint &p) {
        if((x += get_mod() - p.x) >= get_mod()) x -= get_mod();
        return *this;
    }
    ArbitraryModint &operator*=(const ArbitraryModint &p) {
        x = rem((unsigned long long)x * p.x);
        return *this;
    }
    ArbitraryModint &operator/=(const ArbitraryModint &p) {
        *this *= p.inverse();
        return *this;
    }
    ArbitraryModint operator-() const {return ArbitraryModint(-x);};
    ArbitraryModint operator+(const ArbitraryModint &p) const{
        return ArbitraryModint(*this) += p;
    }
    ArbitraryModint operator-(const ArbitraryModint &p) const{
        return ArbitraryModint(*this) -= p;
    }
    ArbitraryModint operator*(const ArbitraryModint &p) const{
        return ArbitraryModint(*this) *= p;
    }
    ArbitraryModint operator/(const ArbitraryModint &p) const {
        return ArbitraryModint(*this) /= p;
    }
    friend ArbitraryModint operator+(long long lhs, const ArbitraryModint &rhs) {
        return ArbitraryModint(lhs) + rhs;
    }
    friend ArbitraryModint operator-(long long lhs, const ArbitraryModint &rhs) {
        return ArbitraryModint(lhs) - rhs;
    }
    friend ArbitraryModint operator*(long long lhs, const ArbitraryModint &rhs) {
        return ArbitraryModint(lhs) * rhs;
    }
    friend ArbitraryModint operator/(long long lhs, const ArbitraryModint &rhs) {
        return ArbitraryModint(lhs) / rhs;
    }
    bool operator==(const ArbitraryModint &p) const {return x == p.x;}
    bool operator!=(const ArbitraryModint &p) const {return x != p.x;}
    ArbitraryModint &operator ++() {if(x == get_mod() - 1) x = 0; else x++; return *this;}
    ArbitraryModint &operator --() {if(x == 0) x = get_mod() - 1; else x--; return *this;} 
    ArbitraryModint operator ++(int) {ArbitraryModint old = *this; ++(*this); return old;}
    ArbitraryModint operator --(int) {ArbitraryModint old = *this; --(*this); return old;}
    ArbitraryModint inverse() const {
        int a = x,b = get_mod(),u = 1,v = 0,t;
        while(b > 0) {
            t = a / b;
            swap(a -= t * b,b);
            swap(u -= t * v,v);
        }
        assert(a == 1);
        return ArbitraryModint(u);
    }
    ArbitraryModint pow(int64_t n) const {
        ArbitraryModint ret(1),mul(x);
        while(n > 0) {
            if(n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }
    friend ostream &operator<<(ostream &os,const ArbitraryModint &p) {
        return os << p.x;
    }
    friend istream &operator>>(istream &is,ArbitraryModint &a) {
        int64_t t;
        is >> t;
        a = ArbitraryModint(t);
        return (is);
    }
    int get() const {return x;}
    static inline unsigned int rem(unsigned long long p) {return barrett().rem(p);};
    static inline Barrett &barrett() {
        static Barrett b;
        return b;
    }
    static inline int &get_mod() {
        static int mod = 0;
        return mod;
    }
    static void set_mod(int md) {
        assert(0 < md && md <= (1LL << 30) - 1);
        get_mod() = md;
        barrett() = Barrett(md);
    }
};