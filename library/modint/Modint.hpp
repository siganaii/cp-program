#pragma once
template <int mod>
struct Modint{
    int x;
    Modint():x(0) {}
    Modint(long long y): x(y % mod) {
        if (x < 0) x += mod;
    }
    Modint &operator += (const Modint &p) {
        // 加算前に判定し、途中の値が int の範囲を超えないようにする。
        if (x >= mod - p.x) x -= mod - p.x;
        else x += p.x;
        return *this;}
    Modint &operator -= (const Modint &p) {
        x -= p.x;
        if (x < 0) x += mod;
        return *this;}
    Modint &operator *= (const Modint &p) {
        x = (int)(1LL * x * p.x % mod);
        return *this;}
    Modint &operator /= (const Modint &p) {
        *this *= p.inverse();
        return *this;}
    Modint operator -() const{return Modint(-x);}
    Modint operator +(const Modint &p) const {return Modint(*this) += p;}
    Modint operator -(const Modint &p) const {return Modint(*this) -= p;}
    Modint operator *(const Modint &p) const {return Modint(*this) *= p;}
    Modint operator /(const Modint &p) const {return Modint(*this) /= p;}
    friend constexpr Modint operator +(long long lhs, const Modint &rhs) {
        return Modint(lhs) + rhs;
    }
    friend constexpr Modint operator -(long long lhs, const Modint &rhs) {
        return Modint(lhs) - rhs;
    }
    friend constexpr Modint operator *(long long lhs, const Modint &rhs) {
        return Modint(lhs) * rhs;
    }
    friend constexpr Modint operator /(long long lhs, const Modint &rhs) {
        return Modint(lhs) / rhs;
    }

    Modint &operator ++() {if(x == mod - 1) x = 0; else x++; return *this;}
    Modint &operator --() {if(x == 0) x = mod - 1; else x--; return *this;} 
    Modint operator ++(int) {Modint old = *this; ++(*this); return old;}
    Modint operator --(int) {Modint old = *this; --(*this); return old;}
    bool operator == (const Modint &p) const {return x == p.x;}
    bool operator != (const Modint &p) const {return x != p.x;}
    Modint inverse() const {
        int a = x, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            swap(a -= t * b, b);
            swap(u -= t * v, v);
        }
        return Modint(u);}
    Modint pow(long long n) const {
        Modint ret(1), mul(x);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;}
    friend ostream &operator<<(ostream &os, const Modint &p) { return os << p.x; }
    friend istream &operator>>(istream &is, Modint &a) {
        long long t;
        is >> t;
        a = Modint<mod>(t);
        return (is);
    }
    int get() const { return x; }
    static constexpr int get_mod() {return mod;}
};
