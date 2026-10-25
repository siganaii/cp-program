#pragma once
template <long long mod>
struct Modint {
    using i64 = long long;
    using i128 = __int128_t;
    using u128 = __uint128_t;

    i64 x;

    constexpr Modint() : x(0) {}
    constexpr Modint(long long y)
        : x(y >= 0 ? y % mod : (mod - (-(__int128_t)y) % mod) % mod) {}

    constexpr Modint& operator+=(const Modint& p) {
        // 加算前に判定し、途中の値が i64 の範囲を超えないようにする。
        if (x >= mod - p.x) x -= mod - p.x;
        else x += p.x;
        return *this;
    }

    constexpr Modint& operator-=(const Modint& p) {
        x -= p.x;
        if (x < 0) x += mod;
        return *this;
    }

    constexpr Modint& operator*=(const Modint& p) {
        x = (i64)((i128)x * p.x % mod);
        return *this;
    }

    constexpr Modint& operator/=(const Modint& p) {
        *this *= p.inverse();
        return *this;
    }

    constexpr Modint operator-() const { return Modint(-x); }
    constexpr Modint operator+(const Modint& p) const { return Modint(*this) += p; }
    constexpr Modint operator-(const Modint& p) const { return Modint(*this) -= p; }
    constexpr Modint operator*(const Modint& p) const { return Modint(*this) *= p; }
    constexpr Modint operator/(const Modint& p) const { return Modint(*this) /= p; }

    constexpr Modint& operator++() {
        if (x == mod - 1) x = 0;
        else ++x;
        return *this;
    }

    constexpr Modint& operator--() {
        if (x == 0) x = mod - 1;
        else --x;
        return *this;
    }

    constexpr bool operator==(const Modint& p) const { return x == p.x; }
    constexpr bool operator!=(const Modint& p) const { return x != p.x; }

    constexpr Modint pow(long long n) const {
        Modint ret(1), mul(*this);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    constexpr Modint inverse() const {
        i64 a = x, b = mod;
        i64 u = 1, v = 0;
        while (b != 0) {
            i64 t = a / b;
            a -= t * b;
            std::swap(a, b);
            u -= t * v;
            std::swap(u, v);
        }
        return Modint(u);
    }

    friend std::ostream& operator<<(std::ostream& os, const Modint& p) {
        return os << p.x;
    }

    friend std::istream& operator>>(std::istream& is, Modint& a) {
        long long t;
        is >> t;
        a = Modint(t);
        return is;
    }

    constexpr long long get() const { return x; }
    static constexpr long long get_mod() { return mod; }
};
