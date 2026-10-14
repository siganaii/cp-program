#pragma once

template <uint32_t mod>
struct LazyMontgomeryModInt {
    using mint = LazyMontgomeryModInt;
    using i32 = int32_t;
    using u32 = uint32_t;
    using u64 = uint64_t;
    static constexpr u32 get_r() {
        u32 ret = mod;
        for (i32 i = 0; i < 4; ++i) {
            ret *= 2 - mod * ret;
        }
        return ret;
    }
    static constexpr u32 r = get_r();
    static constexpr u32 n2 = -u64(mod) % mod;
    static_assert(r * mod == 1);
    static_assert(mod < (1U << 30));
    static_assert((mod & 1) == 1);
    u32 a;
    constexpr LazyMontgomeryModInt() : a(0) {}
    constexpr LazyMontgomeryModInt(int64_t b)
        : a(reduce(u64(b % mod + mod) * n2)) {}
    static constexpr u32 reduce(u64 b) {
        return (b + u64(u32(b) * u32(-r)) * mod) >> 32;
    }
    constexpr mint& operator+=(const mint& b) {
        if (i32(a += b.a - 2 * mod) < 0) {
            a += 2 * mod;
        }
        return *this;
    }
    constexpr mint& operator-=(const mint& b) {
        if (i32(a -= b.a) < 0) {
            a += 2 * mod;
        }
        return *this;
    }
    constexpr mint& operator*=(const mint& b) {
        a = reduce(u64(a) * b.a);
        return *this;
    }
    constexpr mint& operator/=(const mint& b) {
        return *this *= b.inverse();
    }
    friend constexpr mint operator+(mint lhs, const mint& rhs) {
        return lhs += rhs;
    }
    friend constexpr mint operator-(mint lhs, const mint& rhs) {
        return lhs -= rhs;
    }
    friend constexpr mint operator*(mint lhs, const mint& rhs) {
        return lhs *= rhs;
    }
    friend constexpr mint operator/(mint lhs, const mint& rhs) {
        return lhs /= rhs;
    }
    friend constexpr mint operator+(long long lhs, const mint& rhs) {
        return mint(lhs) + rhs;
    }
    friend constexpr mint operator-(long long lhs, const mint& rhs) {
        return mint(lhs) - rhs;
    }
    friend constexpr mint operator*(long long lhs, const mint& rhs) {
        return mint(lhs) * rhs;
    }
    friend constexpr mint operator/(long long lhs, const mint& rhs) {
        return mint(lhs) / rhs;
    }
    constexpr bool operator==(const mint& b) const {
        return (a >= mod ? a - mod : a)
            == (b.a >= mod ? b.a - mod : b.a);
    }
    constexpr bool operator!=(const mint& b) const {
        return !(*this == b);
    }
    constexpr mint& operator++() {
        *this += mint(1);
        return *this;
    }
    constexpr mint& operator--() {
        *this -= mint(1);
        return *this;
    }
    constexpr mint operator++(int) {
        mint old = *this;
        ++*this;
        return old;
    }
    constexpr mint operator--(int) {
        mint old = *this;
        --*this;
        return old;
    }
    constexpr mint operator-() const {
        return mint(0) - *this;
    }
    constexpr mint pow(u64 n) const {
        mint ret(1);
        mint mul(*this);
        while (n > 0) {
            if (n & 1) {
                ret *= mul;
            }
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }
    constexpr mint inverse() const {
        return pow(mod - 2);
    }
    friend std::ostream& operator<<(std::ostream& os, const mint& b) {
        return os << b.get();
    }
    friend std::istream& operator>>(std::istream& is, mint& b) {
        int64_t t;
        is >> t;
        b = mint(t);
        return is;
    }
    constexpr u32 get() const {
        u32 ret = reduce(a);
        return ret >= mod ? ret - mod : ret;
    }
    static constexpr u32 get_mod() {
        return mod;
    }
};