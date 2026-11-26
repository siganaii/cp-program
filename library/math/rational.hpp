#pragma once
#include "library/math/bigint.hpp"

// 既約分数として保持する。分母は正、0 は 0/1。
class rational {
  bigint num, den;
  struct Reduced {};
  static bigint gcd(const bigint& lhs, const bigint& rhs) {
    bigint a = abs(lhs), b = abs(rhs);
    while (!b.is_zero()) {
      a %= b;
      swap(a, b);
    }
    return a;
  }
  void set_reduced(bigint& numerator, bigint& denominator) {
    if (denominator < 0) {
      numerator = -numerator;
      denominator = -denominator;
    }
    swap(num, numerator);
    swap(den, denominator);
  }
  rational(bigint numerator, bigint denominator, Reduced) {
    set_reduced(numerator, denominator);
  }
  void swap_members(rational& other) noexcept {
    swap(num, other.num);
    swap(den, other.den);
  }
  static rational add(const rational& lhs, const rational& rhs, bool subtract) {
    bigint common = gcd(lhs.den, rhs.den);
    bigint left = lhs.num * (rhs.den / common);
    bigint right = rhs.num * (lhs.den / common);
    bigint numerator = subtract ? left - right : left + right;
    if (numerator.is_zero()) return {};
    // 分母の gcd に対してだけ、さらに約分が必要な場合がある。
    bigint divisor = gcd(numerator, common);
    return {numerator / divisor, (lhs.den / common) * (rhs.den / divisor), Reduced{}};
  }

public:
  rational() : den(1) {}
  rational(bigint numerator) : den(1) { swap(num, numerator); }

  template <typename I, enable_if_t<inner::is_broadly_integral_v<I>, int> = 0>
  rational(I numerator) : rational(bigint(numerator)) {}

  // 分母が 0 の場合は domain_error を送出する。
  rational(bigint numerator, bigint denominator) {
    if (denominator.is_zero()) throw domain_error("zero rational denominator");
    if (numerator.is_zero()) {
      den = 1;
      return;
    }
    bigint divisor = gcd(numerator, denominator);
    if (divisor != 1) {
      numerator /= divisor;
      denominator /= divisor;
    }
    set_reduced(numerator, denominator);
  }

  // 十進整数または "分子/分母" を読み取る。
  explicit rational(const string& s) {
    auto slash = s.find('/');
    rational parsed = slash == string::npos
        ? rational(bigint(s))
        : rational(bigint(s.substr(0, slash)), bigint(s.substr(slash + 1)));
    swap_members(parsed);
  }

  rational(const rational&) = default;
  rational(rational&& other) : rational() { swap_members(other); }

  rational& operator=(rational other) noexcept {
    swap_members(other);
    return *this;
  }

  const bigint& numerator() const { return num; }
  const bigint& denominator() const { return den; }
  bool is_zero() const { return num.is_zero(); }

  friend rational operator+(const rational& lhs, const rational& rhs) {
    return add(lhs, rhs, false);
  }
  friend rational operator-(const rational& lhs, const rational& rhs) {
    return add(lhs, rhs, true);
  }

  friend rational operator*(const rational& lhs, const rational& rhs) {
    if (lhs.is_zero() || rhs.is_zero()) return {};
    bigint left = gcd(lhs.num, rhs.den), right = gcd(rhs.num, lhs.den);
    return {(lhs.num / left) * (rhs.num / right),
            (lhs.den / right) * (rhs.den / left), Reduced{}};
  }

  friend rational operator/(const rational& lhs, const rational& rhs) {
    if (rhs.is_zero()) throw domain_error("division by zero rational");
    if (lhs.is_zero()) return {};
    bigint numerators = gcd(lhs.num, rhs.num), denominators = gcd(lhs.den, rhs.den);
    return {(lhs.num / numerators) * (rhs.den / denominators),
            (lhs.den / denominators) * (rhs.num / numerators), Reduced{}};
  }

  rational& operator+=(const rational& rhs) { return *this = *this + rhs; }
  rational& operator-=(const rational& rhs) { return *this = *this - rhs; }
  rational& operator*=(const rational& rhs) { return *this = *this * rhs; }
  rational& operator/=(const rational& rhs) { return *this = *this / rhs; }

  rational operator+() const { return *this; }
  rational operator-() const { return {-num, den, Reduced{}}; }
  friend rational abs(const rational& value) { return {abs(value.num), value.den, Reduced{}}; }

  friend bool operator==(const rational& lhs, const rational& rhs) {
    return lhs.num == rhs.num && lhs.den == rhs.den;
  }
  friend bool operator!=(const rational& lhs, const rational& rhs) { return !(lhs == rhs); }
  friend bool operator<(const rational& lhs, const rational& rhs) {
    if (lhs.den == rhs.den) return lhs.num < rhs.num;
    return lhs.num * rhs.den < rhs.num * lhs.den;
  }
  friend bool operator>(const rational& lhs, const rational& rhs) { return rhs < lhs; }
  friend bool operator<=(const rational& lhs, const rational& rhs) { return !(rhs < lhs); }
  friend bool operator>=(const rational& lhs, const rational& rhs) { return !(lhs < rhs); }

  string to_string() const {
    if (den == 1) return num.to_string();
    return num.to_string() + '/' + den.to_string();
  }

  friend ostream& operator<<(ostream& os, const rational& value) {
    return os << value.to_string();
  }

  friend istream& operator>>(istream& is, rational& value) {
    string token;
    if (!(is >> token)) return is;
    try {
      value = rational(token);
    } catch (const invalid_argument&) {
      is.setstate(ios::failbit);
    } catch (const domain_error&) {
      is.setstate(ios::failbit);
    }
    return is;
  }
};
