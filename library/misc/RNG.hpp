#pragma once
unsigned long long RNG_64() {
  static uint64_t x_
      = uint64_t(chrono::duration_cast<chrono::nanoseconds>(
                     chrono::high_resolution_clock::now().time_since_epoch())
                     .count())
        * 10150724397891781847ULL;
  x_ ^= x_ << 7;
  return x_ ^= x_ >> 9;
}

unsigned long long RNG(unsigned long long lim) { return RNG_64() % lim; }

// l < r を前提として、[l, r) の整数を返す。
long long RNG(long long l, long long r) {
  // 区間幅は符号なしで求め、負の下端への加算は拡張して行う。
  unsigned long long width = (unsigned long long)r - (unsigned long long)l;
  return (long long)((__int128_t)l + RNG_64() % width);
}
