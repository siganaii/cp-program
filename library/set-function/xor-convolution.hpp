#pragma once
#include "library/set-function/walsh-hadamard-transform.hpp"
// 非空入力は同じ長さの2の冪とする。
// modintでは長さの逆元が必要なため、長さ2以上の偶数modには対応しない。
template <typename T>
vector<T> xor_convolution(vector<T> a, vector<T> b) {
  if (a.empty() || b.empty()) return {};
  walsh_hadamard_transform(a);
  walsh_hadamard_transform(b);
  for (int i = 0; i < (int)a.size(); i++) a[i] *= b[i];
  walsh_hadamard_transform(a, true);
  return a;
}
