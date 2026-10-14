#pragma once
#include "formal-power-series.hpp"
// f(x)^k mod g(x)
template <typename mint>
FormalPowerSeries<mint> mod_pow(int64_t k, const FormalPowerSeries<mint>& f,
                                const FormalPowerSeries<mint>& g) {
	using fps = FormalPowerSeries<mint>;
	assert(k >= 0);
	fps modulus = g;
	modulus.shrink();
	assert(!modulus.empty());
	if (modulus.size() == 1) return {};
	fps result{1};
	if (k == 0) return result;
	fps base = f % modulus;
	if (k == 1) return base;
	// 剰余を取った各因子の次数は deg(modulus) 未満になる
	auto inv = modulus.rev().inv((int)modulus.size() - 1);
	auto reduce = [&](fps poly) {
		if (poly.size() >= modulus.size()) {
			int n = poly.size() - modulus.size() + 1;
			auto q = (poly.rev().pre(n) * inv.pre(n)).pre(n).rev();
			poly -= q * modulus;
		}
		poly.shrink();
		return poly;
	};
	while (k) {
		if (k & 1) result = reduce(result * base);
		k >>= 1;
		if (k) base = reduce(base * base);
	}
	return result;
}
