#pragma once
// f(x)^k mod g(x)
template <typename mint>
template <typename mint>
FormalPowerSeries<mint> mod_pow(int64_t k, const FormalPowerSeries<mint>& f,
                                const FormalPowerSeries<mint>& g) {
	using fps = FormalPowerSeries<mint>;
	assert(!g.empty());
	auto inv = g.rev().inv();
	auto quo = [&](const fps& poly) {
		if (poly.size() < g.size()) return fps{};
		int n = poly.size() - g.size() + 1;
		return (poly.rev().pre(n) * inv.pre(n)).pre(n).rev();
	};
	auto res = [&](auto &self, int K) -> fps {
		if(K == 0) return {1};
		if(K == 1) return f;
		fps F = self(self,K/2);
		F *= F;
		F -= quo(F) * g;
		F.shrink();
		if(K & 1) {
			F *= f;
			F -= quo(F) * g;
			F.shrink();
		}
		return F;
	};
	return res(res,k);
}