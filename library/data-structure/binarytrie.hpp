#pragma once
// 格納できる値は 0 <= x < 2^LOG。LOG は T の値ビット数以下。
template <typename T,int LOG = 30,typename D = int>
struct BinaryTrie {
	static_assert(0 <= LOG && LOG <= numeric_limits<T>::digits,
	              "LOG must fit in the value bits of T");
	struct Node {
		D cnt = 0;
		int ch[2] = {-1, -1};
	};
	vector<Node> ns;
	BinaryTrie() : ns(1) {}
	D size() const { return ns[0].cnt; }
	T operator[](D k) const { return kth_element(k, 0); }
	// x ^ xor_addの小さい順のk番目のx ^ xor_addの値を返す
	T kth_element(D k, T xor_add = 0) const {
		assert(0 <= k && k < size());
		int idx = 0;
		T val = 0;
		for(int i = LOG - 1;i >= 0;i--) {
			T c = xor_add >> i & 1;
			int low_ch = ns[idx].ch[c];
			D low_cnt = (low_ch >= 0 ? ns[low_ch].cnt : 0);
			if (k < low_cnt) {
				idx = low_ch;
			} else {
				k -= low_cnt;
				idx = ns[idx].ch[c ^ 1];
				val ^= T(1) << i;
			}
			assert(idx >= 0);
		}
		return val;
  	}
	void add(T val, D cnt = 1) {
		assert(in_range(val));
		int idx = 0;
		for(int i = LOG - 1;i >= 0;i--) {
			ns[idx].cnt += cnt;
			assert(ns[idx].cnt >= 0);
			int &nxt = ns[idx].ch[val >> i & 1];
			if (nxt == -1) {
				idx = nxt = ns.size();
				ns.emplace_back();
			} else {
				idx = nxt;
			}
		}
		ns[idx].cnt += cnt;
		assert(ns[idx].cnt >= 0);
		return;
	}
	void erase(T val, D c = 1) {
		assert(in_range(val));
		assert(c >= 0);
		if (c == 0) return;
		assert(count(val) >= c);

		int idx = 0;
		for (int i = LOG - 1; i >= 0; --i) {
			ns[idx].cnt -= c;
			int b = val >> i & 1;
			idx = ns[idx].ch[b];
			assert(idx >= 0);
		}
		ns[idx].cnt -= c;
	}
	D lower_bound(T val, T xor_add = 0) {
		assert(0 <= val);
		if (!in_range(val)) return size();
		int idx = 0;
		D cnt = 0;
		for(int i = LOG - 1;i >= 0;i--) {
			int b = val >> i & 1, c = xor_add >> i & 1;
			int ch = ns[idx].ch[c];
			cnt += (b & (ch >= 0) ? ns[ch].cnt : 0);
			idx = ns[idx].ch[b ^ c];
			if (idx < 0 or ns[idx].cnt == 0) break;
		}
		return cnt;
	}
	D count(T val) const {
		assert(in_range(val));
		int idx = 0;
		for(int i = LOG - 1;i >= 0;i--) {
			idx = ns[idx].ch[val >> i & 1];
			if (idx < 0 or ns[idx].cnt == 0) return 0;
		}
		return ns[idx].cnt;
	}
	//[L,R)をcount
	D count(T L, T R, T xor_add = 0) {
		assert(0 <= L && L <= R);
		if constexpr (LOG < numeric_limits<T>::digits) {
			assert(R <= (T(1) << LOG));
		}
		return lower_bound(R, xor_add) - lower_bound(L, xor_add);
	}
	T min_element(T xor_add = 0) { return kth_element(0, xor_add); }
	T max_element(T xor_add = 0) { return kth_element(size() - 1, xor_add); }

private:
	static bool in_range(T val) {
		if (val < 0) return false;
		using U = make_unsigned_t<T>;
		if constexpr (LOG == numeric_limits<U>::digits) return true;
		else return U(val) < (U(1) << LOG);
	}
};
