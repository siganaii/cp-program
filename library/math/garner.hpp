#pragma once
#include "library/math/factorize.hpp"
#include "library/modint/ArbitaryModint.hpp"
//存在しない時、-1を返す
template<typename T,typename U>
long long garner(vector<T> rem,vector<U> MOD, long long new_mod = -1, bool coprime = true) {
    assert (rem.size() == MOD.size());
    int n = rem.size();
    using mint = ArbitraryModint;
    if(!coprime) {
        unordered_map<long long,vector<pair<long long,long long>>> mp;
        for(int i = 0;i < n;++i) {
            for(auto [p,e]:prime_factorization(MOD[i])) {
                int m = 1;
                for(int j = 0;j < e;++j) m *= p;
                mp[p].emplace_back(rem[i] % m,m);
            }
        }
        vector<T> xx;
        vector<U> mm;
        for(auto &[p,dat]:mp) {
            long long M = 1;
            long long val = 0;
            for(auto &[x,m]:dat) {
                if(M < m) {
                    M = m;
                    val = x;
                }
            }
            for(auto &[x,m]:dat) {
                if((val - x) % m != 0) return -1;
            }
            xx.emplace_back(val);
            mm.emplace_back(M);
        }
        swap(rem,xx);
        swap(MOD,mm);
        n = (int)rem.size();
    }
    // 呼び出し元のmod設定を、Barrettの状態も含めて保存する。
    int saved_mod = mint::get_mod();
    auto saved_barrett = mint::barrett();
    vector<long long> cfs(n);
    for(int i = 0;i < n;++i) {
        mint::set_mod(MOD[i]);
        long long a = rem[i];
        long long prod = 1;
        for(int j = 0;j < i;++j) {
            a = (mint(cfs[j]) * (-prod) + a).x;
            prod = (mint(prod) * MOD[j]).x;
        }
        cfs[i] = (mint(prod).inverse() * a).x;
    }
    mint::get_mod() = saved_mod;
    mint::barrett() = saved_barrett;
    long long ret = 0;
    long long prod = 1;
    for(int i = 0;i < n;++i) {
        __int128_t next_ret = ret + __int128_t(prod) * cfs[i];
        ret = new_mod == -1 ? next_ret : next_ret % new_mod;
        if(i + 1 < n) {
            __int128_t next_prod = __int128_t(prod) * MOD[i];
            prod = new_mod == -1 ? next_prod : next_prod % new_mod;
        }
    }
    return ret;
}
