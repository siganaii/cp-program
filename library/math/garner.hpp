#pragma once
#include "library/math/factorize.hpp"
#include "library/modint/ArbitaryModint.hpp"
using mint = ArbitraryModint;
//存在しない時、-1を返す
template<typename T,typename U>
long long garner(vector<T> rem,vector<U> MOD, long long new_mod = -1, bool coprime = true) {
    assert (rem.size() == MOD.size());
    int n = rem.size();
    using u64 = unsigned long long;
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
    vector<ll> cfs(n);
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
    long long ret = 0;
    long long prod = 1;
    for(int i = 0;i < n;++i) {
        ret += prod * cfs[i];
        prod *= MOD[i];
        if(new_mod != -1) {
            ret %= new_mod;
            prod %= new_mod;
        }
    }
    return ret;
}