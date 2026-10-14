#pragma once
#include "library/inner/inner-hash.hpp"
// vector<int> a = {0}とb = {0,0}の区別がつかない
template <typename Str, int BASE_NUM = 2>
struct RollingHash {
    using Hash = inner::Hash<BASE_NUM>;
    Str data;
    vector<Hash> hs,pw;
    int s;
    static Hash basis;

    RollingHash(const Str &S = Str()) {build(S);}
    void build(const Str &S) {
        data = S;
        s = S.size();
        hs.resize(s + 1);
        pw.resize(s + 1);
        pw[0] = Hash::set(1);
        hs[0] = Hash::set(0);
        for(int i = 1;i <= s; i++) {
            pw[i] = pw[i-1] * basis;
            hs[i] = pfma(hs[i-1],basis,S[i-1]);
        }
    }

    Hash get(int l,int r = -1) const{
        if(r == -1) r = s;
        return pfma(hs[l],-pw[r-l],hs[r]);
    }

    static Hash get_hash(const Str &T) {
        Hash ret = Hash::set(0);
        for(int i = 0;i < (int)T.size();i++) {
            ret = pfma(ret,basis,T[i]);
        }
        return ret;
    }
    // a + b の hash を返す
    // 引数 : a, b, b の長さ
    static Hash unite(Hash a, Hash b, long long bsize) {
        return pfma(a, basis.pow(bsize), b);
    }
    int find(Str &T, int lower = 0) const {
        auto ths = get_Hash(T);
        for(int i = lower; i <= s -(int)T.size();i++) {
            if(ths == get(i,i + (int)T.size())) return i;
        }  
        return -1;
    }
    friend int LCP(const RollingHash &a,RollingHash &b,int astart,int bstart) {
        int ok = 0,ng = min(a.size() - astart,b.size() - bstart) + 1;
        while (ok + 1 < ng) {
            int med = (ok + ng) / 2;
            (a.get(astart,med + astart) == b.get(bstart,med + bstart) ? ok:ng) = med;
        }
        return ok;
    }
    friend int strcmp(const RollingHash &a,RollingHash &b,int astart,int bstart,int aend = -1,int bend = -1) {
        if (aend == -1) aend = a.size();
        if (bend == -1) bend = b.size();
        if (aend - astart == bend - bstart) {
            if(a.get(astart,aend) == b.get(bstart,bend)) {
                return 0;
            }
        }
        else if(aend - astart < bend - bstart) {
            if(a.get(astart,aend) == b.get(bstart,bstart + aend - astart)) {
                return 1;
            }            
        }
        else {
            if(a.get(astart,astart + bend - bstart) == b.get(bstart,bend)) {
                return -1;
            }               
        }
        if (a.data[astart] < b.data[bstart]) {return 1;}
        else if (a.data[astart] > b.data[bstart]) {return -1;}
        int n = min<int>({LCP(a,b,astart,bstart),aend - astart,bend - bstart});
        return astart + n == aend ? bstart + n == bend ? 0 : 1
            : bstart + n == bend                    ? -1
            : a.data[astart + n] < b.data[bstart + n] ? 1
                                                : -1;
    }
    int size() const{return s;}
};
template <typename Str,int BASE_NUM>
typename RollingHash<Str,BASE_NUM>::Hash RollingHash<Str,BASE_NUM>::basis = inner::Hash<BASE_NUM>::get_basis();
