#pragma once
#include "library/prime/prime_enumerate.hpp"
map<long long,int> Mertens(long long N) {
    map<long long,int> ret;
    ret[0] = 0;
    if(N == 0) return ret;
    int K = pow(N,(double)2 / 3);
    int L = (N + K - 1) / K;
    vector<int> b(K + 1);
    vector<int> primes = prime_enumerate(K);
    b[1] = 1;
    for(int &p:primes) {
        for(int i = 1;i <= K / p;i++) {
            if(b[i] == 0) continue;
            long long n = p * i;
            while(n <= K) {
                b[n] -= b[i];
                n *= p;
            }
        }
    }
    vector<int> B(K+1);
    for(int i = 1;i <= K;i++) B[i] = b[i] + B[i-1];
    vector<long long> B2(L+1,1);
    for(int i = L;i >= 1;i--) {
        long long NI = N / i;
        int M = int(sqrt(NI));
        for(int j = 1;j <= M;j++) {
            long long NIJ = NI / j;
            if(j > 1) {
                if(NIJ <= K) B2[i] -= B[NIJ];
                else B2[i] -= B2[i * j];
            }
            B2[i] -= b[j] * (NIJ - M);
        }
    }
    long long i = 1;
    while(i <= N) {
        long long NI = N / i;
        long long nex = N / NI;
        if(nex > N) nex = N;
        if(NI <= K) ret[NI] = B[NI];
        else ret[NI] = B2[i];
        i = nex + 1;
    }
    return ret;
}
