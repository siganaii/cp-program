#pragma once
vector<pair<long long,int>> prime_factorization(long long n) {
    vector<pair<long long,int>> ret;
    int c = 0;
    while(n % 2 == 0) {
        c++;
        n >>= 1;
    }
    if(c) ret.emplace_back(2,c);
    for(long long i = 3; i * i <= n; i += 2) {
        c = 0;
        while(n % i == 0) {
            n /= i;
            c++;
        }
        if(c) ret.emplace_back(i,c);
    }
    if (n != 1) ret.emplace_back(n,1);
    return ret;
}
vector<long long> divisor(long long n) {
    vector<long long> ret;
    for(long long i = 1; i * i <= n; i++) {
        if (n % i == 0) {
            ret.push_back(i);
            if(i * i != n) {ret.push_back(n / i);}
        }
    }
    sort(ret.begin(),ret.end());
    return ret;
}