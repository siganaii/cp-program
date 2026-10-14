template<typename T>
struct imos_linearfunction {
    int N;
    vector<T> imos1,imos2;
    imos_linearfunction() = default;
    imos_linearfunction(int n) {
        N = n;
        imos1.resize(N + 1);
        imos2.resize(N + 1);
    }
    //[L,R)に(b,a+b,2a+b,...)を加算する
    void add(int L,int R,T a,T b) {
        imos1[L] += b - a * L;
        imos1[R] -= b - a * L;
        imos2[L] += a;
        imos2[R] -= a;
    };
    //buildの後はaddしない
    void build() {
        for(int i = 0;i < N;i++) {
            imos1[i+1] += imos1[i];
            imos2[i+1] += imos2[i];
        }
    };
    T get(int x) {
        return imos1[x] + imos2[x] * x;
    };
};