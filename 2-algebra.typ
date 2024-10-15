
#block( breakable: false,[
= Algebra


== Binary exponentiation

```cpp
ll m_pow(ll base, ll exp, ll mod) {
    base %= mod;
    ll result = 1;
    while (exp > 0) {
        if (exp & 1) result = ((ll)result * base) % mod;
        base = ((ll)base * base) % mod;
        exp >>= 1;
    }
    return result;
}
```
])

#block( breakable: false,[

== Extended euclidean
#label("gcd_ext")


$ a dot x + b dot y = gcd(a, b) $

```cpp
int gcd_ext(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1; y = 0;
        return a;
    }
    int x1, y1;
    int d = gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}
```

== Modular inversion & division

`gcd_ext` defined in #ref(label("gcd_ext"),).

$ exists x (a dot x equiv 1 (mod m)) arrow.l.r.double "gcd"(a,m)=1 $ 

```cpp
int mod_inv(int b, int m) {
    int x, y;
    int g = gcd_ext(b, m, &x, &y);
    if (g != 1) return -1;
    return (x%m + m) % m;
}
int m_divide(ll a, ll b, ll m) {
    int inv = mod_inv(b, m);
    assert(inv != -1);
    return (inv * (a % m)) % m;
}
```
])

#block( breakable: false,[
== Linear Diophantine equation

// `gcd_ext` defined in #ref(label("gcd_ext"),).

$ a dot x + b dot y = c $

$ {x = x_0 + k dot frac(b,g)  ;  y = y_0 - k dot frac(a,g) } $

```cpp
bool find_x0_y0(int a, int b, int c, int &x0, int &y0, int &g) {
    g = gcd_ext(abs(a), abs(b), x0, y0);
    if (c % g) return false;
    x0 *= c / g;
    y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}
```
])


#block( breakable: false,[

== Linear sieve

```cpp
const int N = 10000000;
vector<int> lp(N+1);
vector<int> pr;

for (int i=2; i <= N; ++i) {
    if (lp[i] == 0) {
        lp[i] = i;
        pr.push_back(i);
    }
    for (int j = 0; i * pr[j] <= N; ++j) {
        lp[i * pr[j]] = pr[j];
        if (pr[j] == lp[i]) break;
    }
}
```
])

#block(breakable: false,[

== Matrix multiplication
```cpp
struct Matrix:vector<vector<int>>
{
    // "inherit" vector's constructor
    using vector::vector;
    
    Matrix operator *(Matrix other)
    {
        int rows = size();
        int cols = other[0].size();
        Matrix res(rows, vector<int>(cols));
        for(int i=0;i<rows;i++)
            for(int j=0;j<other.size();j++)
                for(int k=0;k<cols;k++)
                    res[i][k]+=at(i).at(j)*other[j][k];
        return res;
    }
};
```
])


#block(breakable: false,[

== Euler's totient function
```cpp
int phi(int n) {
    int result = n;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if (n > 1)
        result -= result / n;
    return result;
}
void phi_1_to_n(int n) {
    vector<int> phi(n + 1);
    for (int i = 0; i <= n; i++)
        phi[i] = i;

    for (int i = 2; i <= n; i++) {
        if (phi[i] == i) {
            for (int j = i; j <= n; j += i)
                phi[j] -= phi[j] / i;
        }
    }
}
```
])

#block(breakable: false,[

== Gauss method
```cpp
const double EPS = 1e-9;
const int INF = 2; // it doesn't actually have to be infinity or a big number

int gauss (vector < vector<double> > a, vector<double> & ans) {
    int n = (int) a.size();
    int m = (int) a[0].size() - 1;

    vector<int> where (m, -1);
    for (int col=0, row=0; col<m && row<n; ++col) {
        int sel = row;
        for (int i=row; i<n; ++i)
            if (abs (a[i][col]) > abs (a[sel][col]))
                sel = i;
        if (abs (a[sel][col]) < EPS)
            continue;
        for (int i=col; i<=m; ++i)
            swap (a[sel][i], a[row][i]);
        where[col] = row;

        for (int i=0; i<n; ++i)
            if (i != row) {
                double c = a[i][col] / a[row][col];
                for (int j=col; j<=m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        ++row;
    }

    ans.assign (m, 0);
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i]][m] / a[where[i]][i];
    for (int i=0; i<n; ++i) {
        double sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * a[i][j];
        if (abs (sum - a[i][m]) > EPS)
            return 0;
    }

    for (int i=0; i<m; ++i)
        if (where[i] == -1)
            return INF;
    return 1;
}
```
])



#block( breakable: false,[

== FFT

```cpp
using ld = long double;
const int N = 1<<18;
const ld PI = acos(-1.0);
struct T {
  ld x, y;
  T() : x(0), y(0) {}
  T(ld a, ld b=0) : x(a), y(b) {}

  T operator/=(ld k) { x/=k; y/=k; return (*this); }
  T operator*(T a) const { return T(x*a.x - y*a.y, x*a.y + y*a.x); }
  T operator+(T a) const { return T(x+a.x, y+a.y); }
  T operator-(T a) const { return T(x-a.x, y-a.y); }
};

void fft(T* a, int n, int s) {
  for (int i=0, j=0; i<n; i++) {
    if (i>j) swap(a[i], a[j]);
    for (int l=n/2; (j^=l) < l; l>>=1);
  }

  for(int i = 1; (1<<i) <= n; i++){
    int M = 1 << i;
    int K = M >> 1;
    T wn = T(cos(s*2*PI/M), sin(s*2*PI/M));
    for(int j = 0; j < n; j += M) {
      T w = T(1, 0);
      for(int l = j; l < K + j; ++l){
        T t = w*a[l + K];
        a[l + K] = a[l]-t;
        a[l] = a[l] + t;
        w = wn*w;
      }
    }
  }
}

void multiply(T* a, T* b, int n) {
    while (n&(n-1)) n++; // ensure n is a power of two
    fft(a,n,1);
    fft(b,n,1);
    for (int i = 0; i < n; i++) a[i] = a[i]*b[i];
    fft(a,n,-1);
    for (int i = 0; i < n; i++) a[i] /= n;
}

int main() {
  // Example polynomials: (2 + 3x) and (1 - x)
  T a[10] = {T(2), T(3)};
  T b[10] = {T(1), T(-1)};
  multiply(a, b, 4);
  for (int i = 0; i < 10; i++)
    std::cout << int(a[i].x) << " ";
}
```
])

== janY modulo operations

```cpp
int mod;
int mod_f(long long a){
    return ((a%mod)+mod)%mod;
}
int m_add(long long a, long long b){
    return mod_f(mod_f(a)+mod_f(b));
}
int m_mult(long long a, long long b){
    return mod_f(mod_f(a)*mod_f(b));
}
int gcdExt(int a, int b, int *x, int *y){
    if(a==0){*x=0; *y=1; return b;}
    int x1,y1, gcd = gcdExt(b%a,a,&x1,&y1);
    *x = y1 - (b/a)*x1;
    *y = x1;
    return gcd;
}
int modInverse(int b, int m){
    int x,y, g = gcdExt(b,m,&x,&y);
    if(g!=1) return -1;
    return (x%m + m)%m;
}
int m_divide(long long a, long long b){
    a %= mod;
    int inv = modInverse(b, mod);
    return inv == -1 ? -1 : (1LL * inv * a) % mod;
}
int m_pow(long long base, long long exp){
    base %= mod; int res=1;
    while(exp>0){
        if(exp&1) res=(1LL*res*base)%mod;
        base=(1LL*base*base)%mod; exp >>=1;
    }
    return res;
}
```

== janY combinatorics

```cpp
// COMBINATORICS
long long nPr(long long n, long long r){
    if(r>n) return 0;
    if(n-r < r) r = n-r;
    long long res=1; for(long long i=0;i<r;i++) res *= (n--);
    return res;
}
long long nCr(long long n, long long r){
    if(r>n) return 0;
    if(n-r < r) r = n-r;
    long long res=1; for(long long i=0;i<r;i++) res *= (n--);
    for(long long i=1;i<=r;i++) res /= i;
    return res;
}
// Fast nCr with precalc_fact
int MAX_CHOOSE=3e5;
vector<long long> inverse_fact(MAX_CHOOSE+5), fact(MAX_CHOOSE+5);
long long fast_nCr(long long n, long long r){
    if(n<r || r<0) return 0;
    return m_mult(m_mult(fact[n], inverse_fact[r]),
        inverse_fact[n-r]);
}
void precalc_fact(int n){
    fact[0]=fact[1]=1;
    
    for(long long i=2;i<=n;i++) fact[i]=(fact[i-1]*i)%mod;
    inverse_fact[0]=inverse_fact[1]=1;
    
    for(long long i=2;i<=n;i++)
    inverse_fact[i]=(modInverse(i,mod)*inverse_fact[i-1])%mod;
}
```