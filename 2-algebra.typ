
#block( breakable: false,[
= Algebra


== Binary exponentiation

```cpp
ll m_pow(ll base, ll exp, ll mod) {
    base %= mod; ll result = 1;
    while (exp > 0) {
        if (exp & 1) result = (result * base) % mod;
        base = (base * base) % mod; exp >>= 1;
    }
    return result;
}
```
])

#block( breakable: false,[

== Extended euclidean
#label("gcd_ext")

Find integers $x$ and $y$ such that:
$a dot x + b dot y = gcd(a, b)$

```cpp
int gcd_ext(int a, int b, int& x, int& y) {
    if (b == 0) { x = 1; y = 0; return a; }
    int x1, y1; int d = gcd_ext(b, a % b, x1, y1);
    x = y1; y = x1 - y1 * (a / b); return d;
}
```

== Modular inversion & division

Mod inverse exists iff number is coprime with mod.
$ exists x (a dot x equiv 1 (mod m)) arrow.l.r.double "gcd"(a,m)=1 $ 

```cpp
int mod_inv(int b, int m) {
    int x, y; int g = gcd_ext(b, m, &x, &y);
    if (g != 1) return -1;
    return (x%m + m) % m;
}
int m_divide(ll a, ll b, ll m) {
    int inv = mod_inv(b, m); assert(inv != -1);
    return (inv * (a % m)) % m;
}
```
])

#block( breakable: false,[
== Linear Diophantine equation

// `gcd_ext` defined in #ref(label("gcd_ext"),).


$ {(x,y) in ZZ^2 | a dot x + b dot y = c } = {x_0 + k dot (b slash g), y_0 - k dot (a slash g) | k in ZZ } $

```cpp
bool find_x0_y0(int a, int b, int c, int &x0, int &y0, int &g) {
    g = gcd_ext(abs(a), abs(b), x0, y0);
    if (c % g) return false;
    x0 *= c / g; y0 *= c / g;
    if (a < 0) x0 = -x0; if (b < 0) y0 = -y0;
    return true;
}
```
])


#block( breakable: false,[

== Linear sieve

```cpp
const int N=10000000; vector<int> lp(N+1), pr;
for(int i=2;i<=N;i++){
    if(!lp[i]){ lp[i]=i; pr.push_back(i); }
    for(int j=0;j<pr.size() && i*pr[j]<=N;j++){
        lp[i*pr[j]]=pr[j];
        if(pr[j]==lp[i]) break;
    }
}
```
])

#block(breakable: false,[

== Matrix multiplication

Inherit vector's constructor to allow brace initialization.

```cpp
struct Matrix:vector<vector<int>>{
    using vector::vector;
    Matrix operator *(const Matrix& other){
        int rows = size(); int cols = other[0].size();
        Matrix res(rows, vector<int>(cols));
        for(int i=0;i<rows;i++) for(int j=0;j<other.size();j++)
            for(int k=0;k<cols;k++)
                res[i][k]+=at(i).at(j)*other[j][k];
        return res;
    }
};
```
Usage example (prints `403 273 234 442`):

```cpp
    Matrix A = {{19,7},{6, 20}}, B = {{19,7},{6, 20}}, C = A*B;
    for(auto rows: C) for(int cell: rows) cout<<cell<<" ";
```
])


#block(breakable: false,[

== Euler's totient function
```cpp
int phi(int n){
    int res=n; for(int i=2;i*i<=n;i++) if(n%i==0)
        { while(n%i==0)n/=i; res-=res/i; }
    if(n>1) res-=res/n; return res;
}
void phi_1_to_n(int n){
    vector<int> phi(n+1); for(int i=0;i<=n;i++) phi[i]=i;
    for(int i=2;i<=n;i++) if(phi[i]==i)
        for(int j=i;j<=n;j+=i) phi[j]-=phi[j]/i;
}
```
])


#block(breakable: false,[
== Gauss method

System of $n$ linear algebraic equations (SLAE) with $m$ variables.
$ cases(a_(1 1) x_1 + a_(1 2) x_2 + ... + a_(1 m) x_m = b_1,...,a_(n 1) x_1 + a_(n 2) x_2 + ... + a_(n m) x_m = b_n) $

Matrix representation: $A x = b$. Gauss-Jordan elimination impl.:


```cpp
const double EPS=1e-9; const int INF=2;
int gauss(vector<vector<double>> a, vector<double> &ans){
    // last column of matrix a is vector b
    int n=a.size(), m=a[0].size()-1; vector<int> where(m,-1);
    for(int col=0, row=0; col<m && row<n; ++col){
        int sel=row; for(int i=row;i<n;i++)
            if(abs(a[i][col])>abs(a[sel][col])) sel=i;
        if(abs(a[sel][col])<EPS) continue;
        for(int i=col;i<=m;i++) swap(a[sel][i],a[row][i]);
        where[col]=row; for(int i=0;i<n;i++) if(i!=row){
            double c=a[i][col]/a[row][col];
            for(int j=col;j<=m;j++) a[i][j]-=a[row][j]*c;} row++;
    }
    ans.assign(m,0); for(int i=0;i<m;i++)
        if(where[i]!=-1) ans[i]=a[where[i]][m]/a[where[i]][i];
    for(int i=0;i<n;i++){
        double sum=0; for(int j=0;j<m;j++) sum+=ans[j]*a[i][j];
        if(abs(sum-a[i][m])>EPS) return 0;}
    for(int i=0;i<m;i++) if(where[i]==-1) return INF; return 1;
}
```
])



#block( breakable: false,[

== FFT

```cpp
const int N=1<<18;
const ld PI=acos(-1.0);
struct T{
    ld x,y;
    T():x(0),y(0){}
    T(ld a, ld b=0):x(a),y(b){}
    T operator/=(ld k){x/=k;y/=k;return *this;}
    T operator*(const T&a) const {
        return T(x*a.x-y*a.y, x*a.y+y*a.x);}
    T operator+(const T&a) const {return T(x+a.x, y+a.y);}
    T operator-(const T&a) const {return T(x-a.x, y-a.y);}
};
void fft(T*a,int n,int s){
    for(int i=0,j=0;i<n;i++){
        if(i>j) swap(a[i],a[j]);
        for(int l=n/2;(j^=l)<l;l>>=1);
    }
    for(int i=1;(1<<i)<=n;i++){
        int M=1<<i, K=M>>1;
        T wn= T(cos(s*2*PI/M), sin(s*2*PI/M));
        for(int j=0;j<n;j+=M){
            T w=1;
            for(int l=j;l<j+K;l++){
                T t=w*a[l+K];
                a[l+K]=a[l]-t; a[l]=a[l]+t;
                w=wn*w;
            }
        }
    }
}
void multiply(T*a,T*b,int n){
    while(n&(n-1)) n++;
    fft(a,n,1); fft(b,n,1);
    for(int i=0;i<n;i++) a[i]=a[i]*b[i];
    fft(a,n,-1);
    for(int i=0;i<n;i++) a[i]/=n;
}
int main(){
    T a[10]={T(2),T(3)}, b[10]={T(1),T(-1)};
    multiply(a,b,4);
    for(int i=0;i<10;i++) std::cout<<int(a[i].x)<<" ";
}
```
])

== Fast binomial coefficient

```cpp
int MAX_CHOOSE=3e5;
vector<ll> inv_fact(MAX_CHOOSE+5), fact(MAX_CHOOSE+5);
ll fast_nCr(ll n, ll r){
    if(n<r || r<0) return 0;
    return fact[n]*inv_fact[r]%mod*inv_fact[n-r]%mod;
}
void precalc_fact(int n){
    fact[0]=fact[1]=1;
    for(ll i=2;i<=n;i++) fact[i]=(fact[i-1]*i)%mod;
    inv_fact[0]=inv_fact[1]=1;
    for(ll i=2;i<=n;i++)
    inv_fact[i]=(mod_inv(i,mod)*inv_fact[i-1])%mod;
}
```