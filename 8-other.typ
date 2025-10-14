
= I'm running out of time
== Simulated annealing
```cpp

const ld T = (ld)2000;
const ld alpha = 0.999999;
// (new_score - old_score) / (temperature_final) ~ 10 works well

const ld L = (ld)1e6;
ld small_rand(){
	return ((ld)gen(L))/L;
}

ld P(ld old, ld nw, ld temp){
	if(nw > old)
		return 1.0;
	return exp((nw-old)/temp);
}

{
  auto start = chrono::steady_clock::now();
  ld time_limit = 2000; 
  ld temperature = T;
  ld max_score = -1;

  while(elapsed_time < time_limit){
    auto cur = chrono::steady_clock::now();
    elapsed_time = chrono::duration_cast<chrono::milliseconds>(cur - start).count();
    temperature *= alpha;
    
    // try a neighboring state
    // ....
    // ....

    old_score = score(old_state);  
    new_score = score(new_state);
    if(P(old_score, new_score, temperature) >= small_rand()){
      old_state = new_state;
      old_score = new_score;
    }
    if(old_score > max_score){
      max_score = old_score;
      max_state = old_state;
    }
  }
}
```

#columns(2)[
#block( breakable: false,[
== Eulerian Path

#image("./assets/eulerian-path.png", width: 100%)
])

#block(breakable: false)[
== Flows with demands

#image("./assets/flows-with-demands.png", width: 100%)
]

#block( breakable: false,[
== Point in convex polygon $O(log n)$

#image("./assets/point-in-convex-polygon.png", width: 100%)

```cpp
bool pointInTriangle(pt a, pt b, pt c, pt point) {
    long long s1 = abs(a.cross(b, c));
    long long s2 = abs(point.cross(a, b)) + abs(point.cross(b, c)) + abs(point.cross(c, a));
    return s1 == s2;
}

void prepare(vector<pt> &points) {
    n = points.size();
    int pos = 0;
    for (int i = 1; i < n; i++) {
        if (lexComp(points[i], points[pos]))
            pos = i;
    }
    rotate(points.begin(), points.begin() + pos, points.end());

    n--;
    seq.resize(n);
    for (int i = 0; i < n; i++)
        seq[i] = points[i + 1] - points[0];
    translation = points[0];
}
```
])
#block(breakable: false,[
```cpp
bool pointInConvexPolygon(pt point) {
    point = point - translation;
    if (seq[0].cross(point) != 0 &&
            sgn(seq[0].cross(point)) != sgn(seq[0].cross(seq[n - 1])))
        return false;
    if (seq[n - 1].cross(point) != 0 &&
            sgn(seq[n - 1].cross(point)) != sgn(seq[n - 1].cross(seq[0])))
        return false;

    if (seq[0].cross(point) == 0)
        return seq[0].sqrLen() >= point.sqrLen();

    int l = 0, r = n - 1;
    while (r - l > 1) {
        int mid = (l + r) / 2;
        int pos = mid;
        if (seq[pos].cross(point) >= 0)
            l = mid;
        else
            r = mid;
    }
    int pos = l;
    return pointInTriangle(seq[pos], seq[pos + 1], pt(0, 0), point);
}
```
])

#block( breakable: false,[
== Minkowski sum of convex polygons

#image("./assets/minkowski-sum-1.png", width: 100%)

])

#block( breakable: false,[
#image("./assets/minkowski-sum-2.png", width: 100%)

```cpp
void reorder_polygon(vector<pt> & P){
    size_t pos = 0;
    for(size_t i = 1; i < P.size(); i++){
        if(P[i].y < P[pos].y || (P[i].y == P[pos].y && P[i].x < P[pos].x))
            pos = i;
    }
    rotate(P.begin(), P.begin() + pos, P.end());
}

vector<pt> minkowski(vector<pt> P, vector<pt> Q){
    // the first vertex must be the lowest
    reorder_polygon(P); reorder_polygon(Q);
    // we must ensure cyclic indexing
    P.push_back(P[0]); P.push_back(P[1]);
    Q.push_back(Q[0]); Q.push_back(Q[1]);
    vector<pt> result; size_t i = 0, j = 0;
    while(i < P.size() - 2 || j < Q.size() - 2){
        result.push_back(P[i] + Q[j]);
        auto cross = (P[i + 1] - P[i]).cross(Q[j + 1] - Q[j]);
        if(cross >= 0 && i < P.size() - 2) ++i;
        if(cross <= 0 && j < Q.size() - 2) ++j;
    }
    return result
}
```
])

#block(breakable: false,[
== janY template

```cpp
#include<bits/stdc++.h>
using namespace std;
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,popcnt,lzcnt")
#define fo(i,n) for(i=0;i<n;i++)
#define Fo(i,k,n) for(i=k;k<n?i<n:i>n;k<n?i+=1:i-=1)
#define ll long long
#define ld long double
#define all(x) x.begin(),x.end()
#define sortall(x) sort(all(x))
#define rev(x) reverse(x.begin(),x.end())
#define fi first
#define se second
#define pb push_back
#define PI 3.14159265359
typedef pair<int,int> pii;
typedef pair<ll,ll> pl;
typedef vector<int> vi;
typedef vector<ll> vl;
typedef vector<pii> vpii;
typedef vector<pl> vpl;
typedef vector<vi> vvi;
typedef vector<vl> vvl;
bool sortbysec(const pair<int,int> &a,const pair<int,int> &b){return a.second<b.second;}
#define sortpairbysec(x) sort(all(x), sortbysec)
bool sortcond(const pair<int,int> &a,const pair<int,int> &b){
    if(a.fi!=b.fi) return a.fi<b.fi;
    return a.se>b.se;
}
struct myComp {
    constexpr bool operator()(pii const& a, pii const& b) const noexcept{
        if(a.first!=b.first) return a.first<b.first;
        return a.second>b.second;
    }
};
const int mod=1000000007, N=3e5, M=N;
// & - AND; | - OR; ^ - XOR
vl a;
ll n,m,k,q;
void solve(int tc){
    int i,j;
    cin>>n;
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    int t=1; cin>>t; int i;
    fo(i,t) solve(i+1);
    return 0;
}
```
])

]