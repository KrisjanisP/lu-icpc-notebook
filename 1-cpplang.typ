= C++ programming language

== I/O disable sync

```cpp
ios_base::sync_with_stdio(false);
cin.tie(NULL); cout.tie(NULL);
```

== Pragmas

```cpp
// change to O3 to disable fast-math for geometry problems
#pragma GCC optimize("Ofast,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt,tune=native")
```

== Printing structs

```cpp
ostream& operator<<(ostream& os, const pair<int, int>& p) {
    return os << "(" << p.first << ", " << p.second << ")";
}
```

#block(breakable: false,[
== Lambda for sorting

```cpp
using ii = pair<int,int>;
vector<ii> fracs = {{1, 2}, {3, 4}, {1, 3}};
sort(fracs.begin(), fracs.end(),
    [](const ii& a, const ii& b) {
    return a.fi*b.se < b.fi*a.se;
});
```

])

== Abbreviations

```cpp
#define all(x) x.begin(), x.end()
#define fi first
#define se second
#define rep(i,n) for(int i=0;i<n;i++)
#define pb push_back
using namespace std;
using ll = long long;
using vi = vector<int>;
using vvi = vector<vi>;
using pii = pair<int, int>;
using vpii = vector<pii>;
```