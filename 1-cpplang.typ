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

== Some abbreviations

```cpp
#define all(x) x.begin(), x.end()
#define fi first
#define se second
#define pb push_back
using ll = long long;
using vi = vector<ll>;
using vvi = vector<vi>;
using pii = pair<ll, ll>;
using vpii = vector<pii>;
```

== rand() is bad!
// https://codeforces.com/blog/entry/61587

rand() is bad. it may return up to RAND_MAX
which by default may be only 32767

```cpp
int main() {
    mt19937 rng(chrono::steady_clock::now()
        .time_since_epoch().count());
    vector<int> permutation(N);
    for (int i = 0; i < N; i++) permutation[i] = i;
    shuffle(permutation.begin(), permutation.end(), rng);
}
```

for 64-bit, use `mt19937_64`
