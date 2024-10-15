= C++ programming language

== Input/Output disable sync

```cpp
ios_base::sync_with_stdio(false);
cin.tie(NULL); cout.tie(NULL);
```

== Optimization pragmas

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

== Lambda func for sorting

```cpp
using ii = pair<int,int>;
vector<ii> fracs = {{1, 2}, {3, 4}, {1, 3}};
// sort positive rational numbers
sort(fracs.begin(), fracs.end(),
    [](const ii& a, const ii& b) {
    return a.fi*b.se < b.fi*a.se;
});
```