
#block(breakable: false,[
= String Processing
== Knuth-Morris-Pratt (KMP)

```cpp
char s[N], p[N]; int b[N], n, m; // n = strlen(s), m = strlen(p);
void kmppre() {
  b[0] = -1; for (int i = 0, j = -1; i < m; b[++i] = ++j)
    while (j >= 0 and p[i] != p[j]) j = b[j];
}
void kmp() {
  for (int i = 0, j = 0; i < n;) {
    while (j >= 0 and s[i] != p[j]) j=b[j];
    i++, j++; if (j == m) {j = b[j];}
  }
}
```
])

#block( breakable: false,[
== Suffix Array

```cpp
vector<int> suffix_array(string &s){
    int n = s.size(), alph = 256; 
    vector<int> cnt(max(n, alph)), p(n), c(n);
    rep(i,0,n) cnt[s[i]]++;
    rep(i,1,alph) cnt[i]+=cnt[i-1];
    rep(i,0,n) p[--cnt[s[i]]] = i;
    c[p[0]]=0; rep(i,1,n)
    c[p[i]] = (s[p[i]] != s[p[i-1]]) ? c[p[i-1]] + 1 : c[p[i-1]];
    vector<int> p2(n), c2(n);
    for(int k=0; (1<<k)<n; k++){
        rep(i,0,n) p2[i]=(p[i]-(1<<k)+n)%n;
        fill(cnt.begin(), cnt.begin()+c[p[n-1]]+1, 0);
        rep(i,0,n) cnt[c[p2[i]]]++;
        rep(i,1,c[p[n-1]]+1) cnt[i]+=cnt[i-1];
        for(int i=n-1;i>=0;i--) p[--cnt[c[p2[i]]]] = p2[i];
        c2[p[0]]=0; rep(i,1,n) {
            pair<int,int> a1 = {c[p[i]], c[(p[i]+(1<<k))%n]};
            pair<int,int> a2 = {c[p[i-1]], c[(p[i-1]+(1<<k))%n]};
            c2[p[i]] = (a1 != a2) ? c2[p[i-1]] + 1 : c2[p[i-1]];
        }
        c.swap(c2);
    }
    return p;
}
```
])

#block( breakable: false,[
== Longest common prefix (LCP) with SA
```cpp
vector<int> lcp(string &s, vector<int> &p){
    int n = s.size(); vector<int> pi(n), ans(n-1); 
    rep(i,0,n) pi[p[i]] = i;
    int lst = 0;
    rep(i,0,n-1){
        if(pi[i] == n-1){ lst = 0; continue; }
        int j = p[pi[i]+1];
        while(i+lst<n && j+lst<n && s[i+lst] == s[j+lst]) lst++;
        ans[pi[i]] = lst; lst = max(lst-1, 0);
    }
    return ans;
}
```
])

#block( breakable: false,[
== Rabin-Karp pattern match with hashing

```cpp
const int B = 31;
const int MOD = 1e9+7, B = 31;
void rabin(string s, string p){
    int n = s.size(), m = p.size(); if(n<m) return;
    vector<ull> power(max(n, m), 1);
    rep(i,1,power.size()) power[i] = (power[i-1]*B)%MOD;
    ull hp=0, hs=0;
    rep(i,0,m){ hp=(hp*B + p[i])%MOD; hs=(hs*B + s[i])%MOD; }
    if(hs == hp) { /* match at 0 */ }
    rep(i,m,n){
        hs = (hs*B + s[i])%MOD;
        hs = (hs + MOD - (s[i-m]*power[m])%MOD)%MOD;
        if(hs == hp) { /* match at i-m+1 */ }
    }
}
```
])

#block( breakable: false,[
  
== Z-function

The Z-function of a string $s$ is an array $z$ where $z_i$ is the length of the longest substring starting from $s_i$ which is also a prefix of $s$.

Examples:
- "aaaaa": $[0, 4, 3, 2, 1]$
- "aaabaab": $[0,2,1,0,2,1,0]$
- "abacaba": $[0,0,1,0,3,0,1]$

```cpp
vector<int> zfunction(const string& s){
  vector<int> z (s.size());
  for (int i = 1, l = 0, r = 0, n = s.size(); i < n; i++){
    if (i <= r) z[i] = min(z[i-l], r - i + 1);
    while (i + z[i] < n and s[z[i]] == s[z[i] + i]) z[i]++;
    if (i + z[i] - 1 > r) l = i, r = i + z[i] - 1;
  }
  return z;
}
```
])

#block( breakable: false,[
== Manacher's longest palindromic substring

```cpp
int manacher(string s){
    int n = s.size(); string p = "^#";
    rep(i,0,n) p += string(1, s[i]) + "#";
    p += "$"; n = p.size(); vector<int> lps(n, 0);
    int C=0, R=0, m=0;
    rep(i,1,n-1){
        int mirr = 2*C - i;
        if(i < R) lps[i] = min(R-i, lps[mirr]);
        while(p[i + 1 + lps[i]] == p[i - 1 - lps[i]]) lps[i]++;
        if(i + lps[i] > R){ C = i; R = i + lps[i]; }
        m = max(m, lps[i]);
    }
    return m;
}
```
])