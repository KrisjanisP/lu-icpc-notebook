= Data structures

== Policy-based data structures (PBDS)

// from blog https://codeforces.com/blog/entry/60737

preamble
```cpp
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace std;
```

5x faster hash map
```cpp
gp_hash_table<int, int> table;
```

defeating anti-hash tests
```cpp
const int RANDOM = chrono::high_resolution_clock::now().time_since_epoch().count();
struct chash {
    int operator()(int x) const { return x ^ RANDOM; }
};
gp_hash_table<key, int, chash> table;
```

ordered set with support for `find_by_order()` and `order_of_key()`

```cpp
tree<int,null_type,less<int>,
rb_tree_tag,tree_order_statistics_node_update> t;
```

#block(breakable:false,[
== Treap

```cpp
struct Node{
    int value, cnt, pri; Node *left, *right;
    Node(int p) : value(p), cnt(1), pri(gen()),
        left(NULL), right(NULL) {};
};
typedef Node* pnode;
int get(pnode q){if(!q) return 0; return q->cnt;}
void update_cnt(pnode &q){
    if(!q) return; q->cnt=get(q->left)+get(q->right)+1;
}
void merge(pnode &T, pnode lef, pnode rig){
    if(!lef){T=rig;return;} if(!rig){T=lef;return;}
    if(lef->pri>rig->pri){merge(lef->right,lef->right,rig);T=lef;
    }else{merge(rig->left, lef, rig->left); T = rig;}
    update_cnt(T);
}
void split(pnode cur, pnode &lef, pnode &rig, int key){
    if(!cur){lef=rig=NULL;return;} int id=get(cur->left)+1;
    if(id<=key){split(cur->right,cur->right,rig,key-id);lef=cur;}
    else {split(cur->left, lef, cur->left, key); rig = cur;}
    update_cnt(cur);
}
```
])

#block(breakable:false,[
== Lazy segment tree

```cpp
struct SumSegmentTree{
    vector<ll> S, O, L; // S: segment tree, O: original, L: lazy
    void build(ll ti, ll tl, ll tr){
        if(tl==tr){S[ti]=O[tl]; return;}
        build(ti*2,tl,(tl+tr)/2);build((ti*2)+1,(tl+tr)/2+1,tr);
        S[ti]=S[ti*2]+S[(ti*2)+1];
    }
    void push(ll ti, ll tl, ll tr){
        S[ti] += L[ti]*(tr-tl+1); if(tl==tr){L[ti]=0;return;}
        L[ti+ti] += L[ti],L[ti+ti+1] += L[ti]; L[ti] = 0;
    }
    ll query(ll ti, ll tl, ll tr, ll i, ll j){
        push(ti, tl, tr);
        if(i<=tl&&tr<=j) return S[ti]; if(tr<i||tl>j) return 0;
        ll a = query(ti*2, tl, (tl+tr)/2, i, j);
        ll b = query((ti*2)+1, ((tl+tr)/2)+1,tr, i, j);
        return a+b;
    }
    void update(ll ti, ll tl, ll tr, ll i, ll j, ll v){
        if(i<=tl&&tr<=j){L[ti]+=v;return;}
        if(tr<i||tl>j) return; S[ti]+=v*(i-j+1);
        update(ti*2,tl,(tl+tr)/2,i,j,v);
        update((ti*2)+1,(tl+tr)/2+1,tr,i,j,v);
    };
    ST(vector<ll> &V){
        O = V; S.resize(O.size()*4, 0); L.resize(O.size()*4, 0);
        build(1, 0, O.size()-1);
    }
};
```
])

#block(breakable:false,[
== Sparse table

```cpp
const int N, M; //M=log2(N)
int sparse[N][M];
void build() {
  for(int i = 0; i < n; i++) sparse[i][0] = v[i];
  for(int j = 1; j < M; j++) for(int i = 0; i < n; i++)
    sparse[i][j] = i + (1 << j - 1) < n
      ? min(sparse[i][j - 1], sparse[i + (1 << j - 1)][j - 1]) 
      : sparse[i][j - 1];
}
int query(int a, int b){
  int pot = 32 - __builtin_clz(b - a) - 1;
  return min(sparse[a][pot], sparse[b - (1 << pot) + 1][pot]);
}
```
])

#block( breakable: false,[
== Fenwick tree

```cpp
struct FenwickTree {
    int n;vector<ll> bit; // binary indexed tree
    FenwickTree(int n) {this->n=n;bit.assign(n, 0);}
    ll sum(int r) {
        ll ret=0;
        for(;r>=0;r=(r&(r+1))-1) ret+=bit[r];
        return ret;
    }
    ll sum(int l, int r){return sum(r)-sum(l-1);}
    void add(int idx, ll delta){
        for(;idx<n;idx=idx|(idx+1))bit[idx]+=delta;
    }
};
```
])


#block( breakable: false,[

== Trie

```cpp
const int K = 26;

struct Vertex {
    int next[K];
    bool output = false;
    Vertex() {fill(begin(next), end(next), -1);}
};

vector<Vertex> t(1); // trie nodes

void add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if (t[v].next[c] == -1) {
            t[v].next[c] = t.size();
            t.emplace_back(); 
        }
        v = t[v].next[c];
    }
    t[v].output = true;
}
```
])

#block( breakable: false,[

== Aho-Corasick

```cpp
const int K = 26; 
struct V {
    int n[K], go[K], p = -1; // next, go transitions, parent
    char ch;                 // char from parent to this node
    bool out = false;        // is end of a pattern
    int l = -1, d = -1, e = -1; // fail link, depth, exit length
    V(int parent = -1, char c = '$') : p(parent), ch(c) {
        fill(n, n + K, -1);   // initialize transitions
        fill(go, go + K, -1); // initialize go transitions
    }
};
vector<V> t(1);
void add_string(const string& s){ // Add a string to the trie
    int v = 0;
    for(char c : s){
        int ci = c - 'a';
        if(t[v].n[ci] == -1){
            t[v].n[ci] = t.size();
            t.emplace_back(v, c); // create new node
        }
        v = t[v].n[ci];
    }
    t[v].out = true; // mark end of pattern
}
int go_func(int v, char c);
int get_link(int v){ // Get the fail link for node v
    if(t[v].l == -1){
        if(v == 0 || t[v].p == 0) t[v].l = 0;
        else t[v].l = go_func(get_link(t[v].p), t[v].ch);
    } return t[v].l;
}
// Compute the transition for node v with character c
int go_func(int v, char c){
    int ci = c - 'a';
    if(t[v].go[ci] == -1){
        if(t[v].n[ci] != -1) t[v].go[ci] = t[v].n[ci];
        else t[v].go[ci] = (v==0) ? 0 : go_func(get_link(v),c);
    }
    return t[v].go[ci];
}
```

])

#block( breakable: false,[
== Disjoint Set Union

```cpp
struct DSU {
    vector<int> p, r; // p: parent, r: rank
    DSU(int n) {
        p.resize(n); r.resize(n);
        for (int i = 0; i < n; i++) p[i] = i;
    }
    int f(int a){if (p[a] == a) return a; return p[a] = f(p[a]);}
    void unite(int a, int b) {
        a = f(a), b = f(b); if (a == b) return;
        if (r[a] < r[b]) p[a] = b; else if (r[a] > r[b]) p[b] = a;
        else {p[b] = a; r[a]++;}
    }
};
```
])

#block( breakable: false,[

== Merge sort tree

```cpp
struct MergeSortTree{
    int size; vector<vector<ll>> values;
    void init(int n){
        size=1; while(size<n) size*=2;
        values.resize(size*2, vector<ll>());
    }
    void build(vector<ll> &arr, int x, int lx, int rx){
        if(rx-lx==1){
            if(lx<arr.size()) values[x].push_back(arr[lx]);
            else values[x].push_back(-1);
            return;
        }
        int m=(lx+rx)/2;
        build(arr,2*x+1,lx,m);
        build(arr,2*x+2,m,rx);
        int i=0, j=0, asize=values[2*x+1].size(); 
        while(i<asize && j<values[2*x+2].size()){
            if(values[2*x+1][i]<values[2*x+2][j])
              values[x].push_back(values[2*x+1][i++]);
            else values[x].push_back(values[2*x+2][j++]);
        }
        while(i<asize)
          values[x].push_back(values[2*x+1][i++]);
        while(j<values[2*x+2].size())
          values[x].push_back(values[2*x+2][j++]);
    }
    void build(vector<ll> &arr){ build(arr,0,0,size); }
    int calc(int l, int r, int x, int lx, int rx, int k){
        if(lx>=r || rx<=l) return 0;
        if(lx>=l && rx<=r){
            int lft=-1, rght=values[x].size();
            while(rght-lft>1){
                int mid=(lft+rght)/2;
                if(values[x][mid]<k) lft=mid;
                else rght=mid;
            }
            return lft+1;
        }
        int m=(lx+rx)/2;
        return calc(l,r,2*x+1,lx,m,k) + calc(l,r,2*x+2,m,rx,k);
    }
    int calc(int l, int r, int k){ return calc(l,r,0,0,size,k); }
};
```
])

#block(breakable: false, [
== janY mass operations segment tree
```cpp
struct item { ll x; item(ll x=0) : x(x) {} };
struct segtree {
    int size; vector<item> values, ops;
    item NEUTRAL=0, DEFAULT=0, NOOP=0;
    item modify_op(item a, item b, ll len) {
        a.x += b.x*len; return a; }
    void apply_mod_op(item &a, item b, ll len) {
        a = modify_op(a, b, len); }
    item calc_op(item a, item b) { return item(a.x + b.x); }
    void init(int n) {
        size=1; while(size<n) size<<=1;
        values.assign(size<<1, DEFAULT);
        ops.assign(size<<1, NOOP);
    }
    void build(vector<item> &arr, int x=0, int lx=0, int rx=-1) {
        if(rx==-1) rx = size;
        if(rx - lx ==1) {
            values[x] =
                (lx < arr.size()) ? arr[lx] : NEUTRAL; return;
        }
        int m=(lx+rx)/2;
        build(arr,2*x+1,lx,m);
        build(arr,2*x+2,m,rx);
        values[x] = calc_op(values[2*x+1], values[2*x+2]);
    }
    void propagate(int x, int lx, int rx) {
        if(rx - lx ==1) return;
        int m=(lx+rx)/2;
        apply_mod_op(ops[2*x+1], ops[x],1);
        apply_mod_op(values[2*x+1], ops[x],m-lx);
        apply_mod_op(ops[2*x+2], ops[x],1);
        apply_mod_op(values[2*x+2], ops[x],rx-m);
        ops[x] = NOOP;
    }
    void set(int l, int r, ll v, int x=0, int lx=0, int rx=-1) {
        if(rx==-1) rx = size; propagate(x, lx, rx);
        if(lx >= r || rx <= l) return;
        if(lx >= l && rx <= r) {
            apply_mod_op(ops[x], item(v),1);
            apply_mod_op(values[x], item(v), rx-lx); return;
        }
        int m=(lx+rx)/2;
        set(l,r,v,2*x+1,lx,m); set(l,r,v,2*x+2,m,rx);
        values[x] = calc_op(values[2*x+1], values[2*x+2]);
    }
    item calc(int l, int r, int x=0, int lx=0, int rx=-1) {
        if(rx==-1) rx = size; propagate(x, lx, rx);
        if(lx >= r || rx <= l) return NEUTRAL;
        if(lx >= l && rx <= r) return values[x];
        int m=(lx+rx)/2;
        return
            calc_op(calc(l,r,2*x+1,lx,m), calc(l,r,2*x+2,m,rx));
    }
};
```
])

#block( breakable: false,[
== janY fenwick tree range update

```cpp
struct fenwick { // range update
    ll *bit1, *bit2; int fsize;
    void init(int n){
        fsize=n; bit1=new ll[n+1](); bit2=new ll[n+1]();
    }
    ll getSum(ll BIT[], int i){
        ll s=0; i++; while(i>0){ s += BIT[i]; i -= i & -i; }
        return s;
    }
    void updateBIT(ll BIT[], int i, ll v){
        i++; while(i <= fsize){ BIT[i] += v; i += i & -i; }
    }
    ll sum(int x){
        return getSum(bit1,x)*x - getSum(bit2,x);
    }
    void add(int l, int r, ll v){
        updateBIT(bit1,l,v); updateBIT(bit1,r+1,-v);
        updateBIT(bit2,l,v*(l-1)); updateBIT(bit2,r+1,-v*r);
    }
    ll calc(int l, int r){
        return sum(r) - sum(l-1);
    }
};
```
])

#block( breakable: false,[
== Persistent segment tree

```cpp
#define V struct Vertex
struct Vertex { V *l, *r; ll sum; 
    Vertex(ll val){l=r=nullptr; sum=val;} 
    Vertex(V* le, V* ri){l=le;r=ri;sum=(l?l->sum:0)+(r?r->sum:0);} 
}; 
int siz; 
vector<V*> start_nodes; 
V* build(int lx, int rx, vl &a){
    if (lx == rx-1) return new V(a[lx]);
    return new V(build(lx,(lx+rx)/2,a),build((lx+rx)/2,rx,a));
} 
V* build(vl &a){ siz = a.size(); return build(0, siz, a); } 
ll calc(V* v, int lx, int rx, int l, int r){
    if(lx >= r || rx <= l) return 0;
    if(lx >= l && rx <= r) return v->sum;
    int m = (lx + rx) / 2;
    return calc(v->l, lx, m, l, r) + calc(v->r, m, rx, l, r);
} 
ll calc(V* v, int l, int r){
    if (l>r) return 0;
    return calc(v,0,siz,l,r);
} 
V* upd(V* v, int lx, int rx, int i, ll val){
    if(lx == rx-1) return new V(val);
    int m = (lx + rx) / 2;
    if (i < m) return new V(upd(v->l, lx, m, i, val), v->r);
    else return new V(v->l, upd(v->r, m, rx, i, val));
} 
V* upd(V* v, int i, ll val){ return upd(v, 0, siz, i, val); }
```
])
