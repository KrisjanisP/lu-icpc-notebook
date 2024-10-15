= Data structures

#block(breakable:false,[
== Treap

```cpp
// Implicit segment tree implementation

struct Node{
    int value, cnt, priority;
    Node *left, *right;
    Node(int p) : value(p), cnt(1), priority(gen()),
        left(NULL), right(NULL) {};
};
 
typedef Node* pnode;
 
int get(pnode q){
    if(!q) return 0;
    return q->cnt;
}
 
void update_cnt(pnode &q){
    if(!q) return;
    q->cnt = get(q->left) + get(q->right) + 1;
}
 
void merge(pnode &T, pnode lef, pnode rig){
    if(!lef) {T=rig;return;}
    if(!rig){T=lef;return;}
    if(lef->priority > rig->priority){
        merge(lef->right, lef->right, rig);
        T = lef;
    }
    else{
        merge(rig->left, lef, rig->left);
        T = rig;
    }
    update_cnt(T);
}
 
void split(pnode cur, pnode &lef, pnode &rig, int key){
    if(!cur){
        lef = rig = NULL;
        return;
    }
    int id = get(cur->left) + 1;
    if(id <= key){
        split(cur->right, cur->right, rig, key - id);
        lef = cur;
    }
    else{
        split(cur->left, lef, cur->left, key);
        rig = cur;
    }
    update_cnt(cur);
}
```
])

#block(breakable:false,[
== Lazy segment tree

```cpp
struct SumSegmentTree{
    vector<ll> S, O, L;
    void build(ll ti, ll tl, ll tr){
        if(tl==tr){S[ti]=O[tl]; return;}
        build(ti*2, tl, (tl+tr)/2);
        build((ti*2)+1, ((tl+tr)/2)+1, tr);
        S[ti]=S[ti*2]+S[(ti*2)+1];
    }
    void push(ll ti, ll tl, ll tr){
        S[ti] += L[ti]*(tr-tl+1);
        if(tl==tr){L[ti]=0;return;}
        L[ti+ti] += L[ti],L[ti+ti+1] += L[ti];
        L[ti] = 0;
    }
    ll query(ll ti, ll tl, ll tr, ll i, ll j){
        push(ti, tl, tr);
        if(i<=tl&&tr<=j) return S[ti];
        if(tr<i||tl>j) return 0;
        ll a = query(ti*2, tl, (tl+tr)/2, i, j);
        ll b = query((ti*2)+1, ((tl+tr)/2)+1,tr, i, j);
        return a+b;
    }
    void update(ll ti, ll tl, ll tr, ll i, ll j, ll v){
        if(i<=tl&&tr<=j){L[ti]+=v;return;}
        if(tr<i||tl>j) return;
        S[ti]+=v*(i-j+1);
        update(ti*2, tl, (tl+tr)/2, i, j, v);
        update((ti*2)+1, ((tl+tr)/2)+1, tr, i, j, v);
    };
    ST(vector<ll> &V){
        O = V;
        S.resize(O.size()*4, 0);
        L.resize(O.size()*4, 0);
        build(1, 0, O.size()-1);
    }
};
```
])

#block(breakable:false,[
== Sparse table

```cpp
const int N;
const int M; //log2(N)
int sparse[N][M];

void build() {
  for(int i = 0; i < n; i++)
    sparse[i][0] = v[i];

  for(int j = 1; j < M; j++)
    for(int i = 0; i < n; i++)
      sparse[i][j] = 
        i + (1 << j - 1) < n
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
    vector<ll> bit;  // binary indexed tree
    int n;
 
    FenwickTree(int n) {
        this->n = n;
        bit.assign(n, 0);
    }
 
    ll sum(int r) {
        ll ret = 0;
        for (; r >= 0; r = (r & (r + 1)) - 1)
            ret += bit[r];
        return ret;
    }
 
    ll sum(int l, int r) { // l to r of the og array INCLUSIVE
        return sum(r) - sum(l - 1);
    }
 
    void add(int idx, ll delta) {
        for (; idx < n; idx = idx | (idx + 1))
            bit[idx] += delta;
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

struct Vertex {
    int next[K];
    bool output = false;
    int p = -1; // parent node
    char pch; // "transition" character from parent to this node
    int link = -1; // fail link
    int go[K]; // if need more memory can delete this, use “next”

    // additional potentially useful things
    int depth = -1;
    // longest string that has an output from this vertex
    int exitlen = -1; 

    Vertex(int p=-1, char ch='$') : p(p), pch(ch) {
        fill(begin(next), end(next), -1);
        fill(begin(go), end(go), -1);
    }
};
vector<Vertex> t(1);
void add_string(string const& s) {
    int v = 0;
    for (char ch : s) {
        int c = ch - 'a';
        if (t[v].next[c] == -1) {
            t[v].next[c] = t.size();
            t.emplace_back(v, ch); // !!!!! ch not c
        }
        v = t[v].next[c];
    }
    t[v].output = true;
}
int go(int v, char ch);
int get_link(int v) {
    if (t[v].link == -1) {
        if (v == 0 || t[v].p == 0)
            t[v].link = 0;
        else
            t[v].link = go(get_link(t[v].p), t[v].pch);
    }
    return t[v].link;
}
int go(int v, char ch) {
    int c = ch - 'a';
    if (t[v].go[c] == -1) {
        if (t[v].next[c] != -1)
            t[v].go[c] = t[v].next[c];
        else
            // !!!!! ch not c
            t[v].go[c] = v == 0 ? 0 : go(get_link(v), ch); 
    }
    return t[v].go[c];
}
```
])

#block( breakable: false,[
```cpp

// int go(int v, char ch) { // go without the go[K] variable
//     int c = ch - 'a';
//     if (t[v].next[c] == -1) {
//         // !!!!! ch not c
//         t[v].next[c] = v == 0 ? 0 : go(get_link(v), ch); 
//     }
//     return t[v].next[c];
// }

// helper function
int get_depth(int v){
    if (t[v].depth == -1){
        if (v == 0) {
            t[v].depth = 0;
        } else {
            t[v].depth = get_depth(t[v].p)+1;
        }
    }
    return t[v].depth;
}
// helper function
int get_exitlen(int v){
    if (t[v].exitlen == -1){
        if (v == 0){
            t[v].exitlen = 0;
        } else if (t[v].output) {
            t[v].exitlen = get_depth(v);
        } else {
            t[v].exitlen = get_exitlen(get_link(v));
        }
    }
    return t[v].exitlen;
}
```
])

#block( breakable: false,[
== Disjoint Set Union

```cpp
struct DSU {
    vector<int> parent, rank;
    DSU(int n) {
        parent.resize(n); rank.resize(n);
        for (int i = 0; i < n; i++)
            parent[i] = i;
    }
    int root(int a) {
        if (parent[a] == a) return a;
        return parent[a] = find(parent[a]);
    }
    void unite(int a, int b) {
        a = find(a), b = find(b);
        if (a == b) return;
        if (rank[a] < rank[b]) {
            parent[a] = b;
        } else if (rank[a] > rank[b]) {
            parent[b] = a;
        } else {
            parent[b] = a;
            rank[a] = rank[a] + 1;
        }
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
== janY normal segment tree

```cpp
struct item{ll sum;};
struct segtree{
    int size; vector<item> values;
    item merge(item a,item b){return {a.sum+b.sum};}
    item NEUTURAL_ELEMENT={0};
    item single(ll v){return {v};}
    void init(int n){
        size=1; while(size<n) size*=2;
        values.resize(size*2,NEUTURAL_ELEMENT);
    }
    void build(vl &arr,int x,int lx,int rx){
        if(rx-lx==1){
            if(lx<arr.size()) values[x]=single(arr[lx]);
            else values[x]=NEUTURAL_ELEMENT;
            return;
        }
        int m=(lx+rx)/2;
        build(arr,2*x+1,lx,m); build(arr,2*x+2,m,rx);
        values[x]=merge(values[2*x+1],values[2*x+2]);
    }
    void build(vl &arr){build(arr,0,0,size);}
    void set(int i,ll v,int x,int lx,int rx){
        if(rx-lx==1){values[x]=single(v); return;}
        int m=(lx+rx)/2;
        if(i<m) set(i,v,2*x+1,lx,m);
        else set(i,v,2*x+2,m,rx);
        values[x]=merge(values[2*x+1],values[2*x+2]);
    }
    void set(int i,ll v){set(i,v,0,0,size);}
    item calc(int l,int r,int x,int lx,int rx){
        if(lx>=r || rx<=l) return NEUTURAL_ELEMENT;
        if(lx>=l && rx<=r) return values[x];
        int m=(lx+rx)/2;
        return merge(calc(l,r,2*x+1,lx,m),calc(l,r,2*x+2,m,rx));
    }
    item calc(int l,int r){return calc(l,r,0,0,size);}
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
== janY fenwick tree

```cpp
struct fenwick { // point update (delta), range sum
    ll *bit; int fsize;
    void init(int n){
        fsize=n; bit=new ll[n+1]();
    }
    int lsb(int x){ return x & -x; }
    ll query(int v){
        ll s=0; while(v>0){ s += bit[v]; v -= lsb(v); } return s;
    }
    void add(int v, int delta){
        v++; while(v <= fsize){ bit[v] += delta; v += lsb(v); }
    }
    void build(vector<ll> &inp){
        for(int i=1;i<=inp.size();i++) bit[i]=inp[i-1];
        for(int i=1;i<=inp.size();i++){
            int p=i+lsb(i); if(p<=fsize) bit[p]+=bit[i];
        }
    }
    ll calc(int l, int r){
        return query(r+1) - query(l);
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
== janY Aho-Corasick algorithm

```cpp
struct Vertex {
    int next[K], go[K], p=-1; bool output=false;
    char pch; int link=-1, depth=-1, exitlen=-1;
    Vertex(int parent=-1, char ch='$') : p(parent), pch(ch) {
        fill(next, next+K, -1); fill(go, go+K, -1);
    }
};
vector<Vertex> t(1);
void add_string(const string &s){
    int v=0; for(char ch:s){
        int c=ch-'a'; if(t[v].next[c]==-1){
            t[v].next[c]=t.size(); t.emplace_back(v,ch);
        }
        v = t[v].next[c];
    }
    t[v].output=true;
}
int go_func(int v, char ch);
int get_link(int v){
    if(t[v].link==-1){
        if(v==0 || t[v].p==0) t[v].link=0;
        else t[v].link = go_func(get_link(t[v].p), t[v].pch);
    }
    return t[v].link;
}
int go_func(int v, char ch){
    int c=ch-'a';
    if(t[v].go[c]==-1){
        if(t[v].next[c]!=-1) t[v].go[c]=t[v].next[c];
        else t[v].go[c] = (v==0) ? 0 : go_func(get_link(v), ch);
    }
    return t[v].go[c];
}
int get_depth(int v){
    if(t[v].depth==-1){
        t[v].depth = (v==0) ? 0 : get_depth(t[v].p)+1;
    }
    return t[v].depth;
}
int get_exitlen(int v){
    if(t[v].exitlen==-1){
        if(v==0) t[v].exitlen=0;
        else if(t[v].output) t[v].exitlen = get_depth(v);
        else t[v].exitlen = get_exitlen(get_link(v));
    }
    return t[v].exitlen;
}
```
])