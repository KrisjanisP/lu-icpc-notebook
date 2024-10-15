= Data structures

#block(breakable:false,[
== Treap

```cpp
// Implicit segment tree implementation

struct Node{
    int value, cnt, priority;
    Node *left, *right;
    Node(int p) : value(p), cnt(1), priority(gen()), left(NULL), right(NULL) {};
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
    int go[K]; // if need more memory can delete this and use “next”

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
            values[x]=(lx<arr.size())? single(arr[lx]) : NEUTURAL_ELEMENT;
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