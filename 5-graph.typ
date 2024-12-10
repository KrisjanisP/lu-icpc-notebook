
= Graph algorithms
== Bellman-Ford

```cpp
void solve()
{
    vector<int> d(n, INF);
    d[v] = 0;
    for (;;) {
        bool any = false;

        for (Edge e : edges)
            if (d[e.a] < INF)
                if (d[e.b] > d[e.a] + e.cost) {
                    d[e.b] = d[e.a] + e.cost;
                    any = true;
                }

        if (!any)
            break;
    }
    // display d, for example, on the screen
}
```

== Dijkstra
```cpp
vector<int> adj[N], adjw[N];
int dist[N];

memset(dist, 63, sizeof(dist));
priority_queue<pii> pq;
pq.push(mp(0,0));

while (!pq.empty()) {
  int u = pq.top().nd;
  int d = -pq.top().st;
  pq.pop();

  if (d > dist[u]) continue;
  for (int i = 0; i < adj[u].size(); ++i) {
    int v = adj[u][i];
    int w = adjw[u][i];
    if (dist[u] + w < dist[v])
      dist[v] = dist[u]+w, pq.push(mp(-dist[v], v));
  }
}
```

== Floyd-Warshall 

```cpp
int adj[N][N]; // no-edge = INF

for (int k = 0; k < n; ++k)
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      adj[i][j] = min(adj[i][j], adj[i][k]+adj[k][j]);
```

== Bridges & articulations

```cpp
// Articulation points and Bridges O(V+E)
int par[N], art[N], low[N], num[N], ch[N], cnt;

void articulation(int u) {
  low[u] = num[u] = ++cnt;
  for (int v : adj[u]) {
    if (!num[v]) {
      par[v] = u; ch[u]++;
      articulation(v);
      if (low[v] >= num[u]) art[u] = 1;
      if (low[v] >  num[u]) { /* u-v bridge */ }
      low[u] = min(low[u], low[v]);
    }
    else if (v != par[u]) low[u] = min(low[u], num[v]);
  }
}

for (int i = 0; i < n; ++i) if (!num[i])
  articulation(i), art[i] = ch[i]>1;
```

== Dinic's max flow / matching

Time complexity:
- generally: $O(E V^2)$
- small flow: $O(F (V + E))$
- bipartite graph or unit flow: $O(E sqrt(V))$
Usage:
- dinic()
- add_edge(from, to, capacity)   
- recover() (optional) 

#block(breakable: false, [
```cpp
const ll N=1e5+5, INF=1e9;
struct edge{ll v, c, f;};

ll src=0, snk=N-1, h[N], ptr[N];
vector<edge> edgs;

vector<ll> g[N];

void add_edge(ll u, ll v, ll c) {
    edgs.push_back({v,c,0}), edgs.push_back({u,0,0});
    ll k=edgs.size();
    g[u].push_back(k), g[v].push_back(k+1);
}

bool bfs() {
    memset(h, 0, sizeof(h));
    queue<ll> q;
    h[src]=1;
    q.push(src);
    while(!q.empty()){
        ll u=q.front();q.pop();
        for(ll i:g[u]){
            ll v=edgs[i].v;
            if(!h[v]&&edgs[i].f<edgs[i].c)
                q.push(v),h[v]=h[u]+1;
        }
    }
    return h[snk];
}

ll dfs(ll u, ll flow){
    if(!flow or u==snk) return flow;
    for(ll &i=ptr[u];i<g[u].size();i++){
        edge &dir=edgs[g[u][i]],&rev=edgs[g[u][i]^1];
        if(h[dir.v]!=h[u]+1) continue;
        ll inc=min(flow,dir.c-dir.f);
        inc=dfs(dir.v,inc);
        if(inc){ dir.f+=inc,rev.f-=inc; return inc;}
    }
    return 0;
}

ll dinic(){
    ll flow=0;
    while(bfs()){
        memset(ptr,0,sizeof(ptr));
        while(ll inc=dfs(src,INF)) flow += inc;
    }
    return flow;
}

vector<pair<ii,ll>> recover() {
    vector<pair<ii,ll>> res;
    for(ll i=0;i<edgs.size();i+=2){
        if(edgs[i].f>0){
            ll v=edgs[i].v, u=edgs[i^1].v;
            res.push_back({{u,v},edgs[i].f});
        }
    }
    return res;
}
```
])

== Flow with demands

Finding an arbitrary flow
- Assume a network with $[L;R]$ on edges (some may have $L = 0$), let's call it old network.
- Create a New Source and New Sink (this will be the src and snk for Dinic).
- Modelling network:
  + Every edge from the old network will have cost $R - L$
  + Add an edge from New Source to every vertex $v$ with cost:
    - $S(L)$ for every $(u, v)$. (sum all $L$ that LEAVES $v$)
  + Add an edge from every vertex $v$ to New Sink with cost:
    - $S(L)$ for every $(v, w)$. (sum all $L$ that ARRIVES $v$)
  + Add an edge from Old Source to Old Sink with cost INF (circulation problem)
- The Network will be valid if and only if the flow saturates the network (max flow == $S(L)$)

Finding Min Flow
- To find min flow that satisfies just do a binary search in the (Old Sink -> Old Source) edge
- The cost of this edge represents all the flow from old network
- Min flow = $S(L)$ that arrives in Old Sink + flow that leaves (Old Sink -> Old Source)

#block( breakable: false,[
== Kosaraju's SCCs

#text(fill:blue, link("https://judge.yosupo.jp/submission/252308"))

```cpp
namespace kosaraju {
  vi adj[MAX_N], adjt[MAX_N], topo;
  bitset<MAX_N> vis;
  vvi res;
  // first dfs to get topological order
  void visit(int v) {
      vis[v]=1;
      for(int u: adj[v]) if(!vis[u]) visit(u);
      topo.pb(v);
  }
  // second dfs on transpose graph
  void assign(int v, int r) {
      vis[v]=0; res.back().pb(v);
      for(int u: adjt[v]) if(vis[u]) assign(u,r);
  }
  /// @param el 0-indexed edge list
  /// @param n graph node count
  vvi sccs(vpii el, int n){
      rep(i,MAX_N) adj[i].clear(), adjt[i].clear();
      topo.clear(); res.clear();
      for(pii e: el) adj[e.fi].pb(e.se), adjt[e.se].pb(e.fi);
      // topological order
      rep(i,n) if(!vis[i]) visit(i);
      reverse(all(topo));
      rep(i,n) { // sccs on transpose (reversed) graph
          int v=topo[i];
          if(vis[v]) res.push_back({}), assign(v,v);
      }
      return res;
  }
  void test() {
      vpii edges={{1,4},{5,2},{3,0},{5,5},{4,1},{0,3},{4,2}};
      vvi expected={{5},{1,4},{2},{0,3}};
      assert(kosaraju::sccs(edges, 6)==expected);
  }
}
int main() {
    kosaraju::test();

    ios::sync_with_stdio(false);
    int n, m; cin>>n>>m;
    vpii el;
    for(int i=0;i<m;i++){
        int a, b; cin>>a>>b;
        el.push_back({a,b});
    }
    vvi sccs=kosaraju::sccs(el, n);
    cout<<sccs.size()<<"\n";
    for(vi scc: sccs) {
        cout<<scc.size()<<" ";
        for(int v: scc) cout<<v<<" ";
        cout<<endl;
    }
}
```
])

#block( breakable: false,[
== Lowest Common Ancestor

```cpp
const int N = 1e6, M = 25;
int anc[M][N], h[N], rt;

// TODO: Calculate h[u] and set anc[0][u] = parent of node u for each u

// build (sparse table)
anc[0][rt] = rt; // set parent of the root to itself
for (int i = 1; i < M; ++i)
  for (int j = 1; j <= n; ++j)
    anc[i][j] = anc[i-1][anc[i-1][j]];

// query
int lca(int u, int v) {
  if (h[u] < h[v]) swap(u, v);
  for (int i = M-1; i >= 0; --i) if (h[u]-(1<<i) >= h[v])
    u = anc[i][u];

  if (u == v) return u;

  for (int i = M-1; i >= 0; --i) if (anc[i][u] != anc[i][v])
    u = anc[i][u], v = anc[i][v];
  return anc[0][u];
}
```
])

#block( breakable: false,[
== General matching in a graph

```cpp
vector<int> Blossom(vector<vector<int>> graph){
  int n = graph.size();
  int timer = -1;
  vector<int> mate(n, -1), label(n), parent(n), 
              orig(n), aux(n, -1), q;
  auto lca = [&](int x, int y) {
    for (timer++; ; swap(x, y)) {
      if (x == -1) continue;
      if (aux[x] == timer) return x;
      aux[x] = timer;
      x = (mate[x] == -1 ? -1 : orig[parent[mate[x]]]);
    }
  };
  auto blossom = [&](int v, int w, int a) {
    while (orig[v] != a) {
      parent[v] = w; w = mate[v];
      if (label[w] == 1) label[w] = 0, q.push_back(w);
      orig[v] = orig[w] = a; v = parent[w];
    }
  };
  auto augment = [&](int v) {
    while (v != -1) {
      int pv = parent[v], nv = mate[pv];
      mate[v] = pv; mate[pv] = v; v = nv;
    }
  };
  auto bfs = [&](int root) {
    fill(label.begin(), label.end(), -1);
    iota(orig.begin(), orig.end(), 0);
    q.clear(); 
    label[root] = 0; q.push_back(root);
    for (int i = 0; i < (int)q.size(); ++i) {
      int v = q[i];
      for (auto x : graph[v]) {
        if (label[x] == -1) {
          label[x] = 1; parent[x] = v;
          if (mate[x] == -1) 
            return augment(x), 1;
          label[mate[x]] = 0; q.push_back(mate[x]);
        } else if (label[x] == 0 && orig[v] != orig[x]) {
          int a = lca(orig[v], orig[x]);
          blossom(x, v, a); blossom(v, x, a);
        }
      }
    }
    return 0;
  };
  for (int i = 0; i < n; i++) 
    if (mate[i] == -1) 
      bfs(i);
  return mate;
}
```
])