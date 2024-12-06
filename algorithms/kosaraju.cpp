// https://judge.yosupo.jp/submission/252308
#include <bits/stdc++.h>
#define pb push_back
#define fi first
#define se second
#define rep(i,n) for(int i=0;i<n;i++)
#define all(x) x.begin(), x.end()
using namespace std;
using vi = vector<int>;
using vvi = vector<vector<int>>;
using pii = pair<int,int>;
using vpii = vector<pii>;

const int MAX_N=5e5;

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
// second dfs on transpose graph to get sccs.
// pushes reachable nodes to last scc in res
void assign(int v, int r) {
    vis[v]=0; res.back().pb(v);
    for(int u: adjt[v]) if(vis[u]) assign(u,r);
}
// clears previous graph, topo order, found sccs
void clear() {
    rep(i,MAX_N) adj[i].clear(), adjt[i].clear();
    topo.clear(); res.clear();
}
/// @param el 0-indexed edge list
/// @param n graph node count
/// @return strongly connected components
vvi kosaraju(vpii el, int n){
    clear();
    // construct graph
    for(pii e: el) adj[e.fi].pb(e.se), adjt[e.se].pb(e.fi);
    // topological order
    rep(i,n) if(!vis[i]) visit(i); reverse(all(topo));
    // sccs on transpose (reversed) graph
    rep(i,n) {
        int v=topo[i];
        if(vis[v]) res.push_back({}), assign(v,v);
    }
    return res;
}
}


int main() {
    ios::sync_with_stdio(false);
    int n, m; cin>>n>>m;
    vpii el;
    for(int i=0;i<m;i++){
        int a, b; cin>>a>>b;
        el.push_back({a,b});
    }
    vvi sccs=kosaraju::kosaraju(el, n);
    cout<<sccs.size()<<"\n";
    for(vi scc: sccs) {
        cout<<scc.size()<<" ";
        for(int v: scc) cout<<v<<" ";
        cout<<endl;
    }
}