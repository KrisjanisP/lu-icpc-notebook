// https://judge.yosupo.jp/submission/255593
#include <bits/stdc++.h>
#define all(x) x.begin(), x.end()
#define fi first
#define se second
#define rep(i,n) for(int i=0;i<n;i++)
#define pb push_back
using namespace std;
using vi = vector<int>;
using vvi = vector<vi>;
using pii = pair<int, int>;
using vpii = vector<pii>;

const int MAX_N = 5e5;

namespace kosaraju {
    vi adj[MAX_N], adjt[MAX_N], topo;
    bitset<MAX_N> vis;
    void visit(int v) {
        vis[v]=1;
        for(int u: adj[v]) if(!vis[u]) visit(u);
        topo.pb(v);
    }
    vvi res;
    void assign(int v, int r){
        vis[v]=0; res.back().pb(v);
        for(int u: adjt[v]) if(vis[u]) assign(u,r);
    }
    vvi sccs(vpii el, int n){
        rep(i,n) adj[i].clear(), adjt[i].clear();
        topo.clear(); res.clear();
        for(pii e: el) adj[e.fi].pb(e.se), adjt[e.se].pb(e.fi);
        rep(i,n) if(!vis[i]) visit(i);
        reverse(all(topo));
        rep(i, n){
            int v=topo[i];
            if(vis[v]) res.pb({}), assign(v,v);
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
    ios::sync_with_stdio(0); cin.tie(0);
    int n, m; cin>>n>>m;
    vpii edges(m);
    rep(i, m) cin>>edges[i].fi>>edges[i].se;
    vvi sccs = kosaraju::sccs(edges, n);
    cout<<sccs.size()<<"\n";
    for(vi s: sccs) {
        cout<<s.size()<<" ";
        for(int v: s) cout<<v<<" ";
        cout<<"\n";
    }
}