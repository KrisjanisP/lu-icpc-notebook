#include <bits/stdc++.h>
using namespace std;
using ll = long long;
template<typename T>
struct P{T x,y;};
int orientation(P<ll> p1, P<ll> p2, P<ll> p3){
    ll val=(p2.y-p1.y)*(p3.x-p2.x)-(p2.x-p1.x)*(p3.y-p2.y);
    if(val==0) return 0;
    return val>0?1:-1;
}
int quadrant(P<ll> p){
    if(p.x>=0 && p.y>=0) return 1;
    if(p.x<=0 && p.y>=0) return 2;
    if(p.x<=0 && p.y<=0) return 3;
    return 4;
}
int q_o[5]={0,3,4,1,2};
int main(){
    int n; cin>>n;
    P<ll> p[n];
    for(int i=0;i<n;i++){
        int x,y; cin>>x>>y;
        p[i]={x,y};
    }
    sort(p,p+n,[](P<ll> a, P<ll> b){
        if(quadrant(a)==quadrant(b)){
            int o=orientation({0,0},a,b);
            if(o==0) return a.x<b.x;
            return o==1;
        } else {
            return q_o[quadrant(a)]<q_o[quadrant(b)];
        }
    });
    for(int i=0;i<n;i++) cout<<p[i].x<<" "<<p[i].y<<"\n";
}