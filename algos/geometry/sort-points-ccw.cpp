// https://judge.yosupo.jp/submission/320609

#include <iostream>
#include <cassert>
#include <algorithm>
using namespace std;

typedef long long ll;

template<typename T>
struct Point{
    T x,y;
    T dist_sq() {return x*x + y*y;}
};


typedef Point<ll> P;

enum {COL, CW, CCW};
int orientation(P p1, P p2, P p3) {
    ll val = (p2.x-p1.x)*(p3.y-p2.y)-(p3.x-p2.x)*(p2.y-p1.y);
    if(val==0) return COL;
    return val > 0 ? CCW : CW;
}

int quadrant(P p) {
    // origin
    if(p.x==0&&p.y==0) return 0;
    // quadrants
    if(p.x>0&&p.y>0) return 1;
    if(p.x<0&&p.y>0) return 2;
    if(p.x<0&&p.y<0) return 3;
    if(p.x>0&&p.y<0) return 4;
    // x-axes
    if(p.x>0&&p.y==0) return 5;
    if(p.x<0&&p.y==0) return 6;
    // y-axes
    if(p.x==0&&p.y>0) return 7;
    if(p.x==0&&p.y<0) return 8;
    assert(false);
}

//          0 1 2 3 4 5 6 7 8
int q_o[9]={3,5,7,0,2,4,8,6,1};
// (go counterclockwise from third quadrant; end at negative x-axis)

int main() {
    int n; cin>>n;
    P p[n];
    for(int i=0;i<n;i++)
        cin>>p[i].x>>p[i].y;
    sort(p, p+n, [](P a, P b) {
        int q_a = quadrant(a), q_b = quadrant(b);
        if(q_a != q_b) return q_o[q_a] < q_o[q_b];
        int o = orientation(P{0,0}, a, b);
        if(o != COL) return o == CCW;
        else return a.dist_sq() < b.dist_sq();
    });
    for(P p: p) cout<<p.x<<" "<<p.y<<"\n";
}