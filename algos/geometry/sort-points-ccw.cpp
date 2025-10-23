// https://judge.yosupo.jp/submission/322148

#include <iostream>
#include <cassert>
#include <algorithm>
using namespace std;

typedef long long ll;

template<typename T>
struct Point{
    T x,y;
    T dist_sq() {return x*x + y*y;}
    Point<T> operator-(Point<T> p){
        return Point<T>{x-p.x,y-p.y};
    }
};

typedef Point<ll> P;

enum {COL, CW, CCW};
int orientation(P p1, P p2, P p3) {
    ll val = (p2.x-p1.x)*(p3.y-p2.y)-(p3.x-p2.x)*(p2.y-p1.y);
    if(val==0) return COL;
    return val > 0 ? CCW : CW;
}

// go ccw from 3rd quadrant; end at negative x-axis
int quadrant_order(P p) {
    int sx = (p.x > 0) - (p.x < 0);  // -1, 0, or 1
    int sy = (p.y > 0) - (p.y < 0);  // -1, 0, or 1
    static const int map[9] = { 0,1,2,8,3,4,7,6,5 };
    return map[(sy + 1) * 3 + (sx + 1)];
}

void sort_ccw(P points[], int n, P origin){
    sort(points, points+n, [origin](P a, P b) {
        a = a-origin, b = b-origin;
        int q_a = quadrant_order(a), q_b = quadrant_order(b);
        if(q_a != q_b) return q_a < q_b;
        int o = orientation(P{0,0}, a, b);
        if(o != COL) return o == CCW;
        else return a.dist_sq() < b.dist_sq();
    });
}

int main() {
    int n; cin>>n;
    P p[n];
    for(int i=0;i<n;i++)
        cin>>p[i].x>>p[i].y;
    sort_ccw(p,n,P{0,0});
    for(P p: p) cout<<p.x<<" "<<p.y<<"\n";
}
