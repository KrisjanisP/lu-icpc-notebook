// https://codeforces.com/contest/933/submission/345301686

#include<bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long ll;

const ld eps=1e-12;

template<typename T>
struct Point{
T x,y;

typedef Point<T> P;
T dist_sq() {return x*x + y*y;}
P operator-(P p)const{return P{x-p.x,y-p.y};}
P operator+(P p){return P{x+p.x,y+p.y};}
P operator*(T s){return P{x*s,y*s};}
bool operator==(P rhs){
    return fabs(x-rhs.x)<eps && fabs(y-rhs.y)<eps;}
bool operator < (const P &t)const{
    return fabs(x-t.x)<eps ? y<t.y : x<t.x;}
P perp() const {return P{-y, x}; } // +90 degrees
Point<ld> to_ld() const {return Point<ld>{ld(x), ld(y)};}
};

struct Circle{
Point<ll> o; ll r;
// circle-circle intersection
vector<Point<ld>> operator & (const Circle& b){
    if(o==b.o){assert(r!=b.r); return {};}
    Point<ll> vec = b.o - o;
    ll d2 = vec.dist_sq(), sum = r+b.r, dif = r-b.r;
    if (sum*sum < d2 || dif*dif > d2) return {};
    ld p = (d2 + r*r - b.r*b.r)/(ld)(d2*2), h2 = r*r - p*p*d2;
    Point<ld> mid = o.to_ld() + vec.to_ld()*p;
    if(h2<1e-12) return {mid};
    Point<ld> per = vec.perp().to_ld() * sqrt(fmax(0, h2) / d2);
    return {mid + per, mid - per};
}
};

int main() {
    Circle a = {{0, 0}, 5};
    Circle b = {{10, 0}, 5};
    vector<Point<ld>> i = a & b;
    if(!i.empty()) {
        if(i.size() == 1) {
            cout << "1 intersection point: (" << i[0].x << ", " << i[0].y << ")" << endl;
        } else {
            cout << "2 intersection points: (" << i[0].x << ", " << i[0].y << ") and (" << i[1].x << ", " << i[1].y << ")" << endl;
        }
    } else {
        cout << "No intersection" << endl;
    }
    return 0;
}
