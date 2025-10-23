#include <bits/stdc++.h>
using namespace std;
using ld = long double;
using ll = long long;

const ld eps=1e-12;

template<typename T>
struct Point{
    typedef Point<T> P;
    T x,y;
    T dist_sq() {return x*x + y*y;}
    P operator-(P p)const{return P{x-p.x,y-p.y};}
    P operator+(P p){return P{x+p.x,y+p.y};}
    P operator*(T s){return P{x*s,y*s};}
    // be careful with precision issues!
    // consider {return fabs(x-rhs.x)<eps && fabs(y-rhs.y)<eps;}
    bool operator==(P rhs){return x==rhs.x&&y==rhs.y;}
    P perp() const {return P{-y, x}; } //+90 deg

    Point<ld> to_ld() const {return Point<ld>{ld(x),ld(y)};}
    Point<ll> to_ll() const {return Point<ll>{ll(x),ll(y)};}
};

struct C{
    Point<ll> o; ll r;
    // circle-circle intersection
    vector<Point<ld>> operator & (const C& b){
        if(o==b.o){assert(r!=b.r); return {};}
        Point<ll> vec = b.o - o;
        ld d2 = vec.dist_sq(), sum = r+b.r, dif = r-b.r,
        p = (d2 + r*r - b.r*b.r)/(d2*2), h2 = r*r - p*p*d2;
        if (sum*sum < d2 || dif*dif > d2) return {};
        if(h2==0) return {o.to_ld()};
        Point<ld> mid = o.to_ld() + vec.to_ld()*p,
            per = vec.perp().to_ld() * sqrt(fmax(0, h2) / d2);
        return {mid + per, mid - per};
    }
};

int main() {
    C a = {{0, 0}, 5};
    C b = {{10, 0}, 5};
    vector<Point<ld>> intersections = a & b;
    if(!intersections.empty()) {
        if(intersections.size() == 1) {
            std::cout << "One intersection point: (" << intersections[0].x << ", " << intersections[0].y << ")" << std::endl;
        } else {
            std::cout << "Intersection points: (" << intersections[0].x << ", " << intersections[0].y << ") and (" << intersections[1].x << ", " << intersections[1].y << ")" << std::endl;
        }
    } else {
        std::cout << "No intersection" << std::endl;
    }
}
