#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef long double ld;

template<typename T>
struct Point{
    typedef Point<T> P;
    T x,y;
    T dist_sq() {return x*x + y*y;}
    P operator-(P p){return P{x-p.x,y-p.y};}
    P operator+(P p){return P{x+p.x,y+p.y};}
    P operator*(T s){return P{x*s,y*s};}
    bool operator==(P rhs){return x==rhs.x&&y==rhs.y;}
    // rotate +90 degrees
    P perp() const {return P{-y, x}; }
};

typedef Point<ld> P;
struct C{P o; ld r;};

bool circle_circle_i(C a, C b, std::pair<P,P>* out){
    if(a.o==b.o){assert(a.r!=b.r); return false;}
    P vec = b.o - a.o;
    ld d2 = vec.dist_sq(), sum = a.r+b.r, dif = a.r-b.r,
    p = (d2 + a.r*a.r - b.r*b.r)/(d2*2), h2 = a.r*a.r - p*p*d2;
    if (sum*sum < d2 || dif*dif > d2) return false;
    P mid = a.o + vec*p, per = vec.perp() * sqrt(fmax(0, h2) / d2);
    *out = {mid + per, mid - per};

    return true;
}

int main() {
    ll n;
    cin >> n;
    vector<C> circles(n);
    for(int i = 0; i < n; i++) {
        cin >> circles[i].o.x >> circles[i].o.y >> circles[i].r;
    }
    for(int i=0;i<n;i++){
        for(int j=i+1;j<n;j++){
            cout<<i<<" "<<j<<endl;

            std::pair<P,P> intersections;
            bool intersect = circle_circle_i(circles[i], circles[j], &intersections);

            if(intersect) {
                cout << "Yes" << endl;
                cout << "Intersection 1: (" << intersections.first.x << ", " << intersections.first.y << ")" << endl;
                cout << "Intersection 2: (" << intersections.second.x << ", " << intersections.second.y << ")" << endl;
            } else {
                cout << "No" << endl;
            }
        }
    }

}
