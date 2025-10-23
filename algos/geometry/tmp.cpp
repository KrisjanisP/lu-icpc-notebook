#include<bits/stdc++.h>
using namespace std;

typedef long double ld;
typedef long long ll;
const ld eps=1e-12;

template<typename T>
struct Point{
    typedef Point<T> P;
    T x,y;
    Point(){}
    Point(T _x, T _y):x(_x),y(_y){}
    T dist_sq() {return x*x + y*y;}
    P operator-(P p)const{return P{x-p.x,y-p.y};}
    P operator+(P p){return P{x+p.x,y+p.y};}
    P operator*(T s){return P{x*s,y*s};}
    bool operator==(P rhs){return fabs(x-rhs.x)<eps && fabs(y-rhs.y)<eps;}
    bool operator < (const P &t)const
    {
        return fabs(x-t.x)<eps ? y<t.y : x<t.x;
    }
    // rotate +90 degrees
    P perp() const {return P{-y, x}; }
    Point<ll> to_ll() const {return Point<ll>(ll(x), ll(y));}
    Point<ld> to_ld() const {return Point<ld>(ld(x), ld(y));}
};

typedef Point<ld> P;


struct C{
    Point<ll> o; ll r;
    // circle-circle intersection
    vector<Point<ld>> operator & (const C& b){
        if(o==b.o){assert(r!=b.r); return {};}
        Point<ll> vec = b.o - o;
        ll d2 = vec.dist_sq(), sum = r+b.r, dif = r-b.r;
        ld p = (d2 + r*r - b.r*b.r)/(ld)(d2*2), h2 = r*r - p*p*d2;
        if (sum*sum < d2 || dif*dif > d2) return {};
        Point<ld> mid = o.to_ld() + vec.to_ld()*p;
        if(h2<1e-12) return {mid};
        Point<ld> per = vec.perp().to_ld() * sqrt(fmax(0, h2) / d2);
        return {mid + per, mid - per};
    }
}p[5];

struct DSU
{
    int fa[5];
    void init(int n)
    {
        for(int i=1;i<=n;i++)
            fa[i]=i;
    }
    int find(int x)
    {
        return fa[x]==x ? x : fa[x]=find(fa[x]);
    }
    void merge(int x,int y)
    {
        x=find(x),y=find(y);
        if(x!=y)fa[x]=y;
    }
}dsu;

int main()
{
    int n;
    scanf("%d",&n);
    for(int i=1;i<=n;i++)
        cin>>p[i].o.x>>p[i].o.y>>p[i].r;
    for(int i=1;i<=n;i++){
        for(int j=i+1;j<=n;j++){
            vector<Point<ld>> intersections = p[i] & p[j];
            if(!intersections.empty()) {
                cout << "Circle " << i << " and Circle " << j << " intersections: ";
                for(const auto& pt : intersections) {
                    cout << "(" << pt.x << ", " << pt.y << ") ";
                }
                cout << endl;
            }
        }
    }
    vector<P> all;
    dsu.init(n);
    int e=0;
    for(int i=1;i<=n;i++)
    {
        vector<P> pat;
        for(int j=1;j<=n;j++)if(i!=j)
        {
            vector<P> tmp=p[i]&p[j];
            if(!tmp.empty())dsu.merge(i,j);
            for(int k=0;k<(int)tmp.size();k++)
                all.push_back(tmp[k]),pat.push_back(tmp[k]);
        }
        sort(pat.begin(),pat.end());
        e+=unique(pat.begin(),pat.end())-pat.begin();
    }
    sort(all.begin(),all.end());
    int v=unique(all.begin(),all.end())-all.begin(),c=0;
    for(int i=1;i<=n;i++)
        c+=(dsu.find(i)==i);
    cout<<e-v+c+1<<endl;
    return 0;
}
