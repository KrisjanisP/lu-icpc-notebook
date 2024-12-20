
= Geometry

== Dot product

#grid(
    columns: 2,
    gutter: 3em,
[
#align(center)[#image("assets/dot.png", width: 10em)]
],
[
$ a dot b = |a| |b| cos(theta) $
$ a dot b = a_x b_x + a_y b_y $
$ theta = arccos(frac(a_x b_x + a_y b_y , |a| |b|)) $
],
)

Projection of a onto b:
$frac(a dot b, |b|) $

== Cross product

$ a times b = |a| |b| sin(theta) = a_x b_y - a_y b_x $

$theta$ is positive if a is clockwise from b

== Line-point distance

Line given by $a x+b y+c=0$ and point $(x_0, y_0)$.

$ "distance" = frac(|a x_0 + b y_0 +c|, sqrt(a^2+b^2)) $

The coordinates of this point are:

$ x = frac(b(b x_0 - a y_0)- a c,a^2 + b^2) #h(1em) y = frac(a(-b x_0 + a y_0)- b c,a^2 + b^2) $

== Shoelace formula

$ 2A = sum_(i=1)^n space mat(delim:"|", x_i, y_i; x_(i+1), y_(i+1)) , "where" mat(delim: "|", a, b; c, d) = a d - b c $

#block(breakable: false,[
== Circumradius

Let $a$, $b$, $c$ be the sides of a triangle and $A$ the area of the triangle. Then the circumradius $R = a b c slash(4A)$.
Alternatively, using the Law of Sines:

$ R = a/(2 sin(alpha)) = b/(2 sin(beta)) = c/(2 sin(gamma)) $

where $alpha$, $beta$, and $gamma$ are the angles opposite sides $a$, $b$, and $c$ respectively.
])

#block(breakable: false,[
== Law of Sines

In any triangle with sides $a$, $b$, $c$ and opposite angles $alpha$, $beta$, $gamma$ respectively:

$ frac(a,sin(alpha)) = frac(b,sin(beta)) = frac(c,sin(gamma)) = 2R $

where $R$ is the circumradius of the triangle. This can be rearranged to find any side or angle:
$a = 2R sin(alpha)$ and $sin(alpha) = frac(a,2R)$
])

#block(breakable: false,[
== Law of Cosines

In any triangle with sides $a$, $b$, $c$ and opposite angles $alpha$, $beta$, $gamma$ respectively:

$ c^2 = a^2 + b^2 - 2 a b cos(gamma) $

])
    
#block(breakable: false,[
== Median Length Formulas
    
In any triangle with sides $a$, $b$, $c$, the lengths of the medians $m_a$, $m_b$, $m_c$ from the respective vertices are given by:
    
$ m_a = frac(1, 2) sqrt(2b^2 + 2c^2 - a^2) $

These formulas can be derived using the Apollonius's theorem.
])

#block(breakable: false,[
== Segment to line linear equation

Converting segment $((P_x,P_y),(Q_x,Q_y))$ to $A x + B y + C = 0$:

$ ( P_y-Q_y)x + (Q_x-P_x)y + (P_x Q_y - P_y Q_x) = 0 $

])

== Three point orientation

```cpp
int orientation(Point p1, Point p2, Point p3){
    int val = (p2.y-p1.y)*(p3.x-p2.x)-(p2.x-p1.x)*(p3.y-p2.y);
    if (val == 0) return 0; // collinear
    return (val > 0) ? 1 : 2; // clock or counterclock
}```

== Line-line intersection

From system of linear equations derived Cramer's rule:

$ cases(a_1 x + b_1 y = c_1, a_2 x + b_2 y = c_2) => cases(x = (c_1 b_2 - c_2 b_1) slash (a_1 b_2 - a_2 b_1), y = (a_1 c_2 - a_2 c_1) slash (a_1 b_2 - a_2 b_1)) $

If the denominator equals zero, the lines are parallel or coincident.


#block(breakable: false,[
== Check if two segments intersect

```cpp
bool on_seg(Point p, Point q, Point r) { 
    return (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y)) 
} 
bool do_intersect(Point p1, Point q1, Point p2, Point q2){ 
    int o1 = orient(p1, q1, p2), o2 = orient(p1, q1, q2);
    int o3 = orient(p2, q2, p1), o4 = orient(p2, q2, q1);
    if (o1 != o2 && o3 != o4) return true;
    return (o1==0&&on_seg(p1,p2,q1))||(o2==0&&on_seg(p1,q2,q1))||
           (o3==0&&on_seg(p2,p1,q2))||(o4==0&&on_seg(p2,q1,q2));
} 
```
])

#block(breakable: false,[
== Heron's formula
Let $a$, $b$, $c$ - sides of a triangle. Then the area $A$ is:

$ A = 1/4sqrt((a+b+c)(-a+b+c)(a-b+c)(a+b-c)) $

Numerically stable version:

$ a >= b >= c,  A = 1/4sqrt((a+(b+c))(c-(a-b))(c+(a-b))(a+(b-c))) $

])

#block(breakable: false,[
== Graham's scan
Constructs convex hull of a set of points.
```cpp
void convex_hull(vector<pt>&a, bool coll=false){
    pt p = *min_element(all(a), [](const pt&x, const pt&y){
        return make_pair(x.y, x.x) < make_pair(y.y, y.x);
    });
    sort(all(a), [](const pt&x, const pt&y){
        int ori = orientation(p, x, y);
        if(ori == 0)
            return (p.x-x.x)*(p.x-x.x) + (p.y-x.y)*(p.y-x.y) <
                   (p.x-y.x)*(p.x-y.x) + (p.y-y.y)*(p.y-y.y);
        return ori < 0;
    });
    if(coll){
        int i = a.size()-1; while(i>=0 &&
            collinear(a[i], a.back(), p)) i--;
        reverse(a.b()+i+1, a.e());
    }
    vector<pt> s; 
    for(auto &p : a){
        while(s.size() > 1 &&
            !cw(s[s.size()-2], s.back(), p, coll)) s.pop_back();
        s.push_back(p);
    }
    a = s;
}
```
])

#block(breakable: false,[

== Closest pair of points
Finds pair of points with minimum euclidean distance.
```cpp
struct pt { ld x, y; int id; };
ld md; pair<int,int> bp; vector<pt> a, t;
void ua(const pt&a1, const pt&b1){ // update answer
    ld d=sqrtl((a1.x-b1.x)*(a1.x-b1.x)+(a1.y-b1.y)*(a1.y-b1.y));
    if(d<md){md=d; bp={a1.id,b1.id};}
}
void rec(int l, int r){ // recursive function
    if(r - l <= 3){
        rep(i, l, r) rep(j, i+1, r) ua(a[i], a[j]);
        sort(a.b()+l, a.b()+r,
            [](const pt&x, const pt&y){return x.y<y.y;});
        return;}
    int m = (l + r) >> 1, midx = a[m].x; 
    rec(l, m); rec(m, r);
    merge(a.b()+l, a.b()+m, a.b()+m, a.b()+r, t.b(),
        [](const pt&x, const pt&y){return x.y<y.y;});
    copy(t.b(), t.b()+r-l, a.b()+l);
    int ts = 0; rep(i, l, r) if(abs(a[i].x - midx) < md){
        for(int j=ts-1;j>=0&&a[i].y-t[j].y<md;j--)ua(a[i],t[j]);
        t[ts++] = a[i];
    }
}
```
])

