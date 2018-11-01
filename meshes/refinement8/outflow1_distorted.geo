lc = .02378414230005447204;
// lc = 0.16;

r = 0.5;
h1 = 0.1;
l = 4;
h2 = 0.5;
r1 = 0.5;
h3 = 0.00;
h4 = 0.1;

nnorm = Sqrt( (3/16*l)*(3/16*l) + (h2/2)*(h2/2));
n1 = (3/16*l)/nnorm;
n2 = (h2/2)/nnorm;

r1 = 0.5;

t1 = n2;
t2 = -n1;

Point(1) = {3*l/16+l/4, 0, 0 , lc};

Point(15) = {l/2, -(r + h3), 0, lc};
Point(16) = {l/2+l/8, -(r + 3*h3/2), 0, lc};

Point(17) = {13/16*l, -(r + 3*h3/2+h4/2),0, lc};
Point(18) = {l, -(r + 3*h3/2+h4),0, lc};
BSpline(6) = {15, 16, 17, 18};

nnorm = Sqrt( (3/16*l)*(3/16*l) + (h4/2)*(h4/2));
n1 = (3/16*l)/nnorm;
n2 = -(h4/2)/nnorm;

r2 = 0.3;

t1_ = -n2;
t2_ = +n1;

Point(19) = {l + t1_*r2, -(r + 3*h3/2+h4) + t2_*r2,0, lc};
Line(7) = {18, 19};

Point(20) = {13/16*l+ t1_*r2, -(r + 3*h3/2+h4/2)+ t2_*r2,0, lc};
Point(21) = {l/2+l/8+ t1_*r2, -(r + 3*h3/2)+ t2_*r2, 0, lc};

Point(12) = {(l/2+l/8 + t1*r1 + l/2+l/8+ t1_*r2)/2, (r + 3*h1/2 + t2*r1 -(r + 3*h3/2)+ t2_*r2)/2, 0, lc};
//+
BSpline(8) = {19, 20, 21, 12};
//+
Line(9) = {15, 1};
//+
Line(10) = {1, 12};
//+
Line Loop(1) = {9, 10, -8, -7, -6};
//+
Plane Surface(1) = {1};
//+
Physical Line("Wall_down") = {6};
//+
Physical Line("Outflow_1", 2) = {7};
//+
Physical Line("Wall_right", 3) = {8};
//+
Physical Line("Interface_1", 4) = {10};
//+
Physical Line("Interface_2", 5) = {9};
//+
Physical Surface("S", 6) = {1};
