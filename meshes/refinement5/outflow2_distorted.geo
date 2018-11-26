SetFactory("OpenCASCADE");
lc = .018;
// lc = 0.176;

r = 0.5;
h1 = 0.1;
l = 4;

Point(1) = {0, -r, 0, lc};
Point(2) = {0, r, 0, lc};
Point(3) = {l/4, r, 0 , lc};
Point(4) = {l/8+l/4, r+h1/2, 0 , lc};
Point(5) = {l/2, r + h1, 0, lc};

h2 = 0.5;

Point(6) = {l/2+l/8, r + 3*h1/2, 0, lc};

Point(7) = {13/16*l, r + 3*h1/2+h2/2,0, lc};
Point(8) = {l, r + 3*h1/2+h2,0, lc};

nnorm = Sqrt( (3/16*l)*(3/16*l) + (h2/2)*(h2/2));
n1 = (3/16*l)/nnorm;
n2 = (h2/2)/nnorm;

r1 = 0.5;

t1 = n2;
t2 = -n1;

Point(9) = {l + t1*r1, r + 3*h1/2+h2 + t2*r1,0, lc};

Point(10) = {13/16*l + t1*r1, r + 3*h1/2+h2/2 + t2*r1,0, lc};
Point(11) = {l/2+l/8 + t1*r1, r + 3*h1/2 + t2*r1, 0, lc};

h3 = 0.00;
Point(13) = {l/4, -r, 0 , lc};
Point(14) = {l/8+l/4, -(r+h3/2), 0 , lc};
Point(15) = {l/2, -(r + h3), 0, lc};

Point(16) = {l/2+l/8, -(r + 3*h3/2), 0, lc};

h4 = 0.1;
Point(17) = {13/16*l, -(r + 3*h3/2+h4/2),0, lc};
Point(18) = {l, -(r + 3*h3/2+h4),0, lc};

nnorm = Sqrt( (3/16*l)*(3/16*l) + (h4/2)*(h4/2));
n1 = (3/16*l)/nnorm;
n2 = -(h4/2)/nnorm;

r2 = 0.3;

t1_ = -n2;
t2_ = +n1;

Point(19) = {l + t1_*r2, -(r + 3*h3/2+h4) + t2_*r2,0, lc};
//+

Point(20) = {13/16*l+ t1_*r2, -(r + 3*h3/2+h4/2)+ t2_*r2,0, lc};
Point(21) = {l/2+l/8+ t1_*r2, -(r + 3*h3/2)+ t2_*r2, 0, lc};


Point(12) = {(l/2+l/8 + t1*r1 + l/2+l/8+ t1_*r2)/2, (r + 3*h1/2 + t2*r1 -(r + 3*h3/2)+ t2_*r2)/2, 0, lc};
Point(22) = {3*l/16+l/4, 0, 0 , lc};
//+
Line(1) = {22, 5};
//+
Line(2) = {22, 12};
//+
BSpline(3) = {5, 6, 7, 8};
//+
BSpline(4) = {12, 11, 10, 9};
//+
Line(5) = {9, 8};
//+
Line Loop(1) = {2, 4, 5, -3, -1};
//+
Plane Surface(1) = {1};
//+
Physical Line("Interface_1", 1) = {2};
//+
Physical Line("Wall_down", 2) = {4};
//+
Physical Line("Outflow_2", 3) = {5};
//+
Physical Line("Wall_up", 4) = {3};
//+
Physical Line("Interface_2", 5) = {1};
//+
Physical Surface("S",6) = {1};
