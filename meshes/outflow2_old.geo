lc = 0.05;

r = 0.3;
h = 3.0/7.0*r;
l1 = 0.8;
l2 = l1/10;
l3 = l1/2;
l4 = 2;

Point(1) = {l1+l2, -h, 0, lc};
Point(2) = {l1+l2+l3, -h/2, 0, lc};
Point(3) = {l1+l2+l3, 2*h, 0, lc};
Point(4) = {l4-(3*h+h/3), 3*h, 0, lc};
Point(5) = {l4, 3*h, 0, lc};
Point(6) = {l4, 5*h, 0, lc};
Point(7) = {l4-(3*h+h/3), 5*h, 0, lc};
Point(8) = {l1+2*h, r+2*h, 0, lc};
Point(9) = {l1, r, 0, lc};

Line(1) = {1,2};
BSpline(2) = {2, 3, 4, 5};

Line(3) = {5,6};
BSpline(4) = {6, 7, 8, 9};

Line(5) = {9,1};
//+
Physical Line(1) = {1};
//+
Physical Line(2) = {2};
//+
Physical Line(3) = {3};
//+
Physical Line(4) = {4};
//+
Physical Line(5) = {5};
//+
Line Loop(1) = {4, 5, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("S") = {1};
