lc = 0.05;

r = 0.3;
h = 3.0/7.0*r;
l1 = 0.8;
l2 = l1/10;

Point(1) = {0, -r, 0, lc};
Point(2) = {0, r, 0, lc};
Point(3) = {l1, r, 0, lc};
Point(4) = {l1+l2, -h, 0, lc};
Point(5) = {l1, -r, 0 ,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,1};

Line Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
//+
Physical Line(1) = {5};
//+
Physical Line(2) = {4};
//+
Physical Line(3) = {3};
//+
Physical Line(4) = {2};
//+
Physical Line(5) = {1};
//+
Physical Surface("S", 6) = {1};
