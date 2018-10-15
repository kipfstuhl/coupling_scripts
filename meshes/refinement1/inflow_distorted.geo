lc = .09600000000000000000;

r = 0.5;
h1 = 0.1;
l = 4;

Point(1) = {0, -r, 0, lc};
Point(2) = {0, r, 0, lc};
Line(1) = {1,2};
Point(3) = {l/4, r, 0 , lc};
Point(4) = {l/8+l/4, r+h1/2, 0 , lc};
Point(5) = {l/2, r + h1, 0, lc};

BSpline(2) = {2, 3, 4, 5};

h3 = 0.00;
Point(13) = {l/4, -r, 0 , lc};
Point(14) = {l/8+l/4, -(r+h3/2), 0 , lc};
Point(15) = {l/2, -(r + h3), 0, lc};

BSpline(5) = {1, 13, 14, 15};

Point(16) = {3*l/16+l/4, 0, 0 , lc};
//+
Line(6) = {15, 16};
//+
Line(7) = {16, 5};
//+
Physical Line("Wall_down", 1) = {5};
//+
Physical Line("Interface_1", 2) = {6};
//+
Physical Line("Interface_2", 3) = {7};
//+
Physical Line("Wall_up", 4) = {2};
//+
Physical Line("Inflow", 5) = {1};
//+
Line Loop(1) = {2, -7, -6, -5, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("S") = {1};
