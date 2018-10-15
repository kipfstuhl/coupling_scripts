lc = 0.06;

Point(4) = {0, 0, 0, lc};
Point(5) = {1, 0, 0, lc};
Point(6) = {1, 1, 0, lc};
Point(7) = {0, 1, 0, lc};
Line(3) = {7, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line Loop(7) = {3, 4, 5, 6};
Plane Surface(9) = {7};
//+
Physical Line("Bottom", 1) = {4};
//+
Physical Line("Right", 2) = {5};
//+
Physical Line("Top", 3) = {6};
//+
Physical Line("Left", 4) = {3};
//+
Physical Surface("Internal", 5) = {9};
