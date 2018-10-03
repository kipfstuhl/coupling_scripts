lc = .05232511306011979183;

Point(3) = {2, 0.5, 0, lc};
Point(5) = {2.3, 0.6, 0, lc};
Point(6) = {2.8, 0.7, 0, lc};
Point(7) = {3.5, 0.7, 0, lc};
Point(8) = {3.5, 0.3, 0, lc};
Point(9) = {2.7, 0.3, 0, lc};
Point(10) = {2.3, 0.2, 0, lc};
Point(11) = {2.1, 0.1, 0, lc};
Point(18) = {1.5, -0, 0, lc};
//+
Line(1) = {18, 3};
//+
Line(2) = {3, 5};
//+
Line(3) = {5, 6};
//+
Line(4) = {6, 7};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 9};
//+
Line(7) = {9, 10};
//+
Line(8) = {10, 11};
//+
Line(9) = {11, 18};
//+
Physical Line("Interface_2", 1) = {9};
//+
Physical Line("Wall_bottom", 2) = {8, 7, 6};
//+
Physical Line("Outflow", 3) = {5};
//+
Physical Line("Wall_up", 4) = {2, 3, 4};
//+
Physical Line("Interface_3", 5) = {1};
//+
Line Loop(1) = {9, 1, 2, 3, 4, 5, 6, 7, 8};
//+
Plane Surface(1) = {1};
//+
Physical Surface("S", 10) = {1};
