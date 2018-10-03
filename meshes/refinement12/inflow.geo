lc = .01427048538003268959;

Point(1) = {0, -0.5, 0, lc};
Point(2) = {0, 0.5, 0, lc};
Point(3) = {2, 0.5, 0, lc};
Point(4) = {1.7, -0.5, 0, lc};
Point(5) = {1.5, -0, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 5};
Line(4) = {5, 4};
Line(5) = {4, 1};

Physical Line("Wall_bottom", 1) = {5};
Physical Line("Interface_1", 2) = {4};
Physical Line("Interface_3", 3) = {3};
Physical Line("Wall_up", 4) = {2};
Physical Line("Left", 5) = {1};
//+
Line Loop(1) = {2, 3, 4, 5, 1};
//+
Plane Surface(1) = {1};
//+
Physical Surface("S", 6) = {1};
