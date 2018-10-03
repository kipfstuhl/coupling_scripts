lc = .03363585661014864017;

Point(4) = {1.7, -0.5, 0, lc};
Point(11) = {2.1, 0.1, 0, lc};
Point(12) = {2.3, -0.1, 0, lc};
Point(13) = {2.9, -0.4, 0, lc};
Point(14) = {3.5, -0.7, 0, lc};
n1 = -0.3 / Sqrt(0.3*0.3 + 0.6*0.6);
n2 = -0.6 / Sqrt(0.3*0.3 + 0.6*0.6);

r2 = 0.3;

Point(15) = {3.5 + n1*r2, -0.7 + n2*r2, 0, lc};
Point(16) = {2.9 + n1*r2, -0.4 + n2*r2, 0, lc};
Point(17) = {(1.7 + 2.9 + n1*r2)/2, (-0.5 -0.4+n2*r2)/2+0.08, 0, lc};
Point(18) = {1.5, -0, 0, lc};
//+
Line(1) = {18, 4};
//+
Line(2) = {4, 17};
//+
Line(3) = {17, 16};
//+
Line(4) = {16, 15};
//+
Line(5) = {15, 14};
//+
Line(6) = {14, 13};
//+
Line(7) = {13, 12};
//+
Line(8) = {12, 11};
//+
Line(9) = {11, 18};
//+
Line Loop(1) = {2, 3, 4, 5, 6, 7, 8, 9, 1};
//+
Plane Surface(1) = {1};
//+
Physical Line("Wall_bottom", 1) = {2, 3, 4};
//+
Physical Line("Outflow", 2) = {5};
//+
Physical Line("Wall_up", 3) = {6, 7, 8};
//+
Physical Line("Interface_2", 4) = {9};
//+
Physical Line("Interface_1", 5) = {1};
//+
Physical Surface("S") = {1};
