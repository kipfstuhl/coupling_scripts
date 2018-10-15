lc = 0.006;
// lc = 0.02378414230005447204;

//+
Point(1) = {0, -0.5, 0, lc};
//+
Point(2) = {0, 0.5, 0, lc};
//+
Point(3) = {1, 0.5, 0, lc};
//+
Coherence;
//+
Recursive Delete {
  Point{3};
}
//+
Point(3) = {2, 0.5, 0, lc};
//+
Point(4) = {1.7, -0.5, 0, lc};
//+
Point(5) = {2.3, 0.6, 0, lc};
//+
Point(6) = {2.8, 0.7, 0, lc};
//+
Point(7) = {3.5, 0.7, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 7};
//+
Point(8) = {3.5, 0.3, 0, lc};
//+
Point(9) = {2.7, 0.3, 0, lc};
//+
Point(10) = {2.3, 0.2, 0, lc};
//+
Point(11) = {2.1, 0.1, 0, lc};
//+
Point(12) = {2.3, -0.1, 0, lc};
//+
Point(13) = {2.9, -0.4, 0, lc};
//+
Point(14) = {3.5, -0.7, 0, lc};
//+
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 10};
//+
Line(9) = {10, 11};
//+
Line(10) = {11, 12};
//+
Line(11) = {12, 13};
//+
Line(12) = {13, 14};

n1 = -0.3 / Sqrt(0.3*0.3 + 0.6*0.6);
n2 = -0.6 / Sqrt(0.3*0.3 + 0.6*0.6);

r2 = 0.3;

Point(15) = {3.5 + n1*r2, -0.7 + n2*r2, 0, lc};
Point(16) = {2.9 + n1*r2, -0.4 + n2*r2, 0, lc};
Point(17) = {(1.7 + 2.9 + n1*r2)/2, (-0.5 -0.4+n2*r2)/2+0.08, 0, lc};
//+
Line(13) = {14, 15};
//+
Line(14) = {15, 16};
//+
Line(15) = {16, 17};
//+
Line(16) = {4, 1};
//+
Line(17) = {17, 4};
//+
Point(18) = {1.5, -0, 0, lc};
//+
Line Loop(1) = {16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17};
//+
Plane Surface(1) = {1};
//+
Physical Line("Wall_bottom", 1) = {16, 17, 15, 14};
//+
Physical Line("Outflow_1", 2) = {13};
//+
Physical Line("Wall_right", 3) = {12, 11, 10, 9, 8, 7};
//+
Physical Line("Outflow_2", 4) = {6};
//+
Physical Line("Wall_top", 5) = {5, 4, 3, 2};
//+
Physical Line("Inflow", 6) = {1};
//+
Physical Surface("S", 20) = {1};
//+
Hide "*";
//+
Show {
Point{14,15};
Line{13};
}
