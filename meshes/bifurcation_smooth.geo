lc = 0.2;

r = 0.5;
h = 0.8;
l = 4;

Point(1) = {0, -r, 0, lc};
Point(2) = {0, r, 0, lc};
Line(1) = {1,2};
Point(3) = {l/2, r, 0 , lc};
Point(4) = {3*l/4, r+h/2, 0 , lc};
Point(5) = {l, r + h, 0, lc};

dx = 1*l/4;
dy = h/2;
n1 = -dy/Sqrt(dx*dx + dy*dy);
n2 = dx/Sqrt(dx*dx + dy*dy);

BSpline(2) = {2, 3, 4, 5};

// do the same for the bottom

Point(6) = {l/2, -r, 0 , lc};
Point(7) = {3*l/4, -r-h/2, 0 , lc};
Point(8) = {l, -r - h, 0, lc};

BSpline(3) = {1, 6, 7, 8};

r2 = 0.3;
Point(9) = {l - n1*2*r2, r + h - n2*2*r2, 0, lc};
Point(10) = {l - n1*2*r2, -r - h + n2*2*r2, 0, lc};

Line(4) = {5, 9};
Line(5) = {10, 8};

Point(11) = {3*l/4 - n1*2*r2, r+h/2 - n2*2*r2, 0 , lc};
Point(12) = {3*l/4 - n1*2*r2 - l/4, 0, 0 , lc};
// Point(12) = {1*l/3, 0, 0 ,lc};
Point(13) = {3*l/4 - n1*2*r2, -r-h/2 + n2*2*r2, 0 , lc};
//+
BSpline(6) = {10, 13, 12, 11, 9};
//+
Physical Line("Inlet", 1) = {1};
//+
Physical Line("Wall_1", 2) = {2};
//+
Physical Line("Outlet_1", 3) = {4};
//+
Physical Line("Wall_2", 4) = {6};
//+
Physical Line("Outlet_2", 5) = {5};
//+
Physical Line("Wall_3", 6) = {3};
//+
Line Loop(1) = {2, 4, -6, 5, -3, 1};
Plane Surface(1) = {1};
//+
Physical Surface("Internal", 7) = {1};
