cl__1 = 0.005;
Point(1) = {0, 0, 0, 0.005};
Point(2) = {0, 0.5, 0, 0.005};
Point(3) = {0.5, 0.5, 0, 0.005};
Point(4) = {0.5, 0, 0, 0.005};
Point(5) = {0.24, 0, 0, 0.005};
Point(6) = {0.26, 0, 0, 0.005};
Point(7) = {0.26, 0.16, 0, 0.005};
Point(8) = {0.24, 0.16, 0, 0.005};
Point(9) = {0.24, 0.34, 0, 0.005};
Point(10) = {0.26, 0.34, 0, 0.005};
Point(11) = {0.26, 0.5, 0, 0.005};
Point(12) = {0.24, 0.5, 0, 0.005};
Line(1) = {2, 1};
Line(2) = {1, 5};
Line(3) = {5, 8};
Line(4) = {8, 7};
Line(5) = {7, 6};
Line(6) = {6, 4};
Line(7) = {4, 3};
Line(8) = {3, 11};
Line(9) = {11, 10};
Line(10) = {10, 9};
Line(11) = {9, 12};
Line(12) = {12, 2};
Line Loop(14) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(14) = {14};
Line Loop(15) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(15) = {15};
Line Loop(16) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(16) = {16};
Line Loop(17) = {12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
Plane Surface(17) = {17};
Physical Line(18) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Physical Surface(19) = {14};

Mesh.Algorithm = 2
