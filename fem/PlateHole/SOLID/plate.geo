Include "input.geo";
R = (dH*W)/2;
SS = Round(Max(W,H)/10);

// MESH SIZE
lc = 0.1;
ni = size*SS;
elsize = SS/ni;
nH = Round((H-2*SS)/(4*elsize));
nW = Round((W-2*SS)/(4*elsize));
delta = 0.9;

// POINTS
Point (1) = {0, 0, 0, lc};
Point (2) = {.5*R, .5*R, 0, lc};
Point (3) = {.5*R, -.5*R, 0, lc};
Point (4) = {-.5*R, -.5*R, 0, lc};
Point (5) = {-.5*R, .5*R, 0, lc};

Point (6) = {SS, SS, 0, lc};
Point (7) = {SS, -SS, 0, lc};
Point (8) = {-SS, -SS, 0, lc};
Point (9) = {-SS, SS, 0, lc};

Point (10) = {W/2, H/2, 0, lc};
Point (11) = {W/2, SS, 0, lc};
Point (12) = {W/2, -SS, 0, lc};
Point (13) = {W/2, -H/2, 0, lc};
Point (14) = {SS, -H/2, 0, lc};
Point (15) = {-SS, -H/2, 0, lc};
Point (16) = {-W/2, -H/2, 0, lc};
Point (17) = {-W/2, -SS, 0, lc};
Point (18) = {-W/2, SS, 0, lc};
Point (19) = {-W/2, H/2, 0, lc};
Point (20) = {-SS, H/2, 0, lc};
Point (21) = {SS, H/2, 0, lc};

// LINES
Line (1) = {6, 2};
Line (2) = {7, 3};
Line (3) = {8, 4};
Line (4) = {9, 5};

// CIRCLES
Circle (5) = {2, 1, 3};
Circle (6) = {3, 1, 4};
Circle (7) = {4, 1, 5};
Circle (8) = {5, 1, 2};

// INTERNAL SQUARE
Line (9) = {6, 7};
Line (10) = {7, 8};
Line (11) = {8, 9};
Line (12) = {9, 6};

// EXTERNAL SQUARE
Line (13) = {10, 11};
Line (14) = {11, 12};
Line (15) = {12, 13};
Line (16) = {13, 14};
Line (17) = {14, 15};
Line (18) = {15, 16};
Line (19) = {16, 17};
Line (20) = {17, 18};
Line (21) = {18, 19};
Line (22) = {19, 20};
Line (23) = {20, 21};
Line (24) = {21, 10};

// INNER LINES
Line (25) = {6, 11};
Line (26) = {7, 12};
Line (27) = {7, 14};
Line (28) = {8, 15};
Line (29) = {8, 17};
Line (30) = {9, 18};
Line (31) = {9, 20};
Line (32) = {6, 21};

// TRANSFINITE LINES
Transfinite Curve {1:4} = ni Using Progression delta;
Transfinite Curve {5:14,17,20,23} = ni Using Progression 1;
Transfinite Curve {16,18,22,24,25,26,29,30} = nW Using Progression 1;
Transfinite Curve {13,15,19,21,27,28,31,32} = nH Using Progression 1;


// CURVE LOOPS
Curve Loop (1) = {1, 5, -2, -9};
Curve Loop (2) = {2, 6, -3, -10};
Curve Loop (3) = {3, 7, -4, -11};
Curve Loop (4) = {-1, -12, 4, 8};

Curve Loop (11) = {-13, -24, -32, 25};
Curve Loop (12) = {-14, -25, 9, 26};
Curve Loop (13) = {-15, -26, 27, -16};
Curve Loop (14) = {-27, 10, 28, -17};
Curve Loop (15) = {-28, 29, -19, -18};
Curve Loop (16) = {11, 30, -20, -29};
Curve Loop (17) = {31, -22, -21, -30};
Curve Loop (18) = {-23, -31, 12, 32};

// SURFACES
Plane Surface (1) = {1};
Plane Surface (2) = {2};
Plane Surface (3) = {3};
Plane Surface (4) = {4};
Plane Surface (11) = {11};
Plane Surface (12) = {12};
Plane Surface (13) = {13};
Plane Surface (14) = {14};
Plane Surface (15) = {15};
Plane Surface (16) = {16};
Plane Surface (17) = {17};
Plane Surface (18) = {18};

// TRANSFINITE SURFACES
Transfinite Surface {1} = {2,6,7,3};
Transfinite Surface {2} = {3,7,8,4};
Transfinite Surface {3} = {4,8,9,5};
Transfinite Surface {4} = {5,9,6,2};

Transfinite Surface {11} = {10,11,6,21};
Transfinite Surface {12} = {11,12,7,6};
Transfinite Surface {13} = {12,13,14,7};
Transfinite Surface {14} = {14,15,8,7};

Transfinite Surface {15} = {15,16,17,8};
Transfinite Surface {16} = {17,18,9,8};
Transfinite Surface {17} = {18,19,20,9};
Transfinite Surface {18} = {20,21,6,9};

Recombine Surface {1:4,11:18};

// PHYSICAL ENTITIES
Physical Line("fix") = {16:18};
Physical Line("load") = {22:24};
Physical Line("hole") = {5:8};
Physical Surface('plate') = {1:5,11:18};

// MESH
Mesh 2;
Mesh.SaveGroupsOfNodes=1;
Save "plate.inp";
