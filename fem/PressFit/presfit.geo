Geometry.CopyMeshingMethod=1;
d=25;
D=50;
r=5;
L1=2*D;
L2=D;
//+
Point(1) = {0, r+d/2-(2.5*r*Sin(Pi/4)), 0, 1.0};
//+
Point(2) = {L1+L2, r+d/2-(2.5*r*Sin(Pi/4)), 0, 1.0};
//+
Point(3) = {L1+L2, d/2, 0, 1.0};
//+
Point(4) = {L1+r, d/2, 0, 1.0};
//+
Point(5) = {L1, r+d/2, 0, 1.0};
//+
Point(6) = {L1, d, 0, 1.0};
//+
Point(7) = {0, d, 0, 1.0};
//+
Point(8) = {L1+r, r+d/2, 0, 1.0};
//+
Point(9) = {L1+r-(r*Cos(Pi/4)), r+d/2-(r*Sin(Pi/4)), 0, 1.0};
//+
Point(10) = {L1+r-(2.5*r*Cos(Pi/4)), r+d/2-(2.5*r*Sin(Pi/4)), 0, 1.0};
//+
Point(11) = {L1+r, r+d/2-(2.5*r*Sin(Pi/4)), 0, 1.0};
//+
Point(12) = {L1+r-(2.5*r*Cos(Pi/4)), r+d/2, 0, 1.0};
//+
Point(13) = {L1+r-(2.5*r*Cos(Pi/4)), d, 0, 1.0};
//+
Point(14) = {0, r+d/2, 0, 1.0};
//+
Point(15)={0,d/5,0,1};
//+
Point(16)={0,0,d/5,1};
//+
Point(17)={0,-d/5,0,1};
//+
Point(18)={0,0,-d/5,1};
//+
Line(1) = {1, 10};
//+
Line(2) = {10, 11};
//+
Line(3) = {11, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {3, 4};
//+
Circle(6) = {4, 8, 9};
//+
Circle(7) = {9, 8, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {6, 13};
//+
Line(10) = {13, 7};
//+
Line(11) = {7, 14};
//+
Line(12) = {14, 1};
//+
Line(13) = {4, 11};
//+
Line(14) = {9, 10};
//+
Line(15) = {5, 12};
//+
Line(16) = {12, 14};
//+
Line(17) = {13, 12};
//+
Line(18) = {12, 10};
//+
Transfinite Curve {-1, 16, 10} = 19 Using Progression 1;
//+
Transfinite Curve {3, -5} = 19 Using Progression 1;
//+
Transfinite Curve {12, 18, 7} = 11 Using Progression 1;
//+
Transfinite Curve {11, 17, 8} = 6 Using Progression 1;
//+
Transfinite Curve {13, -4, 14, 15, 9} = 6 Using Progression 1.25;
//+
Transfinite Curve {6, 2} = 11 Using Progression 1;
//+
Curve Loop(1) = {1, -18, 16, 12};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, -13, 6, 14};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {14, -18, -15, -7};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {3, 4, 5, 13};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {15, -17, -9, -8};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {16, -11, -10, 17};
//+
Surface(6) = {6};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};
//+
Transfinite Surface {6};
//+
Transfinite Surface {5};
//+
Recombine Surface {1, 3, 2, 4, 5, 6};

//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Surface{4}; Surface{2}; Surface{3}; Surface{5}; Surface{1}; Surface{6}; Layers{10}; Recombine;
}


//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Surface{128}; Surface{150}; Surface{106}; Surface{84}; Surface{62}; Surface{40}; Layers{10}; Recombine;
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Surface{172}; Surface{194}; Surface{216}; Surface{238}; Surface{260}; Surface{282}; Layers{10}; Recombine;
}
//+
Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Surface{304}; Surface{326}; Surface{348}; Surface{370}; Surface{392}; Surface{414}; Layers{10}; Recombine;
}
//+
Line(528) = {16, 15};
//+
Line(529) = {15, 18};
//+
Line(530) = {18, 17};
//+
Line(531) = {17, 16};
//+
Curve Loop(7) = {529, 530, 531, 528};
//+
Surface(533) = {7};
//+
Transfinite Curve {530, 529, 528, 531} = 11 Using Progression 1;
//+
Transfinite Surface {533};
//+
Recombine Surface {533};
//+
Line(532) = {18, 220};
//+
Line(533) = {17, 164};
//+
Line(534) = {16, 136};
//+
Line(535) = {15, 1};
//+
Curve Loop(8) = {535, -421, -532, -529};
//+
Surface(534) = {8};
//+
Curve Loop(9) = {532, -289, -533, -530};
//+
Surface(535) = {9};
//+
Curve Loop(10) = {533, -157, -534, -531};
//+
Surface(536) = {10};
//+
Curve Loop(11) = {534, -113, -535, -528};
//+
Surface(537) = {11};
//+
Transfinite Curve {535, 532, 533, 534} = 4 Using Progression 1;
//+
Transfinite Surface {534, 535, 536, 537};
//+
Recombine Surface {534, 535, 536, 537};

//+
Extrude {L1+r-(2.5*r*Cos(Pi/4)), 0, 0} {
  Surface{533}; Surface{537}; Surface{534}; Surface{535}; Surface{536}; Layers{18}; Recombine;
}


//+
Extrude {(2.5*r*Cos(Pi/4)), 0, 0} {
  Surface{647}; Surface{581}; Surface{603}; Surface{625}; Surface{559}; Layers{10}; Recombine;
}

//+
Extrude {L2-r, 0, 0} {
  Surface{691}; Surface{757}; Surface{713}; Surface{735}; Surface{669}; Layers{18}; Recombine;
}

//+
Physical Surface("Fix") = {185, 171, 141, 127, 448, 435, 317, 303, 535, 537, 534, 533, 536};
//+
Physical Surface("Load") = {31, 273, 405, 528, 779, 867, 845, 823, 801};
//+
Physical Volume("Shaft") = {14, 13, 28, 25, 26, 5, 6, 20, 19, 27, 29, 7, 8, 15, 16, 17, 33, 34, 31, 2, 3, 4, 21, 22, 23, 32, 30, 11, 10, 9, 18, 38, 36, 35, 1, 39, 12, 37, 24};
