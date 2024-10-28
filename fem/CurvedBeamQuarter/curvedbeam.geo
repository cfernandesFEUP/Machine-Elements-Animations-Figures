//+
Include "var.geo";
//+
b = 19.05;
re = 152.4;
ri = 50.8;
deg = -Pi/2;
//+
Point(1) = {0, 0, 0};
//+
Point(2) = {ri, 0, 0};
//+
Point(3) = {re, 0, 0};
//+
Point(4) = {re*Cos(deg), re*Sin(deg), 0};
//+
Point(5) = {ri*Cos(deg), ri*Sin(deg), 0};
//+
Line(1) = {2, 3};
//+
Circle(2) = {3, 1, 4};
//+
Line(3) = {4, 5};
//+
Circle(4) = {5, 1, 2};
//+
Transfinite Curve {1, 3} = nr Using Progression 1;
//+
Transfinite Curve {4, 2} = nt Using Progression 1;
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
Recombine Surface{1};
//+
Transfinite Surface {1} = {2, 3, 4, 5};
//+
Extrude {0, 0, b} {
  Surface{1}; Layers{na}; Recombine;
}
//+
Physical Curve("LOAD") = {20};
//+
Physical Surface("LEFT") = {21};
//+
Physical Surface("RIGHT") = {13};
//+
Physical Surface("SHAFT") = {25};
//+
Physical Volume("BODY") = {1};
// MESH
Mesh 3;
Mesh.SaveGroupsOfNodes=1;
Save "curvedbeam.inp";
