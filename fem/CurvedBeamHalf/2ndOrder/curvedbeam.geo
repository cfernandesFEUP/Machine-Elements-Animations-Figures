//+
Include "var.geo";
//+
b = 19.05;
zo = -b/2;
re = 152.4;
ri = 50.8;
deg = Pi/2;
//+
Point(1) = {0, 0, zo};
//+
Point(2) = {0, -ri, zo};
//+
Point(3) = {0, -re, zo};
//+
Point(4) = {re*Cos(deg), re*Sin(deg), zo};
//+
Point(5) = {ri*Cos(deg), ri*Sin(deg), zo};
//+
Line(1) = {2, 3};
//+
Circle(2) = {3, 1, 4};
//+
Line(3) = {4, 5};
//+
Circle(4) = {2, 1, 5};
//+
Transfinite Curve {1, 3} = nr Using Progression 1;
//+
Transfinite Curve {4, 2} = nt Using Progression 1;
//+
Curve Loop(1) = {1, 2, 3, -4};
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
Physical Curve("LOAD") = {11};
//+
Physical Surface("TOP") = {21};
//+
Physical Surface("BOTTOM") = {13};
//+
Physical Volume("BODY") = {1};
// MESH
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh 3;
Mesh.SaveGroupsOfNodes=1;
Save "curvedbeam.inp";
