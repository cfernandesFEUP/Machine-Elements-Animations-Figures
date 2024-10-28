//+
nz = 51;
nt = 31;
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
//+
Transfinite Curve {1} = nt Using Progression 1;
//+
Extrude {0, 0, 100} {
  Curve{1}; Layers{nz}; Recombine;
}
//+
Physical Curve("fix") = {1};
Physical Curve("load") = {3};
Physical Surface("cylinder") = {1};
//+
Mesh 2;
Mesh.SaveGroupsOfNodes=1;
Save "cylinder.inp";
