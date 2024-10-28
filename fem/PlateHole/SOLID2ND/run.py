import os

width = input('Plate width: ')
height = input('Plate height: ')
thickness = input('Plate thickness: ')
ratio = input('Hole with ratio: ')
size = input('Mesh Size: ')

with open('input.geo', 'w') as f:
	f.write('dH = '+ ratio +';\n')
	f.write('W = '+ width +';\n')
	f.write('H = '+ height +';\n')
	f.write('size = '+ size +';')

with open('thickness.inp', 'w') as f:
	f.write(thickness)

## RUN Gmsh
os.system("gmsh plate.geo")

## REMOVE LINES
linelist = open("plate.inp").readlines()
newfile = open('mesh.inp', 'w')
flag = 1

for line in linelist:
    if line.startswith("*ELEMENT, type=T3D3"):
        flag = 0
    if line.startswith("*ELEMENT, type=CPS8"):
        flag = 1
    if flag and not line.startswith("EndModuleData"):
       newfile.writelines(line)

## LOAD TO HAVE 100 MPa
t = float(thickness)
db = float(ratio)
W = float(width)
P = 100*(W-db*W)*t

with open('load.inp', 'w') as file:
    file.write('DCload,2,'+str(P))

## RUN CalculiX
os.system("ccx stress")

## Convert to Paraview
os.system("ccx2paraview stress.frd vtk")

