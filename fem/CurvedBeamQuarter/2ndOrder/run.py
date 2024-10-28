import os

R = input('Number of nodes in radial direction: ')
C = input('Number of nodes in circunferential direction: ')
Z = input('Number of nodes in axial direction: ')

with open('var.geo', 'w') as f:
	f.write('nr = '+ R +';\n')
	f.write('nt = '+ C +';')
	f.write('na = '+ Z +';')

os.system("gmsh curvedbeam.geo")

## REMOVE LINES
linelist = open("curvedbeam.inp").readlines()
newfile = open('mesh.inp', 'w')
flag = 1

for line in linelist:
    if line.startswith("*ELEMENT, type=T3D3"):
        flag = 0
    if line.startswith("*ELEMENT, type=C3D27"):
        flag = 1
    if flag and not line.startswith("EndModuleData"):
       newfile.writelines(line)

with open('mesh.inp', 'r') as file :
  filedata = file.read()
# Replace the target string
filedata = filedata.replace('C3D27', 'C3D20')

with open('mesh.inp', 'w') as file:
  file.write(filedata)

## RUN CalculiX
os.system("ccx solve")

## CONVERT TO PARAVIEW
os.system("ccx2paraview solve.frd vtu") 
