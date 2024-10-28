import os

os.system("gmsh disc.geo")

## REMOVE LINES
linelist = open("disc.inp").readlines()
newfile = open('mesh.inp', 'w')
flag = 1

for line in linelist:
    if line.startswith("*ELEMENT, type=CPS4"):
        flag = 0
    if line.startswith("*ELEMENT, type=C3D8"):
        flag = 1
    if flag and not line.startswith("EndModuleData"):
       newfile.writelines(line)

## RUN CalculiX
os.system("ccx solve")

## CONVERT TO PARAVIEW
os.system("ccx2paraview solve.frd vtu") 
