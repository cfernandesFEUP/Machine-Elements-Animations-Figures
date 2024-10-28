import os

## RUN Gmsh
os.system("gmsh cylinder.geo")

## REMOVE LINES
linelist = open("cylinder.inp").readlines()
newfile = open('mesh.inp', 'w')
flag = 1

for line in linelist:
    if line.startswith("*ELEMENT, type=T3D2"):
        flag = 0
    if line.startswith("*ELEMENT, type=CPS4"):
        flag = 1
    if flag and not line.startswith("EndModuleData"):
       newfile.writelines(line)

with open('mesh.inp', 'r') as file :
  filedata = file.read()
# Replace the target string
filedata = filedata.replace('CPS4', 'S4')

with open('mesh.inp', 'w') as file:
  file.write(filedata)

## RUN CalculiX
os.system("ccx load")

## CONVERT TO PARAVIEW
os.system("ccx2paraview load.frd vtu") 
