import os

os.system("gmsh disc.geo")

## REMOVE LINES
linelist = open("disc.inp").readlines()
newfile = open('mesh.inp', 'w')
flag = 1

for line in linelist:
    if line.startswith("*ELEMENT, type=M3D9"):
        flag = 0
    if line.startswith("*ELEMENT, type=C3D27"):
        flag = 1
    if flag and not line.startswith("EndModuleData"):
       newfile.writelines(line)

with open('mesh.inp', 'r') as file :
  filedata = file.read()
# Replace the target string
filedata = filedata.replace('C3D27', 'C3D20R')

with open('mesh.inp', 'w') as file:
  file.write(filedata)

## RUN CalculiX
os.system("ccx solve")

## CONVERT TO PARAVIEW
os.system("ccx2paraview solve.frd vtu") 
