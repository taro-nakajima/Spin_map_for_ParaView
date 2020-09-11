import math
import numpy as np

#system size
points_x = 80
points_y = 80
num_points=points_x*points_y

#modulation period
Lambda_m=10

#scale factor for the length of arrows
Arrow_scale=0.2

#120 degree rotation matrix about the z axis
rad120=math.radians(120.0)
Rot120 = np.matrix([
[np.cos(rad120),-np.sin(rad120),0.0],
[np.sin(rad120),np.cos(rad120),0.0],
[0.0,0.0,1.0]
])

z_ofst=0.001

#definitions of Q-vectors
q_len = 2.0*math.pi/Lambda_m
Q1=np.matrix([[0.0],[q_len],[0.0]])
Q2=Rot120*Q1
Q3=Rot120*Rot120*Q1

print(Q2)
print(Q3)

FH = open("SkL_02.vtk","w")

FH.write("# vtk DataFile Version 2.0\n")
FH.write("test\n")
FH.write("ASCII\n")
FH.write("DATASET UNSTRUCTURED_GRID\n")
FH.write("POINTS %d float\n"%(num_points))

for ii in range(points_x):
    x=float(ii)
    for jj in range(points_y):
        y=float(jj)
        z=0.0
        FH.write("{0} {1} {2}\n".format(x,y,z))

FH.write("POINT_DATA %d\n"%(num_points))
FH.write("VECTORS vector float\n")

FHT=open("temporary.txt","w")

for ii in range(points_x):
    for jj in range(points_y):
        r=np.matrix([[float(ii),float(jj),0.0]])
        phase1=r*Q1
        phase2=r*Q2
        phase3=r*Q3
        S1=np.matrix([[np.sin(phase1)],[0.0],[-np.cos(phase1)]])

        S2=np.matrix([[np.sin(phase2)],[0.0],[-np.cos(phase2)]])
        S2_rotated=Rot120*S2

        S3=np.matrix([[np.sin(phase3)],[0.0],[-np.cos(phase3)]])
        S3_rotated=Rot120*Rot120*S3

        S_all=S1+S2_rotated+S3_rotated
        S_all=S_all/np.sqrt(np.dot(S_all.A1,S_all.A1))*Arrow_scale
        FH.write("{0} {1} {2}\n".format(float(S_all[0,0]),float(S_all[1,0]),float(S_all[2,0])))
        FHT.write("{0}\n".format(float(S_all[2,0])))

FH.write("SCALARS z float\n")
FH.write("LOOKUP_TABLE defalut\n")

FHT.close()
FHT=open("temporary.txt","r")
for line in FHT:
    FH.write(line)

FHT.close()
FH.close()
#for ii in range(points_x):
#    x=math.cos(math.radians(float(ii)*5))
#    for jj in range(points_y):
#        y=math.sin(math.radians(float(ii)*5))
#        z=0.0
#        FH.write("{0}\n".format(x))
