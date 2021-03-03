import math
import numpy as np

#system size
points_x = 80
points_y = 80
num_points=points_x*points_y

#modulation period
Lambda_m=44.3

#scale factor for the length of arrows
Arrow_scale=0.2

#180 degree rotation matrix about the z axis
rad180=math.radians(180.0)
Rot180 = np.matrix([
[np.cos(rad180),-np.sin(rad180),0.0],
[np.sin(rad180),np.cos(rad180),0.0],
[0.0,0.0,1.0]
])

#180 degree rotation matrix about the z axis
rad120=math.radians(120.0)
Rot120 = np.matrix([
[np.cos(rad120),-np.sin(rad120),0.0],
[np.sin(rad120),np.cos(rad120),0.0],
[0.0,0.0,1.0]
])

#30 degree rotation matrix about the z axis
rad30=math.radians(30.0)
Rot30 = np.matrix([
[np.cos(rad30),-np.sin(rad30),0.0],
[np.sin(rad30),np.cos(rad30),0.0],
[0.0,0.0,1.0]
])

#-30 degree rotation matrix about the z axis
radm30=math.radians(-30.0)
Rotm30 = np.matrix([
[np.cos(radm30),-np.sin(radm30),0.0],
[np.sin(radm30),np.cos(radm30),0.0],
[0.0,0.0,1.0]
])

#60 degree rotation matrix about the z axis
rad60=math.radians(60.0)
Rot60 = np.matrix([
[np.cos(rad60),-np.sin(rad60),0.0],
[np.sin(rad60),np.cos(rad60),0.0],
[0.0,0.0,1.0]
])

z_ofst=0.001
#definitions of lattice constants
a_x=4.075
b_x=4.075*np.cos(rad120)
b_y=4.075*np.sin(rad120)



#definitions of Q-vectors
q_len = 0.14*2*math.pi/a_x
Q1=np.matrix([[q_len],[0.0],[0.0]])
Q2=np.matrix([[0.0],[-q_len],[0.0]])
Q3=np.matrix([[q_len],[-q_len],[0.0]])


FH = open("SkL_02_1016.vtk","w")

FH.write("# vtk DataFile Version 2.0\n")
FH.write("test\n")
FH.write("ASCII\n")
FH.write("DATASET UNSTRUCTURED_GRID\n")
FH.write("POINTS %d float\n"%(num_points))

for ii in range(points_x):
    for jj in range(points_y):
        x=float(ii*a_x+jj*b_x)
        y=float(jj*b_y)
        z=0.0
        FH.write("{0} {1} {2}\n".format(x,y,z))
FH.write("POINT_DATA %d\n"%(num_points))
FH.write("VECTORS vector float\n")

FHT=open("temporary.txt","w")
#formatrix, S스핀 돌려주고, 원하는 방향으로 돌리기
for ii in range(points_x):
    for jj in range(points_y):
        r=np.matrix([[float(ii),float(jj),0.0]])
        phase1=r*Q1
        phase2=r*Q2
        phase3=r*Q3
        S1=np.matrix([[0.0],[np.sin(phase1)],[-np.cos(phase1)]])#a direction? or a* direction?
        S1_rotated=Rot30*S1 #a*direction?

        S2=np.matrix([[np.sin(phase2)],[0.0],[-np.cos(phase2)]])#b* direction????

        S3=np.matrix([0.0],[[np.sin(phase3)],[-np.cos(phase3)]])
        S3_rotated=Rotm30*S3

        S_all=S1_rotated+S2+S3_rotated
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
