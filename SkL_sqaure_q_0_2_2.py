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
rad90=math.radians(90.0)

#definitions of lattice constants
a_len=4.234
c_len=9.901

a_vec=np.array([a_len,0.0,0.0])
b_vec=np.array([0.0,a_len,0.0])
c_vec=np.array([0.0,0.0,c_len])


#definitions of Q-vectors
#q-vectorはreciprocal lattice unit (a*, b*, c*のunit)で書くようにしましょう。a*を掛け算する必要はないです。
#q_len = 0.14*2*math.pi/a_x

q_len = 0.2968 #qa*=0.2*1.484
Q1=np.array([q_len,0.0,0.0])
Q2=np.array([0.0,q_len,0.0])

FH = open("SkL_sqaure_q_0_2_3.vtk","w")   #ファイル名にはスペースを使わないようにしましょう。

#ここから、位置座標(x, y, z)を出力します。出力する際には直交座標系に直す必要があります。

FH.write("# vtk DataFile Version 2.0\n")
FH.write("test\n")
FH.write("ASCII\n")
FH.write("DATASET UNSTRUCTURED_GRID\n")
FH.write("POINTS %d float\n"%(num_points))

for ii in range(points_x):
    for jj in range(points_y):
        r_vec=float(ii)*a_vec+float(jj)*b_vec #결정구조?
        FH.write("{0} {1} {2}\n".format(r_vec[0],r_vec[1],r_vec[2])) #r=(a,b,c)플로트?

#位置座標の出力はここまで。

#ここからスピンベクトルの計算

SA_vec = np.array([0,0,-1]) #이게 뭐냐돌릴 스핀 기본?
SB_vec = np.array([1,0,0])

FH.write("POINT_DATA %d\n"%(num_points))
FH.write("VECTORS vector float\n")

FHT=open("temporary.txt","w")   #これはz成分だけ書いておくためのバッファー

#formatrix, S스핀 돌려주고, 원하는 방향으로 돌리기
for ii in range(points_x):
    for jj in range(points_y):
        r_vec2=np.array([float(ii),float(jj),0.0])
        phase1=(2.0*np.pi)*np.dot(r_vec2,Q1)        # r=(na,nb,0)とQ1=(q,0,0)の内積(dot積)。rはa,b,c unit, Qはa*,b*,c* unit.
        phase2=(2.0*np.pi)*np.dot(r_vec2,Q2)

        S1=np.cos(phase1)*SA_vec+np.sin(phase1)*SB_vec+np.cos(phase2)*SA_vec+np.sin(phase2)*SB_vec #a direction? or a* direction?
        #b* direction????

        S_all=S1
        S_all=S_all/np.linalg.norm(S_all)*Arrow_scale
        FH.write("{0} {1} {2}\n".format(float(S_all[0]),float(S_all[1]),float(S_all[2])))
        FHT.write("{0}\n".format(float(S_all[2])))

FH.write("SCALARS z float\n")
FH.write("LOOKUP_TABLE defalut\n")

FHT.close()

#スピンベクトル、Sx, Sy, Szはここまでで終わり。

#最後に色付けのためにSz成分だけ書き出す。

FHT=open("temporary.txt","r")
for line in FHT:
    FH.write(line)

FHT.close()
'''
FH.close()
#for ii in range(points_x):
#    x=math.cos(math.radians(float(ii)*5))
#    for jj in range(points_y):
#        y=math.sin(math.radians(float(ii)*5))
#        z=0.0
#        FH.write("{0}\n".format(x))
'''