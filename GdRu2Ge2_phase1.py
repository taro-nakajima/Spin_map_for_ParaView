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

#90 degree rotation matrix about the z axis
rad90=math.radians(90.0)

#definitions of lattice constants
a_len=4.165
c_len=9.610

a_vec=np.array([a_len,0.0,0.0])
b_vec=np.array([0.0,a_len,0.0])
c_vec=np.array([0.0,0.0,c_len])


#definitions of Q-vectors
#q-vectorはreciprocal lattice unit (a*, b*, c*のunit)で書くようにしましょう。a*を掛け算する必要はないです。
q_len = 0.343 #=qa*=0.226*1.51812
Q1=np.array([q_len,0.0,0.0])
Q2=np.array([0.0,q_len,0.0])

FH = open("GdRu2Ge2_phase1.vtk","w")   #ファイル名にはスペースを使わないようにしましょう。

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

        SA1_vec=np.array([0,0,-1])
        SB1_vec=np.array([0,1,0])
        SA2_vec=np.array([0,0,-1])
        SB2_vec=np.array([-1,0,0])
        Mz=0#2가 되면 스킬미온이 아니게됨. 스핀적분이 4pi가 되도록하려면..? 계산할 수도 있지만, 1(?)~1.9정도쯤 중앙부가 나오고 메론-안티메론상태가 없어지도록 해야함. 0일땐 제로자기장이라서 스킬미온이 아닌 메론-안티메론 상태가 됨.(IC-1)
        Sz_vec=np.array([0,0,Mz])
        S1=np.cos(phase1)*SA1_vec+np.sin(phase1)*SB1_vec #a direction? or a* direction?
        S2=np.cos(phase2)*SA2_vec+np.sin(phase2)*SB2_vec

        S_all=S1+S2+Sz_vec
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