import numpy as np
import abdbeam as ab
from matplotlib import pyplot as plt
from itertools import chain
import scipy.stats
from sympy import Curve, line_integrate, E, ln
from sympy import *
from sympy.abc import x, y, t
from sympy.vector import vector_integrate
import pandas as pd


def laminaload2(qmfx,abd_c,TT,QQ):
    Epsilon12fx = []
    Sigma12fx = []
    ekfx = np.matmul(abd_c, qmfx.transpose())#přetvoření
    exyfx = np.delete(ekfx, [3, 4, 5])
    e12fx = []
    sig12fx = []
    for j in range(len(TT)):
        ee12fx = np.matmul(np.linalg.inv(TT[j].transpose()), exyfx)
        e12fx.append(ee12fx)
        sigfx = np.round(np.matmul(QQ[j], ee12fx), 7)
        sig12fx.append(sigfx)
    Epsilon12fx.append(e12fx)
    Sigma12fx.append(sig12fx)

    eps=np.zeros([np.size(Epsilon12fx,axis=1),3])
    sigma=np.zeros([np.size(Sigma12fx,axis=1),3])
    for i in range(np.size(Epsilon12fx,axis=1)):
        eps[:][i]=Epsilon12fx[0][i]
    for i in range(np.size(Sigma12fx,axis=1)):
        sigma[:][i]=Sigma12fx[0][i]
    Epsilon12fx=eps
    Sigma12fx=sigma
    return(ekfx,Epsilon12fx,Sigma12fx)


profile = np.loadtxt("31mm_critical_crosssection.txt")
points = len(profile)

x = []
y = []

# Delka tetivy

chord = 276.82

for i in range(0,points):

    bod = profile[i]

    x.append((bod[0])*chord)
    y.append((bod[1])*chord)


#                           Stanoveni odstredive sily

#    Otacky [s-1]

n = 100

#    V1 - Vnejsi cast - objem [mm3]
#    V2 - Prostredni cast - objem [mm3]
#    V3 - Vnitrni cast - objem [mm3]

V1 = 9970
V2 = 25521.7
V3 = 26554.4


#    A1 - Vnejsi cast - povrch [mm2]
#    A2 - Prostredni cast - povrch [mm2]
#    A3 - Vnitrni cast - povrch [mm2]

A1 = 10194.2
A2 = 15524.4
A3 = 10271.7


#    hustota jadra g/mm3

jadro_hustota = 0.00017

#    Vrstveni

#    Vnejsi cast listu - plosne hmotnosti vrstev - vrstvy_1

vrstvy_1 = [4*0.000045, 4*0.000095, 1*0.0002]
vrstvy_2 = [4*0.000045, 4*0.000095, 2*0.0002]
vrstvy_3 = [4*0.000045, 4*0.000095, 3*0.0002]


m_j_1 = jadro_hustota * V1
m_j_2 = jadro_hustota * V2
m_j_3 = jadro_hustota * V3


m1 = 0
m2 = 0
m3 = 0

#    Hmotnost vnejsi cast listu - vrstvy

for i in range(0, len(vrstvy_1)):
    m1 += A1 * vrstvy_1[i]

#    Hmotnost prostredni cast listu - vrstvy

for i in range(0, len(vrstvy_2)):
    m2 += A2 * vrstvy_2[i]

#    Hmotnost vnitrni cast listu - vrstvy

for i in range(0, len(vrstvy_3)):
    m3 += A3 * vrstvy_3[i]

#                                   Vnejsi tretina listu
#    Stred segmentu 1 [m]

r1 = 0.3175
Fd1 = 4 * (np.pi**2) * n**2 * r1 * ((m_j_1 + m1)/1000)


#                                   Prostredni tretina listu
#    Stred segmentu 2 [m]

r2 = 0.1905
Fd2 = 4 * (np.pi**2) * n**2 * r2 * ((m_j_2 + m2)/1000)

#                                   Vnitrni tretina listu
#    Stred segmentu 3 [m]

r3 = 0.079
Fd3 = 4 * (np.pi**2) * n**2 * r3 * ((m_j_3 + m3)/1000)

#    Celkova odstrediva sila pusobici v miste rezu

Fd = Fd1 + Fd2 + Fd3

print("Odstrediva sila: ", Fd)

print("Pocet definovanych bodu:", points)

# abdbeam

sc = ab.Section()

# materials



mts = dict()
mts[1] = ab.Laminate()
mts[2] = ab.Laminate()

# Unidirectional carbon fiber - SIGRATEX U200-0/SO
# t = 0.30 mm
# [1]

c_ud = ab.PlyMaterial(0.3, 86460, 3400, 2581, 0.28)
mts[1].ply_materials[1] = c_ud
mts[2].ply_materials[1] = c_ud

# Carbon fabric SIGRATEX® C W95-PL1/1 (KDL 8023), 95g/m2
# t = 0.15 mm
# [2] # stojina +45 -45

c_fabric2 = ab.PlyMaterial(0.075, 60543, 60543, 2581, 0.037)
mts[1].ply_materials[2] = c_fabric2
mts[2].ply_materials[2] = c_fabric2

# Glass fabric Interglas 02034/106 24,5 g/m2
# t = 0.05 mm
# [3]

g_fabric1 = ab.PlyMaterial(0.025, 16600, 16600, 2800, 0.037)
mts[1].ply_materials[3] = g_fabric1
mts[2].ply_materials[3] = g_fabric1

# Layering

mts[1].plies = [[+45,3], [-45,3], [+45,2], [-45,2], [0,1], [0,1], [0,1]]
mts[2].plies = [[+45,2], [-45,2]]

# No symmetry
mts[1].symmetry = "T"
mts[1].calculate_properties()
mts[2].symmetry = "T"
mts[2].calculate_properties()


# points

pts = {}
for i in range (0, len(x)):
   pts[i+1] = ab.Point(x[i], y[i])


# segments

sgs = {}

for i in range(0, len(x)-1):
    sgs[i+1] = ab.Segment(i+1, i+2, 1)


# Propojeni na odtokove hrane

sgs[len(x)] = ab.Segment(len(x), 1, 1)


# Tloustka

t1 = 0.05 + 0.15 + 3*0.3
t_st = 0.15
print("Tloustka profil", t1)
print("Tloustka stojina", t_st)

#                                             Stojina

#                                             Definovani polohy stojiny

poloha_stojiny = 0.4

#               Vzdalenost od nabezneho bodu

l_stojiny = chord * poloha_stojiny


#               Nalezeni nejblizsich bodu

min_distance = 1000
m_0 = 1000

for i in range(0, points):
    if m_0 < abs(x[i] - l_stojiny):
        break

    else:

        min_distance_i = abs(x[i] - l_stojiny)

        if min_distance > min_distance_i:
            min_distance = min_distance_i
            point_i = i


min_distance_2 = 1000


n_0 = 1000

for n in reversed(range(0, points)):
    if n_0 < abs(x[n] - l_stojiny):
        break

    else:
        min_distance_i_2 = abs(x[n] - l_stojiny)

        if min_distance_2 > min_distance_i_2:

            min_distance_2 = min_distance_i_2
            point_i_2 = n

            n_0 = min_distance_i_2


# Segment - stojina, material 2


sgs[len(x)+1] = ab.Segment(point_i + 1, point_i_2+1, 2)


sc.materials = mts
sc.points = pts
sc.segments = sgs

#Vypocet

sc.calculate_properties()
sc.summary()

ab.plot_section(sc, segment_coord=True, figsize=(6.4*0.8, 4.8*0.8))

# Zatizeni

sc.loads[1] = ab.Load(Px_c=Fd, My=-38.55*1000, Mz=-0.48*1000, Vz_s=167, Vy_s=2.06, Tx=0.94*1000, yv=0.25*chord, zv=0)
#Nmm, N


sc.calculate_internal_loads()

# Vykresleni vnitrnich ucinku


ab.plot_section_loads(sc, 1, int_load_list=['Nx'],
                      title_list=['Abdbeam - Nx (N/mm)',
                      'Abdbeam - Nx (N/mm)'], figsize=(6.4*1.2, 4.8*1.2))

ab.plot_section_loads(sc, 1, int_load_list=['Nxy'],
                      title_list=['Abdbeam - Nxy (N/mm)',
                      'Abdbeam - Nxy (N/mm)'], figsize=(6.4*1.2, 4.8*1.2))



ab.plot_section_loads(sc, 1, int_load_list=['Mx'],
                      title_list=['Abdbeam - Mx (N/mm)',
                      'Abdbeam - Mx (N/mm)'], figsize=(6.4*1.2, 4.8*1.2))

ab.plot_section_loads(sc, 1, int_load_list=['My'],
                      title_list=['Abdbeam - My (N/mm)',
                      'Abdbeam - My (N/mm)'], figsize=(6.4*1.2, 4.8*1.2))

ab.plot_section_loads(sc, 1, int_load_list=['Mxy'],
                      title_list=['Abdbeam - Mxy (N/mm)',
                      'Abdbeam - Mxy (N/mm)'], figsize=(6.4*1.2, 4.8*1.2))

#ab.plot_section_loads(sc, 1, int_load_list=['Nx', 'Ny', 'Nxy', 'Mx', 'My', 'Mxy'],
#                      title_list=['Abdbeam - Nx (N/mm)', 'Abdbeam - Ny (N/mm)','Abdbeam - Nxy (N/mm)',
#'Abdbeam - Mx (N/mm)','Abdbeam - My (N/mm)','Abdbeam - Mxy (N/mm)'], figsize=(6.4*0.8, 4.8*0.8))

# Vypocet deformace v obecnem s.s. ke strednicove rovine

# Prepocet deformace na deformace jednotlivych vrstev

# Prepocet deformace vrstvy na napeti v materialovem s.s.

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

Nxy = sc.sgs_int_lds_df.get("Nxy")
Nxy_max = Nxy.get("Max")
Nx = sc.sgs_int_lds_df.get("Nx")
Nx_max = Nx.get("Max")


sigma_red = []

# Redukovane napeti

red_max = 0


for i in range(0, points-1):
    if (((Nx_max[i]/t1)**2)+(4*((Nxy_max[i]/t1)**2)))**(0.5) > red_max:
        red_max = (((Nx_max[i]/t1)**2)+(4*((Nxy_max[i]/t1)**2)))**(0.5)
        red_max_pozice = i


print("V profilu max. red. napeti:", red_max, "v segmentu", red_max_pozice)


red_max_stojina = ((((Nx_max[points]/t_st)**2)+(4*((Nxy_max[points]/t_st)**2)))**(0.5))
red_max_pozice_stojina = points
print("Maximalni redukovane napeti ve stojine", red_max_stojina)

Nx = sc.sgs_int_lds_df.get("Nx")
Nx_max_profil = Nx.get("Max")[red_max_pozice]
Nxy = sc.sgs_int_lds_df.get("Nxy")
Nxy_max_profil = Nxy.get("Max")[red_max_pozice]
Mx = sc.sgs_int_lds_df.get("Mx")
Mx_max_profil = Mx.get("Max")[red_max_pozice]
My = sc.sgs_int_lds_df.get("My")
My_max_profil = My.get("Max")[red_max_pozice]
Mxy = sc.sgs_int_lds_df.get("Mxy")
Mxy_max_profil= Mxy.get("Max")[red_max_pozice]



Nx = sc.sgs_int_lds_df.get("Nx")
Nx_max_stojina = Nx.get("Max")[red_max_pozice_stojina]
Nxy = sc.sgs_int_lds_df.get("Nxy")
Nxy_max_stojina = Nxy.get("Max")[red_max_pozice_stojina]
Mx = sc.sgs_int_lds_df.get("Mx")
Mx_max_stojina = Mx.get("Max")[red_max_pozice_stojina]
My = sc.sgs_int_lds_df.get("My")
My_max_stojina = My.get("Max")[red_max_pozice_stojina]
Mxy = sc.sgs_int_lds_df.get("Mxy")
Mxy_max_stojina= Mxy.get("Max")[red_max_pozice_stojina]


zatizeni_profil = np.zeros(6)

zatizeni_profil[0] = Nx_max_profil
zatizeni_profil[1] = 0
zatizeni_profil[2] = Nxy_max_profil
zatizeni_profil[3] = Mx_max_profil
zatizeni_profil[4] = My_max_profil
zatizeni_profil[5] = Mxy_max_profil

zatizeni_stojina = np.zeros(6)

zatizeni_stojina[0] = Nx_max_stojina
zatizeni_stojina[1] = 0
zatizeni_stojina[2] = Nxy_max_stojina
zatizeni_stojina[3] = Mx_max_stojina
zatizeni_stojina[4] = My_max_stojina
zatizeni_stojina[5] = Mxy_max_stojina


napeti_profil = laminaload2(zatizeni_profil,mts[1].abd_c,mts[1].TT,mts[1].QQ)[2]
napeti_stojina = laminaload2(zatizeni_stojina,mts[2].abd_c,mts[2].TT,mts[2].QQ)[2]


print("Napeti ve vrstvach laminatu v kritickem segmentu na profilu", napeti_profil)

print("Napeti ve vrstvach laminatu ve stojine", napeti_stojina)



for i in (1,len(mts[1].plies)):
    if napeti_profil[i-1,0] > 50:
        print("Napeti ve smeru vlaken v profilu vetsi, nez maximalni napeti ve vrstve", i, "Velikost napeti = ", napeti_profil[i,0])
    if napeti_profil[i-1,1] > 250:
        print("Napeti proti smeru vlaken v profilu vetsi, nez maximalni napeti ve vrstve", i, "Velikost napeti = ", napeti_profil[i,1])
    if napeti_profil[i-1,2] > 70:
        print("Smykove napeti v profilu vetsi, nez maximalni napeti ve vrstve", i, "Velikost napeti = ", napeti_profil[i,2])

for i in (1,len(mts[2].plies)):
    if napeti_profil[i-1,0] > 50:
        print("Napeti ve smeru vlaken ve stojine vetsi, nez maximalni napeti ve vrstve", i, "Velikost napeti = ", napeti_profil[i,0])
    if napeti_profil[i-1,1] > 250:
        print("Napeti proti smeru vlaken ve stojine vetsi, nez maximalni napeti ve vrstve", i, "Velikost napeti = ", napeti_profil[i,1])
    if napeti_profil[i-1,2] > 70:
        print("Smykove napeti ve stojine vetsi, nez maximalni napeti ve vrstve", i, "Velikost napeti = ", napeti_profil[i,2])
