import abdbeam as ab
import numpy as np
def laminaload2(alfa,qmfx,abd_c,TT,QQ):
  """
  alfa - vektor úhlů orientace vrstev
  qmfx - matice zatížení v jednotkách [N/m] - ve tvaru:
  [liniové zatížení na segment - osa x]
  [liniové zatížení na segment - osa y]
  [liniové zatížení na segment - osa z]
  [smykový tok osa - x]
  [smykový tok osa - y]
  [smykový tok osa - z] 
  abd_c - inverzní matice ABD
  TT - transformační matice 
  """
  Epsilon12fx = []
  Sigma12fx = []
  ekfx = np.matmul(abd_c, qmfx.transpose())#přetvoření
  exyfx = np.delete(ekfx, [3, 4, 5])
  e12fx = []
  sig12fx = []
  for j in range(len(alfa)):
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



sc = ab.Section()
mts = dict()
mts[1] = ab.Laminate()
ply_mat = ab.PlyMaterial(0.166666, 148000, 9650, 4550, 0.3)
mts[1].ply_materials[1] = ply_mat
mts[1].plies = [[45,1], [-45,1], [0,1], [90,1]]
mts[1].calculate_properties()
# print(mts[1].TT[2])#transformační matice
# print(mts[1].QQ[2])#matice tuhosti

alfa=45,-45,0,90#zkus dostat tohle z plies, neměl jsem na to čas, ale nemělo by to být složité
qmfx=np.zeros(6)
qmfx[0]=5#[N/m]vektor zatížení, které dostaneš z abdbeamu
print(laminaload2(alfa,qmfx,mts[1].abd_c,mts[1].TT,mts[1].QQ)[2])#sloupce jsou směry napětí X,Y,Z, řádky jsou vrstvy