import numpy as np
from sklearn.cluster import KMeans,MiniBatchKMeans,AgglomerativeClustering,OPTICS
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist,squareform
import math
import sys
import random


random.seed(42)

f=sys.argv[1]
nrep=int(sys.argv[2])
nc=int(sys.argv[3])
mode=int(sys.argv[4])

#readin positions from xyz
pos=[]
xyz = open(f, "r")
xyz.readline()
xyz.readline()
x=0
y=0
z=0
N=0
for line in xyz:
  atom, x1, y1, z1 = line.split()
#  print(atom,x1,y1,z1)
  x=x+float(x1)
  y=y+float(y1)
  z=z+float(z1)
  N=N+1
  if N==3:
    N=0
    pos.append([x/3.0,y/3.0,z/3.0])
    x=0
    y=0
    z=0

xyz.close()
print(len(pos))
print(pos[0])

X=pos

DMAT=np.zeros((len(pos),len(pos)))

basis=(9.8528*nrep,9.8528*nrep,9.8528*nrep)

if mode>=0:
  if mode==0:
    for i in range(0,len(pos)):
      for j in range(0,i+1):
        d=(pos[i][0]-pos[j][0]-0*basis[0])**2+(pos[i][1]-pos[j][1]-0*basis[1])**2+(pos[i][2]-pos[j][2]-0*basis[2])**2
        DMAT[i,j]=math.sqrt(d)
        DMAT[j,i]=DMAT[i,j]

  elif mode==1:
    for i in range(0,len(pos)):
      for j in range(0,i+1):
        d=1000000000.0
        a=range(0,2)
        for i1 in a:
          for i2 in a:
            for i3 in a:
              d=min(d,(pos[i][0]-pos[j][0]-i1*basis[0])**2+(pos[i][1]-pos[j][1]-i2*basis[1])**2+(pos[i][2]-pos[j][2]-i3*basis[2])**2)

        DMAT[i,j]=math.sqrt(d)
        DMAT[j,i]=DMAT[i,j]

  #db = KMeans(n_clusters=nc,init='k-means++',verbose=1).fit(DMAT)
  db = AgglomerativeClustering(n_clusters=nc,affinity='precomputed',linkage='single').fit(DMAT)
else:
  db = KMeans(n_clusters=nc,init='k-means++',verbose=1).fit(X)
  #db = MiniBatchKMeans(n_clusters=nc,init='k-means++',verbose=1).fit(X)
labels = db.labels_

for i in range(0,len(pos)):
  print("result",i+1,pos[i][0],pos[i][1],pos[i][2],labels[i])

