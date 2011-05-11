import numpy as np
from test import *
from graphLib6 import *

g = Graph(3)

g.setVecf(g.addNode(0),Vecf(1,0,0,0))
g.setVecf(g.addNode(0),Vecf(1,1,0,0))
g.setVecf(g.addNode(0),Vecf(1,0,1,0))
g.setVecf(g.addNode(0),Vecf(1,0,0,1))
a = 1/SQRT(2)
g.setVecf(g.addNode(0),Vecf(1,0,a,a))
g.setVecf(g.addNode(0),Vecf(1,a,0,a))
g.setVecf(g.addNode(0),Vecf(1,a,a,0))
e12 = g.addNode(1); g.addArch(1,e12); g.addArch(2,e12);
e13 = g.addNode(1); g.addArch(1,e13); g.addArch(3,e13);
e14 = g.addNode(1); g.addArch(1,e14); g.addArch(4,e14);

e53 = g.addNode(1); g.addArch(5,e53); g.addArch(3,e53);
e54 = g.addNode(1); g.addArch(5,e54); g.addArch(4,e54);
e62 = g.addNode(1); g.addArch(6,e62); g.addArch(2,e62);
e64 = g.addNode(1); g.addArch(6,e64); g.addArch(4,e64);
e72 = g.addNode(1); g.addArch(7,e72); g.addArch(2,e72);
e73 = g.addNode(1); g.addArch(7,e73); g.addArch(3,e73);

mat456 = np.array([[1,1,1,1],[0,0,1,1],[0,a,a,1],[a,0,a,1]])
mat267 = np.array([[1,1,1,1],[1,0,0,1],[a,0,a,1],[a,a,0,1]])
mat357 = np.array([[1,1,1,1],[0,1,0,1],[0,a,a,1],[a,a,0,1]])


def minor(arr,i,j):
    # ith row, jth column removed
    return arr[np.array(range(i)+range(i+1,arr.shape[0]))[:,np.newaxis],
               np.array(range(j)+range(j+1,arr.shape[1]))]


def plane(theMat):
	a = np.linalg.det(minor(theMat,0,0))
	b = -np.linalg.det(minor(theMat,0,1))
	c = np.linalg.det(minor(theMat,0,2))
	d = -np.linalg.det(minor(theMat,0,3))
	return a,b,c,d

	
def intersect(plane1,plane2,plane3):
	A = np.array([plane1[:-1],plane2[:-1],plane3[:-1]])
	b = np.array([[-plane1[-1]],[-plane2[-1]],[-plane3[-1]]])
	Ainv = np.linalg.inv(A)
	return CAT(np.dot(Ainv,b).tolist())

x,y,z = intersect(plane(mat456),plane(mat267),plane(mat357))

v = g.addNode(0); g.setVecf(v,Vecf(1,x,y,z))

DRAW(g)()