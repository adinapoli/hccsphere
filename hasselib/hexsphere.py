import numpy as np
from test import *
from graphLib6 import *

def minor(arr,i,j):
    # ith row, jth column removed
    return arr[np.array(range(i)+range(i+1,arr.shape[0]))[:,np.newaxis],
               np.array(range(j)+range(j+1,arr.shape[1]))]


def plane(theMat):
	d = np.linalg.det(minor(theMat,0,0))
	a = -np.linalg.det(minor(theMat,0,1))
	b = np.linalg.det(minor(theMat,0,2))
	c = -np.linalg.det(minor(theMat,0,3))
	return d,a,b,c

	
def intersect(plane1,plane2,plane3):
	A = np.array([plane1[1:],plane2[1:],plane3[1:]])
	b = np.array([[-plane1[0]],[-plane2[0]],[-plane3[0]]])
	A_inv = np.linalg.inv(A)
	return CAT([[1.0]] + np.dot(A_inv,b).tolist())


def reflect(points,i):	
	ret,d = [],len(points[0])
	for point in points:
		p=[]
		for k in range(d):
			if k==i: p += [-point[k]]
			else: p += [point[k]]
		ret += [p]
	return ret

def cell(args):
	i,j,k = args
	return [index+1 for index,v in enumerate(verts) if AND([ 
		(GE(0) if i==1 else LE(0))(v[1]), 
		(GE(0) if j==1 else LE(0))(v[2]), 
		(GE(0) if k==1 else LE(0))(v[3]) ])  ]

g = Graph(3)

a = 1/SQRT(2)
i,j,k = [[1,1,0,0],[1,0,1,0],[1,0,0,1]]
ij,jk,ki = [[1,a,a,0],[1,0,a,a],[1,a,0,a]]
zero = [1,0,0,0]
mat_ij = np.array([zero,ki,i,ij])
mat_jk = np.array([zero,ij,j,jk])
mat_ki = np.array([zero,jk,k,ki])
ijk = intersect(plane(mat_ij),plane(mat_jk),plane(mat_ki))
verts = [zero,i,j,k,ij,jk,ki,ijk]
for k in range(1,4):
	verts += reflect(verts,k)
myprint("g",g)
myprint("len(verts)",len(verts))
verts = AA(eval)(list(set(AA(repr)(verts))))
myprint("len(verts)",len(verts))
for vert in verts:
	v = g.addNode(0); g.setVecf(v,Vecf(vert))
	
myprint("g",g)
DRAW(g)()

		
vertices = [v[1:] for v in verts]
pols = AA(cell)([[1,1,1],[-1,1,1],[1,-1,1],[1,1,-1],[-1,-1,1],[-1,1,-1],[1,-1,-1],[-1,-1,-1]])
out = MKPOL([vertices, pols, None])

myprint("pols",pols)

for cell in pols:
	node = g.addNode(3); 
	for vert in cell: g.addArch(node,vert)
	
	
edges = [[7,16],[21,16],[5,16],[14,16],[17,16],[25,16],
[7,11],[11,21],[21,15],[15,5],[5,0],[0,7], [14,9],[9,17],[17,6],[6,25],[25,4],[4,14],
[21,2],[2,25],[5,26],[26,14],[7,18],[18,17],[17,8],[8,5],[14,19],[19,21],[25,23],[23,7],
[1,6],[1,18],[1,23],[3,0],[3,11],[3,15],[10,2],[10,11],[10,23],[12,15],[12,19],[12,26],
[13,8],[13,9],[13,26],[20,0],[20,8],[20,18],[22,4],[22,6],[22,9],[24,2],[24,4],[24,19]]
for edge in edges:
	node = g.addNode(1); g.addArch(edge[0]+1,node);g.addArch(edge[1]+1,node)

#	CELLSPERLEVEL(g)(1) = [36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89]
#							1	2	3	4	5	6	7	8	9	10	11	12

#DRAW(g)([36, 42, 43, 37, 40, 41, 42, 43, 44, 45, 46, 47])
DRAW(g)([1, 2, 3, 4, 8, 12, 16, 20])

faces = [[36,42,43,37],[37,44,45,38],[36,38,46,47]]
for face in faces:
	node = g.addNode(2) 
	for k in range(4): g.addArch(face[k],node)

DRAW(g,[1.5,1.5,1.5])()

VIEW(out)
VIEW(SKELETON(1)(out))

myprint("CELLSPERLEVEL(g)(0)",CELLSPERLEVEL(g)(0))
myprint("CELLSPERLEVEL(g)(1)",CELLSPERLEVEL(g)(1))
myprint("CELLSPERLEVEL(g)(2)",CELLSPERLEVEL(g)(2))
myprint("CELLSPERLEVEL(g)(3)",CELLSPERLEVEL(g)(3))

# DRAW(g)(CELLSPERLEVEL(g)(0))
# DRAW(g)(CELLSPERLEVEL(g)(1))
# DRAW(g)(CELLSPERLEVEL(g)(2))
# DRAW(g)(CELLSPERLEVEL(g)(3))

myprint("g",g)
