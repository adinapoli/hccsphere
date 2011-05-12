""" To generate the initial approximation of sphere with hexahedra """

import numpy as np
from hasselib.graphLib6 import *

#/////////////////////////////////////////////////////////////////////
# Utility functions 
#/////////////////////////////////////////////////////////////////////
def minor(arr,i,j):
    # ith row, jth column removed
    return arr[np.array(range(i)+range(i+1,arr.shape[0]))[:,np.newaxis],
               np.array(range(j)+range(j+1,arr.shape[1]))]

def plane(theMat):
	""" Compute the coefficients of the plane equations in 3D """
	d = np.linalg.det(minor(theMat,0,0))
	a = -np.linalg.det(minor(theMat,0,1))
	b = np.linalg.det(minor(theMat,0,2))
	c = -np.linalg.det(minor(theMat,0,3))
	return d,a,b,c
	
def intersect(plane1,plane2,plane3):
	""" Intersection point of 3 planes in 3D -- first coord homogeneous """
	A = np.array([plane1[1:],plane2[1:],plane3[1:]])
	b = np.array([[-plane1[0]],[-plane2[0]],[-plane3[0]]])
	A_inv = np.linalg.inv(A)
	x0, x1,x2,x3 = CAT([[1.0]] + np.dot(A_inv,b).tolist())
	return [x0, x1,x2,x3]

def reflect(points,i):	
	ret,d = [],len(points[0])
	for point in points:
		p=[]
		for k in range(d):
			if k==i: p += [-point[k]]
			else: p += [point[k]]
		ret += [p]
	return ret
		
def classifyVerts(g):
	verts,vclass = CELLSPERLEVEL(g)(0),[[],[],[],[]]
	for v in verts:
		vert = eval(str(g.getVecf(v)))[1:]
		for k in range(4):
			if len([x for x in vert if x==0.0])==k:
				vclass[k] += [v]
				continue				
	return vclass

def signature(vert):
	point = eval(str(g.getVecf(vert)))[1:]
	ret = []
	for x in point: ret += [SIGN(x) if x!=0 else 0]
	return ret

def test(v,w): 
	return INNERPROD([signature(v),signature(w)])

def matching(vclass1,signtr):		
	return [v for v in vclass1 if len([k for k in range(3) 
		if signature(v)[k]==signtr[k]]) == 2]

def cyclicPairs(seq):
    return TRANS([seq,seq[1:]+[seq[0]]])

def initialSphere(g):

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
	verts = AA(eval)(list(set(AA(repr)(verts))))
	for vert in verts:
		v = g.addNode(0); g.setVecf(v,Vecf(vert))

	def cell(args):
		i,j,k = args
		return [index+1 for index,v in enumerate(verts) if AND([ 
			(GE(0) if i==1 else LE(0))(v[1]), 
			(GE(0) if j==1 else LE(0))(v[2]), 
			(GE(0) if k==1 else LE(0))(v[3]) ])  ]
			
	vertices = [v[1:] for v in verts]
	pols = AA(cell)([[1,1,1],[-1,1,1],[1,-1,1],[1,1,-1],[-1,-1,1],[-1,1,-1],[1,-1,-1],[-1,-1,-1]])
	out = MKPOL([vertices, pols, None])

	for pol in pols:
		node = g.addNode(3); 
		for vert in pol: g.addArch(node,vert)
	
	vclass = classifyVerts(g)

	edges = [[v,vclass[3][0]] for v in vclass[2]]
	edges += CAT([DISTL([v,[w for w in vclass[2] if test(v,w)==1]]) for v in vclass[1]])
	edges += CAT([DISTL([v,[w for w in vclass[1] if test(v,w)==2]]) for v in vclass[0]])

	for edge in edges:
		 node = g.addNode(1); g.addArch(edge[0],node);g.addArch(edge[1],node)
			
	faces = CAT([[[v,u,w,vclass[3][0]] for [u,w] in CART([vclass[2],vclass[2]]) if u<w and VECTSUM([
		signature(u),signature(w) ])==signature(v)] for v in vclass[1]])
	
	faces = []
	o = vclass[3][0]
	for v in vclass[1]:
		for [u,w] in CART([vclass[2],vclass[2]]):
			if u<w and VECTSUM([ signature(u),signature(w) ])==signature(v):
				faces += [CAT(AA(GETINTERSECTION(g))([[v,u],[v,w],[o,u],[o,w]]))]

	for vert in vclass[0]:
		u,v,w = matching(vclass[1],signature(vert))
		vertTriples = [[u,v,vert],[v,w,vert],[w,u,vert]]
		for k in range(3):
			signatures = AA(signature)(vertTriples[k])
			common = COMP([AA(PROD),TRANS])(signatures)
			w = [v for v in vclass[2] if signature(v) == common]+vertTriples[k]
			faces += [CAT(AA(GETINTERSECTION(g))([[w[0],w[1]],[w[0],w[2]],[w[1],w[3]],[w[2],w[3]]]))]
		
	for face in faces:
		node = g.addNode(2) 
		for k in range(4): g.addArch(face[k],node)
		
	return g

g = Graph(3)
g = initialSphere(g)
DRAW(g,[1.5,1.5,1.5])()	
