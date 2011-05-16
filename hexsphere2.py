""" To generate the initial approximation of the sphere with hexahedra """

import numpy as np
from hasselib.graphLib6 import *
from utils import *

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

def minus(list1):
    def minus0(list2):
        return list(set(list2).difference(list1))
    return minus0

def intrsz(list1):
    def intrsz0(list2):
        return list(set(list2).intersection(list1))
    return intrsz0


######################################################################
#/////////////////////////////////////////////////////////////////////
# Sphere generation 
#/////////////////////////////////////////////////////////////////////
def initialSphere(g):
    """ Generation of an approximation of the standard unit sphere with 8 hexahedra.
        Return a Graph instance.
    """
    #//////////////////////////////////////////////////////////
    # generation of vertex nodes
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

    vclass = classifyVerts(g)


    #//////////////////////////////////////////////////////////
    # generation of edge nodes

    edges = [[v,vclass[3][0]] for v in vclass[2]]
    edges += CAT([DISTL([v,[w for w in vclass[2] if test(v,w)==1]]) for v in vclass[1]])
    edges += CAT([DISTL([v,[w for w in vclass[1] if test(v,w)==2]]) for v in vclass[0]])

    for edge in edges:
         node = g.addNode(1); g.addArch(edge[0],node);g.addArch(edge[1],node)

    #//////////////////////////////////////////////////////////
    # generation of face nodes

    faces,o = [],vclass[3][0]

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
            faceEdgeVertices = [[w[0],w[1]],[w[0],w[2]],[w[1],w[3]],[w[2],w[3]]]
            faceEdges = CAT(AA(GETINTERSECTION(g))(faceEdgeVertices))
            faces += [faceEdges]

    for face in faces:
        node = g.addNode(2)
        for k in range(4): g.addArch(face[k],node)

    #//////////////////////////////////////////////////////////
    # generation of down arcs of 3D cell nodes

    for vert in vclass[0]:
        edges = UPCELLS(g)(vert)
        faces = list(set(CAT(AA(UPCELLS(g))(edges))))
        edgepairs = AA(DOWNCELLS(g))(faces)
        edgepairs = CAT(AA(minus(edges))(edgepairs))

        def int_test(pair):
            out = intrsz(DOWNCELLS(g)(pair[0]))(DOWNCELLS(g)(pair[1]))
            if len(out) != 1: return 0
            else: return out[0]

        epairs = [pair for pair in CART([edgepairs,edgepairs]) if int_test(pair) in vclass[1]]
        faces += CAT(AA(GETINTERSECTION(g))(epairs))

        node = g.addNode(3)
        verts = list(set(CAT(AA(DOWNCELLS(g))(edges+edgepairs))))
        for vert in verts: g.addArch(node,vert)
        for face in faces: g.addArch(face,node)

    return g


#/////////////////////////////////////////////////////////////////////
# Layered sphere 
#/////////////////////////////////////////////////////////////////////


def   hexSphere(g,scaling=2):

    #/////////////////////////////////////////////////////////////////////
    #   STEP 1:  d-1  boundary-complex extraction
    #/////////////////////////////////////////////////////////////////////
    bComplex = boundaryComplex(g)
    n = g.getPointDim()
    DRAW(g,[1.5,1.5,1.5])(CAT(bComplex))
    

    #/////////////////////////////////////////////////////////////////////
    # STEP 2: d-2 wall-complex extraction 
    #/////////////////////////////////////////////////////////////////////
    wallComplex = bComplex[:-1]
    DRAW(g,[1.5,1.5,1.5])(CAT(wallComplex))

    # STEP 3: basket separators construction \\\\\\\\\\\\\\\\
    # (d-1)-complex construction
    mapping = {}

    #Compute the parallel edges
    planes = [[1/SQRT(2), 1/SQRT(2), 0], [-1/SQRT(2), 1/SQRT(2), 0],
              [1, 0, 0], [0, 1, 0]]

    meridian_vtx = [meridians(g)(plane)[0:10] for plane in planes]
    meridian_edges = CAT([meridians_edges(g)(cells) for cells in meridian_vtx])

    DRAW(g, [1.5, 1.5, 1.5])(meridian_edges)

    #duplicate wall-complex
    for skeleton in wallComplex:
        
        #add 0-cells of mapped wall-complex cells
        for cell in skeleton:

            if cell in meridian_edges or cell in wallComplex[0]:
                newCell = g.addNode(0)
                mapping.update({cell:newCell})
                point = [CENTROID(g)(cell).get(i) for i in range(1,n+1)]
                point = SCALARVECTPROD([UNITVECT(point),scaling])
                g.setVecf(newCell,Vecf([1.0]+point))

    DRAW(g, [1.5, 1.5, 1.5])()

    #add 1-cell extensions to top corners
    for node in wallComplex[0]:
        newNode = mapping[node]
        newArc = g.addNode(1)
        g.addArch(node,newArc)
        g.addArch(newNode,newArc)


    #add 1-cells of top-lateral boundary of added polytopes
    para = meridians(g)([0.707, 0.707, 0])
    print para
    DRAW(g)(para[0:8])
    newBoundary0 = zip(para, para[1:] + [para[0]])
    for pair in newBoundary0:
        newArc = g.addNode(1)
        g.addArch(pair[0],newArc)
        g.addArch(pair[1],newArc)
    DRAW(g)()


#/////////////////////////////////////////////////////////////////////
# Local testing 
#/////////////////////////////////////////////////////////////////////

if __name__=="__main__":

    g = Graph(3)
    g = initialSphere(g)
    DRAW(g,[1.5,1.5,1.5])()
    hexSphere(g)
