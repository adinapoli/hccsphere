# -*- coding: utf-8 -*-
from graphLib6 import *
""" Thihedral decomposition of the Schönhardt polyhedron. """

#/////////////////////////////////////////////////
# Graph of a Hpc object (simple, 3D)
#/////////////////////////////////////////////////

def GRAPH(hpcGraph):
    g1 = Graph(hpcGraph.getPointDim())
    d = hpcGraph.getMaxDimCells()
    nodes = []

    for k in range(d+1):
        nodes += [{}]
        for node in CELLSPERLEVEL(hpcGraph)(k):
            newNode = g1.addNode(k)
            nodes[k].update({node:newNode})            
            if k==0 or k==2:
                g1.setVecf(newNode,hpcGraph.getVecf(node))
        if k>0:
            for node in CELLSPERLEVEL(hpcGraph)(k-1):
                for father in UPCELLS(hpcGraph)(node):                
                    g1.addArch(nodes[k-1][node],nodes[k][father])
        for node in CELLSPERLEVEL(g1)(d):
            verts = list(set(CAT(AA(DOWNCELLS(g1))(
                list(set(CAT(AA(DOWNCELLS(g1))(DOWNCELLS(g1)(node)))))))))
            [g1.addArch(node,vert) for vert in verts]
    return g1


#/////////////////////////////////////////////////
# Some utilities
#/////////////////////////////////////////////////

def cyclicPairs(seq):
    return TRANS([seq,seq[1:]+[seq[0]]])

def combinations(Set):
    def combinations0(n):
        return [(i[1],j[1]) for i in enumerate(Set)
                for j in enumerate(Set) if j[0]>i[0]]
    return combinations0


#/////////////////////////////////////////////////
# Schönhardt's kernel construction
#/////////////////////////////////////////////////

def SchonhardtKernel(vertices):
    
    g=Graph(3)
    
    V2,V10,V20,V01,V11,V21 = vertices
    
    v2=g.addNode(0);g.setVecf(v2,V2)
    v10=g.addNode(0);g.setVecf(v10,V10)
    v20=g.addNode(0);g.setVecf(v20,V20)
    v01=g.addNode(0);g.setVecf(v01,V01)
    v11=g.addNode(0);g.setVecf(v11,V11)
    v21=g.addNode(0);g.setVecf(v21,V21)

    triangles = CELLSPERLEVEL(g)(2)
    triaEdges = [DOWNCELLS(g)(triangle) for triangle in triangles]
    triaVerts = [list(set(CAT(AA(DOWNCELLS(g))(triangle))))
                 for triangle in triaEdges]
    cellVerts = [[1, 2, 5,6], [1, 4, 5,3], [2, 3, 6,4],
                 [2, 5, 6,1], [1, 3, 4,5], [3, 4, 6,2]]
    for cell in cellVerts:
        c = g.addNode(3)
        for vert in cell: g.addArch(c,vert)

    top = Vecf(1,0,0,0.2679491)
    bottom = Vecf(1,0,0,-0.2679491)

    p,v = range(8),range(8)
    for i in range(1,7):
        p[i] = Graph(3); v[i]=p[i].addNode(0); p[i].setVecf(v[i],g.getVecf(i))
    p[0]=Graph(3); v[0]=p[0].addNode(0); p[0].setVecf(v[0],top)
    p[7]=Graph(3); v[7]=p[7].addNode(0); p[7].setVecf(v[7],bottom)

    Top = JOIN(AA(Hpc)([p[0],p[1],p[2],p[3]]))
    Bottom = JOIN(AA(Hpc)([p[4],p[5],p[6],p[7]]))
    boundary = SKELETON(1)(Hpc(g))

    solidcells = (AA)(COMP([JOIN,AA(Hpc)]))(
        [[p[1],p[2],p[5],p[6]],[p[1],p[4],p[5],p[3]],[p[2],p[3],p[6],p[4]],
         [p[2],p[5],p[6],p[1]],[p[1],p[3],p[4],p[5]],[p[3],p[4],p[6],p[2]] ])

    kernel = INTERSECTION(solidcells)
    #VIEW(STRUCT([Top,kernel,Bottom,boundary]))
    return kernel


#/////////////////////////////////////////////////
# Graph of a polytopal decomposition of Schönhardt's polyhedron
#/////////////////////////////////////////////////

##def Schonhardt(angle):

angle = PI/6

V2 = Vecf(1,COS(0),SIN(0),1)
V10 = Vecf(1,COS(2*PI/3),SIN(2*PI/3),1)
V20 = Vecf(1,COS(4*PI/3),SIN(4*PI/3),1)

V01 = Vecf(1,COS(0+PI/6),SIN(0+PI/6),-1)
V11 = Vecf(1,COS(2*PI/3 + angle),SIN(2*PI/3 + angle),-1)
V21 = Vecf(1,COS(4*PI/3 + angle),SIN(4*PI/3 + angle),-1)

vertices = [V2,V10,V20,V01,V11,V21]
kernel = SchonhardtKernel(vertices)
g = GRAPH(kernel.g)
kernelEdges = CELLSPERLEVEL(g)(1)

def vertsPerFace(g):
    def vertsPerFace0(faceChain):
        edgesPerFace = AA(DOWNCELLS(g))(faceChain)
        vertsPerEdgePerFace = AA(AA(DOWNCELLS(g)))(edgesPerFace)
        vertsPerFace = AA(COMP([list,set,CAT]))(vertsPerEdgePerFace)
        return vertsPerFace,edgesPerFace,vertsPerEdgePerFace
    return vertsPerFace0

def cap(verts,v):
        return CAT([GETINTERSECTION(g)(getEdges(pair,v))
                for pair in cyclicPairs(verts)])

def topBottom(verts):
    faces = cap(verts[:-1],verts[-1])
    edges = CAT([GETINTERSECTION(g)([u,w])
                 for u,w in cyclicPairs(verts[:-1])])
    face = g.addNode(2);
    [g.addArch(e,face) for e in edges]
    solid = g.addNode(3);
    [g.addArch(f,solid) for f in faces+[face]]
    links = [g.addArch(solid,v) for v in verts]

def getEdges(pair,v):
    u,w = pair
    e1 = GETINTERSECTION(g)([u,w])
    e2 = GETINTERSECTION(g)([u,v])
    e3 = GETINTERSECTION(g)([w,v])
    return CAT([e1,e2,e3])


##    verts = []
##    for k in range(1,9):
##        verts += [k]
##        DRAW(g)(CELLSPERLEVEL(g)(1)+verts)



# nodes 6,7 are kernel's top and bottom, respectively
""" top face vertices """
v2=g.addNode(0);g.setVecf(v2,V2)
v10=g.addNode(0);g.setVecf(v10,V10)
v20=g.addNode(0);g.setVecf(v20,V20)
v5 = g.addNode(0); g.setVecf(v5,(V2+V10)/2)
v120 = g.addNode(0); g.setVecf(v120,(V10+V20)/2)
v3 = g.addNode(0); g.setVecf(v3,(V20+V2)/2)
v6 = g.addNode(0); g.setVecf(v6,(V20+V10+V2)/3)
""" bottom face vertices """
v01=g.addNode(0);g.setVecf(v01,V01)
v11=g.addNode(0);g.setVecf(v11,V11)
v21=g.addNode(0);g.setVecf(v21,V21)
v011 = g.addNode(0); g.setVecf(v011,(V01+V11)/2)
v121 = g.addNode(0); g.setVecf(v121,(V11+V21)/2)
v201 = g.addNode(0); g.setVecf(v201,(V21+V01)/2)
v0121 = g.addNode(0); g.setVecf(v0121,(V21+V11+V01)/3)

e = g.addNode(1); g.addArch(4,e); g.addArch(v20,e);
##    e = g.addNode(1); g.addArch(4,e); g.addArch(v10,e);
##    e = g.addNode(1); g.addArch(2,e); g.addArch(v20,e);

DRAW(g)(CELLSPERLEVEL(g)(1) 
##        +[6,3,2,5,v2,v5,v3,v6]+[v10,v20]
            +[6,1,4,3,v20]+[v2,v3,v6]
##            +[6,5,8,1,v10]
##            +[7,2,3,4,v21]
##            +[7,4,1,8,v11]
##            +[7,8,2,5,v01]
        )

##    return g


##g = Schonhardt(PI/6)
