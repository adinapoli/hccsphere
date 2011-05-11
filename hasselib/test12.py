# -*- coding: utf-8 -*-
from graphLib6 import *
""" Construct the Schönhardt polyhedron and some decompositions of it. """

#/////////////////////////////////////////////////
# Schönhardt's kernel construction
#/////////////////////////////////////////////////
g=Graph(3)

V00=Vecf(1,COS(0),SIN(0),1)
V10=Vecf(1,COS(2*PI/3),SIN(2*PI/3),1)
V20=Vecf(1,COS(4*PI/3),SIN(4*PI/3),1)
V01=Vecf(1,COS(0+PI/6),SIN(0+PI/6),-1)
V11=Vecf(1,COS(2*PI/3+PI/6),SIN(2*PI/3+PI/6),-1)
V21=Vecf(1,COS(4*PI/3+PI/6),SIN(4*PI/3+PI/6),-1)

v00=g.addNode(0);g.setVecf(v00,V00)
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
VIEW(STRUCT([Top,kernel,Bottom,boundary]))

#/////////////////////////////////////////////////
# Graph of an HPC (simple, 3D) value
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
            myprint("node",node)
            verts = list(set(CAT(AA(DOWNCELLS(g1))(
                list(set(CAT(AA(DOWNCELLS(g1))(DOWNCELLS(g1)(node)))))))))
            [g1.addArch(node,vert) for vert in verts]
    return g1

h = GRAPH(kernel.g)

#/////////////////////////////////////////////////
# Graph of a polytopal decomposition
#/////////////////////////////////////////////////

# nodes 6,7 are kernel's top and bottom, respectively
""" top face vertices """
v00=h.addNode(0);h.setVecf(v00,V00)
v10=h.addNode(0);h.setVecf(v10,V10)
v20=h.addNode(0);h.setVecf(v20,V20)
""" bottom face vertices """
v01=h.addNode(0);h.setVecf(v01,V01)
v11=h.addNode(0);h.setVecf(v11,V11)
v21=h.addNode(0);h.setVecf(v21,V21)

""" top face edges """
e010=h.addNode(1);h.addArch(v00,e010);h.addArch(v10,e010)
e120=h.addNode(1);h.addArch(v10,e120);h.addArch(v20,e120)
e200=h.addNode(1);h.addArch(v20,e200);h.addArch(v00,e200)
""" top lateral faces' edges """
e016=h.addNode(1);h.addArch(6,e016);h.addArch(v10,e016)
e126=h.addNode(1);h.addArch(6,e126);h.addArch(v20,e126)
e206=h.addNode(1);h.addArch(6,e206);h.addArch(v00,e206)
myprint("CELLSPERLEVEL(h)(1)",CELLSPERLEVEL(h)(1))

""" bottom face edges """
e011=h.addNode(1);h.addArch(v01,e011);h.addArch(v11,e011)
e121=h.addNode(1);h.addArch(v11,e121);h.addArch(v21,e121)
e201=h.addNode(1);h.addArch(v21,e201);h.addArch(v01,e201)
""" bottom lateral faces' edges """
e017=h.addNode(1);h.addArch(7,e017);h.addArch(v11,e017)
e127=h.addNode(1);h.addArch(7,e127);h.addArch(v21,e127)
e207=h.addNode(1);h.addArch(7,e207);h.addArch(v01,e207)
myprint("CELLSPERLEVEL(h)(1)",CELLSPERLEVEL(h)(1))

""" top faces """
f0 = h.addNode(2);h.addArch(34,f0);h.addArch(35,f0);h.addArch(36,f0)
f01 = h.addNode(2);h.addArch(37,f01);h.addArch(38,f01);h.addArch(34,f01)
f02 = h.addNode(2);h.addArch(38,f02);h.addArch(39,f02);h.addArch(35,f02)
f03 = h.addNode(2);h.addArch(39,f03);h.addArch(37,f03);h.addArch(36,f03)
myprint("CELLSPERLEVEL(h)(2)",CELLSPERLEVEL(h)(2))

""" bottom faces """
f1 = h.addNode(2);h.addArch(40,f1);h.addArch(41,f1);h.addArch(42,f1)
f11 = h.addNode(2);h.addArch(43,f11);h.addArch(44,f11);h.addArch(40,f11)
f12 = h.addNode(2);h.addArch(44,f12);h.addArch(45,f12);h.addArch(41,f12)
f13 = h.addNode(2);h.addArch(45,f13);h.addArch(43,f13);h.addArch(42,f13)
myprint("CELLSPERLEVEL(h)(2)",CELLSPERLEVEL(h)(2))

""" top solid """
c0 = h.addNode(3);
h.addArch(f0,c0);h.addArch(f01,c0);h.addArch(f02,c0);h.addArch(f03,c0)
h.addArch(c0,v00);h.addArch(c0,v10);h.addArch(c0,v20);h.addArch(c0,6)
""" bottom solid """
c1 = h.addNode(3);
h.addArch(f1,c1);h.addArch(f11,c1);h.addArch(f12,c1);h.addArch(f13,c1)
h.addArch(c1,v01);h.addArch(c1,v11);h.addArch(c1,v21);h.addArch(c1,7)

""" boundary "diagonal" lines """
e012=h.addNode(1);h.addArch(v00,e012);h.addArch(v11,e012)
e122=h.addNode(1);h.addArch(v10,e122);h.addArch(v21,e122)
e202=h.addNode(1);h.addArch(v20,e202);h.addArch(v01,e202)
""" boundary "vertical" lines """
e013=h.addNode(1);h.addArch(v00,e013);h.addArch(v01,e013)
e123=h.addNode(1);h.addArch(v10,e123);h.addArch(v11,e123)
e203=h.addNode(1);h.addArch(v20,e203);h.addArch(v21,e203)

def cyclicPairs(seq):
    return TRANS([seq,seq[1:]+[seq[0]]])
def combinations(Set):
    def combinations0(n):
        return [(i[1],j[1]) for i in enumerate(Set)
                for j in enumerate(Set) if j[0]>i[0]]
    return combinations0

edges = []
for v,w in combinations([1,8,v00,v10])(2):
    e = h.addNode(1); h.addArch(v,e); h.addArch(w,e)
    edges += [e]
for ei,ej,ek in [[edges[k],edges[k+1],edges[k+2]] for k in range(len(edges)-2)]:
    f = h.addNode(2); h.addArch(ei,f);h.addArch(ej,f);h.addArch(ek,f);


DRAW(h)()

