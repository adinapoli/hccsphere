from graphLib6 import *


pointdim=3
g=Graph(pointdim)

V00=Vecf(1,COS(0),SIN(0),1)
V10=Vecf(1,COS(2*PI/3),SIN(2*PI/3),1)
V20=Vecf(1,COS(4*PI/3),SIN(4*PI/3),1)
V01=Vecf(1,COS(0+PI/6),SIN(0+PI/6),-1)
V11=Vecf(1,COS(2*PI/3+PI/6),SIN(2*PI/3+PI/6),-1)
V21=Vecf(1,COS(4*PI/3+PI/6),SIN(4*PI/3+PI/6),-1)

V02=(V00+V01)/2
V12=(V10+V11)/2
V22=(V20+V21)/2
V03=(V00+V11)/2
V13=(V10+V21)/2
V23=(V20+V01)/2


v00=g.addNode(0);g.setVecf(v00,V00)
v10=g.addNode(0);g.setVecf(v10,V10)
v20=g.addNode(0);g.setVecf(v20,V20)
v01=g.addNode(0);g.setVecf(v01,V01)
v11=g.addNode(0);g.setVecf(v11,V11)
v21=g.addNode(0);g.setVecf(v21,V21)
    
##    v02=g.addNode(0);g.setVecf(v02,V02)
##    v12=g.addNode(0);g.setVecf(v12,V12)
##    v22=g.addNode(0);g.setVecf(v22,V22)
##    v03=g.addNode(0);g.setVecf(v03,V03)
##    v13=g.addNode(0);g.setVecf(v13,V13)
##    v23=g.addNode(0);g.setVecf(v23,V23)

e010=g.addNode(1);g.addArch(v00,e010);g.addArch(v10,e010)
e120=g.addNode(1);g.addArch(v10,e120);g.addArch(v20,e120)
e200=g.addNode(1);g.addArch(v20,e200);g.addArch(v00,e200)

e011=g.addNode(1);g.addArch(v01,e011);g.addArch(v11,e011)
e121=g.addNode(1);g.addArch(v11,e121);g.addArch(v21,e121)
e201=g.addNode(1);g.addArch(v21,e201);g.addArch(v01,e201)

e012=g.addNode(1);g.addArch(v00,e012);g.addArch(v11,e012)
e122=g.addNode(1);g.addArch(v10,e122);g.addArch(v21,e122)
e202=g.addNode(1);g.addArch(v20,e202);g.addArch(v01,e202)

e013=g.addNode(1);g.addArch(v00,e013);g.addArch(v01,e013)
e123=g.addNode(1);g.addArch(v10,e123);g.addArch(v11,e123)
e203=g.addNode(1);g.addArch(v20,e203);g.addArch(v21,e203)

##    e0140=g.addNode(1);g.addArch(v02,e0140);g.addArch(v03,e0140)
##    e0141=g.addNode(1);g.addArch(v03,e0141);g.addArch(v12,e0141)
##    e0142=g.addNode(1);g.addArch(v12,e0142);g.addArch(v13,e0142)
##    e0143=g.addNode(1);g.addArch(v13,e0143);g.addArch(v22,e0143)
##    e0144=g.addNode(1);g.addArch(v22,e0144);g.addArch(v23,e0144)
##    e0145=g.addNode(1);g.addArch(v23,e0145);g.addArch(v02,e0145)
    
##    f010=g.addNode(2);
##    g.addArch(e010,f010);g.addArch(e120,f010);g.addArch(e200,f010)
##    f011=g.addNode(2);
##    g.addArch(e011,f011);g.addArch(e121,f011);g.addArch(e201,f011)

f012=g.addNode(2);
g.addArch(e010,f012);g.addArch(e012,f012);g.addArch(e123,f012)
f013=g.addNode(2);
g.addArch(e011,f013);g.addArch(e012,f013);g.addArch(e013,f013)
f014=g.addNode(2);
g.addArch(e120,f014);g.addArch(e122,f014);g.addArch(e203,f014)
f015=g.addNode(2);
g.addArch(e121,f015);g.addArch(e122,f015);g.addArch(e123,f015)
f016=g.addNode(2);
g.addArch(e200,f016);g.addArch(e202,f016);g.addArch(e013,f016)
f017=g.addNode(2);
g.addArch(e201,f017);g.addArch(e202,f017);g.addArch(e203,f017)

##    f018=g.addNode(2);
##    g.addArch(e0140,f018);g.addArch(e0141,f018);g.addArch(e0142,f018)
##    g.addArch(e0143,f018);g.addArch(e0144,f018);g.addArch(e0145,f018)

##f = g.addNode(2)
##for e in [1, 2, 5]: g.addArch(e,f)
##f = g.addNode(2)
##for e in [2, 3, 6]: g.addArch(e,f)
##f = g.addNode(2)
##for e in [1, 3, 4]: g.addArch(e,f)
##
##f = g.addNode(2)
##for e in [1, 4, 5]: g.addArch(e,f)
##f = g.addNode(2)
##for e in [2, 5, 6]: g.addArch(e,f)
##f = g.addNode(2)
##for e in [3, 4, 6]: g.addArch(e,f)

if True:
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


    p1=Graph(3); v1=p1.addNode(0); p1.setVecf(v1,g.getVecf(1))
    p2=Graph(3); v2=p2.addNode(0); p2.setVecf(v2,g.getVecf(2))
    p3=Graph(3); v3=p3.addNode(0); p3.setVecf(v3,g.getVecf(3))
    p4=Graph(3); v4=p4.addNode(0); p4.setVecf(v4,g.getVecf(4))
    p5=Graph(3); v5=p5.addNode(0); p5.setVecf(v5,g.getVecf(5))
    p6=Graph(3); v6=p6.addNode(0); p6.setVecf(v6,g.getVecf(6))

    p0=Graph(3); v0=p0.addNode(0); p0.setVecf(v0,top)
    p7=Graph(3); v7=p7.addNode(0); p7.setVecf(v7,bottom)
    
    Top = JOIN(AA(Hpc)([p0,p1,p2,p3]))
    Bottom = JOIN(AA(Hpc)([p4,p5,p6,p7]))
    boundary = SKELETON(1)(Hpc(g))

    solidcells = (AA)(COMP([JOIN,AA(Hpc)]))(
        [[p1,p2,p5,p6],[p1,p4,p5,p3],[p2,p3,p6,p4],
         [p2,p5,p6,p1],[p1,p3,p4,p5],[p3,p4,p6,p2]])

    kernel = INTERSECTION(solidcells)
    other = DIFFERENCE([UNION(solidcells),INTERSECTION(solidcells)])
    VIEW(STRUCT([Top,kernel,Bottom,boundary]))
    VIEW(other)

    top = Vecf(1,0,0,0.2679491)
    bottom = Vecf(1,0,0,-0.2679491)
    


DRAW(g)()
DRAW(g,[1.5,1.5,1.5])()


[[0.0, 0.26794919, -0.071796879],
 [0.0, 0.0, -0.26794919],
 [0.0, 0.0, 0.26794919],
 [0.1339747, 0.232050687, 0.07179671],
 [-0.2320508, -0.13397458, -0.0717967599],
 [0.1339747, -0.23205079, 0.07179671],
 [0.2320508, -0.1339747, -0.07179671],
 [-0.26794919, 0.0, 0.07179671]]


hpcGraph=kernel.g

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

e012=h.addNode(1);h.addArch(v00,e012);h.addArch(v11,e012)
e122=h.addNode(1);h.addArch(v10,e122);h.addArch(v21,e122)
e202=h.addNode(1);h.addArch(v20,e202);h.addArch(v01,e202)

e013=h.addNode(1);h.addArch(v00,e013);h.addArch(v01,e013)
e123=h.addNode(1);h.addArch(v10,e123);h.addArch(v11,e123)
e203=h.addNode(1);h.addArch(v20,e203);h.addArch(v21,e203)


##t1 = h.addNode(3);
##for v in [1,8,v00,v10]: h.addArch(f11,c1)


DRAW(h)()

