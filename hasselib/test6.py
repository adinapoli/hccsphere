from test5 import *
from numpy import arange
from math import pi,sin,cos

def ngon(n):
    points = [[cos(u), sin(u)] for u in arange(0,2*pi,2*pi/n)]    
    g = Graph(2)
    for i in range(n):          # add 0-cells
        vert = g.addNode(0);
        g.setVecf(vert,Vecf([1.0]+points[i]))        
    for i in range(n):          # add 1-cells
        edge = g.addNode(1);
        g.addArch(i+1,edge);
        g.addArch((i+2 if i<n-1 else 1),edge);        
    face = g.addNode(2)         # add 2-cells
    for i in range(1,n+1):
        g.addArch(n+i,face)
        g.addArch(face,i)
    return g


def boundaryChain(g):
    chain = CELLSPERLEVEL(g)(g.getMaxDimCells())
    return [k for k in CAT(AA(DOWNCELLS(g))(chain)) if len(UPCELLS(g)(k))==1]

def boundaryComplex(g):
    chain = boundaryChain(g)
    out = [chain]
    for d in range(g.getMaxDimCells()-1):
        chain = list(set(CAT(AA(DOWNCELLS(g))(chain))))
        out.append(chain)
    return out[::-1]

def extrude(g):
    d = g.getMaxDimCells()
    h = Quote(Graph(1),[1])
    g = Graph.power(Matf(1),Matf(1),  g,None,None,  h,None,None)
    return g

def sweep(g,scaling=2.0):
    # compute boundary complex
    bComplex = boundaryComplex(g)
    d = g.getMaxDimCells()
    n = g.getPointDim()
    maps = dict([[cell,g.addNode(0)] for cell in bComplex[0]]) 
    # add information to 0-cells
    for cell in maps.iterkeys():
        vert = maps[cell]
        point = [g.getVecf(cell).get(i) for i in range(1,n+1)]
        point = SCALARVECTPROD([UNITVECT(point),scaling])
        g.setVecf(vert,Vecf([1]+point))
    # add information to 1-cells
    for node in maps.iterkeys():
        newNode = maps[node]
        newArc = g.addNode(1)
        g.addArch(node,newArc)
        g.addArch(newNode,newArc)
    # compute pivotal vertices
    for cell in bComplex[-1]:        
        newcell = g.addNode(d); g.addArch(cell,newcell)
        pivotPoint = list(CENTROID(g)(cell).get(i)*scaling for i in range(1,n+1))
        pivotPoint = SCALARVECTPROD([UNITVECT(pivotPoint),scaling])
        pivotCell = g.addNode(0);g.setVecf(pivotCell,Vecf([1]+pivotPoint))
        # compute pivotal cells
        for c in DOWNCELLS(g)(cell):
            newNode = g.addNode(1)
            g.addArch(maps[c],newNode)
            g.addArch(pivotCell,newNode)
            g.addArch(newNode,newcell)
            g.addArch(GETINTERSECTION(g)([c,maps[c]])[0],newcell)    
    return g


if __name__ == "__main__":
    
    
    h = ngon(3)
    DRAW(h)()
    h = sweep(h)
    DRAW(h)()
    h = sweep(h,3)
    DRAW(h)()
    h = sweep(h,4)
    DRAW(h)()
    h = hccmesh(h)
    DRAW(h)()
    h = extrude(h)
    DRAW(h,[2,2,2])()
