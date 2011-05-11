from graphLib6 import *
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


def addLevelSet(g,dr):

    def boundary(g):
        return [k for k in CELLSPERLEVEL(g)(1) if len(UPCELLS(g)(k))==1]

    bndry = boundary(g)
    n = len(bndry)
    # get the first vertex on the boundary manifold
    firstBndrVert = UPCELLS(g)(UPCELLS(g)(bndry[0])[0])[0]
    r = VECTNORM([g.getVecf(firstBndrVert).get(k) for k in range(2)])
    
    points = [[(r+dr)*cos(u), (r+dr)*sin(u)] for u in arange(0,2*pi,2*pi/(2*n))]

    m = g.getNumNode()
    for i in range(2*n):          # add 0-cells
        vert = g.addNode(0);
        g.setVecf(vert,Vecf([1.0]+points[i]))
        
    for i in range(n):          # add 1-cells
        edge = g.addNode(1);
        g.addArch(m+2*i+1,edge);
        g.addArch(edge,m+2*i+1);
        g.addArch((m+2*i+2 if 2*i<2*n-1 else m+1),edge);        
        g.addArch(edge, (m+2*i+2 if 2*i<2*n-1 else m+1));        
        edge = g.addNode(1);
        g.addArch(m+(2*i+1)+1,edge);
        g.addArch(edge, m+(2*i+1)+1);
        g.addArch((m+(2*i+1)+2 if (2*i+1)<2*n-1 else m+1),edge);        
        g.addArch(edge, (m+(2*i+1)+2 if (2*i+1)<2*n-1 else m+1));        
        edge = g.addNode(1);
        g.addArch(i+1,edge);
        g.addArch(edge,i+1);
        g.addArch(m+(2*i)+1,edge);        
        g.addArch(edge, m+(2*i)+1);        
        
    for i in range(n):      # add 2-cells
        face = g.addNode(2)
        g.addArch(n+1+i,face)
        g.addArch(2*m+3*i,face)
        g.addArch(2*m+3*i+1,face)
        g.addArch(2*m+3*i+2,face)
        g.addArch((2*m+3*i+5 if (3*i+1)<(3*n-2) else 2*m+2), face)

        faceAsEdges = DOWNCELLS(g)(face)        
        facesAsVerts = list(set(CAT(AA(DOWNCELLS(g))(faceAsEdges))))                        
        for i in facesAsVerts: g.addArch(face,i)

    return g


if __name__=="__main__":

    g = ngon(4)
    DRAW(ngon(4))()

    g = addLevelSet(g,1)
    g = hccmesh(g)
    DRAW(g)()
    g = hccmesh(g)
    DRAW(g,[1.5,1.5,1.5])()
