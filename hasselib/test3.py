from graphLib6 import *

    
def crown2D(points1,points2):
    n = len(points1)
    g = Graph(2)
    # insert 0-nodes ////////////////////////////////
    for i in range(n):
        vert = g.addNode(0); g.setVecf(vert,Vecf([1.0]+points1[i]))

    for i in range(n):
        vert = g.addNode(0); g.setVecf(vert,Vecf([1.0]+points2[i]))
        
    # insert 1-nodes ////////////////////////////////
    for i in range(1,n):
        edge = g.addNode(1); g.addArch(i,edge); g.addArch(i+1,edge);
        g.setVecf(edge,g.getFittingPlane(edge))
    edge = g.addNode(1); g.addArch(n,edge); g.addArch(1,edge);

    for i in range(1,n):
        edge = g.addNode(1);
        g.addArch(n+i,edge); g.addArch(n+i+1,edge);
        g.setVecf(edge,g.getFittingPlane(edge))
    edge = g.addNode(1); g.addArch(2*n,edge); g.addArch(n+1,edge);
    m = edge

    for i in range(1,n+1):
        edge = g.addNode(1); g.addArch(i,edge); g.addArch(n+i,edge);
        g.setVecf(edge,g.getFittingPlane(edge))

    # insert 2-nodes ////////////////////////////////
    for i in range(1,n):
        face = g.addNode(2);
        g.addArch(2*n+i,face);g.addArch(3*n+i,face);
        g.addArch(4*n+i,face);g.addArch(4*n+i+1,face);
        g.addArch(face,i);g.addArch(face,i+1);
        g.addArch(face,n+i);g.addArch(face,n+i+1);
        
    face = g.addNode(2);
    g.addArch(3*n,face);g.addArch(4*n,face);
    g.addArch(5*n,face);g.addArch(4*n+1,face);
    g.addArch(face,n);g.addArch(face,2*n);
    g.addArch(face,1);g.addArch(face,1+n);
    

    face = g.addNode(2);
    for edge in range(2*n+1,3*n+1):
        g.addArch(edge,face)

    return g


if False:
    points = [[COS(i*(2*PI)/6),SIN(i*(2*PI)/6)] for i in range(6)]
    points1 = [[p[0]*1, p[1]*1] for p in points]
    points2 = [[p[0]*2, p[1]*2] for p in points]
    
    g = crown2D(points1,points2)
    VIEW(Hpc(g)) # ??? 

    g = mymesh(g) # OK
    SHOW(g)(CELLSPERLEVEL(g)(0) + CELLSPERLEVEL(g)(1) + CELLSPERLEVEL(g)(2))

    g = mymesh(g) # OK
    SHOW(g)(CELLSPERLEVEL(g)(0) + CELLSPERLEVEL(g)(1) + CELLSPERLEVEL(g)(2))
    
    VIEW(Hpc(g)) # ??? 

if True:
    points = [[COS(i*(2*PI)/6),SIN(i*(2*PI)/6)] for i in range(6)]
    points1 = [[p[0]*1, p[1]*1] for p in points]
    points2 = [[p[0]*2, p[1]*2] for p in points]
    
    g = crown2D(points1,points2)
    g1 = Graph.cuboid(1)
    vmat = Matf(1); hmat = Matf(1)
    gg = Graph.power(vmat,hmat,  g,None,None,  g1,None,None)
    VIEW(SKELETON(1)(Hpc(gg)))
     

    g = hccmesh(gg)
    VIEW(SKELETON(1)(Hpc(g)))
    #DRAW(g,[1.5,1.5,1.5])(range(1,g.getNumNode()+1))
    g = hccmesh(g)
    SHOW(g,[1.5,1.5,1.5])
    #DRAW(g,[1.5,1.5,1.5])()

