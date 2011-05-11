from graphLib6 import *


def ngon(n): return CIRCUMFERENCE(1)(n)

def tetrahedron():
    return JOIN([ COMP([T(3)(-1.0/3),EMBED(1), ngon ])(3), MK([0,0,1]) ])

def mapping(g):
    for vert in CELLSPERLEVEL(g)(0):
        point = [g.getVecf(vert)[k] for k in range(1,4)]
        size = VECTNORM(point)
        newpoint = SCALARVECTPROD([1./size,point])
        g.setVecf(vert,Vecf([1.0] + newpoint))


# /////////////////////////////////////////////
# creo un tetraedro
# /////////////////////////////////////////////

def tetra():
    
    pointdim=3
    g=Graph(pointdim)
    nav=GraphNavigator()

    verts = S1(UKPOL(tetrahedron()))

    v0=g.addNode(0);g.setVecf(v0,Vecf(*([1]+verts[0])))
    v1=g.addNode(0);g.setVecf(v1,Vecf(*([1]+verts[1])))
    v2=g.addNode(0);g.setVecf(v2,Vecf(*([1]+verts[2])))
    v3=g.addNode(0);g.setVecf(v3,Vecf(*([1]+verts[3])))

    e01=g.addNode(1);g.addArch(v0,e01);g.addArch(v1,e01)
    e12=g.addNode(1);g.addArch(v1,e12);g.addArch(v2,e12)
    e20=g.addNode(1);g.addArch(v2,e20);g.addArch(v0,e20)
    e03=g.addNode(1);g.addArch(v0,e03);g.addArch(v3,e03)
    e13=g.addNode(1);g.addArch(v1,e13);g.addArch(v3,e13)
    e23=g.addNode(1);g.addArch(v2,e23);g.addArch(v3,e23)

    # attenzione: non sono correttamente orientate! bisognerebbe crearle coerentemente orientate
    f012=g.addNode(2);g.addArch(e01,f012);g.addArch(e12,f012);g.addArch(e20,f012)
    f013=g.addNode(2);g.addArch(e01,f013);g.addArch(e13,f013);g.addArch(e03,f013)
    f123=g.addNode(2);g.addArch(e12,f123);g.addArch(e23,f123);g.addArch(e13,f123)
    f203=g.addNode(2);g.addArch(e20,f203);g.addArch(e03,f203);g.addArch(e23,f203)

    # per pointdim>=3, le (pointdim-1)-celle devono avere l'informazione circa l'iperpiano passante
    g.setVecf(f012,g.getFittingPlane(f012))
    g.setVecf(f013,g.getFittingPlane(f013))
    g.setVecf(f123,g.getFittingPlane(f123))
    g.setVecf(f203,g.getFittingPlane(f203))

    tet=g.addNode(3)
    g.addArch(f012,tet)
    g.addArch(f013,tet)
    g.addArch(f123,tet)
    g.addArch(f203,tet)

    g.addArch(tet,v0)
    g.addArch(tet,v1)
    g.addArch(tet,v2)
    g.addArch(tet,v3)

    return g

def pyramid():
    g = Graph(3)
    
    v0=g.addNode(0); g.setVecf(v0,Vecf(1,0,0,1))
    v1=g.addNode(0); g.setVecf(v1,Vecf(1,COS(0*PI/2),SIN(0*PI/2),-1./3))
    v2=g.addNode(0); g.setVecf(v2,Vecf(1,COS(1*PI/2),SIN(1*PI/2),-1./3))
    v3=g.addNode(0); g.setVecf(v3,Vecf(1,COS(2*PI/2),SIN(2*PI/2),-1./3))
    v4=g.addNode(0); g.setVecf(v4,Vecf(1,COS(3*PI/2),SIN(3*PI/2),-1./3))

    e01=g.addNode(1);g.addArch(v0,e01);g.addArch(v1,e01)
    e02=g.addNode(1);g.addArch(v0,e02);g.addArch(v2,e02)
    e03=g.addNode(1);g.addArch(v0,e03);g.addArch(v3,e03)
    e04=g.addNode(1);g.addArch(v0,e04);g.addArch(v4,e04)
    e12=g.addNode(1);g.addArch(v1,e12);g.addArch(v2,e12)
    e23=g.addNode(1);g.addArch(v2,e23);g.addArch(v3,e23)
    e34=g.addNode(1);g.addArch(v3,e34);g.addArch(v4,e34)
    e41=g.addNode(1);g.addArch(v4,e41);g.addArch(v1,e41)

    f012=g.addNode(2);g.addArch(e01,f012);g.addArch(e12,f012);g.addArch(e02,f012)
    f023=g.addNode(2);g.addArch(e02,f023);g.addArch(e23,f023);g.addArch(e03,f023)
    f034=g.addNode(2);g.addArch(e03,f034);g.addArch(e34,f034);g.addArch(e04,f034)
    f041=g.addNode(2);g.addArch(e04,f041);g.addArch(e41,f041);g.addArch(e01,f041)
    f1234=g.addNode(2);g.addArch(e12,f1234);g.addArch(e23,f1234)
    g.addArch(e34,f1234);g.addArch(e41,f1234);

    pyr=g.addNode(3)
    g.addArch(f012,pyr);g.addArch(f023,pyr);g.addArch(f034,pyr);
    g.addArch(f041,pyr);g.addArch(f1234,pyr)
    g.addArch(pyr,v0);g.addArch(pyr,v1);g.addArch(pyr,v2);
    g.addArch(pyr,v3);g.addArch(pyr,v4)



    return g

if __name__ == "__main__":


    DRAW(pyramid(),[1.5,1.5,1.5])()


if False:    

    VIEW(STRUCT([tetrahedron(), SKELETON(1)(SPHERE(1)([12,18]))]))

    myprint("UKPOL(tetrahedron())",UKPOL(tetrahedron()))

    g = tetra()
    g.Print()

    Plasm.View(Hpc(g))

    DRAW(g,[1.5,1.5,1.5])()

    g = hccmesh(g)
    DRAW(g,[1.5,1.5,1.5])()
    g = hccmesh(g)
    DRAW(g,[1.5,1.5,1.5])()
    g = hccmesh(g)
    DRAW(g,[1.5,1.5,1.5])()
    myprint("CELLSPERLEVEL(g)(2)",CELLSPERLEVEL(g)(2))


    
##    g = mymesh(g)
##    #mapping(g)
##    DRAW(g,[2,2,2])()
##    #g = mymesh(g)
##    DRAW(g,[2,2,2])()

#  CAT ~ [  DISTR ~ [[DOWNCELLS(g) ~ S1], S1],  TAIL  ]

