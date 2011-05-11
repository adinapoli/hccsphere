from graphLib6 import *

def schonhardt():

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


    return g


if True:
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

    pointdim=3

    # //////////////////////////////////////////////////////////////////////////////
    # correct orientation 
    planes,centroid={} , (V02+V12+V22+V03+V13+V23)/6.0

    g = schonhardt()
                     
    # esempio di uso di BSP
    nav,tolerance,max_try=GraphNavigator(),1e-6,10 
    GBsp=Graph.cuboid(pointdim,-10.0,+10.0)
    cell=GBsp.each(3).getNode()
                     
    for f2d in CELLSPERLEVEL(g)(2):
        h=g.getFittingPlane(f2d).forceBelow(centroid)
        [cell_below,cell_equal,cell_above]=GBsp.split(nav,cell,h,tolerance,max_try)
        cell=cell_below
        if cell_above!=0: GBsp.remNode(cell_above)
                                       
    Plasm.View(Hpc(GBsp))
    # //////////////////////////////////////////////////////////////////////////////

if __name__ == "__main__":

    DRAW(schonhardt())()
    DRAW(schonhardt(),[1.5,1.5,1.5])()
