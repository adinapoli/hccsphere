from pyplasm import *

# convex combination operator
CCOMB = COMP([ SCALARVECTPROD,
              CONS([ COMP([ DIV, CONS([K(1),LEN]) ]), VECTSUM ]) ])

def EXPLODE (sx,sy,sz):
    def explode0 (scene):
        centers = [CCOMB(S1(UKPOL(obj))) for obj in scene]
        scalings = len(centers) * [S([1,2,3])([sx,sy,sz])]        
        scaledCenters = [UK(APPLY(pair)) for pair in
            zip(scalings, [MK(p) for p in centers])]
        translVectors = [ VECTDIFF((p,q)) for (p,q) in zip(scaledCenters, centers) ]
        translations = [ T([1,2,3])(v) for v in translVectors ]
        return STRUCT([ APPLY((t,obj)) for (t,obj) in zip(translations,scene) ])
    return explode0


def tetra2hexs (V0):
    # Original points
    v0,v1,v2,v3 = V0
    # Steiner points
    V1 = AA(CCOMB)([[v0,v1],[v0,v2],[v0,v3],[v1,v2],[v1,v3],[v2,v3]])
    V2 = AA(CCOMB)([[v0,v1,v2],[v0,v1,v3],[v0,v2,v3],[v1,v2,v3]])
    V3 = AA(CCOMB)([[v0,v1,v2,v3]])
    # Generated cells
    c0 = JOIN(AA(MK)([V0[0],V3[0],V1[0],V1[1],V1[2],V2[0],V2[1],V2[2]]))
    c1 = JOIN(AA(MK)([V0[1],V3[0],V1[0],V1[3],V1[4],V2[0],V2[1],V2[3]]))
    c2 = JOIN(AA(MK)([V0[2],V3[0],V1[1],V1[3],V1[5],V2[0],V2[2],V2[3]]))
    c3 = JOIN(AA(MK)([V0[3],V3[0],V1[2],V1[4],V1[5],V2[1],V2[2],V2[3]]))
    return [c0,c1,c2,c3]


if __name__ == "__main__":
    V0 = [0,0,0],[1,0,0],[0,1,0],[0,0,1]
    VIEW(EXPLODE(1.2,1.2,1.2)(tetra2hexs(V0)))
