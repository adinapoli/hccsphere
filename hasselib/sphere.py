from graphLib6 import *

def mycube():
    fun = COMP([R([1,2])(PI/4),ROTN([ATAN(SQRT(2)),[1,-1,0]]),
                STRUCT,CONS([ID,K(POLYLINE([[-0.5,-0.5,0.5],[0.5,0.5,-0.5]]))]),
                SKELETON(1),T([1,2,3])([-0.5,-0.5,-0.5]),CUBE])
    return fun(1)

VIEW(mycube())
