#Test con Cython
from hccsphere import hexsphere
from utils import DRAW
import time
from pyplasm import VIEW

start = time.clock()
g = hexsphere(2, 1.2)
DRAW(g, [2.0, 2.0, 2.0])()
end = time.clock()

print "Sphere builded in ", end - start, " seconds."