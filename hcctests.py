def sphere_kernel():
    
    g = Graph(3)
    origin = g.addNode(0); g.setVecf(origin, Vecf([1.0, 0, 0,0]))
    x = g.addNode(0); g.setVecf(x, Vecf([1.0, 1, 0,0]))
    y = g.addNode(0); g.setVecf(y, Vecf([1.0, 0, 1,0]))
    z = g.addNode(0); g.setVecf(z, Vecf([1.0, 0, 0,1]))
    xn = g.addNode(0); g.setVecf(xn, Vecf([1.0, -1, 0, 0]))
    yn = g.addNode(0); g.setVecf(yn, Vecf([1.0, 0, -1, 0]))
    zn = g.addNode(0); g.setVecf(zn, Vecf([1.0, 0, 0, -1]))
    
    points = [[cos(u), sin(u), 0] for u in arange(pi/4, 2*pi, pi/4)]
         
    for i in range(len(points)):
        vert = g.addNode(0)
        g.setVecf(vert, Vecf([1.0] + points[i]))
        
    points = [[0, cos(u), sin(u)] for u in arange(pi/4, 2*pi, pi/2)]

    for i in range(len(points)):
        vert = g.addNode(0)
        g.setVecf(vert, Vecf([1.0] + points[i]))
        
    points = [[cos(u), 0, sin(u)] for u in arange(pi/4, 2*pi, pi/2)]

    for i in range(len(points)):
        vert = g.addNode(0)
        g.setVecf(vert, Vecf([1.0] + points[i]))
        
    return g



#Sembrava stampare (quasi) i paralleli
norm_zero = []
for id, vect in points.iteritems():
    if INNERPROD([vect, [0,0,1]]) > 0:
        norm_zero.append(id)