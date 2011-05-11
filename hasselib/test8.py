from test6 import *
from test4 import tetra
from OrderedDict import *

def unitSimplex(d):
    for k in range(d+1):
        v = g.addNode(0); g.setVecf(v,Vecf(*([1]+p)))
    

def combinations(Set):
    def combinations0(n):
        return [(i[1],j[1]) for i in enumerate(Set)
                for j in enumerate(Set) if j[0]>i[0]]
    return combinations0


def pcsphere(g,scaling=2):    

    # STEP 1: d-1 boundary-complex extraction \\\\\\\\\\\\\\\
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    bComplex = boundaryComplex(g)
    'myprint("bComplex",bComplex)'
    for level in bComplex:
        'myprint("level",(level,AA(g.Level)(level)))'
    d = g.getMaxDimCells()
    n = g.getPointDim()
    DRAW(g,[1.5,1.5,1.5])(CAT(bComplex))

    # STEP 2: d-2 wall-complex extraction \\\\\\\\\\\\\\\\\\\
    wallComplex = bComplex[:-1]
    'myprint("wallComplex",wallComplex)'
    DRAW(g,[1.5,1.5,1.5])(CAT(wallComplex))

    # STEP 3: basket separators construction \\\\\\\\\\\\\\\\
    """ (d-1)-complex construction """
    mapping = {}
    """ duplicate wall-complex """
    for skeleton in wallComplex:
        """ add 0-cells of mapped wall-complex cells """
        for cell in skeleton:
            'myprint("cell",cell)'
            newCell = g.addNode(0)
            mapping.update({cell:newCell})
            point = [CENTROID(g)(cell).get(i) for i in range(1,n+1)]
            point = SCALARVECTPROD([UNITVECT(point),scaling])
            g.setVecf(mapping[cell],Vecf([1]+point))

    inverted = dict([v,k] for k,v in mapping.iteritems())
                
    """ add 1-cell extensions to top corners """
    for node in wallComplex[0]:
        newNode = mapping[node]
        newArc = g.addNode(1)
        g.addArch(node,newArc)
        g.addArch(newNode,newArc)

    DRAW(g)()
        
    """ add 1-cells of top-lateral boundary of added polytopes """
    for node in wallComplex[0]:
        for arc in UPCELLS(g)(node):
            if arc in wallComplex[1]:
                newArc = g.addNode(1)
                g.addArch(mapping[node],newArc)
                g.addArch(mapping[arc],newArc)

    DRAW(g)(CELLSPERLEVEL(g)(0)+CELLSPERLEVEL(g)(1)+CELLSPERLEVEL(g)(3))

    """ construction of 2-faces of wall-complex """
    for edge in wallComplex[1]:
        newFace = g.addNode(2)
        g.addArch(edge,newFace)
        for node in DOWNCELLS(g)(edge):
            arc = GETINTERSECTION(g)([node,mapping[node]])[0]
            g.addArch(arc,newFace)
            faceTop = GETINTERSECTION(g)([mapping[edge],mapping[node]])[0]
            g.addArch(faceTop,newFace)

    # STEP 4: basket content construction \\\\\\\\\\\\\\\\\\\

    """ for each boundary facet do: """
    
    for cell in bComplex[-1]:
        faces = []
        'myprint("cell",cell)'
        cellEdges = DOWNCELLS(g)(cell)
        
        """ add 0-tips of added polytopes """
        topCell = g.addNode(0)
        mapping.update({cell:topCell})
        point = [CENTROID(g)(cell).get(i) for i in range(1,n+1)]
        point = SCALARVECTPROD([UNITVECT(point),scaling])
        g.setVecf(topCell,Vecf([1]+point))
        
        """ add 1-edges from tips to top nodes of lateral basket"""
        for node in mapping.iterkeys():
            'myprint("cell,node",(cell,node))'
            if (cell != node) and (node in DOWNCELLS(g)(cell)):
                newArc = g.addNode(1)
                g.addArch(mapping[node],newArc)
                g.addArch(topCell,newArc)

        def cyclicPairs(seq):
            return TRANS([seq,seq[1:]+[seq[0]]])        

        """ add 1-edges from tips to mapped 0-nodes of bottom faces"""
        for node in list(set(CAT(AA(DOWNCELLS(g))(DOWNCELLS(g)(cell))))):
            'myprint("node",node)'
            newArc = g.addNode(1)
            g.addArch(mapping[node],newArc)
            g.addArch(topCell,newArc)

        """ build (basket's) top cycle of edges """
        vertPairs = CAT([DISTR([DOWNCELLS(g)(edge),edge])
              for edge in DOWNCELLS(g)(cell)])
        mappedVertPairs = [[mapping[i],mapping[j]] for (i,j) in vertPairs]
        edges = [GETINTERSECTION(g)(pair)[0] for pair in mappedVertPairs]

        inverted = dict([v,k] for k,v in mapping.iteritems())

        """ build basket's cap of triangles """
        def test(cell,edge,face):
            edges = set(DOWNCELLS(g)(face)).difference([edge])
            if len(edges)>3: return True
            else:
                vertpairs = AA(DOWNCELLS(g))(edges)
                commonVert = list(set(vertpairs[0]).intersection(vertpairs[1]))[0]
                return inverted[commonVert]==cell

        for edge in edges:
            newFace = g.addNode(2)
            g.addArch(edge,newFace)
            g.addArch(GETINTERSECTION(g)([DOWNCELLS(g)(edge)[0],topCell])[0],newFace)
            g.addArch(GETINTERSECTION(g)([DOWNCELLS(g)(edge)[1],topCell])[0],newFace)
            faces += [f for f in UPCELLS(g)(edge) if test(cell,edge,f)]

            
        """ build solid basket """
        newSolid = g.addNode(3)
        faces = list(set(faces + [cell]))
        for face in faces: g.addArch(face,newSolid)
        nodes = DOWNCELLS(g)(newSolid)
        for i in range(2):
            nodes = list(set(CAT(AA(DOWNCELLS(g))(nodes))))
            
        for node in nodes:
            g.addArch(newSolid,node)

    'myprint("CELLSPERLEVEL(g)(0)",CELLSPERLEVEL(g)(0))'
    'myprint("CELLSPERLEVEL(g)(1)",CELLSPERLEVEL(g)(1))'
    'myprint("CELLSPERLEVEL(g)(2)",CELLSPERLEVEL(g)(2))'
    'myprint("CELLSPERLEVEL(g)(3)",CELLSPERLEVEL(g)(3))'
    'myprint("mapping",mapping)'
            
    DRAW(g,[2.5,2.5,2.5])(CELLSPERLEVEL(g)(3))


    return g


h = ngon(6)
h = pcsphere(tetra())
h = pcsphere(h,3)
g = hccmesh(h)
DRAW(g,[2.5,2.5,2.5])(CELLSPERLEVEL(g)(3))
h = pcsphere(h,4)
