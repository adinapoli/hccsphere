from test6 import *
from test4 import tetra

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
    myprint("bComplex",bComplex)
    for level in bComplex:
        myprint("level",(level,AA(g.Level)(level)))
    d = g.getMaxDimCells()
    n = g.getPointDim()
    DRAW(g)(CAT(bComplex))

    # STEP 2: d-2 wall-complex extraction \\\\\\\\\\\\\\\\\\\
    wallComplex = bComplex[:-1]
    'myprint("wallComplex",wallComplex)'
    DRAW(g)(CAT(wallComplex))

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
                
    """ add 1-cell extensions to top corners """
    for node in wallComplex[0]:
        newNode = mapping[node]
        newArc = g.addNode(1)
        g.addArch(node,newArc)
        g.addArch(newNode,newArc)
        
        """ add 1-cells of top-lateral boundary of added polytopes """
        for arc in UPCELLS(g)(node)[:-1]:
            newArc = g.addNode(1)
            g.addArch(mapping[node],newArc)
            g.addArch(mapping[arc],newArc)

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
        'myprint("cell",cell)'
        cellEdges = DOWNCELLS(g)(cell)
        
        """ add 0-tips of added polytopes """
        newCell = g.addNode(0)
        mapping.update({cell:newCell})
        point = [CENTROID(g)(cell).get(i) for i in range(1,n+1)]
        point = SCALARVECTPROD([UNITVECT(point),scaling])
        g.setVecf(mapping[cell],Vecf([1]+point))
        
        """ add 1-edges from tips to top nodes of lateral basket"""
        tips = []
        for node in mapping.iterkeys():
            'myprint("cell,node",(cell,node))'
            if (cell != node) and (node in DOWNCELLS(g)(cell)):
                tips += [mapping[node]]
                newArc = g.addNode(1)
                g.addArch(mapping[node],newArc)
                g.addArch(newCell,newArc)

        def cyclicPairs(seq):
            return TRANS([seq,seq[1:]+[seq[0]]])

        
##        myprint("cyclicPairs(cellEdges)",cyclicPairs(cellEdges))
##        for pair in cyclicPairs(cellEdges):
##            commonNode = AA(DOWNCELLS(g))(pair)
##            commonNode = list(set(commonNode[0]).intersection(commonNode[1]))
##            nodeTriple = [mapping[k] for k in CAT([commonNode,pair])]
##            edge1 = GETINTERSECTION(g)([nodeTriple[0],nodeTriple[1]])[0]
##            edge2 = GETINTERSECTION(g)([nodeTriple[0],nodeTriple[2]])[0]
##            myprint("nodeTriple",nodeTriple)
##            myprint("edge1,edge2",(edge1,edge2))
##            face = g.addNode(2)
##            g.addArch(edge1,face);g.addArch(edge2,face)
##            g.setVecf(face,g.getFittingPlane(face))
##            myprint("g.getVecf(face)",g.getVecf(face))
##        
##
##        roofTriples = AA(AL)(DISTL([newCell,cyclicPairs(UPCELLS(g)(newCell))]))
##        myprint("roofTriples",roofTriples)
##        for triple in roofTriples:
##            face = g.addNode(2)
##            for node in triple: g.addArch(node,face)
##            g.setVecf(face,g.getFittingPlane(face))
        

        """ add 1-edges from tips to mapped 0-nodes of bottom faces"""
        for node in list(set(CAT(AA(DOWNCELLS(g))(DOWNCELLS(g)(cell))))):
            'myprint("node",(node))'
            newArc = g.addNode(1)
            g.addArch(mapping[node],newArc)
            g.addArch(newCell,newArc)

    myprint("CELLSPERLEVEL(g)(0)",CELLSPERLEVEL(g)(0))
    myprint("CELLSPERLEVEL(g)(1)",CELLSPERLEVEL(g)(1))
    myprint("CELLSPERLEVEL(g)(2)",CELLSPERLEVEL(g)(2))
    myprint("mapping",mapping)
            
    DRAW(g)()
            

    # STEP 5: (d-1)-chain pivoting \\\\\\\\\\\\\\\\\\\\\\\\\\

    # STEP 6: (d-1)-chain roofing \\\\\\\\\\\\\\\\\\\\\\\\\\\

    # STEP 7: (d)-chain construction \\\\\\\\\\\\\\\\\\\\\\\\

    return g


h = ngon(6)
#h = pcsphere(h)
h = pcsphere(tetra())


