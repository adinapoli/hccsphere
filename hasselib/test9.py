from test6 import *
from test4 import tetra,pyramid
""" Geodetic sphere generation by superimposition of polytopal layers. """

def unitSimplex(d):
    for k in range(d+1):
        v = g.addNode(0); g.setVecf(v,Vecf(*([1]+p)))
    

def combinations(Set):
    def combinations0(n):
        return [(i[1],j[1]) for i in enumerate(Set)
                for j in enumerate(Set) if j[0]>i[0]]
    return combinations0


def pcsphere(g,scaling=2,trihedrons=False):    

    # STEP 1: d-1 boundary-complex extraction \\\\\\\\\\\\\\\
    # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    bComplex = boundaryComplex(g)
    'myprint("bComplex",bComplex)'
    for level in bComplex:
        'myprint("level",(level,AA(g.Level)(level)))'
    d = g.getMaxDimCells()
    n = g.getPointDim()
    #DRAW(g,[1.5,1.5,1.5])(CAT(bComplex))

    # STEP 2: d-2 wall-complex extraction \\\\\\\\\\\\\\\\\\\
    wallComplex = bComplex[:-1]
    'myprint("wallComplex",wallComplex)'
    #DRAW(g,[1.5,1.5,1.5])(CAT(wallComplex))

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

        
    """ add 1-cells of top-lateral boundary of added polytopes """
    for node in wallComplex[0]:
        for arc in UPCELLS(g)(node):
            if arc in wallComplex[1]:
                newArc = g.addNode(1)
                g.addArch(mapping[node],newArc)
                g.addArch(mapping[arc],newArc)
                
##
##    DRAW(g)(CELLSPERLEVEL(g)(0)+CELLSPERLEVEL(g)(1)+CELLSPERLEVEL(g)(3))
##
    """ construction of 2-faces of wall-complex """
    for edge in wallComplex[1]:
        
        if trihedrons == False:
            # build pentagonal faces
            newFace = g.addNode(2)
            g.addArch(edge,newFace)
            for node in DOWNCELLS(g)(edge):
                arc = GETINTERSECTION(g)([node,mapping[node]])[0]
                g.addArch(arc,newFace)
                faceTop = GETINTERSECTION(g)([mapping[edge],mapping[node]])[0]
                g.addArch(faceTop,newFace)

        # build missing edges (iff trihedrons=True)
        if trihedrons == True:
            verts = DOWNCELLS(g)(edge)
            hverts = [mapping[vert] for vert in verts]
            newEdge = g.addNode(1)
            g.addArch(hverts[0],newEdge)
            g.addArch(hverts[1],newEdge)
            'myprint("newEdge",newEdge)'

            """ build 1 quad + 1 triangle faces """
            # triangle face
            triangleFace = g.addNode(2)
            g.addArch(newEdge,triangleFace)
            g.addArch(GETINTERSECTION(g)([hverts[0],mapping[edge]])[0],triangleFace)
            g.addArch(GETINTERSECTION(g)([hverts[1],mapping[edge]])[0],triangleFace)            
            # quad face
            quadFace = g.addNode(2)
            g.addArch(edge,quadFace)
            g.addArch(newEdge,quadFace)
            g.addArch(GETINTERSECTION(g)([verts[0],hverts[0]])[0],quadFace)
            g.addArch(GETINTERSECTION(g)([verts[1],hverts[1]])[0],quadFace)

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

        """ build solid cells """
        edges = DOWNCELLS(g)(cell)
        verts = list(set(CAT(AA(DOWNCELLS(g))(edges))))
        top = mapping[cell]
        hverts = [mapping[vert] for vert in verts]
        wedgeFaces = [cell]
        
        if trihedrons == False:
            
            """ build solid basket """
            newSolid = g.addNode(3)
            faces = list(set(faces + [cell]))
            for face in faces: g.addArch(face,newSolid)
            nodes = DOWNCELLS(g)(newSolid)
            for i in range(2):
                nodes = list(set(CAT(AA(DOWNCELLS(g))(nodes))))                
            for node in nodes:
                g.addArch(newSolid,node)
        else:
            """ build basket's trihedral decomposition """
            hedges = [ CAT(GETINTERSECTION(g)([hpair[0],hpair[1]]))
                    for hpair in cyclicPairs(hverts)]
            
            # one internal tetrahedron
            vertices = hverts + [top]
            tetraEdges = CAT(AA(GETINTERSECTION(g))(cyclicPairs(hverts)))
            'myprint("tetraEdges",tetraEdges)'
            tetraCell = g.addNode(3)
            for vert in vertices:
                g.addArch(tetraCell,vert);
            bottomFace = g.addNode(2)
            wedgeFaces += [bottomFace]
            g.addArch(bottomFace,tetraCell);
            for edge in tetraEdges:
                g.addArch(edge,bottomFace);
                edge0 = CAT(GETINTERSECTION(g)([DOWNCELLS(g)(edge)[0],top]))
                edge1 = CAT(GETINTERSECTION(g)([DOWNCELLS(g)(edge)[1],top]))
                newFace = g.addNode(2)
                g.addArch(edge0,newFace); g.addArch(edge1,newFace);g.addArch(edge,newFace);
                g.addArch(newFace,tetraCell);
            'myprint("bottomFace",bottomFace)'

            # three external tetrahedra
            for k in range(3):
                hvert2 = cyclicPairs(hverts)[k]
                nod = mapping[GETINTERSECTION(g)(
                    [inverted[v] for v in hvert2])[0]]
                tetraEdges = GETINTERSECTION(g)(hvert2)
                tetraEdges += GETINTERSECTION(g)([hvert2[0],top])
                tetraEdges += GETINTERSECTION(g)([hvert2[1],top])
                tetraEdges += GETINTERSECTION(g)([nod,top])
                tetraEdges += GETINTERSECTION(g)([hvert2[0],nod])
                tetraEdges += GETINTERSECTION(g)([hvert2[1],nod])
                face0 = GETINTERSECTION(g)([tetraEdges[0],tetraEdges[1],tetraEdges[2]])
                face1 = GETINTERSECTION(g)([tetraEdges[1],tetraEdges[3],tetraEdges[4]])
                face2 = GETINTERSECTION(g)([tetraEdges[2],tetraEdges[3],tetraEdges[5]])
                face3 = GETINTERSECTION(g)([tetraEdges[0],tetraEdges[4],tetraEdges[5]])
                #DRAW(g)(CELLSPERLEVEL(g)(0)+tetraEdges)
                newCell = g.addNode(3)
                for face in CAT([face0,face1,face2,face3]):
                    g.addArch(face,newCell)
                vertices = [nod,top]+hvert2
                for vertex in vertices:
                    g.addArch(newCell,vertex)
            
            # one wedge solid attached to triangle in bComplex
            wedgeCell = g.addNode(3)
            for vert in verts+hverts: g.addArch(wedgeCell,vert)

            sideFaces = list(set(CAT(AA(UPCELLS(g))(edges))).intersection(
                CAT(AA(UPCELLS(g))(hedges))))
                
            wedgeFaces += sideFaces
            for face in wedgeFaces: g.addArch(face,wedgeCell)

    'myprint("CELLSPERLEVEL(g)(0)",CELLSPERLEVEL(g)(0))'
    'myprint("CELLSPERLEVEL(g)(1)",CELLSPERLEVEL(g)(1))'
    'myprint("CELLSPERLEVEL(g)(2)",CELLSPERLEVEL(g)(2))'
    'myprint("CELLSPERLEVEL(g)(3)",CELLSPERLEVEL(g)(3))'
    'myprint("mapping",mapping)'
            
    #DRAW(g,[2.5,2.5,2.5])(CELLSPERLEVEL(g)(3))

    return g

h = pcsphere(tetra(),trihedrons=True)
DRAW(h)()
##g = hccmesh(h)
##DRAW(g)()
##h = pcsphere(h,2,True)
##DRAW(h)()
##g = hccmesh(h)
##DRAW(g)()
##h = pcsphere(h,3,True)
##DRAW(h)()
##g = hccmesh(h)
##DRAW(g)()


##h = pyramid()
##DRAW(h)(CELLSPERLEVEL(h)(0)+CELLSPERLEVEL(h)(1))
##h = hccmesh(h)
##DRAW(h,[1.5,1.5,1.5])(CELLSPERLEVEL(h)(3))
##h = hccmesh(h)
##DRAW(h,[1.5,1.5,1.5])(CELLSPERLEVEL(h)(3))
##h = hccmesh(h)
##DRAW(h,[1.5,1.5,1.5])(CELLSPERLEVEL(h)(3))

##h = pcsphere(pyramid(),trihedrons=True)
##DRAW(h,[2.0,2.0,2.0])(CELLSPERLEVEL(h)(3))
##h = hccmesh(pcsphere(pyramid(),trihedrons=True))
##DRAW(h,[2.0,2.0,2.0])(CELLSPERLEVEL(h)(3))
##
##
##
###h = ngon(6)
##h = hccmesh(hccmesh(pcsphere(tetra(),trihedrons=False)))
##DRAW(h,[2.5,2.5,2.5])(CELLSPERLEVEL(h)(3))
##h = pcsphere(tetra(),trihedrons=False)
##DRAW(h,[1.5,1.5,1.5])(CELLSPERLEVEL(h)(3))
##h = pcsphere(tetra(),trihedrons=True)
##DRAW(h,[1.5,1.5,1.5])(CELLSPERLEVEL(h)(3))
##h = hccmesh(pcsphere(tetra(),trihedrons=True))
##DRAW(h,[1.5,1.5,1.5])(CELLSPERLEVEL(h)(3))
##h = pcsphere(h,3,trihedrons=True)
##g = boundaryComplex(h)         #TODO => deleting nodes in h ...
##
##DRAW(h,[1.5,1.5,1.5])([78,79,80,100,101,102,122,123,124,144,145,146]+CELLSPERLEVEL(h)(1))
##
###h = pcsphere(h,3)
##h = hccmesh(h)
##
##
##h = pcsphere(h,4)
