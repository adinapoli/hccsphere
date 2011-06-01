""" To generate the initial approximation of the sphere with hexahedra """

import numpy as np
from hasselib.graphLib6 import *
from utils import *

#/////////////////////////////////////////////////////////////////////
# Utility functions 
#/////////////////////////////////////////////////////////////////////
def minor(arr,i,j):
    # ith row, jth column removed
    return arr[np.array(range(i)+range(i+1,arr.shape[0]))[:,np.newaxis],
               np.array(range(j)+range(j+1,arr.shape[1]))]

def plane(theMat):
    """ Compute the coefficients of the plane equations in 3D """
    d = np.linalg.det(minor(theMat,0,0))
    a = -np.linalg.det(minor(theMat,0,1))
    b = np.linalg.det(minor(theMat,0,2))
    c = -np.linalg.det(minor(theMat,0,3))
    return d,a,b,c

def intersect(plane1,plane2,plane3):
    """ Intersection point of 3 planes in 3D -- first coord homogeneous """
    A = np.array([plane1[1:],plane2[1:],plane3[1:]])
    b = np.array([[-plane1[0]],[-plane2[0]],[-plane3[0]]])
    A_inv = np.linalg.inv(A)
    x0, x1,x2,x3 = CAT([[1.0]] + np.dot(A_inv,b).tolist())
    return [x0, x1,x2,x3]

def reflect(points,i):	
    ret,d = [],len(points[0])
    for point in points:
        p=[]
        for k in range(d):
            if k==i: p += [-point[k]]
            else: p += [point[k]]
        ret += [p]
    return ret

def classifyVerts(g):
    verts,vclass = CELLSPERLEVEL(g)(0),[[],[],[],[]]
    for v in verts:
        vert = eval(str(g.getVecf(v)))[1:]
        for k in range(4):
            if len([x for x in vert if x==0.0])==k:
                vclass[k] += [v]
                continue
    return vclass

def signature(vert):
    point = eval(str(g.getVecf(vert)))[1:]
    ret = []
    for x in point: ret += [SIGN(x) if x!=0 else 0]
    return ret

def test(v,w): 
    return INNERPROD([signature(v),signature(w)])

def matching(vclass1,signtr):		
    return [v for v in vclass1 if len([k for k in range(3)
        if signature(v)[k]==signtr[k]]) == 2]

def cyclicPairs(seq):
    return TRANS([seq,seq[1:]+[seq[0]]])

def minus(list1):
    def minus0(list2):
        return list(set(list2).difference(list1))
    return minus0

def intrsz(list1):
    def intrsz0(list2):
        return list(set(list2).intersection(list1))
    return intrsz0


######################################################################
#/////////////////////////////////////////////////////////////////////
# Sphere generation 
#/////////////////////////////////////////////////////////////////////
def initialSphere(g):
    """ Generation of an approximation of the standard unit sphere with 8 hexahedra.
        Return a Graph instance.
    """
    #//////////////////////////////////////////////////////////
    # generation of vertex nodes
    a = 1/SQRT(2)
    i,j,k = [[1,1,0,0],[1,0,1,0],[1,0,0,1]]
    ij,jk,ki = [[1,a,a,0],[1,0,a,a],[1,a,0,a]]
    zero = [1,0,0,0]
    mat_ij = np.array([zero,ki,i,ij])
    mat_jk = np.array([zero,ij,j,jk])
    mat_ki = np.array([zero,jk,k,ki])
    ijk = intersect(plane(mat_ij),plane(mat_jk),plane(mat_ki))
    verts = [zero,i,j,k,ij,jk,ki,ijk]
    for k in range(1,4):
        verts += reflect(verts,k)
    verts = AA(eval)(list(set(AA(repr)(verts))))
    for vert in verts:
        v = g.addNode(0); g.setVecf(v,Vecf(vert))

    vclass = classifyVerts(g)


    #//////////////////////////////////////////////////////////
    # generation of edge nodes

    edges = [[v,vclass[3][0]] for v in vclass[2]]
    edges += CAT([DISTL([v,[w for w in vclass[2] if test(v,w)==1]]) for v in vclass[1]])
    edges += CAT([DISTL([v,[w for w in vclass[1] if test(v,w)==2]]) for v in vclass[0]])

    for edge in edges:
         node = g.addNode(1); g.addArch(edge[0],node);g.addArch(edge[1],node)

    #//////////////////////////////////////////////////////////
    # generation of face nodes

    faces,o = [],vclass[3][0]

    for v in vclass[1]:
        for [u,w] in CART([vclass[2],vclass[2]]):
            if u<w and VECTSUM([ signature(u),signature(w) ])==signature(v):
                faces += [CAT(AA(GETINTERSECTION(g))([[v,u],[v,w],[o,u],[o,w]]))]

    for vert in vclass[0]:
        u,v,w = matching(vclass[1],signature(vert))
        vertTriples = [[u,v,vert],[v,w,vert],[w,u,vert]]
        for k in range(3):
            signatures = AA(signature)(vertTriples[k])
            common = COMP([AA(PROD),TRANS])(signatures)
            w = [v for v in vclass[2] if signature(v) == common]+vertTriples[k]
            faceEdgeVertices = [[w[0],w[1]],[w[0],w[2]],[w[1],w[3]],[w[2],w[3]]]
            faceEdges = CAT(AA(GETINTERSECTION(g))(faceEdgeVertices))
            faces += [faceEdges]

    for face in faces:
        node = g.addNode(2)
        for k in range(4): g.addArch(face[k],node)

    #//////////////////////////////////////////////////////////
    # generation of down arcs of 3D cell nodes

    for vert in vclass[0]:
        edges = UPCELLS(g)(vert)
        faces = list(set(CAT(AA(UPCELLS(g))(edges))))
        edgepairs = AA(DOWNCELLS(g))(faces)
        edgepairs = CAT(AA(minus(edges))(edgepairs))

        def int_test(pair):
            out = intrsz(DOWNCELLS(g)(pair[0]))(DOWNCELLS(g)(pair[1]))
            if len(out) != 1: return 0
            else: return out[0]

        epairs = [pair for pair in CART([edgepairs,edgepairs]) if int_test(pair) in vclass[1]]
        faces += CAT(AA(GETINTERSECTION(g))(epairs))

        node = g.addNode(3)
        verts = list(set(CAT(AA(DOWNCELLS(g))(edges+edgepairs))))
        for vert in verts: g.addArch(node,vert)
        for face in faces: g.addArch(face,node)

    return g


#/////////////////////////////////////////////////////////////////////
# Layered sphere 
#/////////////////////////////////////////////////////////////////////


def hexSphere(g,scaling=1.2):

    #/////////////////////////////////////////////////////////////////////
    #   STEP 1:  d-1  boundary-complex extraction
    #/////////////////////////////////////////////////////////////////////
    bComplex = boundaryComplex(g)
    n = g.getPointDim()
    DRAW(g,[1.5,1.5,1.5])(CAT(bComplex))

    #/////////////////////////////////////////////////////////////////////
    # STEP 2: d-2 wall-complex extraction
    #/////////////////////////////////////////////////////////////////////
    wallComplex = bComplex[:-1]
    DRAW(g,[1.5,1.5,1.5])(CAT(wallComplex))

    # STEP 3: basket separators construction \\\\\\\\\\\\\\\\
    # (d-1)-complex construction
    mapping = {}
    centroidMap = {}
    edges2downcells = {}

    #duplicate wall-complex
    for skeleton in wallComplex:

        #add 0-cells of mapped wall-complex cells
        for cell in skeleton:

            #Conserva un mapping fra gli spigoli e le downcells
            edges2downcells.update({cell: DOWNCELLS(g)(cell)})

            #conserva un mapping tra i centroidi degli spigoli ai livelli n-1 -> n
            upperCell = g.addNode(0)
            mapping.update({cell:upperCell})
            point = [CENTROID(g)(cell).get(i) for i in range(1,n+1)]

            g.setVecf(upperCell,Vecf([1.0]+SCALARVECTPROD([point,scaling])))
            lowerCell = g.addNode(0)
            g.setVecf(lowerCell,Vecf([1.0]+point))
            centroidMap.update({upperCell:lowerCell})

    DRAW(g, [1.5, 1.5, 1.5])()

    #add 1-cell extensions to top corners
    for node in wallComplex[0]:
        newNode = mapping[node]
        newArc = g.addNode(1)
        g.addArch(node,newArc)
        g.addArch(newNode,newArc)

    #add 1-cells of top-lateral boundary of added polytopes
    for cell in wallComplex[1]:

        #Get the top vertex the centroid is contained between
        vtx = [mapping[c] for c in edges2downcells[cell]]

        #Get the vertices to connect, something like this:
        #[vtx1, roof-centroid], [roof-centroid, vtx2]
        #roof-centroid is taken from the dict mapping, given
        #an edge(cell)
        pairs = zip(vtx, [mapping[cell]]*2)

        for pair in pairs:
            newArc = g.addNode(1)
            g.addArch(pair[0],newArc)
            g.addArch(pair[1],newArc)

    DRAW(g)()

    # In order to recognize the right pairs of centroids,
    # we use this algorithm:
    # 1. Iterate on every 2d cell (a facet)
    # 1.1 Compute the centroid of the facet
    # 2. Get the downcells of every cells, i.e. 4 edges
    # 3. Connect their centroids.
    
    wallComplex = bComplex
    facet2centroid = {}
    upper2lower = {}
    for facet in wallComplex[2]:

        #EXPERIMENTAL: Connette fra loro corner e lower-entroid. Deve
        #gestire gli spigoli duplicati!!
        l_corners_vtx = corners(g)(facet)
        u_corner_vtx = [mapping[vtx] for vtx in l_corners_vtx]
        upper_centroids = [mapping[edge] for edge in DOWNCELLS(g)(facet)]
        lower_centroids = [centroidMap[edge] for edge in upper_centroids]

        #Gestisce il mapping fra lo spigolo superiore <corner,centroid>
        #e lo spigolo inferiore <corner,centroid>

        for corner in u_corner_vtx:
            for centroid in upper_centroids:
                inters = GETINTERSECTION(g)([corner, centroid])

                if inters and inters[0] not in upper2lower.keys():
                    newArc = g.addNode(1)
                    corner_idx = u_corner_vtx.index(corner)
                    centroid_idx = upper_centroids.index(centroid)
                    g.addArch(l_corners_vtx[corner_idx], newArc)
                    g.addArch(lower_centroids[centroid_idx], newArc)
                    upper2lower.update({inters[0]:newArc})

        #END EXPERIMENTAL
        
        #aggiunta alla centroidMap i centroidi delle facce
        point = [CENTROID(g)(facet).get(i) for i in range(1,n+1)]
        upperCentroid = g.addNode(0)
        g.setVecf(upperCentroid, Vecf([1.0]+SCALARVECTPROD([point,scaling])))
        facet2centroid.update({facet: upperCentroid})
        lowerCentroid = g.addNode(0)
        g.setVecf(lowerCentroid, Vecf([1.0]+point))
        centroidMap.update({upperCentroid:lowerCentroid})
        edges = DOWNCELLS(g)(facet)

        for edge in edges:
            newArc = g.addNode(1)
            g.addArch(mapping[edge], newArc)
            g.addArch(upperCentroid, newArc)

        #Trying to connect the lower edges centroids to the central
        #lower face centroid. We need first to get, given an upper
        #edge its centroids, then go lower getting from centroidMap
        #the lower centroids and we go.
        for edge in edges:
            newArc = g.addNode(1)
            g.addArch(centroidMap[mapping[edge]], newArc)
            g.addArch(lowerCentroid, newArc)

    #connette i centroidi degli spigoli e delle facce
    for upper,lower in centroidMap.iteritems():
        newArc = g.addNode(1)
        g.addArch(lower, newArc)
        g.addArch(upper, newArc)

    DRAW(g)()

    #/////////////////////////////////////////////////////////////////////
    # STEP 3: d-2 and d-3 levels construction
    #/////////////////////////////////////////////////////////////////////

    # Per disegnare una faccia 2D ho bisogno di quattro spigoli.
    # Chiamiamo lower-corners i vertici "angolo" della generica faccia 2D
    # del livello n-1 su cui stiamo iterando.
    # Chiamiamo lower-centroids i centroidi degli spigoli del livello
    # n che stiamo costruendo
    # Chiamiamo upper{corners, centroids} i vertici del sottolivello superiore
    # ma sempre nel livello n che stiamo costruendo.
    # Per ogni coppia (upper corner, upper centroid) ci chiediamo
    # 1. Esiste uno ed un solo spigolo che li unisce?
    # 2. Se si andiamo a trovare gli altri spigoli con dei giochi sui
    # dizionari e connettiamo tutto.
    
    #Lista per gestire i duplicati
    not_dup_edges = []
    facet_count = 0
    for facet in wallComplex[2]:

        l_corners_vtx = corners(g)(facet)
        u_corners_vtx = [mapping[vtx] for vtx in l_corners_vtx]
        facet_edges = DOWNCELLS(g)(facet)
        u_edges_centroids = [mapping[edge] for edge in facet_edges]
        l_edges_centroids = [centroidMap[vtx] for vtx in u_edges_centroids]

        for corner in u_corners_vtx:
            for centroid in u_edges_centroids:

                e1 = GETINTERSECTION(g)([corner, centroid])
                corner_idx = u_corners_vtx.index(corner)
                e2 = GETINTERSECTION(g)([l_corners_vtx[corner_idx],
                                             centroidMap[centroid]])
                e3 = GETINTERSECTION(g)([corner, centroidMap[corner]])
                e4 = GETINTERSECTION(g)([centroid, centroidMap[centroid]])

                if len(e1) == len(e2) == len(e3) == len(e4) == 1:
                    if e1 not in not_dup_edges:
                        new_facet = g.addNode(2)
                        g.addArch(e1[0], new_facet); g.addArch(e2[0], new_facet)
                        g.addArch(e3[0], new_facet); g.addArch(e4[0], new_facet)
                        not_dup_edges.append(e1)

        #Adesso dobbiamo creare le facce "a croce", le 4 facce interne
        #al basket.
        upper_face_centroid = facet2centroid[facet]
        lower_face_centroid = centroidMap[upper_face_centroid]

        for centroid in u_edges_centroids:
            e1 = GETINTERSECTION(g)([upper_face_centroid, centroid])
            e2 = GETINTERSECTION(g)([centroid, centroidMap[centroid]])
            e3 = GETINTERSECTION(g)([centroidMap[centroid],
                                     lower_face_centroid])
            e4 = GETINTERSECTION(g)([upper_face_centroid,
                                     lower_face_centroid])

            if len(e1) == len(e2) == len(e3) == len(e4) == 1: 
                new_facet = g.addNode(2)
                g.addArch(e1[0], new_facet); g.addArch(e2[0], new_facet)
                g.addArch(e3[0], new_facet); g.addArch(e4[0], new_facet)


        #Ora e' il turno delle quattro lower facets e degli upper
        #Mi serve il prodotto cartesiano fra i lower centroid per poter
        #iterare su coppie di centroid.

        lower_centroids_cart = CART([l_edges_centroids, l_edges_centroids])
        lower_centroids_cart = filter(lambda x: x[0] < x[1] and x[0] != x[1],
                                      lower_centroids_cart)

        for corner in l_corners_vtx:
            for centroid_pair in lower_centroids_cart:

                e1 = GETINTERSECTION(g)([corner, centroid_pair[0]])
                e2 = GETINTERSECTION(g)([centroid_pair[0], lower_face_centroid])
                e3 = GETINTERSECTION(g)([lower_face_centroid, centroid_pair[1]])
                e4 = GETINTERSECTION(g)([centroid_pair[1], corner])

                if len(e1) == len(e2) == len(e3) == len(e4) == 1:

                    #Aggiungo una lower facet
                    new_facet = g.addNode(2)
                    g.addArch(e1[0], new_facet); g.addArch(e2[0], new_facet)
                    g.addArch(e3[0], new_facet); g.addArch(e4[0], new_facet)

                    #Se abbiamo trovato 4 spigoli buoni per le lower facet
                    #allora andiamo a prendere pure i punti superiori per
                    #disegnare direttamente le upper facet
                    facet_count += 1
                    upper_centroids_idx = l_edges_centroids.index(centroid_pair[0])
                    e1 = GETINTERSECTION(g)([mapping[corner],
                                             u_edges_centroids[upper_centroids_idx]])
                    e2 = GETINTERSECTION(g)([u_edges_centroids[upper_centroids_idx],
                                             upper_face_centroid])
                    upper_centroids_idx = l_edges_centroids.index(centroid_pair[1])
                    e3 = GETINTERSECTION(g)([u_edges_centroids[upper_centroids_idx],
                                             upper_face_centroid])
                    e4 = GETINTERSECTION(g)([mapping[corner],
                                             u_edges_centroids[upper_centroids_idx]])

                    new_facet = g.addNode(2)
                    g.addArch(e1[0], new_facet); g.addArch(e2[0], new_facet)
                    g.addArch(e3[0], new_facet); g.addArch(e4[0], new_facet)


        #Ora e' il turno dei cubi gialli, per cui abbiamo bisogno di
        # 8 vertici
        # uf = Upper face
        # lf = Lower face
        # l1 = lateral 1
        # l2 = lateral 2
        # l3 = lateral 3
        # l4 = lateral 4
        upper_centroids_cart = CART([u_edges_centroids, u_edges_centroids])
        upper_centroids_cart = filter(lambda x: x[0] < x[1] and x[0] != x[1],
                                      upper_centroids_cart)
        
        for centroid_pair in upper_centroids_cart:
            common_corner = get_corner_from(g)(centroid_pair)

            if common_corner[0] in u_corners_vtx:
                corner_idx = u_corners_vtx.index(common_corner[0])
                uf = [upper_face_centroid, common_corner[0], centroid_pair[0],
                      centroid_pair[1]]
                lf = [lower_face_centroid, l_corners_vtx[corner_idx],
                      centroidMap[centroid_pair[0]], centroidMap[centroid_pair[1]]]
                l1 = [lower_face_centroid, upper_face_centroid, centroid_pair[0],
                      centroidMap[centroid_pair[0]]]
                l2 = [lower_face_centroid, upper_face_centroid, centroid_pair[1],
                      centroidMap[centroid_pair[1]]]
                l3 = [centroid_pair[0], centroidMap[centroid_pair[0]],
                      common_corner[0], l_corners_vtx[corner_idx]]
                l4 = [centroid_pair[1], centroidMap[centroid_pair[1]],
                      common_corner[0], l_corners_vtx[corner_idx]]

                list_set = [uf, lf, l1, l2, l3, l4]
                faces = [get_facet_from(g)(l) for l in list_set]
                verts = CAT(list_set)
                node = g.addNode(3)
                for vert in verts: g.addArch(node,vert)
                for face in faces: g.addArch(face,node)

    assert(facet_count == len(wallComplex[2])*4)
    return g


#/////////////////////////////////////////////////////////////////////
# Local testing
#/////////////////////////////////////////////////////////////////////

if __name__=="__main__":

    g = Graph(3)
    g = initialSphere(g)
    DRAW(g,[1.5,1.5,1.5])()
    g = hexSphere(g)
    DRAW(g, [2.0, 2.0, 2.0])()


    #Attenzione! per estrarre correttamente il bordo al passo
    #successivo e' necessario disegnare PRIMA spigoli e vertici LOWER,
    #e poi gli UPPER
    bComplex = boundaryComplex(g)
    wallComplex = bComplex
    DRAW(g,[1.5,1.5,1.5])(CAT(wallComplex))
