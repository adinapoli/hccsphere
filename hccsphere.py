""" To generate the initial approximation of the sphere with hexahedra """

import numpy as np
from utils import *
import time

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

def signature(g):
    def signature0(vert):
        point = eval(str(g.getVecf(vert)))[1:]
        ret = []
        for x in point: ret += [SIGN(x) if x!=0 else 0]
        return ret
    return signature0

def test(g,v,w):
    return INNERPROD([signature(g)(v),signature(g)(w)])

def matching(g,vclass1,signtr):
    return [v for v in vclass1 if len([k for k in range(3)
        if signature(g)(v)[k]==signtr[k]]) == 2]

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
    edges += CAT([DISTL([v,[w for w in vclass[2] if test(g,v,w)==1]]) for v in vclass[1]])
    edges += CAT([DISTL([v,[w for w in vclass[1] if test(g,v,w)==2]]) for v in vclass[0]])

    for edge in edges:
         node = g.addNode(1); g.addArch(edge[0],node);g.addArch(edge[1],node)

    #//////////////////////////////////////////////////////////
    # generation of face nodes

    faces,o = [],vclass[3][0]

    for v in vclass[1]:
        for [u,w] in CART([vclass[2],vclass[2]]):
            if u<w and VECTSUM([ signature(g)(u),signature(g)(w) ])==signature(g)(v):
                faces += [CAT(AA(GETINTERSECTION(g))([[v,u],[v,w],[o,u],[o,w]]))]

    for vert in vclass[0]:
        u,v,w = matching(g,vclass[1],signature(g)(vert))
        vertTriples = [[u,v,vert],[v,w,vert],[w,u,vert]]
        for k in range(3):
            signatures = AA(signature(g))(vertTriples[k])
            common = COMP([AA(PROD),TRANS])(signatures)
            w = [v for v in vclass[2] if signature(g)(v) == common]+vertTriples[k]
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

CURRENT_EXTERNAL_FACETS = []

def impose_layer(g, layer = 1, scaling=1.2):

    offset = (scaling*layer+1)/float(scaling)

    global CURRENT_EXTERNAL_FACETS
    #/////////////////////////////////////////////////////////////////////
    #   STEP 1: boundary-complex extraction
    #/////////////////////////////////////////////////////////////////////
    wallComplex = boundaryComplex(g, CURRENT_EXTERNAL_FACETS)
    n = g.getPointDim()

    #/////////////////////////////////////////////////////////////////////
    # STEP 2: iterations over 2d facets
    #/////////////////////////////////////////////////////////////////////

    edge2UpCentroid = {}
    edge2LwCentroid = {}
    #centroidMap = {}
    edges2downcells = {}
    cornersLow2Up = {}
    facet2UpCentroid = {}
    facet2LwCentroid = {}
    upper2lower = {}

    for facet in wallComplex[2]:

        #/////////////////////////////////////////////////////////////////////
        # STEP 2.1: vertices creation
        #/////////////////////////////////////////////////////////////////////

        #Per prima cosa estrudo i corner
        corners_vtx = corners(g)(facet)

        for vtx in corners_vtx:
            if vtx not in cornersLow2Up.keys():

                upperCorner = g.addNode(0)
                cornersLow2Up.update({vtx:upperCorner})
                point = get_coords_from(g)(vtx)[1:]
                g.setVecf(upperCorner,Vecf([1.0]+SCALARVECTPROD([UNITVECT(point),offset])))

                upper2lower.update({upperCorner:vtx})
                #Connessione tra i corner
                newArc = g.addNode(1)
                g.addArch(vtx,newArc)
                g.addArch(upperCorner,newArc)


        #Ora creo i centroidi
        facet_edges = DOWNCELLS(g)(facet)

        for edge in facet_edges:

            #Gestione duplicati
            if edge not in edges2downcells.keys():
                #Conserva un mapping fra gli spigoli e le downcells
                edges2downcells.update({edge: DOWNCELLS(g)(edge)})

                #Conserva un mapping tra i centroidi degli spigoli
                #ai livelli n-1 -> n
                upperCell = g.addNode(0)
                edge2UpCentroid.update({edge: upperCell})
                point = [CENTROID(g)(edge).get(i) for i in range(1,n+1)]
                g.setVecf(upperCell,Vecf([1.0]+SCALARVECTPROD([UNITVECT(point),offset])))

                lowerCell = g.addNode(0)
                edge2LwCentroid.update({edge: lowerCell})
                g.setVecf(lowerCell,Vecf([1.0]+point))

                upper2lower.update({upperCell:lowerCell})
                #connessione fra upper e lower
                newArc = g.addNode(1)
                g.addArch(upperCell,newArc)
                g.addArch(lowerCell,newArc)

                #Grazie al mapping edge-downcells so quali sono le downcells
                #e quindi i corner! Posso creare subito gli edge che collegano
                #i corner ai centroidi!
                corn = edges2downcells[edge]

                #Creo gli spigoli che collegano i lower corner con il
                #lower centroid
                for corner in corn:
                    newArc = g.addNode(1)
                    g.addArch(lowerCell, newArc)
                    g.addArch(corner, newArc)

                for corner in corn:
                    newArc = g.addNode(1)
                    g.addArch(cornersLow2Up[corner], newArc)
                    g.addArch(upperCell, newArc)


        #Infine e' il turno dei centroidi delle facce 2d
        point = [CENTROID(g)(facet).get(i) for i in range(1,n+1)]
        upperCentroid = g.addNode(0)

        g.setVecf(upperCentroid, Vecf([1.0]+SCALARVECTPROD([UNITVECT(point),offset])))
        facet2UpCentroid.update({facet: upperCentroid})

        lowerCentroid = g.addNode(0)
        g.setVecf(lowerCentroid, Vecf([1.0]+point))
        facet2LwCentroid.update({facet: lowerCentroid})

        upper2lower.update({upperCentroid:lowerCentroid})
        #connessione upper-lower
        newArc = g.addNode(1)
        g.addArch(upperCentroid,newArc)
        g.addArch(lowerCentroid,newArc)

        #Connetto subito i centroidi delle facce 2d ai centroidi, sia
        #lower che upper!

        u_centroid_list = [edge2UpCentroid[e] for e in facet_edges]
        for c in u_centroid_list:
            newArc = g.addNode(1)
            g.addArch(upperCentroid, newArc)
            g.addArch(c, newArc)

        l_centroid_list = [edge2LwCentroid[e] for e in facet_edges]
        for c in l_centroid_list:
            newArc = g.addNode(1)
            g.addArch(lowerCentroid, newArc)
            g.addArch(c, newArc)


    CURRENT_EXTERNAL_FACETS = []
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
        u_corners_vtx = [cornersLow2Up[vtx] for vtx in l_corners_vtx]
        facet_edges = DOWNCELLS(g)(facet)
        u_edges_centroids = [edge2UpCentroid[edge] for edge in facet_edges]
        l_edges_centroids = [upper2lower[vtx] for vtx in u_edges_centroids]

        for corner in u_corners_vtx:
            for centroid in u_edges_centroids:

                e1 = GETINTERSECTION(g)([corner, centroid])
                corner_idx = u_corners_vtx.index(corner)
                e2 = GETINTERSECTION(g)([l_corners_vtx[corner_idx],
                                             upper2lower[centroid]])
                e3 = GETINTERSECTION(g)([corner, upper2lower[corner]])
                e4 = GETINTERSECTION(g)([centroid, upper2lower[centroid]])

                if len(e1) == len(e2) == len(e3) == len(e4) == 1:
                    if e1 not in not_dup_edges:
                        new_facet = g.addNode(2)
                        g.addArch(e1[0], new_facet); g.addArch(e2[0], new_facet)
                        g.addArch(e3[0], new_facet); g.addArch(e4[0], new_facet)
                        not_dup_edges.append(e1)

        #Adesso dobbiamo creare le facce "a croce", le 4 facce interne
        #al basket.
        upper_face_centroid = facet2UpCentroid[facet]
        lower_face_centroid = upper2lower[upper_face_centroid]

        for centroid in u_edges_centroids:
            e1 = GETINTERSECTION(g)([upper_face_centroid, centroid])
            e2 = GETINTERSECTION(g)([centroid, upper2lower[centroid]])
            e3 = GETINTERSECTION(g)([upper2lower[centroid],
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
                    e1 = GETINTERSECTION(g)([cornersLow2Up[corner],
                                             u_edges_centroids[upper_centroids_idx]])
                    e2 = GETINTERSECTION(g)([u_edges_centroids[upper_centroids_idx],
                                             upper_face_centroid])
                    upper_centroids_idx = l_edges_centroids.index(centroid_pair[1])
                    e3 = GETINTERSECTION(g)([u_edges_centroids[upper_centroids_idx],
                                             upper_face_centroid])
                    e4 = GETINTERSECTION(g)([cornersLow2Up[corner],
                                             u_edges_centroids[upper_centroids_idx]])

                    new_facet = g.addNode(2)
                    CURRENT_EXTERNAL_FACETS.append(new_facet)
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
                      upper2lower[centroid_pair[0]], upper2lower[centroid_pair[1]]]
                l1 = [lower_face_centroid, upper_face_centroid, centroid_pair[0],
                      upper2lower[centroid_pair[0]]]
                l2 = [lower_face_centroid, upper_face_centroid, centroid_pair[1],
                      upper2lower[centroid_pair[1]]]
                l3 = [centroid_pair[0], upper2lower[centroid_pair[0]],
                      common_corner[0], l_corners_vtx[corner_idx]]
                l4 = [centroid_pair[1], upper2lower[centroid_pair[1]],
                      common_corner[0], l_corners_vtx[corner_idx]]

                list_set = [uf, lf, l1, l2, l3, l4]
                faces = [get_facet_from(g)(l) for l in list_set]
                verts = CAT(list_set)
                node = g.addNode(3)
                for vert in verts: g.addArch(node,vert)
                for face in faces: g.addArch(face,node)

    assert(facet_count == len(wallComplex[2])*4)
    return g


def hexsphere(layers = 2, scaling = 1.2):
    g = Graph(3)
    g = initialSphere(g)

    for i in range(1,layers+1):
        g = impose_layer(g, i, scaling)

    return g


#/////////////////////////////////////////////////////////////////////
# Local testing
#/////////////////////////////////////////////////////////////////////

if __name__=="__main__":

    start = time.clock()
    g = hexsphere(3)
    end = time.clock()

    print "Sphere builded in ", end - start, " seconds."
    VOLVIEW(g)