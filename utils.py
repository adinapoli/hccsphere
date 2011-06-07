#Utilities functions
from hasselib.graphLib6 import *


def boundaryChain(g):
    myprint("g.getMaxDimCells()",g.getMaxDimCells())
    chain = CELLSPERLEVEL(g)(g.getMaxDimCells())
    out = [k for k in CAT(AA(DOWNCELLS(g))(chain)) if len(UPCELLS(g)(k))==1]
    return out


def boundaryComplex(g, chain = []):

    if not chain:
        chain = boundaryChain(g)
        
    out = [chain]
    for d in range(g.getMaxDimCells()-1):
        chain = list(set(CAT(AA(DOWNCELLS(g))(chain))))
        out.append(chain)
    return out[::-1]


def get_coords_from(g):
    def get_coords_from0(id):
        vect = g.getVecf(id)
        return [vect[0], vect[1], vect[2], vect[3]]

    return get_coords_from0


def id2vect(g):
    "Returns a dictionary like this {node_id: [x,y,z]}"

    def id2vect0(id_list):
        mapping = {}
        for id in id_list:
            mapping.update({id: get_coords_from(g)(id)[1:]})
        return mapping

    return id2vect0


def unique(alist):
    """Removes duplicates from a list."""

    no_dupes = []
    [no_dupes.append(i) for i in alist if not no_dupes.count(i)]
    return no_dupes


def angle_between(v1, v2):
    return ACOS(INNERPROD([v1, v2]) / (VECTNORM(v1)*VECTNORM(v2)))

def meridians(g):
    def meridians0(plane):
    
        boundary = boundaryComplex(g)[0]
        points = id2vect(g)(boundary)

        projected = [[k, ORTHOPROJECT(plane)(v)] for k,v in points.iteritems()]
        sorted_by_distance = sorted(projected,
                                    key= lambda x: VECTNORM(x[1]),
                                    reverse=True)

        #Punti a zero norma, ovvero quelli che proiettati su un piano
        #hanno distanza 0 dall'origine.
        norm_zero = []
        for id, vect in points.iteritems():
            if angle_between(vect, [0,0,1]) == 0 or \
               angle_between(vect, [0,0,1]) == PI:
                norm_zero.append(id)

        return norm_zero + [el[0] for el in sorted_by_distance]
    return meridians0


def meridians_edges(g):
    def meridians_edges0(id_list):
        edges = [GETINTERSECTION(g)(el) for el in CART([id_list, id_list])]
        edges = unique(CAT(filter(lambda x: len(x) == 1, edges)))

        return edges
    return meridians_edges0


def polar_edges(g):

    """Retrieves the edges that indice over the polar vertices"""
    result = []
    for cell in CELLSPERLEVEL(g)(1):
        vertices = [get_coords_from(g)(v) for v in DOWNCELLS(g)(cell)]
        for vtx in vertices:
            if vtx[1] == 0 and vtx[2] == 0:
                result.append(cell)

    return result

def corners(g):
    """Given a 2d-facet id, retrieves the four corner vertexes."""
    def corners0(facet_id):
        downcells = DOWNCELLS(g)(facet_id)

        corners = []
        for edge in downcells:
            corners.extend(DOWNCELLS(g)(edge))

        return unique(corners)
    return corners0


def get_facet_from(g):
    def get_facet_from0(id_list):
        pairs = CART([id_list, id_list])
        pairs = filter(lambda x: x[0] < x[1] and x[0] != x[1], pairs)
        inters = CAT(unique([GETINTERSECTION(g)(p) for p in pairs]))
        upcells = CAT([UPCELLS(g)(i) for i in inters])
        result = filter(lambda x: upcells.count(x) >= 3, upcells)[0]
        return result
    return get_facet_from0


#Attenzione, deve ritornare una lista di 1 solo elemento!
def get_corner_from(g):
    def get_corner_from0(id_list):
        edges = CAT([UPCELLS(g)(v) for v in id_list])
        edges = filter(lambda x: g.Level(x) == 1, edges)
        vtxes = CAT([DOWNCELLS(g)(e) for e in edges])
        result = filter(lambda x: vtxes.count(x) == 2, vtxes)
        return [result[0]]
    return get_corner_from0


def cmp_by_vecf(g):
    def cmp_by_vecf0(id):
        centroid = g.getVecf(id)
        return abs(centroid[1]) + abs(centroid[2]) + abs(centroid[3])
    return cmp_by_vecf0

