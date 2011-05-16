#Utilities functions
from hasselib.graphLib6 import *
from pyplasm import *


def boundaryChain(g):
    myprint("g.getMaxDimCells()",g.getMaxDimCells())
    chain = CELLSPERLEVEL(g)(g.getMaxDimCells())
    myprint("chain",chain)
    out = [k for k in CAT(AA(DOWNCELLS(g))(chain)) if len(UPCELLS(g)(k))==1]
    myprint("out",out)
    return out


def boundaryComplex(g):
    chain = boundaryChain(g)
    myprint("chain",chain)
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

def parallels(g):
    def parallels0(plane):
    
        boundary = boundaryComplex(g)[0]
        points = id2vect(g)(boundary)

        projected = [[k, DIRPROJECT(plane)(v)] for k,v in points.iteritems()]
        sorted_by_distance = sorted(projected,
                                    key= lambda x: VECTNORM(x[1]),
                                    reverse=True)

        #Punti a zero norma, ovvero quelli che proiettati su un piano
        #hanno distanza 0 dall'origine.
        norm_zero = []
        for id, vect in points.iteritems():
            if INNERPROD([vect, plane]) == 0:
                norm_zero.append(id)

        return norm_zero + [el[0] for el in sorted_by_distance]
    return parallels0


def parallel_edges(g):
    def parallel_edges0(id_list):
        edges = [GETINTERSECTION(g)(el) for el in CART([id_list, id_list])]
        edges = unique(CAT(filter(lambda x: len(x) == 1, edges)))

        return edges
    return parallel_edges0
  