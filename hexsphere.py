""" To generate the initial approximation of the sphere with hexahedra """

import numpy as np
from hasselib.graphLib6 import *

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

    if True:
        assert(len(CELLSPERLEVEL(g)(0)) == 27)


    vclass = classifyVerts(g)

    #//////////////////////////////////////////////////////////
    # generation of nodes of 3D cells
    def cell(args):
        i,j,k = args
        return [index+1 for index,v in enumerate(verts) if AND([
            (GE(0) if i==1 else LE(0))(v[1]),
            (GE(0) if j==1 else LE(0))(v[2]),
            (GE(0) if k==1 else LE(0))(v[3]) ])  ]

    vertices = [v[1:] for v in verts]
    pols = AA(cell)(CART([[1,-1],[1,-1],[1,-1]]))

    for pol in pols:
        node = g.addNode(3);
        for vert in pol: g.addArch(node,vert)

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

    #Con questo for vengono create le facce interne, sono 12.
    for v in vclass[1]:
        for [u,w] in CART([vclass[2],vclass[2]]):
            if u<w and VECTSUM([ signature(u),signature(w) ])==signature(v):
                faces += [CAT(AA(GETINTERSECTION(g))([[v,u],[v,w],[o,u],[o,w]]))]

    
    #Con questo for vengono create le facce esterne, 24
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

    return g


def graph2batches(g,expl=[1,1,1]):

    n = g.getMaxDimCells()
    m = g.getPointDim()

    def offset(point,expl=[1,1,1]):
        scaledpoint = [point[k]*expl[k] for k in range(3)]
        vect = VECTDIFF([scaledpoint,point])
        return vect


    def spheres(points,expl=[1,1,1]):
        batches = []
        sx = 0.05

        for point in points:
            batchSphere = get_sphere_batch().next()
            vect = offset(point,expl)
            if len(point) == 2:
                point = point + [0.0]
                vect = vect + [0.0]
            batchSphere.matrix =  Mat4f.translate(*vect) * \
                Mat4f.translate(*point)*Mat4f.scale(sx,sx,sx)
            batchSphere.diffuse=CYAN
            batches += [batchSphere]
        return batches


    def transfCylr(batchCylinder,pointpair,expl=[1,1,1]):
        vect,point = VECTDIFF(REVERSE(pointpair)),pointpair[0]
        sx = 0.025

        def vectTransform(vect):
            qz = UNITVECT(vect)
            qx = UNITVECT(VECTPROD([ vect,[0,0,1] ]))
            qy = VECTPROD([ qz,qx ])
            Rot = TRANS([qx,qy,qz])
            Rot = CAT([ Rot[0]+[0.], Rot[1]+[0.], Rot[2]+[0.], [0.,0.,0.,1.] ])
            h = VECTNORM(vect)

            def isclose (a,b,filter_threshold=0.5):
                if abs(a-b)<filter_threshold: return True
                else: return False

            if isclose (Mat4f.determinant(Mat4f(*Rot)),
                        0.0, 1E-5):
                return h,Mat4f.scale(1,SIGN(vect[1]),SIGN(vect[2]))
            else: return h,Mat4f(*Rot)

        h,rot = vectTransform(vect)
        center = [c/2.0 for c in VECTSUM(pointpair)]
        vect = offset(center,expl)
        batchCylinder.matrix = Mat4f.translate(*vect) *\
            Mat4f.translate(*point) * rot * Mat4f.scale(sx,sx,h)
        batchCylinder.diffuse = MAGENTA
        return batchCylinder

    def cylinders(batches,edgepoints,expl=[1,1,1]):

        vects = [VECTDIFF(edge) for edge in edgepoints]
        for pointpair in edgepoints:
            batchCyl = get_cylinder_batch().next()
            batchCyl = transfCylr(batchCyl,pointpair,expl)
            batches += [batchCyl]
        return batches

    def planecells(batches,facepoints,expl=[1,1,1]):
        for points in facepoints:
            n = len(points)
            center = [coord/float(n) for coord in VECTSUM(points)]
            vect = offset(center,expl)
            points = [[point[k]+vect[k] for k in range(3)] for point in points]
            def sign(points):
                return SIGN(VECTPROD(AA(C(VECTDIFF)(center))(points[2:0:-1])))
            face = MKPOL([points,[range(1,n+1)],None])
            faceBatch = Plasm.getBatches(face)
            faceBatch[0].diffuse = WHITE
            batches += faceBatch
        return batches

    def cells(batches,cellpoints,expl=[1,1,1]):
        for points in cellpoints:
            n = len(points)
            center = [coord/float(n) for coord in VECTSUM(points)]
            vect = offset(center,expl)
            points = [[point[k]+vect[k] for k in range(3)] for point in points]
            cell = MKPOL([points,[range(1,n+1)],None])
            cellBatch = Plasm.getBatches(cell)
            cellBatch[0].diffuse = YELLOW
            batches += cellBatch
            # view rotation
            rot = ROTN([ ACOS(INNERPROD([ [1,1,1],[0,0,1] ])), VECTPROD([ [1,1,1],[0,0,1] ]) ])
            batches += Plasm.getBatches(STRUCT([rot, MK([1,1,1])]))
        return batches

    def DRAW0(chain=range(1,g.getNumNode()+1)):

        m = g.getPointDim()
        d = g.getMaxDimCells()
        mapping = [dict(zip(CELLSPERLEVEL(g)(k),range(len(CELLSPERLEVEL(g)(k)))))
                   for k in range(d+1)]

        chains = [[],[],[],[]]
        [chains[g.Level(node)].append(node) for node in chain[::-1]]

        nodepoints = [[g.getVecf(node)[i] for i in range(1,m+1)]
                      if m>2 else
                      [g.getVecf(node)[i] for i in range(1,m+1)]+[0.0]
                      for node in CELLSPERLEVEL(g)(0)]


        def translate(pointIds):

            return [nodepoints[mapping[0][vert[1]]] for vert in
                          enumerate(pointIds)]

        if m==2: m+=1
        if chains[0]: vertpoints = translate(chains[0])
        if chains[1]:
            edges = [DOWNCELLS(g)(edge) for edge in chains[1]]

            edgepoints = AA(translate)(edges)

        if chains[2]:
            facesAsEdges = [DOWNCELLS(g)(face) for face in chains[2]]
            facesAsVerts = [list(set(CAT(AA(DOWNCELLS(g))(face))))
                            for face in facesAsEdges]
            facepoints = AA(translate)(facesAsVerts)
        if d == 3 and chains[3]:
            solidsAsVerts = [UPCELLS(g)(cell) for cell in chains[3]]
            cellpoints = AA(translate)(solidsAsVerts)


        batches = []
        if chains[0]:
            batches = spheres(vertpoints,expl)
        if chains[1]:
            batches = cylinders(batches,edgepoints,expl)
        if chains[2]:
            batches = planecells(batches,facepoints,expl)
        if n == 3 and chains[3]:
            batches = cells(batches,cellpoints,expl)

        return batches

    return DRAW0


def PEEK(batches):
    """Like DRAW and VIEW, but more general.

    Display batches.
    """
    if not isinstance(batches,list):
        batches = [batches]

    octree = Octree(Batch.Optimize(batches))
    viewer = Viewer(octree)
    viewer.Run()

#/////////////////////////////////////////////////////////////////////
# Local testing 
#/////////////////////////////////////////////////////////////////////

if __name__=="__main__":

    g = Graph(3)
    hull = SKELETON(1)(SPHERE(1)([8,8]))
    kernel = initialSphere(g)
    batches = graph2batches(kernel, [1.5, 1.5, 1.5])(range(36,126))
    batches += Plasm.getBatches(hull)
    #DRAW(kernel,[1.5,1.5,1.5])([2,1,6])
    #out = MKPOL([vertices, pols, None])
    #VIEW(SKELETON(1)(out))
    PEEK(batches)
    print classifyVerts(kernel)

