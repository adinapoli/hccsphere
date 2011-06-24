from config import *

#//////////////////////////////////////////////
# Usando questa si velocizza il caricamento dei 
# batch
#//////////////////////////////////////////////
def get_sphere_batch():
    yield Batch(SPHERE_BATCH)
    
def get_cylinder_batch():
    yield Batch(CYLINDER_BATCH)


# /////////////////////////////////////////////
# celle dw di una cella
# /////////////////////////////////////////////
def DOWNCELLS(g):
    def DOWNCELLS0(cell):
        it=g.goDw(cell);ret=[]
        while not it.end(): ret+=[it.getNode()];it.goForward()
        return ret
    return DOWNCELLS0

# /////////////////////////////////////////////
# celle up di una cella
# /////////////////////////////////////////////
def UPCELLS(g):
    def UPCELLS0(cell):
        it=g.goUp(cell);ret=[]
        while not it.end(): ret+=[it.getNode()];it.goForward()
        return ret
    return UPCELLS0

# /////////////////////////////////////////////
#  celle ad un certo livello del grafo
# /////////////////////////////////////////////
def CELLSPERLEVEL(g):
    def CELLSPERLEVEL0(level):
        it=g.each(level) 
        ret=[]
        while not it.end():
            ret+=[it.getNode()];it.goForward()
        return REVERSE(ret)  
    return CELLSPERLEVEL0


# /////////////////////////////////////////////
# ad esempio per trovare tutte le 0-celle 
# /////////////////////////////////////////////
def FINDCELLS(g,nav):
    def FINDCELLS0(level,cell):
        num=g.findCells(level,cell,nav)
        return [nav.getCell(level,N) for N in range(0,num)] 
    return FINDCELLS0

# ////////////////////////////////////////////////
def GETINTERSECTION(g):
    def GETINTERSECTION0(from_cells):

        result = []
        getNode = g.getNode
        num = len(from_cells)

        # use the Tmp info to reset the "reachable" info
        for c in from_cells:
            it = g.goUp(c)
            while not it.end(): 
                reached = getNode(it.getNode())
                reached.Tmp = 0
                it.goForward()
            
        # increment the "reachable" Tmp, if all nodes reaches one common father
        # then the father is added to the list
        for c in from_cells:
            it = g.goUp(c)
            while not it.end(): 
                reached = getNode(it.getNode())
                reached.Tmp += 1
                if reached.Tmp==num: 
                    result+=[it.getNode()]
                it.goForward()
            
        return result
    return GETINTERSECTION0

###################################################################

# /////////////////////////////////////////////
# Generation of 1D polyhedron
# /////////////////////////////////////////////

def Quote(g,numList):
    """ To create the graph of a 1D cell complex """
    sizes = [abs(num) for num in numList]
    points = [Vecf([1.0, x]) for x in AL([0,PROGRESSIVESUM(sizes)])]
    for point in points: node = g.addNode(0); g.setVecf(node,point)
    nodes = CELLSPERLEVEL(g)(0)
    edges = [[nodes[k],nodes[k+1]] for k in range(len(nodes)-1)]
    for edge in edges: g.addNode(1)
    nodes = CELLSPERLEVEL(g)(1)
    for k in range(len(edges)):
        if numList[k] > 0:
            g.addArch(edges[k][0],nodes[k]);g.addArch(edges[k][1],nodes[k])
            #aggiungi la doppia connettivita' dei nodi a livello top
            g.addArch(nodes[k],edges[k][0]);g.addArch(nodes[k],edges[k][1])
        else: g.remNode(nodes[k])
    return g
    

# /////////////////////////////////////////////
# Generation of nD grids
# /////////////////////////////////////////////

def Grid(listOfListsOfNum):
    numList = listOfListsOfNum[0]
    g = Quote(Graph(1),numList)
    for numList in listOfListsOfNum[1:]:
        g1 = Quote(Graph(1),numList)
        g = Graph.power(Matf(1),Matf(1),  g,None,None,  g1,None,None)
    return g


# /////////////////////////////////////////////
# Centroid operator
# /////////////////////////////////////////////


def CENTROID(g):
    def CENTROID0(cell):
        """ Compute the centroid of a cell.
            Returns a value of class pyplasm.xge.xgepy.Vecf
        """
        if g.Level(cell) == 0: return g.getVecf(cell)
        out = DOWNCELLS(g)(cell)
        while g.Level(out[0]) > 0:
            out = list(set(CAT(AA(DOWNCELLS(g))(out))))
        points = [g.getVecf(cell) for cell in out]
        centroid = VECTSUM([[point.get(i) for i in range(g.getPointDim()+1)]
                            for point in points])
        centroid = SCALARVECTPROD([1/centroid[0],centroid])
        return Vecf(centroid)
    return CENTROID0


###################################################################

# /////////////////////////////////////////////
# Drawing 3-cells of Graph instances with Hpc
# /////////////////////////////////////////////

def graph2hexs(g):
    def graph2hexs0(chain):
        cells = [UPCELLS(g)(cell) for cell in chain]
        cellverts = [[[g.getVecf(node).get(i) for i in range(1,4)]
                 for node in cell] for cell in cells]
        return [MKPOL([verts,[range(1,len(verts)+1)],None]) for verts in cellverts]
    return graph2hexs0

def SHOW(g,expl=[1,1,1]):
    VIEW(EXPLODE(*expl)(graph2hexs(g)(CELLSPERLEVEL(g)(3))))


###################################################################

# /////////////////////////////////////////////
# Drawing Graph instances with spheres and cylinders
# /////////////////////////////////////////////


def DRAW(g,expl=[1,1,1]):

    n = g.getMaxDimCells()
    m = g.getPointDim()

    def offset(point,expl=[1,1,1]):
        scaledpoint = [point[k]*expl[k] for k in range(3)]
        vect = VECTDIFF([scaledpoint,point])
        return vect
        

    def spheres(points,expl=[1,1,1]):
        batches = []
        sx = 0.05
        #unitSphere = get_sphere_obj()
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
        #unitCylinder = get_cylinder_obj()      
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
        #myprint("chains",chains)
        nodepoints = [[g.getVecf(node)[i] for i in range(1,m+1)]
                      if m>2 else
                      [g.getVecf(node)[i] for i in range(1,m+1)]+[0.0]
                      for node in CELLSPERLEVEL(g)(0)]
        #myprint("nodepoints",nodepoints)

        def translate(pointIds):
            #myprint("pointIds",pointIds)
            return [nodepoints[mapping[0][vert[1]]] for vert in 
                          enumerate(pointIds)]
            
        if m==2: m+=1
        if chains[0] != []: vertpoints = translate(chains[0])
        if chains[1] != []:
            edges = [DOWNCELLS(g)(edge) for edge in chains[1]]
            #myprint("edges",edges)
            edgepoints = AA(translate)(edges)
            #myprint("edgepoints",edgepoints)
        if chains[2] != []:
            facesAsEdges = [DOWNCELLS(g)(face) for face in chains[2]]
            facesAsVerts = [list(set(CAT(AA(DOWNCELLS(g))(face))))
                            for face in facesAsEdges]
            facepoints = AA(translate)(facesAsVerts)
        if d == 3:
            if chains[3] != []:
                solidsAsVerts = [UPCELLS(g)(cell) for cell in chains[3]]
                cellpoints = AA(translate)(solidsAsVerts)

        #this is the list of batches you want to display
        batches = []
        if chains[0] != []:
            batches = spheres(vertpoints,expl)
        if chains[1] != []:
            batches = cylinders(batches,edgepoints,expl)
        if chains[2] != []:
            batches = planecells(batches,facepoints,expl)
        if n == 3:
            if chains[3] != []: batches = cells(batches,cellpoints,expl)

        # organize the batch in a loose octree
        octree=Octree(batches)

        # create the viewer and run it
        viewer=Viewer(octree)
        viewer.Run()

    return DRAW0

###################################################################

# ////////////////////////////////////////////////////////
# Meshing with d-hypercubes (HCC - Hyper-Cuboidal-Complex)
# ///////////////////////////////////////////////////////

def hccmesh(g):
    
    d = g.getMaxDimCells()
    g1 = Graph(d)


    def getfathers(tuple):
        return [GETINTERSECTION(g1)([tuple[i],tuple[i+1]])
                for i in range(len(tuple)-1)]
    
    def grouping(tuples):
        out, n = [], len(tuples)
        groups = [[tuples[0]]]
        first,last = tuples[0][0],tuples[0][-1]
        for k in range(1,len(tuples)):
            if (tuples[k][0],tuples[k][-1]) == (first,last):
                groups[-1].append(tuples[k])
            else:
                groups.append([tuples[k]])
                first,last = tuples[k][0],tuples[k][-1]
        return groups

    def DOWNTRAVERSE(g,nrecursion,cell,up_direction=False):
        def multiTraverse(g,nrecursion,cell,up_direction):
            if nrecursion == 0: return [[cell]]
            ret = []
            for Down in DOWNCELLS(g)(cell) if not up_direction else UPCELLS(g)(cell):
                ret += [[cell] + L for L in
                        multiTraverse(g,nrecursion-1,Down,up_direction)]
            return ret
        return grouping(sorted(AA(REVERSE)(multiTraverse(g,nrecursion,cell,up_direction))))   

    def UPTRAVERSE(subgraph,k):            
        if k == 1: return CAT(subgraph)
        else: return COMP([list,set,CAT,AA(COMP([CAT,getfathers]*(k-1)))])(subgraph)
    
    d = g.getMaxDimCells()
    g1 = Graph(d)
    # create the 0-layer of HCC
    for node in range(1,g.getNumNode()+1):
        newnode = g1.addNode(0)
        g1.setVecf(newnode,CENTROID(g)(node))
    # create higher level layers of HCC
    for k in range(1,d):
        # create Nk and Ak
        for h in range(k,d+1):
            for root in CELLSPERLEVEL(g)(h):
                # forward search in g for the isomorphic subgraphs
                subgraphs = DOWNTRAVERSE(g,k,root)
                # backtrack upon g1: looking for (k-1)-faces of each newnode
                faces = [UPTRAVERSE(sg,k) for sg in subgraphs]
                # build the k-layer of HCC
                for face in faces:
                    newnode = g1.addNode(k)
                    for node in face:
                        g1.addArch(node,newnode)                            
    # create the last layer of HCC
    for root in CELLSPERLEVEL(g)(d):
        subgraphs = DOWNTRAVERSE(g,d,root)
        facets = [UPTRAVERSE(sg,d) for sg in subgraphs]
        vertices = [list(set(CAT(sg))) for sg in subgraphs]
        # build the d-layer of HCC
        for facet,verts in zip(facets,vertices):
            newnode = g1.addNode(d)
            for node in facet:
                g1.addArch(node,newnode)
            for vert in verts:
                g1.addArch(newnode,vert)                                                                
    # return the output HCC graph
    return g1

def hccMeshSize(d):
    g = Grid([[1]]*d)
    g1 = hccmesh(g)
    return [len(CELLSPERLEVEL(g1)(k)) for k in range(d+1)]

if __name__ == "__main__":

    g = Grid([[1],[1],[1]])
    DRAW(g,[2.0,2.0,2.0])()
    g = hccmesh(g)
    DRAW(g,[2.0,2.0,2.0])()
    g = hccmesh(g)
    DRAW(g,[2.0,2.0,2.0])()


    
    g = Grid([[1,-1,1],[1,-1,1]])
    g = Grid([[1,-1,1],[1,-1,1],[1,-1,1]])

    g = hccmesh(g)
    VIEW(SKELETON(1)(Hpc(g)))
    SHOW(g,[1.5,1.5,1.5])
    DRAW(g,[1.5,1.5,1.5])()
    
##    myprint("hccMeshSize(1)",hccMeshSize(1))
##    myprint("hccMeshSize(2)",hccMeshSize(2))
##    myprint("hccMeshSize(3)",hccMeshSize(3))
##    myprint("hccMeshSize(4)",hccMeshSize(4))
##    myprint("hccMeshSize(5)",hccMeshSize(5))
    
    g = hccmesh(g)
    SHOW(g,[1.5,1.5,1.5])

