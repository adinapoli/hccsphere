from ballorganizer import *
from Bio.PDB import *
from atomic_radius import *

def myprint (string):
    print "\n" + string + " ->", eval(string)

parser=PDBParser()
structure=parser.get_structure('atp', 'ATP.pdb')
myprint("structure")
model = structure[0]
myprint("model")
chain = model[' ']
myprint("chain")
residue = chain[0] 
myprint("residue")

def getPDBconnect (filename):
    myfile = open(filename,'r')
    arcs = []
    for record in myfile:
        terms = record.split()
        if terms[0] == "CONECT":
            pairs = DISTL([ eval(terms[1]), AA(eval)(terms[2:]) ])
            arcs += [arc for arc in pairs if arc[0] < arc[1]]
    myfile.close()
    return arcs

#myprint("getPDBconnect('BTC.pdb')")


def graph (filename):
    parser=PDBParser()
    structure=parser.get_structure('molecule', filename)
    model = structure[0]
    chain = model[' ']
    residue = chain[0]
    
    nodes = [atom.get_coord().tolist() for atom in residue]
    edges = getPDBconnect(filename)
    atoms = list([atom.get_id()[0] for atom in residue])
    return nodes,edges,atoms

#myprint("graph('BTC.pdb')")
#myprint("list([atom.get_id()[0] for atom in residue])")


molecule = graph('ATP.pdb')
for node,atom in zip(molecule[0],molecule[2]):
	myprint("atomic_radius[atom][RADIUS_TYPE]/100,Vec3f(node)")
	
	
moleculeBallSet = BallSet()
for node,atom in zip(molecule[0],molecule[2]):
	# ball=Ball3f(atomic_radius[atom][RADIUS_TYPE]/100.,Vec3f(node))
	ball=Ball3f(atomic_radius[atom][RADIUS_TYPE]/(100.*2),Vec3f(node))
	moleculeBallSet.addBall(ball)
moleculeBallSet.createOctree() 

#myprint("moleculeBallSet")

balls_ = []
circles = 0
for ball in moleculeBallSet.balls:
	balls_ += [ball]
	myprint("(ball,ball.ID)")
	myprint("moleculeBallSet.findIntersectingBalls(ball)")
	circles += len(moleculeBallSet.findIntersectingBalls(ball))
	
	
if __name__ == "__main__":
		
	myprint("circles")
	g=buildGraphFromBalls(balls_)
	Plasm.View(SKELETON(1)(Hpc(g)))
	
# for edge in myGraph[1]:
# 	myprint("VECTNORM(VECTDIFF([myGraph[0][edge[1]-1], myGraph[0][edge[0]-1]]))")