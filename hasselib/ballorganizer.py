from pyplasm import *
import numpy
import math

# //////////////////////////////////////////////////////////////////
class BallSet:

	#__________________________________________________________________________
	# costruttore
	def __init__(self,max_depth=8):
		self.max_depth = max_depth
		self.ballmap={}         # ogni bounding ball verra' inserita in un python dictionary 
		self.world_box=Box3f()	# questo e' il bounding box globale di tutte le bounding ball
		self.balls=[]           # mi tengo una copia delle balls per l'inserimento successivo

		
	#__________________________________________________________________________
	# aggiunge una bounding ball, call addBall(None) when all balls have been added
	def addBall(self,ball):
	
		# finished adding all the balls
		if ball is None:
		
			print "Creating loose octree with world_box",str(self.world_box)," max_depth",self.max_depth	
			self.octree=Octree(self.world_box,self.max_depth)
			 
			for ball in self.balls:
				R,C=ball.radius,ball.center
				box=Box3f(C-Vec3f([R,R,R]),C+Vec3f([R,R,R]))
				ID=self.octree.getNode(box).getId()
				# assign each ball to it's slot of the loose octree
				if not self.ballmap.has_key(ID): self.ballmap[ID]=[]
				self.ballmap[ID]+=[ball]
				
			return
			
		#assign an unique id
		ball.ID=len(self.balls)
		
		# is a bounding ball to add
		R,C=ball.radius,ball.center
		# note: the inner loose octree works with box instead of balls 
		self.world_box.add(C-Vec3f([R,R,R]))
		self.world_box.add(C+Vec3f([R,R,R]))
		self.balls+=[ball]

	#__________________________________________________________________________
	# test if two balls intersect
	def intersect(self,ball1,ball2):
		return (ball1.center - ball2.center).module()<=(ball1.radius+ball2.radius)	

	#__________________________________________________________________________
	# find all balls intersecting with the given ball
	def intersection(self,ball,node=None):

		# start with the root node of the octree	
		if node is None:
			node=self.octree.root
			
		R,C=ball.radius,ball.center
		box=Box3f(C-Vec3f([R,R,R]),C+Vec3f([R,R,R]))
			
		ret=[]
	
		# il nodo attuale ha un bounding box che non ha intersezione
		if not node.box.intersection(box).isValid(): 
			return ret
			
		# vedi se qualche ball ha intersezione
		id=node.getId()
		
		if self.ballmap.has_key(id):
			for B in self.ballmap[id]:
				if self.intersect(ball,B): ret+=[B]
				
		# vedi anche nei figli (se esistono)
		for I in range(0,8):
			child=node.getChild(I)
			if child: ret+=self.intersection(ball,child)

		return ret
		
		
		
		
# //////////////////////////////////////////////////////////////////
# esempio di creazione di un pol complex a partire dalle balls
# //////////////////////////////////////////////////////////////////
def MkPolOfBallSet(ballset,STEP=16,TOLERANCE=1e-6):

	
	# this are the points of a ball centered in the origin with radius 1
	POINTS=[Vec3f(math.sin(alpha)*math.cos(beta),math.sin(alpha)*math.sin(beta),math.cos(alpha))
                for alpha in numpy.arange(math.pi/24,math.pi*(23./24),math.pi*(22./24)/STEP)
                for beta in numpy.arange(0,2*math.pi,2*math.pi/STEP)]
	
	# _________________ STEP 1 per ogni (bi,bj) intersecting find cutting plane h
	for ball in ballset.balls: ball.planes=[]
	for bi in ballset.balls:
		for bj in ballset.intersection(bi):
			if bj.ID<=bi.ID: continue
			ci,cj,ri,rj=bi.center,bj.center,bi.radius,bj.radius
			n=(cj-ci).normalize()
			p=(ci+cj+n*(ri-rj))*0.5
			bi.planes+=[Plane4f(n*+1.0,p)]
			bj.planes+=[Plane4f(n*-1.0,p)]
			
	# _________________ STEP 2 create points for all spheres projecting inside ball.planes
	for ball in ballset.balls:
		ball.points=[]
		for i in range(len(POINTS)):
			point=(POINTS[i]*ball.radius+ball.center)
			for h in ball.planes:	
				if h.getDistance(point)>=0: 
					point=h.projectPoint(point)
				
			ball.points+=[point.x,point.y,point.z]
				
	# _________________ STEP 3 create mkpol with points
	hpc=Hpc()
	for ball in ballset.balls:
		vmat,hmat=Matf(3),Matf(3)
		npoints=len(ball.points)/3
		g=Graph.mkpol(vmat,hmat,3,npoints,ball.points,TOLERANCE)
		hpc.add(Hpc(g,vmat,hmat))
		
	return hpc
	
		
# //////////////////////////////////////////////////////////////////
# esempio di utilizzo
# //////////////////////////////////////////////////////////////////

if __name__ == "__main__":

	ball1=Ball3f(1.0,Vec3f([2,2,0]))
	ball2=Ball3f(1.5,Vec3f([4,2,0]))

	ballset=BallSet()
	ballset.addBall(ball1)
	ballset.addBall(ball2)
	ballset.addBall(None)
	
	#print "Balls intersecting",str(ball1),ballset.intersection(ball1)
	#print "Balls intersecting",str(ball2),ballset.intersection(ball2)
	
	#other=Ball3f(1,Vec3f([4,4,4]))
	#print "Balls intersecting",str(other),ballset.intersection(other)
	
	# view the hpc
	hpc=MkPolOfBallSet(ballset)
	Plasm.View(
            #SKELETON(1)
            (hpc))


