# Questo file ha lo scopo di creare semplici molecole PDB per mezzo
# della sfera creata con la decomposizione ad esaedri.

from graphics import *
from hccsphere import *
import pstats
import cProfile
from math import log


class Atom:
    def __init__(self, name, coords, radius = 3):
        self.name = name.lower().capitalize()
        self.radius = atomic_radius[self.name[0]]
        self.scale = self.radius[radius]/300.0
        self.color = atom_color[self.name[0]]
        self.graph = None
        self.coords = coords
        self.batches = None
        self.volume = None

        
    def build(self):
        self.graph = hexsphere()
        self.batches = smart_hccbatches(self)

        
    def clone(self, graph):
        self.graph = Graph(graph)
        self.batches = smart_hccbatches(self)

        
    def view(self):
        MOLVIEW(self.batches)

        
    def volview(self):
        if not self.volume:
            self.volume = sectionize(self)

        print self.volume
        MOLVIEW(self.volume)

    def sectionize(self):
        self.volume = sectionize(self)


class Molecule:
    def __init__(self, filename):
        self.filename = filename
        self.atoms = []
        self.batches = []
        self.volumes = []


    def build(self):
        """Il trucco e' copiare il grafo invece di
        computarne uno nuovo. All'inizio, tutte le palle
        sono uguali!"""
        graph = get_trigraph(self.filename)
        it = izip(graph[2], graph[0])
        first_pair = it.next()
        first_atom = Atom(first_pair[0], first_pair[1])
        first_atom.build()

        self.atoms = [Atom(i[0], i[1]) for i in it]
        for a in self.atoms:
            a.clone(first_atom.graph)
        self.atoms.insert(0, first_atom)

        add = self.batches.extend
        for a in self.atoms:
            add(a.batches)

    def sectionize(self):
        for a in self.atoms:
            a.sectionize()
            self.volumes.extend(a.volume)

    def view(self):
        MOLVIEW(self.batches)


    def volview(self):
        if not self.volumes:
            for a in self.atoms:
                a.sectionize()
                self.volumes.extend(a.volume)

        MOLVIEW(self.volumes)


if __name__ == "__main__":

    start = time.clock()
    glu = Molecule("3NK2.pdb")
    glu.build()
    end = time.clock()

    print "Molecule builded in ", end - start, " seconds."
    glu.view()
    glu.volview()
