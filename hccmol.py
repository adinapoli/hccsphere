# Questo file ha lo scopo di creare semplici molecole PDB per mezzo
# della sfera creata con la decomposizione ad esaedri.

from graphics import *
from hccsphere import *


class Atom:
    def __init__(self, name, coords):
        self.name = name.lower().capitalize()
        self.radius = atomic_radius[self.name[0]]
        self.scale = self.radius[3]/200.0
        self.color = atom_color[self.name[0]]
        self.graph = hexsphere(1)
        self.coords = coords


    def view(self):
        MOLVIEW(hccbatches(self))


class Molecule:
    def __init__(self, filename):
        self.filename = filename
        self.atoms = []
        self.batches = []


    def build(self):
        graph = get_trigraph(self.filename)
        self.atoms = [Atom(i[0], i[1]) for i in izip(graph[2], graph[0])]

        for a in self.atoms:
            self.batches.extend(hccbatches(a))

    def view(self):
        MOLVIEW(self.batches)


if __name__ == "__main__":

    glu = Molecule("GLU.pdb")
    glu.build()
    glu.view()
