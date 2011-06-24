# Questo file ha lo scopo di creare semplici molecole PDB per mezzo
# della sfera creata con la decomposizione ad esaedri.

from atomic_radius import atomic_radius
from graphics import atom_color, MOLVIEW
from hccsphere import hexsphere

class Atom:
    def __init__(self, name):
        self.name = name
        self.radius = atomic_radius[self.name]
        self.color = atom_color[self.name]
        self.graph = hexsphere()


    def view(self):
        MOLVIEW(self.graph, self.color)()


if __name__ == "__main__":
    h = Atom("O")
    print h.radius
    print h.color
    h.view()