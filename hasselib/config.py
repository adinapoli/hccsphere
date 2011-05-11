from explode import *
from os import getcwd, path


def myprint(arg0,arg1): print "\n"+arg0+" =", arg1


objs_path = path.dirname(__file__)

#Oggetto caricato da disco che modella una sphera batch e che viene usato passato ogni
#volta come parametro nel costruttore Batch(), impiegato nel generatore get_sphere().
#L'idea e' avere un unico batch caricato in RAM e prendere copie identiche clonate da esso.
SPHERE_BATCH = Batch.openObj("".join([objs_path, "/objs/sphere18x27.obj"]))[0]
CYLINDER_BATCH = Batch.openObj("".join([objs_path, "/objs/cylinder4x27.obj"]))[0]


COMPOUNDS_DIR = "".join([getcwd(), "/compounds/"])