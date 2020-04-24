import sys
import math
from clipper_python import _clipper as clipper
from clipper_tools.metrics import MetricsModel

#set a dictionary of residues as a global variable
#allows multiple functions to draw from same dictionary
residues = dict(ALA = "A", ARG = "R", HIS = "H", LYS = "K", ASP = "D",
GLU = "E", SER = "S", THR = "T", ASN = "N", GLN = "Q", PRO = "P",
VAL = "V", ILE = "I", LEU = "L", MET = "M", PHE = "F", TYR = "Y",
TRP = "W", CYS = "C", GLY = "G", UNK = "X") #dictionary of natural AAs.

def clipper_readin(filename): #from Hackers tutorial
    model_file = clipper.MMDBfile()
    model_file.read_file(filename) #supports .pdb and .cif by default
    molecule = clipper.MiniMol()
    model_file.import_minimol(molecule)
    return molecule

def find_cys(molecule):
    model = molecule.model()
    cysteines = []
    for polymer in model:
        for monomer in polymer:
            if monomer.type().trim() == "CYS":
                cysteines.append(monomer.id().trim())
    return cysteines

def ss_bridge(molecule=None):
    model = molecule.model()
    cysteines = []
    bridge_pairs = {} #set up lists to add outputs to
    i = 0
    for polymer in model:
        for monomer in polymer:
            if monomer.type().trim() == "CYS":
                cysteines.append(monomer) #assemble a list of all CYS residues
    for cys in cysteines:
        cys_1_S = cys.find(clipper.String("SG")).coord_orth #set location of first cys
        for next in cysteines:
            #test first cys against all subsequent, looking for closeness (hence a bridge)
            cys_2_S = next.find(clipper.String("SG")).coord_orth
            if clipper.Coord_orth.length(cys_1_S, cys_2_S) < 2.6 and clipper.Coord_orth.length(cys_1_S, cys_2_S) != 0:
                i += 1#gives a numbered key for use as our bridge annotation
                bridge_pairs[i] = (cys.id().trim(), next.id().trim())
                cysteines.remove(cys) #prevents repetition

    return bridge_pairs




molecule = clipper_readin(sys.argv[1])
print find_cys(molecule)
print ss_bridge(molecule)
