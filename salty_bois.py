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

def find_charges(molecule):
    positives = []
    negatives = []
    model = molecule.model()
    for polymer in model:
        for monomer in polymer:
        	if monomer.type().trim() == "ARG" or monomer.type().trim() == "LYS":
        		positives.append (monomer.id().trim())
        	elif monomer.type().trim() == "ASP" or monomer.type().trim() == "GLU":
        		negatives.append(monomer.id().trim())
    print "Positives:"
    print positives
    print "Negatives:"
    print negatives

def salt_bridges (molecule=None):
    model = molecule.model()
    metrics = MetricsModel(molecule)
    index = 0
    complete_sides = []
    list_of_dicts = []#set up lists to add outputs to
    x = 0

    for chain in metrics.chains:
        for residue in chain:
            if residue.is_sidechain_complete == True:
                complete_sides.append(True)
            else:
                complete_sides.append(False)
#many files are missing side chains for some residues
#these cause errors if not excepeted!

    for polymer in model:
        print polymer.id()
        positives = []
        negatives = []
        bridge_pairs = {}
        for monomer in polymer:
            if complete_sides[index]: #checks if current mresiude is complete
                if monomer.type().trim() in ("ARG", "LYS"):
                    positives.append(monomer)
                elif monomer.type().trim() in ("ASP", "GLU"):
        		    negatives.append(monomer)
            index +=1 #increments index for next residue

        for positive in positives:
            if positive.type().trim() == "ARG":
                plus_coord = positive.find(clipper.String("CZ")).coord_orth
            elif positive.type().trim() == "LYS":
                plus_coord = positive.find(clipper.String("NZ")).coord_orth
            closest_range = 999

    	    for negative in negatives:
                if negative.type().trim() == "ASP":
                    neg_coord = negative.find(clipper.String("CG")).coord_orth
                elif negative.type().trim() == "GLU":
                    neg_coord = negative.find(clipper.String("CD")).coord_orth

                if clipper.Coord_orth.length(plus_coord, neg_coord) < closest_range and clipper.Coord_orth.length(plus_coord, neg_coord) < 4.6:
                    closest_range = clipper.Coord_orth.length(plus_coord, neg_coord)
                    paired_aa = negative
            if closest_range != 999:
                x += 1
                bridge_pairs[x] = (positive.id().trim(), positive.type().trim(), closest_range, paired_aa.id().trim(), paired_aa.type().trim())

        for negative in negatives:
            rep_counter = 0
            distances = [0]
            for key, value in bridge_pairs.items():
                if negative.id().trim() in value:
                    rep_counter += 1
                    distances.append(value[2])
                if rep_counter == 2:
                    for key, value in bridge_pairs.items():
                        if max(distances) in value:
                            del bridge_pairs[key]
                if rep_counter > 2:
                    for key, value in bridge_pairs.items():
                        if negative.id().trim() in value:
                            del bridge_pairs[key]

        list_of_dicts.append(bridge_pairs)

    return list_of_dicts


molecule = clipper_readin(sys.argv[1])
salt_bridges(molecule)
