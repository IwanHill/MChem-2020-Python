import sys
import math
from clipper_python import _clipper as clipper

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

def fasta(polymer):
    sequence = "" #need something to add residues to
    for monomer in polymer:
        if monomer.type().trim() in residues:
            sequence += residues[monomer.type().trim()]
        else:
            pass #cuts out *all* non standard amino acids, heteroatoms, worth looking at later
    return sequence

def get_colour(residue):
    colours = dict(D = "#D71313", E = "#D71313", R = "#0027FF",
    H = "#0027FF", K = "#0027FF", S = "#15E4EA", T = "#15E4EA",
    N = "#15E4EA", Q = "#15E4EA", C = "#FF00AB", U = "#FF00AB",
    G = "#FF00AB", P = "#FF00AB", A = "#C6BE2A", V = "#C6BE2A",
    I = "#C6BE2A", L = "#C6BE2A", M = "#C6BE2A", F = "#C6BE2A",
    Y = "#C6BE2A", W = "#C6BE2A") # in order:
    #negative chains in red #positive side chains in dark blue
    #polar side chains in sky blue #special side chains in hot pink
    #hydrophobic side chains in oily yellow #won't style unknown or unnatrual AAs
    if residue in colours:
        return colours[residue]
    else:
        return "#000000"

def calc_avg_B(monomer):
    atom_Bs = []
    for atom in monomer:
        u = atom.u_iso #draws u value out - don't use brackets!
        bfac = u * 8 * (math.pi)**2 #convert U factor to B factor for atom.
        atom_Bs.append(bfac)
    avg_B = round(sum(atom_Bs)/len(atom_Bs), 1) #brute force calulation of mean B factor
    return avg_B #need to call function within loops.

def background_B_colour(list_of_Bs):
    bkgrnd_colours = []
    for avg_B in list_of_Bs:
        if avg_B <= (0.05 * max(list_of_Bs)):
            bkgrnd_colours.append("#FFFFFF")
        elif avg_B >= (0.05 * max(list_of_Bs)) and avg_B <= (0.15 * max(list_of_Bs)):
            bkgrnd_colours.append("#E6E6E6")
        elif avg_B >= (0.15 * max(list_of_Bs)) and avg_B <= (0.25 * max(list_of_Bs)):
            bkgrnd_colours.append("#CCCCCC")
        elif avg_B >= (0.25 * max(list_of_Bs)) and avg_B <= (0.35 * max(list_of_Bs)):
            bkgrnd_colours.append("#B3B3B3")
        elif avg_B >= (0.35 * max(list_of_Bs)) and avg_B <= (0.45 * max(list_of_Bs)):
            bkgrnd_colours.append("#999999")
        elif avg_B >= (0.45 * max(list_of_Bs)) and avg_B <= (0.55 * max(list_of_Bs)):
            bkgrnd_colours.append("#808080")
        elif avg_B >= (0.55 * max(list_of_Bs)) and avg_B <= (0.65 * max(list_of_Bs)):
            bkgrnd_colours.append("#666666")
        elif avg_B >= (0.65 * max(list_of_Bs)) and avg_B <= (0.75 * max(list_of_Bs)):
            bkgrnd_colours.append("#4D4D4D")
        elif avg_B >= (0.75 * max(list_of_Bs)) and avg_B <= (0.85 * max(list_of_Bs)):
            bkgrnd_colours.append("#333333")
        elif avg_B >= (0.85 * max(list_of_Bs)) and avg_B <= (0.95 * max(list_of_Bs)):
            bkgrnd_colours.append("#1A1A1A")
        elif avg_B >= (0.95 * max(list_of_Bs)):
            bkgrnd_colours.append("#000000")
    return(bkgrnd_colours)

def calc_sec():
    pass

def calc_tert():
    pass

def write_boilerplate_start(filename = "2DREP_Output.html"):
    with open(filename, "a") as output: #use "a" to append, prevents overwriting.
        #below: HTML boilerplate with inline CSS
        output.write("""
            <html lang="en" dir="ltr">
            <head>
                <meta charset="utf-8" name="viewport" content="width=device-width, initial-scale=1.0">
                <title>"2DREP"</title>
                </head>
                <body>
                <img src="logo_and_legend.png" alt="2DREP Logo" style="width:100%;height:40%;">
                """
                )

def write_boilerplate_end(filename = "2DREP_Output.html"):
    with open(filename, "a") as output: #"a" to append to end of sequence output.
        output.write("""
                </div>
            </p>
        </body>
    </html>""") #close HTML

def write_data(filename = "2DREP_Output.html", molecule = None):
    model = molecule.model()
    with open(filename, "a") as output:
        for polymer in model:
            character_number = 78
            counter = 0 #allows us to integrate line_breaker functionality into our writer
            residue_count = 0 #tracked separately, won't reset with each new line.
            output.write("""<div style="font-family:courier; font-size:2vw; text-decoration:underline;"><br> Chain """)
            output.write(str(polymer.id().trim()))
            output.write(": </div>")
            bfacs_list = [] #we'll pull a B factor for display later
            for monomer in polymer:
                if monomer.type().trim() in residues: #need to act before our residues are a string
                    bfacs_list.append(calc_avg_B(monomer)) #adds the b factor for each residue in the ATOM list - will match our sequence!
            background_shades = background_B_colour(bfacs_list)
            sequence = fasta(polymer) #we'll write out this string styled
            for residue in sequence:
                residue_count += 1
                if counter == 0:
                    secondary = calc_sec()
                    output.write("<div style=\"font-family:courier;font-size:2vw;\">&nbsp;&nbsp;&nbsp;&nbsp;" + str(secondary))
                    output.write("</div><div style=\"font-family:courier;font-size:2vw;\">")
                    #call some function that will write the  annotations we want! Probs need to check in with Jon for this...
                    output.write(str(residue_count))
                    if residue_count < 100:
                        output.write("&nbsp;")
                    if residue_count < 10:
                        output.write("&nbsp;")
                    output.write("&nbsp")
                colour = get_colour(residue)
                background = background_shades[residue_count - 1] #residue count starts at 1, indices start at 0.
                output.write("<span style=\"color:" + colour + "; background-color:" + background + "\">" + residue + "</span>") #colours with the appropriate hexcode for our residue
                counter += 1
                if len(str(residue_count)) != 1 and len(str(residue_count)) != 2:
                    character_number = 79 #ends up with 74 residues to a line, so new lines are on multiples of 75, for easier counting.
                if counter == (character_number - len(str(residue_count))):
                    output.write("</div>")
                    counter = 0

molecule = clipper_readin(sys.argv[1]) #reads in file from command line
write_boilerplate_start(sys.argv[2])
write_data(sys.argv[2], molecule)
write_boilerplate_end(sys.argv[2])
