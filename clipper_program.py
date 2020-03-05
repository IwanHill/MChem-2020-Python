import sys
from clipper_python import _clipper as clipper

def clipper_readin(filename): #from Hackers tutorial
    model_file = clipper.MMDBfile()
    model_file.read_file(filename) #supports .pdb and .cif by default
    molecule = clipper.MiniMol()
    model_file.import_minimol(molecule)
    return molecule

def fasta(polymer):
    sequence = "" #need something to add residues to.
    residues = dict(ALA = "A", ARG = "R", HIS = "H", LYS = "K", ASP = "D",
    GLU = "E", SER = "S", THR = "T", ASN = "N", GLN = "Q", PRO = "P",
    VAL = "V", ILE = "I", LEU = "L", MET = "M", PHE = "F", TYR = "Y",
    TRP = "W", CYS = "C", GLY = "G", UNK = "X") #dictionary of natural AAs.
    for monomer in polymer:
        if monomer.type().trim() in residues:
            sequence += residues[monomer.type().trim()]
        else:
            pass #cuts out *all* non standard amino acids, heteroatoms, worth looking at later,
    return sequence

def get_colour(residue):
    if residue in ("D", "E"):
        return("#D71313") #negative chains in red
    elif residue in ("R", "H", "K"):
        return("#0027FF") #positive side chains in dark blue
    elif residue in ("S", "T", "N", "Q"):
        return("#15E4EA") #polar side chains in sky blue
    elif residue in ("C", "U", "G", "P"):
        return("#FF00AB") #special side chains in hot pink
    elif residue in ("A", "V", "I", "L", "M", "F", "Y", "W"):
        return("#C6BE2A") #hydrophobic side chains in oily yellow
    else:
        return False #won't throw style at any unnatural amino acid or formatting characters.

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
                <img src="2DREP_logo.png" alt="2DREP Logo" style="width:40%;height:25%;">
                """
                ) #Working title is 2DREP
                #logo styling is not working inline but other bits are...

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
            residue_count = 0
            output.write("""<div style="font-family:courier; font-size:2vw; text-decoration:underline;"><br> Chain """)
            output.write(str(polymer.id().trim()))
            output.write(": </div>")
            sequence = fasta(polymer)
            for residue in sequence:
                residue_count += 1
                if counter == 0:
                    output.write("""
                    <div style="font-family:courier; font-size:2vw;">
                    &nbsp;&nbsp;&nbsp;&nbsp;S A M P L E L E T T E R S ~ALPHA1~ <img src="sample_helix.png" alt="sample_helix" style="width:10%;heigth:10%">
                    _BETA1_ <img src="sample_sheet.png" alt="sample_sheet" style="width:10%;heigth:10%"> </div><div style="font-family:courier;
                    font-size:2vw;"> """)
                    #call some function that will write the  annotations we want! Probs need to check in with Jon for this...
                    output.write(str(residue_count))
                    if residue_count < 100:
                        output.write("&nbsp;")
                    if residue_count < 10:
                        output.write("&nbsp;")
                    output.write("&nbsp")
                colour = get_colour(residue)
                output.write("<span style=\"color:" + colour + "\">" + residue + "</span>") #colours with the appropriate hexcode for our residue
                counter += 1
                if len(str(residue_count)) != 1 and len(str(residue_count)) != 2:
                    character_number = 79
                if counter == (character_number - len(str(residue_count))):
                    output.write("</div>")
                    counter = 0

molecule = clipper_readin(sys.argv[1]) #reads in file from command line
write_boilerplate_start(sys.argv[2])
write_data(sys.argv[2], molecule)
write_boilerplate_end(sys.argv[2])
