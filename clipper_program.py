import sys
from clipper_python import _clipper as clipper

def clipper_readin(filename): #from Hackers tutorial
    model_file = clipper.MMDBfile()
    model_file.read_file(filename) #supports .pdb and .cif by default
    molecule = clipper.MiniMol()
    model_file.import_minimol(molecule)
    return molecule

def fasta(molecule):
    model = molecule.model()
    sequence = "" #need something to add residues to
    #counter = 0
    residues = dict(ALA = "A", ARG = "R", HIS = "H", LYS = "K", ASP = "D",
    GLU = "E", SER = "S", THR = "T", ASN = "N", GLN = "Q", PRO = "P",
    VAL = "V", ILE = "I", LEU = "L", MET = "M", PHE = "F", TYR = "Y",
    TRP = "W", CYS = "C", GLY = "G", UNK = "X") #dictionary of natural AAs.
    for polymer in model:
        #sequence += "<br>"
        #counter += 1
        #sequence += "chain " + str(counter) + ": <br>"
        for monomer in polymer:
            if monomer.type().trim() in residues:
                sequence += residues[monomer.type().trim()]
            else:
                pass #cuts out *all* non standard amino acids, heteroatoms, worth looking at later,
    return sequence

def line_breaker(input): #slices and rejoins every 80 characters with 2 HTML line breaks
    string_record = str(input)
    string_record = '</div><br><div style="font-family:courier; font-size:2vw;">'.join(string_record[i:i + 80] for i in range(0, len(string_record), 80))
    return string_record

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
                    <div style="font-family:courier; font-size:2vw;">
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

def write_data(filename = "2DREP_Output.html", sequence = None):
    counter = 0 #allows us to integrate line_breaker functionality into our
    with open(filename, "a") as output:
        for residue in sequence:
            if counter == 0:
                output.write("""
                S A M P L E L E T T E R S ~ALPHA1~ <img src="sample_helix.png" alt="sample_helix" style="width:10%;heigth:10%">
                _BETA1_ <img src="sample_sheet.png" alt="sample_sheet" style="width:10%;heigth:10%"> </div><div style="font-family:courier;
                font-size:2vw;"> """)
            #call some function that will write the  annotations we want! Probs need to check in with Jon for this...
            counter += 1
            colour = get_colour(residue)
            if colour:
                output.write("<span style=\"color:" + colour + "\">" + residue + "</span>") #colours with the appropriate hexcode for our residue
            else:
                output.write(residue)
            if counter == 80:
                output.write("</div><div style=\"font-family:courier; font-size:2vw;\">")
                counter = 0

molecule = clipper_readin(sys.argv[1]) #reads in file from command line
sequence = fasta(molecule)
write_boilerplate_start(sys.argv[2])
write_data(sys.argv[2], sequence)
write_boilerplate_end(sys.argv[2])
