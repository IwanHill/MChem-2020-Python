import sys
from clipper_python import _clipper as clipper

def clipper_readin(filename): #from Hackers tutorial
    model_file = clipper.MMDBfile()
    model_file.read_file(filename) #supports .pdb and .cif by default
    molecule = clipper.MiniMol()
    model_file.import_minimol(molecule)
    return molecule

def three_letter_rep(monomer):
    for polymer in model:
        for monomer in polymer:
            if monomer.type().trim() in ("ALA", "ARG", "ASN", "ASP", "CYS",
            "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
            "PRO", "SER", "THR", "TRP", "TYR", "VAL"):
#need to isolate amino acids from monosaccahrides, water, ligands etc
                print monomer.type().trim()

def fasta(molecule):
    model = molecule.model()
    for polymer in model:
        sequence = "" #starts empty
        for monomer in polymer: #fasta covnersion by brute force
        #do we wanna use a dictonary here?
            if monomer.type().trim() == "ALA":
                sequence += "A"
            elif monomer.type().trim() == "ARG":
                sequence += "R"
            elif monomer.type().trim() == "HIS":
                sequence += "H"
            elif monomer.type().trim() == "LYS":
                sequence += "K"
            elif monomer.type().trim() == "ASP":
                sequence += "D"
            elif monomer.type().trim() == "GLU":
                sequence += "E"
            elif monomer.type().trim() == "SER":
                sequence += "S"
            elif monomer.type().trim() == "THR":
                sequence += "T"
            elif monomer.type().trim() == "ASP":
                sequence += "N"
            elif monomer.type().trim() == "GLU":
                sequence += "Q"
            elif monomer.type().trim() == "CYS":
                sequence += "C"
            elif monomer.type().trim() == "GLY":
                sequence += "G"
            elif monomer.type().trim() == "PRO":
                sequence += "P"
            elif monomer.type().trim() == "VAL":
                sequence += "V"
            elif monomer.type().trim() == "ILE":
                sequence += "I"
            elif monomer.type().trim() == "LEU":
                sequence += "L"
            elif monomer.type().trim() == "MET":
                sequence += "M"
            elif monomer.type().trim() == "PHE":
                sequence += "F"
            elif monomer.type().trim() == "TYR":
                sequence += "Y"
            elif monomer.type().trim() == "TRP":
                sequence += "W"
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
        return False #won't throw style at any unnatural amino acid or

def write_out(filename = "2DREP_Output.html", sequence = None):
    with open(filename, "w") as output: #below: HTML boilerplate with inline CSS
        output.write("""
            <html lang="en" dir="ltr">
            <head>
                <meta charset="utf-8" name="viewport" content="width=device-width, initial-scale=1.0">
                <style>
                img {
                width: 50%;
                height: 25%;
                }
                </style>
                <title>"2DREP"</title>
                </head>
                <body>
                <img src="2DREP_logo.png" alt="2DREP Logo">
                    <div style="font-family:courier; font-size:2vw;">
                """
                ) #Working title is 2DREP
        for residue in sequence: #sequence is seen as a string here for printing.
            colour = get_colour(residue) #chooses the hexcode for the appropriate class of residue
            if colour: #non-residue letters such as line break characters return False instead of a hexcode in get_colour
                output.write("<span style=\"color:" + colour + "\">" + residue +"</span>")
                #writes a span tag for each character with its appropriate hexcode colour
            else:
                output.write(residue)
                #writes your unknown residue or line break charater without style
        output.write("""
                </div>
            </p>
        </body>
    </html>""") #close HTML


molecule = clipper_readin(sys.argv[1]) #reads in file from command line
sequence = fasta(molecule)
br_sequence = line_breaker(sequence)
write_out(sys.argv[2], br_sequence)
