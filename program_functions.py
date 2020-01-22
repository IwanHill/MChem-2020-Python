from Bio import SeqIO #needed for read_in

def read_in(filename):
    if ".pdb" in filename:
        mode = "pdb-seqres" #full pdb structure
    elif ".cif" in filename:
        mode = "cif-seqres" #full mmcif structure
    else:
        raise ValueError("Please give a .pdb or .cif file!")
        #aborts program if pdb or mmcif not given
    for record in SeqIO.parse(filename, mode):
        return record.seq #returns a Sequence object (iterable)

def line_breaker(input): #slices and rejoins every 80 characters with 2 HTML line breaks
    string_record = str(input)
    string_record = '<br><br>'.join(string_record[i:i + 80] for i in range(0, len(string_record), 80))
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
    elif residue in ("<", "b", "r", ">"):
        return False #won't throw style at the line breaks!
    else:
        return("#000000") #unknown residues will stay black

def write_out(filename = "2DREP_Output.html", sequence = None):
    with open(filename, "w") as output: #below: HTML boilerplate with inline CSS
        output.write("""
            <html lang="en" dir="ltr">
            <head>
                <meta charset="utf-8">
                <title>"2DREP"</title>
                </head>
                <body>
                    <p style="font-family:courier; font-size:2vw;">
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
            </p>
        </body>
    </html>""") #close HTML
