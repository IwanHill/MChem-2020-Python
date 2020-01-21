from Bio import SeqIO
import sys
from get_colour import get_colour

def read_in(filename):
    for record in SeqIO.parse(filename, "pdb-seqres"):
        return record.seq

def line_breaker(input):
    string_record = str(input)
    string_record = '<br><br>'.join(string_record[i:i + 80] for i in range(0, len(string_record), 80))
    return string_record

def write_out(filename = "html_out.html", sequence = None):
    with open(filename, "w") as output:
        output.write("""
            <html lang="en" dir="ltr">
            <head>
                <meta charset="utf-8">
                <title>"2DREP"</title>
                </head>
                <body>
                    <p style="font-family:courier; font-size:2vw;">
                """
                )
        for residue in sequence:
            colour = get_colour(residue)
            if colour:
                output.write("<span style=\"color:" + colour + "\">" + residue +"</span>")
            else:
                output.write(residue)
        output.write("""
            </p>
        </body>
    </html>""")


sequence = read_in(sys.argv[1])
molecule = line_breaker(sequence)
write_out(sys.argv[2], molecule)
