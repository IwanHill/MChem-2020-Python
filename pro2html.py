from Bio import SeqIO
import sys

def read_in(filename = None):
    if filename == None:
        print("Please specify a file to read in!")
    for record in SeqIO.parse(filename, "pdb-seqres"):
        return record.seq

def write_out(filename = "html_out.html", input = None):
    html_wrapper = """
    <html lang="en" dir="ltr">
      <head>
        <meta charset="utf-8">
        <title>"2DREP"</title>
      </head>
      <body>
      %s
      </body>
    </html>"""
    page = html_wrapper % input
    with open(filename, "w") as output:
        output.write(page)

molecule = read_in(sys.argv[1])
write_out(sys.argv[2], molecule)
