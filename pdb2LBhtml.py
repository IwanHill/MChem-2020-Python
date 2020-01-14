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
      <ul>
      <li>%s</li>
      <li>%s</li>
      <li>%s</li>
      </ul>
      </body>
    </html>"""
    page = html_wrapper % (input[0:10], input[10:20], input[20:30])
    with open(filename, "w") as output:
        output.write(page)

molecule = read_in(sys.argv[1])
write_out(sys.argv[2], molecule)
