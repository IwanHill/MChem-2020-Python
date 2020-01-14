from Bio import SeqIO

for record in SeqIO.parse("3qvp.pdb", "pdb-seqres"):
    print(record.seq)
