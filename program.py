import sys
import program_functions as pf

sequence = pf.read_in(sys.argv[1])
molecule = pf.line_breaker(sequence)
pf.write_out(sys.argv[2], molecule)
