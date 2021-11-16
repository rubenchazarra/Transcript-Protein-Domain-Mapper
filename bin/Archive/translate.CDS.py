

## TODO in this script
# 1) Devide code in functions - At leats one for translateing the sequence

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' #This is related to an OPENBLAS error we where facing

import argparse
import sys
from Bio import SeqIO

def main(): 

	global args
	global fasta
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-fasta', metavar='<fasta>', dest="fasta", help="Fasta file with CDS with a single record")
	args = parser.parse_args()
	
	## Prep global objects
	fasta = str(args.fasta)
	
	## Read FASTA file
	record = SeqIO.read(fasta, "fasta")
	# get fasta sequence	
	sequence = record._seq
	# translate sequence
	prot_seq = sequence.translate(to_stop=True)# this argument prevents the introdution of the STOP codon symbol '*'
	print(prot_seq)

if __name__ == "__main__":
   try:
      main()
   except KeyboardInterrupt:
      # do nothing here
      main()
