import argparse
import requests
import sys

def main():


	global args
	global transcript_id
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-transcript_id', metavar='<transcript_id>', dest="transcript_id", help="Transcript ID for CDS retrieval by ENSEMBL REST API")
	parser.add_argument('-seq_type', metavar='<seq_type>', dest="seq_type", help="Sequence type to retrieve from Ensembl")
	
	args = parser.parse_args()
	
	## Prep global objects
	transcript_id = str(args.transcript_id)
	seq_type = str(args.seq_type)

	server = "https://rest.ensembl.org"
	ext = "/sequence/id/"
	 
	r = requests.get(server+ext+transcript_id+"?type="+seq_type, headers={ "Content-Type" : "text/x-fasta"})
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
	 
	 
	print(r.text)

if __name__ == "__main__":
    main()
