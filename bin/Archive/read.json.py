
import argparse
import sys
import json

def main(): 

	global args
	global json_file
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-json_file', metavar='<json_file>', dest="json_file", help="JSON file output by PFAM database")
	args = parser.parse_args()
	
	## Prep global objects
	json_file = str(args.json_file)
	
	#Read JSON file
	with open(json_file, 'r') as json_file: 
		data = json.load(json_file)[0]
	#Print
	print("Model length ----------------------")
	print(data['model_length'])
	
	print("Alignment ----------------------")
	print(data['align'])
	
	print("Env ----------------------")
	print(data['env'])
	
	print("Name ----------------------")
	print(data['name'])
	
	print("Accession ----------------------")
	print(data['acc'])
	
	print("Significative ----------------------")
	print(data['sig'])
	
	print("Expectation value (E-value) ----------------------")
	print(data['evalue'])
	
	print("Description ----------------------")
	print(data['desc'])
	
	print("Env ----------------------")
	print(data['env'])
	
	print("HMM ----------------------")
	print(data['hmm'])
	
	print("Active site ----------------------")
	print(data['act_site'])
	
	print("Type ----------------------")
	print(data['type'])
	
	print("Bits ----------------------")
	print(data['bits'])
	
	print("Clan ----------------------")
	print(data['clan'])
	
	print("Sequence ----------------------")
	print(data['seq'])

if __name__ == "__main__":
      main()
