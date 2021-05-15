import requests, sys
 
server = "https://rest.ensembl.org"
ext = "/sequence/id/ENST00000288602?type=cds"
 
r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
 
print(r.text)
