#!/usr/bin/env python

# based on https://www.biostars.org/p/77526/

import requests
import csv

SRA_id = 'SRR11092056'
url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + SRA_id

with requests.Session() as request:
    r = requests.post(url)
    csv = csv.DictReader(r.content.decode('utf-8').splitlines(), delimiter=',')
    for row in csv: print row['Platform']

