#!/usr/bin/python

from __future__ import print_function
import sys
import yaml
import operator
import datetime
from titlecase import titlecase

include_fields = ['AB', 'AID', 'AU', 'DEP', 'DP', 'FAU', 'JT', 'LID', 'PMID', 'TA', 'TI']
# include_fields = ['DEP', 'DP', 'TI']
pubs = {}

field = ''
contents = ''
pmid = ''
with open(sys.argv[1],'r') as f:
  for line in f:
    if not line.strip():
      continue
    
    if line[0:4].strip():
      # New field
      # Store last complete field
      if pmid:        
        if field in pubs[pmid]:
          pubs[pmid][field].append(contents)
        else:
          pubs[pmid][field] = [contents]
        
        # Add consortia to authors
        if field == 'CN':
          if 'AU' in pubs[pmid]:
            pubs[pmid]['AU'].append(contents)
          else:
            pubs[pmid]['AU'] = [contents]
          if 'FAU' in pubs[pmid]:
            pubs[pmid]['FAU'].append(contents)
          else:
            pubs[pmid]['FAU'] = [contents]
          
      # Process line
      field = line[0:4].strip()
      contents = line[6:].strip()
      
      if field == 'PMID':
        # New publication
        pmid = contents
        pubs[pmid] = {}
      
    else:
      # Continued field
      contents = contents+' '+line[6:].strip()

# Store last field
if field in pubs[pmid]:
  pubs[pmid][field].append(contents)
else:
  pubs[pmid][field] = [contents]

pubs_list = []
for pmid, pub in pubs.items():
  pub_reformat = {}
  
  pub_reformat['title'] = pub['TI'][0]
  
  pub_reformat['authors'] = pub['AU']
  pub_reformat['authors_full'] = pub['FAU']
  
  pub_reformat['abstract'] = pub['AB'][0]
  pub_reformat['journal'] = pub['TA'][0]
  pub_reformat['journal_full'] = titlecase(pub['JT'][0])
  pub_reformat['pmid'] = pub['PMID'][0]
  
  try:
    date = datetime.datetime.strptime(pub['DP'][0], '%Y %b %d')
  except ValueError:
    try:
      date = datetime.datetime.strptime(pub['DP'][0], '%Y %b')
    except ValueError:
      date = datetime.datetime.strptime(pub['DEP'][0], '%Y%m%d')
  pub_reformat['date'] = date
  pub_reformat['year'] = date.year
  
  ids = []
  if 'AID' in pub:
    ids = ids+pub['AID']
  if 'LID' in pub:
    ids = ids+pub['LID']
  doi_ids = [x for x in ids if ' [doi]' in x]
  doi = doi_ids[0].replace(' [doi]','')
  pub_reformat['doi'] = doi
  
#   AID LID
  
  pubs_list.append(pub_reformat)


# Sort list of publications by date
pubs_list.sort(key=lambda p: p['date'], reverse=True)
# l_pubs = [x[1] for x in sorted(pubs.items(), key=lambda p: p[1]['date'], reverse=True)]
print(yaml.dump(pubs_list, default_flow_style=False))


# print(yaml.dump(pubs, default_flow_style=False))