#!/usr/bin/python

from __future__ import print_function
import sys
import yaml
import operator
import datetime
from titlecase import titlecase

def readCitations(fn, type):
  pubs_list = []
  pub = {}
  field = ''
  contents = ''
  with open(fn,'r') as f:
    for line in f:
      if not line.strip():
        # Add last publication
        if pub:
          pubs_list.append(pub)
        pub = {}
        continue
      
      if line[0:4].strip():
        # Process line
        field = line[0:4].strip()
        contents = line[6:].strip()
        if field in pub:
          pub[field].append(contents)
        else:
          pub[field] = [contents]
        
        # Add consortia to authors
        if field == 'CN':
          if 'AU' in pub:
            pub['AU'].append(contents)
          else:
            pub['AU'] = [contents]
          if 'FAU' in pub:
            pub['FAU'].append(contents)
          else:
            pub['FAU'] = [contents]
         
      else:
        # Continued field
        pub[field][-1] = pub[field][-1]+' '+line[6:].strip()

  pubs_list_reformatted = []
  for pub in pubs_list:
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
  
    pubs_list_reformatted.append(pub_reformat)


  # Sort list of publications by date
  pubs_list_reformatted.sort(key=lambda p: p['date'], reverse=True)
  # l_pubs = [x[1] for x in sorted(pubs.items(), key=lambda p: p[1]['date'], reverse=True)]
  print(yaml.dump(pubs_list_reformatted, default_flow_style=False))


if __name__ == "__main__":
  readCitations(sys.argv[1], 'nbib')