from __future__ import division
from __future__ import print_function
import re
import gzip
import itertools

variants = {} # key = chrom:pos, value = allele
sampleid = []
scores = open('./scores.txt','w')

# Create list of all 50 sample ID #'s
with open('./vcfs.txt','r') as vcf_paths:
  for line in vcf_paths:
    m = re.search(r'Sample_(\d+)',line)
    sampleid.append(m.group(1))
sampleid.sort()

# Iterate through the sample files. For each iteration, create dict of
# current file then compare to each remaining file. Avoid repeating
# comparisons and comparing one to itself.
for ref in sampleid:
  fin = gzip.open('/projects/big/BIG-Pilot-Projects/021115-Statin-Induced-Myopathy/'
                'Project_NOW_01348_WGS_2015-02-03/Sample_{0}/analysis/'
                '{0}.recalibrated.haplotypeCalls.vcf.gz'.format(ref),'r')
  for i,line in enumerate(fin):
    m = re.search(r'^(\d{1,2})\s+(\d+)\s+\S+\s+[ATCG]+\s+([,ATCG]+)',line)
    if m:
      vkey = m.group(1) + ':' + m.group(2)
      variants[vkey] = m.group(3)
  fin.close()
  
  for sid in sampleid:
    if sid < ref: # abc = already been compared
      print('{0} <> {1} = abc'.format(ref,sid),file = scores)
    if sid == ref: # skip compare, similarity to self is 1.0
      print('{0} <> {1} = 1.00'.format(ref,sid),file = scores) 
    if sid > ref:
      count = 0
      matches = 0
      fin = gzip.open('/projects/big/BIG-Pilot-Projects/021115-Statin-Induced-Myopathy/'
               'Project_NOW_01348_WGS_2015-02-03/Sample_{0}/analysis/'
               '{0}.recalibrated.haplotypeCalls.vcf.gz'.format(sid),'r')
      for line in (fin):
        m = re.search(r'(^\d{1,2})\s+(\d+)\s+\S+\s+[ATCG]+\s+([,ATCG]+)',line)
        if m:
          count += 1
          vkey = m.group(1) + ':' + m.group(2)
          if vkey in variants and variants[vkey] == m.group(3)::
            matches += 1
      fin.close()
      uniquehits = len(variants) + count - matches
      sscore = matches / uniquehits
      print('{0} <> {1} = {2:.3f} '.format(ref,sid,sscore),file = scores)
scores.close()
