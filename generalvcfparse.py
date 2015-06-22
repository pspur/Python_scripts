from __future__ import division
from __future__ import print_function
import re
import gzip
import itertools

variants = {} # key = chrom:pos, value = allele
sampleid = []

# Create list of all 50 sample ID #'s
with open('/gpfs/user/paulspur/vcfs.txt','r') as vcf_paths:
  for line in vcf_paths:
    m = re.search(r'Sample_(\d+)',line)
    sampleid.append(m.group(1))
sampleid.sort()

scores = open('/gpfs/user/paulspur/scores.txt','w')
#print('sID  ', end = '', file = scores)
#print('\n', file = scores)
 
for ref in sampleid:
  fin = gzip.open('/projects/big/BIG-Pilot-Projects/021115-Statin-Induced-Myopathy/'
                'Project_NOW_01348_WGS_2015-02-03/Sample_{0}/analysis/'
                '{0}.recalibrated.haplotypeCalls.vcf.gz'.format(ref),'r')
  for i,line in enumerate(fin):
#    if i % 500000 == 0:
#      print('building dict...')
    m = re.search(r'^(\d{1,2})\s+(\d+)\s+\S+\s+[ATCG]+\s+([,ATCG]+)',line)
    if m:
      vkey = m.group(1) + ':' + m.group(2)
      variants[vkey] = m.group(3)
  fin.close()
  
  for sid in sampleid:
    if sid < ref:
      print('{0} <> {1} = abc'.format(ref,sid),file = scores)
    if sid == ref:
#      print('{0} <> {1} -- similarity score: 1.000'.format(ref,sid))
      print('{0} <> {1} = 1.00'.format(ref,sid),file = scores) 
    if sid > ref:
#      print('Comparing {0} to {1}'.format(ref,sid))
      count = 0
      matches = 0
      fin = gzip.open('/projects/big/BIG-Pilot-Projects/021115-Statin-Induced-Myopathy/'
               'Project_NOW_01348_WGS_2015-02-03/Sample_{0}/analysis/'
               '{0}.recalibrated.haplotypeCalls.vcf.gz'.format(sid),'r')
      for line in (fin):
        m = re.search(r'(^\d{1,2})\s+(\d+)\s+\S+\s+[ATCG]+\s+([,ATCG]+)',line)
        if m:
          count += 1
 #         if count % 1000000 == 0:
 #           print('processing hits {0}++..'.format(count))
          vkey = m.group(1) + ':' + m.group(2)
          if vkey in variants:
            if variants[vkey] == m.group(3):
              matches += 1
      fin.close()
      uniquehits = len(variants) + (count - matches)
      sscore = matches / uniquehits
#      print('matches: {0}, unique hits: {1}'.format(matches,uniquehits))
#      print('{0} <> {1} -- similarity score: {2:.3f}'.format(ref,sid,sscore))
      print('{0} <> {1} = {2:.3f} '.format(ref,sid,sscore),file = scores)
scores.close()
