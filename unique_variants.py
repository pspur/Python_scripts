'''
When executing this script, follow this format:
python unique_variants.py [pipeline #] [sample #]
Entering a pipeline other than '1' or '2' fails gracefully, entering an
invalid sample # splats.

Overall flow of code: 
1. Open vcf file from one pipeline. For each variant call, if it passes logic
   test store in dict; key = chrom:pos, value = alt.
2. Open equivalent vcf file from second pipeline. For each variant call, if 
   chrom:pos is not in dict or if chrom:pos are present but alt doesn't match
   then it's unique. 
3. Output chrom:pos, ref, alt from second pipeline

Group ids: 1 = chrom, 2 = pos, 3 = ref, 4 = alt, 5 = filter, 6 = GT value
'''

from __future__ import print_function
from __future__ import division
import re
import gzip
import sys

def builddict(fin):
	vdict = {}	
	for i,line in enumerate(fin):
		m = re.search(r'^([0-9XY]+)\s+(\d+)\s+\S+\s+([,ATCG]+)\s+([,ATCG]+)\s+'
					 '[.\d]+\s+([.\w]+)\s+\S+\s+\S+\s+([\d/]+)',line)  		
		if i % 500000 == 0:	print('building dict...')
		if m and 'PASS' in m.group(5) and '0/0' not in m.group(6):
			vkey = m.group(1) + ':' + m.group(2)
			vdict[vkey] = m.group(4)
	return(vdict)

def plcompare(vardict, fin):
	count = 0	
	rarecount = 0	
	var_diffs = open('./p'+pline+'_'+sam+'_diffs.test','w')	
#	rejects = open('./p'+pline+'_'+sam+'_diffs.trash','w')
	print('Pipeline {0} {1}.vcf unique variants:'.format(pline,sam),
		 file=var_diffs)
	for i,line in enumerate(fin):
		m = re.search(r'^([0-9XY]+)\s+(\d+)\s+\S+\s+([,ATCG]+)\s+([,ATCG]+)\s+'
					 '[.\d]+\s+([.\w]+)\s+\S+\s+\S+\s+([\d/]+)',line)
		if (i % 500000 == 0):	print('parsing file...')
		if m and 'PASS' in m.group(5) and '0/0' not in m.group(6):
			vkey = m.group(1) + ':' + m.group(2)
			if vkey not in vardict:
				count += 1
				print('{0}\t{1}\t{2}\t{3}\t{4}'.format(vkey,m.group(3),
					 m.group(4),m.group(5),m.group(6)),file=var_diffs)
			elif vkey in vardict:
				if vardict[vkey] != m.group(4):
					count += 1
					rarecount += 1
					print('**{0}\t{1}\t{2}\t{3}\t{4}'.format(vkey,m.group(3),
						 m.group(4),m.group(5),m.group(6)),file=var_diffs)
#		elif m and 'PASS' in m.group(5):
#			print('{0}'.format(m.groups()),file=rejects)
	rarefrac = rarecount / count
	print('{0} total rare cases, {1:.2f} of total'.format(rarecount,rarefrac),
		 file=var_diffs)
	var_diffs.close()
#	rejects.close()		

pline = sys.argv[1]
sam = sys.argv[2]
finp1 = gzip.open('/projects/big/BIG-Pilot-Projects/021115-Statin-Induced-'
					'Myopathy/Project_NOW_01348_WGS_2015-02-03/Sample_'+sam+
					'/analysis/'+sam+'.recalibrated.haplotypeCalls.vcf.gz','r')
finp2 = open('/gpfs/scratch/jw24/gatk/Debug/phaseIII_3.3.0/'
                  'split_vcf/p'+sam+'.vcf','r')

if '1' in pline:
	variants = builddict(finp2)
	print(len(variants))
	plcompare(variants, finp1)
elif '2' in pline:
	variants = builddict(finp1)
	print(len(variants)) 
	plcompare(variants, finp2)
else:
	print('Usage: python unique_variants.py [pipeline #] [sample #]')
	sys.exit(1)

finp1.close()
finp2.close()
