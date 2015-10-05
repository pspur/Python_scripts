# This script takes in 2 command-line arguments, the names of the input and output files.

from __future__ import print_function
import sys
import re

# Taken from http://annovar.openbioinformatics.org/en/latest/user-guide/gene/
# with 'synonymous' replaced by '.' to account for presence of 'splicing' Func
severity = {'frameshift insertion':1,
            'frameshift deletion':2,
            'frameshift block substition':3,
            'stopgain':4,
            'stoploss':5,
            'nonframeshift insertion':6,
            'nonframeshift deletion':7,
            'nonframeshift block substitution':8,
            'nonsynonymous SNV':9,
            '.':10,
            'unknown':11}

# Taken from rush:/projects/big/SIM/WGS/vcf/SIM_snp_recalibrated_indel_recalibrated.vcf
samples = ('1344','1394','1406','1409','2003','2017','2021','2023','2024','2042',
          '2048','2052','2059','2063','2066','2072','2074','2077','2078','2087',
          '2102','2134','2142','2145','2152','2155','2157','2161','2162','2166',
          '2168','2170','2171','2173','2179','2181','2184','2185','2193','2194',
          '2196','2198','2208','2214','2221','2239','2257','2261','2269','2517')

def collapse(f):
  samdict = {}
  with open(f,'r') as fin:
    for line in fin:
      cols = line.split('\t')
      Func_refGene = cols[5]
      Gene_refGene = cols[6]
      ExonicFunc_refGene = cols[8]
      Ref = cols[3]
      Alt = cols[4]
      if ((Func_refGene.startswith('exonic')) or (Func_refGene.startswith('splicing'))) and \
          (not ExonicFunc_refGene.startswith('synonymous')) and \
          (len(Ref) == 1) and (len(Alt) == 1): # only SNPs
        sev = severity[cols[8]]
        if (',' in Gene_refGene): # if multiple genes in field
          genes = Gene_refGene.split(',')
        else:
          genes = []
          genes.append(Gene_refGene)
        for gene in genes:
          samdict.setdefault(gene,{})
          for i,sam in enumerate(cols[-50:]): # the 50 genotype fields  
            if (('0/0' in sam) or ('./.' in sam)):
              samdict[gene].setdefault(samples[i],99)
            elif (samples[i] in samdict[gene]) and (sev < samdict[gene][samples[i]]):
              samdict[gene][samples[i]] = sev
            else:
              samdict[gene].setdefault(samples[i],sev) 
  return(samdict)


if __name__ == '__main__':
  filein = sys.argv[1]
  fileout = sys.argv[2]
  samdict = collapse(filein)
  with open(fileout,'w') as fout:
    fout.write('\t')
    for s in samples:
      fout.write('{}\t'.format(s))
    fout.write('\n')
    for gene,value in samdict.iteritems():
      fout.write('{}\t'.format(gene))
      for k,v in sorted(value.iteritems()):
        if v == 99:
          fout.write('\t')
        else:
          fout.write('{}\t'.format(v))
      fout.write('\n')
