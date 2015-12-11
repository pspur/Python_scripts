# This script takes in 2 command-line arguments, the names of the input and 
# output files. 
# Input file is expected to contain fields with format:
# Chr Start End Ref Alt Func.refGene Gene.refGene GeneDetail.refGene ExonicFunc.refGene, etc
# Output is 2 separate files, one with .gt appended and one with .sev appended 
# to provided output name.

from __future__ import print_function
import sys

# Taken from http://annovar.openbioinformatics.org/en/latest/user-guide/gene/
# with 'synonymous' replaced by '.' to account for presence of 'splicing' Func
# and removal of any synonymous variant
severity = {'frameshift insertion':11,
            'frameshift deletion':10,
            'frameshift block substition':9,
            'stopgain':8,
            'stoploss':7,
            'nonframeshift insertion':6,
            'nonframeshift deletion':5,
            'nonframeshift block substitution':4,
            'nonsynonymous SNV':3,
            '.':2,
            'unknown':1}

# Find header line which lists sample names. Parse them out and store them in list. 
def getSampleIDs(f):
  discards = ('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
  with open(f,'r') as fin:
    for line in fin:
      line = line.strip()
      if line.startswith('#CHROM'):
        sams = line.split('\t')
        sams = [item for item in sams if item not in discards]
        break
  return(sams)

# General flow: for each line in file: split by tabs, process only lines with
# desired attributes. For each one of those, convert ExonicFunc.refGene (severity)
# field to numerical value and for each gene in the Gene.refGene field add the
# sample:severity to the gene's severity dict, or update it if present and 
# current severity > stored severity.
# At the same time, perform the same process on each line for the genotype dict,
# updating it only when currenty severity >= stored severity and gt count > stored
# count.
def collapse_to_matrix(f,num):
  nt = ['A','C','T','G']
  sevdict = {}
  gtdict = {}
  with open(f,'r') as fin:
    for line in fin:
      line = line.strip()
      if not line.startswith('#'):
        cols = line.split('\t')
        Func_refGene = cols[5]
        Gene_refGene = cols[6]
        ExonicFunc_refGene = cols[8]
        Ref = cols[3]
        Alt = cols[4]
        if ((Func_refGene.startswith('exonic')) or (Func_refGene.startswith('splicing'))) and \
            (not ExonicFunc_refGene.startswith('synonymous')) and \
            (Ref in nt) and (Alt in nt): # only SNPs
          sev = severity[cols[8]]
          if (',' in Gene_refGene): # if multiple genes in field
            genes = Gene_refGene.split(',')
          else:
            genes = []
            genes.append(Gene_refGene)
          for gene in genes:
            sevdict.setdefault(gene,{})
            gtdict.setdefault(gene,{})
            for i,sam in enumerate(cols[-num:]): # the genotype fields  
              gt = sam.split(':')[0].split('/')
              gtsum = 0
              if '.' in gt:
                gtsum = 0
              else:
                gtsum = int(gt[0]) + int(gt[1])
              if gtsum == 0:
                sevdict[gene].setdefault(samples[i],0)
                gtdict[gene].setdefault(samples[i],0)
              elif (samples[i] in sevdict[gene]) and (sev >= sevdict[gene][samples[i]]):
                sevdict[gene][samples[i]] = sev
                if (gtsum > gtdict[gene][samples[i]]): 
                  gtdict[gene][samples[i]] = gtsum
              else:
                sevdict[gene].setdefault(samples[i],sev)
       	    gtdict[gene].setdefault(samples[i],gtsum)
  return(sevdict,gtdict)

def print_severity_matrix(sam,sev_matrix,outfile):
  with open(outfile+'.sev','w') as fout:
    fout.write('\t')
    for s in sam:
      fout.write('{}\t'.format(s))
    fout.write('\n')
    for gene,value in sev_matrix.iteritems():
      fout.write('{}\t'.format(gene))
      for k,v in sorted(value.iteritems()):
        fout.write('{}\t'.format(v))
      fout.write('\n')

def print_gt_matrix(sam,gt_matrix,outfile):
  with open(outfile+'.gt','w') as fout:
    fout.write('\t')
    for s in sam:
      fout.write('{}\t'.format(s))
    fout.write('\n')
    for gene,value in gt_matrix.iteritems():
      fout.write('{}\t'.format(gene))
      for k,v in sorted(value.iteritems()):
        fout.write('{}\t'.format(v))
      fout.write('\n')

if __name__ == '__main__':
  filein = sys.argv[1]
  fileout = sys.argv[2]
  samples = getSampleIDs(filein)
  number_of_samples = len(samples)
  sev_dict,gt_dict = collapse_to_matrix(filein,number_of_samples)
  print_severity_matrix(samples,sev_dict,fileout)
  print_gt_matrix(samples,gt_dict,fileout)
