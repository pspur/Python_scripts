# command line ex: python sample_concordance_vcf.py [vcf to parse] [type of correlation] > [fileout] 
# correlation type: c for concordance, p for pearsons



from __future__ import division
from __future__ import print_function
import sys
import numpy

# Find header line which lists sample names. Parse them out and store them in list. 
def getSampleInfo(f):
    discards = ('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')
    with open(f,'r') as fin:
        for line in fin:
            line = line.strip()
            if (line.startswith('#') and not line.startswith('##')):
                header_line = line.split('\t')
                sams = [item for item in header_line if item not in discards]
                break
    return(sams)


# For every sample, if variant is present in sample append 1, else append 0
def buildSampleVarMatrix(f,sams):
    nt = ['A','C','T','G']
    var_total = 0
    sam_count = len(sams)
    var_counts = [[] for x in range(sam_count)]
    with open(f,'r') as fin:
        for line in fin:
            line = line.strip()
            if (line.startswith('#') and not line.startswith('##')):
                headers = line.split('\t')
            elif (not line.startswith('#')):
                cols = line.split('\t')
                gatk_filter = cols[headers.index('FILTER')]
                ref = cols[headers.index('REF')]
                alt = cols[headers.index('ALT')]
                alleles = cols[-sam_count:]
                if ((gatk_filter == 'PASS') and (ref in nt) and (alt in nt)) and ('./.' not in alleles):
                    var_total += 1
                    for sam_parent,allele_parent in enumerate(alleles):
                        if (not '0/0' in allele_parent):
                            var_counts[sam_parent].append(1)
                        else:
                            var_counts[sam_parent].append(0)
    return(var_counts, var_total)


# Compare each sample to every sample. Skip comparisons that are the reverse of ones
# that have already been done then fill in after. For each comparison, get sum of variant logical ANDS,
# calculate overall concordance. 
def calcConcordance(vm, total):
    concordMatrix = [[0 for i in range(len(vm))] for j in range(len(vm))]
    for s1idx,sam1 in enumerate(vm):
        for vidx,var1 in enumerate(sam1):
            for s2idx in range(s1idx,len(vm)):
                if (var1 == vm[s2idx][vidx]): 
                    concordMatrix[s1idx][s2idx] += 1
    # fill in bottom half of matrix by copying from top half
    for xidx in range(len(vm)):
        for yidx in range(xidx,len(vm)):
            concordMatrix[yidx][xidx] = concordMatrix[xidx][yidx] 
    # calculate concordances
    for i in range(len(concordMatrix)):
        for j in range(len(concordMatrix)):
            concordMatrix[i][j] = concordMatrix[i][j]/total
    return(concordMatrix)


def writeMatrix(m,sams):
    for s in sams:
        print('\ts{}'.format(s),end='')
    for i,row in enumerate(m):
        print()
        print('s{}\t'.format(sams[i]),end='')
        for value in row:
            print('{:.3f}\t'.format(value),end='')
    print()
     

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python sample_concordance_vcf.py [vcf to parse] [type of correlation] > [fileout]\n'
              'correlation type: c for concordance, p for pearsons',file=sys.stderr)
        sys.exit(1)
    if sys.argv[2].lower() not in ('c','p'):
        print('Usage: python sample_concordance_vcf.py [vcf to parse] [type of correlation] > [fileout]\n'
              'correlation type: c for concordance, p for pearsons',file=sys.stderr)
        sys.exit(1)
    filein = sys.argv[1]
    samples = getSampleInfo(filein)
    var_matrix, total_vars = buildSampleVarMatrix(filein, samples)
    if (sys.argv[2].lower() == 'c'):
        ccords = calcConcordance(var_matrix, total_vars)
    elif (sys.argv[2].lower() == 'p'):
        ccords = numpy.corrcoef(var_matrix)
    writeMatrix(ccords, samples)
