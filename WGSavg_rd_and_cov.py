#Running script requires one command-line argument: the sample # you want to process.

from __future__ import division
from __future__ import print_function
import subprocess
import re
import sys

chroms = [n for n in range(1,23)]
chroms.append('X')

# Takes in a sample #.
# Returns a list of chromosome lengths. Chromosomes included are: 1-22,X.
def chr_lengths(sample):
	chrlens = []
	for chrom in chroms:
		head = subprocess.Popen(['head', '-n2', '/gpfs/scratch/paulspur/gatk/'
				'sample/{0}/{0}chr{1}_depth.out'.format(sample,chrom)],
				 stdout=subprocess.PIPE)
		tail = subprocess.Popen(['tail', '-1', '/gpfs/scratch/paulspur/gatk/'
				'sample/{0}/{0}chr{1}_depth.out'.format(sample,chrom)], 
				 stdout=subprocess.PIPE)
		headline, errhead = head.communicate()
		tailline, errtail = tail.communicate()
		chrstart = re.search('\w+:(\d+)',headline)
		chrend = re.search('\w+:(\d+)',tailline)
		chrlen = int(chrend.group(1)) - int(chrstart.group(1)) + 1
		chrlens.append(chrlen)
	return(chrlens)


# Takes in a sample # and a list of chromosome lengths.
# Returns nothing, outputs a file containing average read depth and coverages
# for each chromosome and average read depth and coverages for entire sample.
def avg_rd_and_coverage(sample,chrlengths):
	fout = open('./{0}_RD_avg.out'.format(sam),'w')
	sampledepth = 0
	samplelen = sum(chrlengths)
	sam10x = 0
	sam15x = 0
	sam20x = 0
	sam40x = 0
	sam50x = 0
	for i, chrom in enumerate(chroms):
		with open('/gpfs/scratch/paulspur/gatk/sample/{0}/'
			  '{0}chr{1}_depth.out'.format(sample,chrom),'r') as fin:
			chrdepth = 0
			chr10x = 0
			chr15x = 0
			chr20x = 0
			chr40x = 0
			chr50x = 0			
			fin.next()
			for line in fin:
				depth = int(re.search('\w+:\d+\s(\d+)',line).group(1))
				chrdepth += depth
				sampledepth += depth
				if depth >= 50:		
					sam50x += 1
					chr50x += 1
				if depth >= 40:		
					sam40x += 1
					chr40x += 1
				if depth >= 20:		
					sam20x += 1
					chr20x += 1
				if depth >= 15:		
					sam15x += 1
					chr15x += 1
				if depth >= 10:		
					sam10x += 1
					chr10x += 1
		chr50x_cov = chr50x / chrlengths[i]
		chr40x_cov = chr40x / chrlengths[i]
		chr20x_cov = chr20x / chrlengths[i]
		chr15x_cov = chr15x / chrlengths[i]
		chr10x_cov = chr10x / chrlengths[i]
		chr_rd_avg = chrdepth / chrlengths[i]
		fout.write('chr{0} read depth avg: {1:.2f}\n'.format(chrom,chr_rd_avg))
		fout.write('chr{0} coverages:\n10x\t15x\t20x\t40x\t50x\n'.format(chrom))
		fout.write('{0:.2%}\t{1:.2%}\t{2:.2%}\t{3:.2%}\t{4:.2%}\n\n'.format
					(chr10x_cov,chr15x_cov,chr20x_cov,chr40x_cov,chr50x_cov))
	sam50x_cov = sam50x / samplelen
	sam40x_cov = sam40x / samplelen
	sam20x_cov = sam20x / samplelen
	sam15x_cov = sam15x / samplelen
	sam10x_cov = sam10x / samplelen
	sample_rd_avg = sampledepth / samplelen
	fout.write('sample {0} total read depth avg: {1:.2f}\n'.format
				(sample,sample_rd_avg))
	fout.write('sample {0} coverages:\n10x\t15x\t20x\t40x\t50x\n'.format(sample))
	fout.write('{0:.2%}\t{1:.2%}\t{2:.2%}\t{3:.2%}\t{4:.2%}\n'.format
				(sam10x_cov,sam15x_cov,sam20x_cov,sam40x_cov,sam50x_cov))
	fout.close()

if len(sys.argv) != 2: 
	print('Usage: python avg_readdepth.py [sample #]')
	sys.exit(1)

sam = sys.argv[1]
chrlengths = chr_lengths(sam)
avg_rd_and_coverage(sam,chrlengths)
