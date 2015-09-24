#Running script requires one command-line argument: the sample # you want to process.

from __future__ import division
from __future__ import print_function
import subprocess
import re
import sys


# Takes in a sample #.
# Returns length of sample.
def sam_length(sample):
    wc = subprocess.Popen(['wc', '-l', '/gpfs/scratch/paulspur/seqdataoutput/'
	            		'read_depth_coverage/output/{0}_depth.out'.format(sample)],
				         stdout=subprocess.PIPE)
    lenline, errlen = wc.communicate()
    len = int(re.search('(\d+)',lenline).group(1))
    return(len)


# Takes in a sample # and a sample length.
# Returns nothing, outputs a file containing average read depth and coverages
# for sample.
def avg_rd_and_coverage(sample,samlength):
	fout = open('./{0}_RD_avg.out'.format(sam),'w')
	sampledepth = 0
	sam10x = 0
	sam15x = 0
	sam20x = 0
	sam40x = 0
	sam50x = 0
	with open('/gpfs/scratch/paulspur/seqdataoutput/read_depth_coverage/'
              'output/{0}_depth.out'.format(sample),'r') as fin:
		fin.next()
		for line in fin:
			depth = int(re.search('\w+:\d+\s(\d+)',line).group(1))
			sampledepth += depth
			if depth >= 50:		
				sam50x += 1
			if depth >= 40:		
				sam40x += 1
			if depth >= 20:		
				sam20x += 1
			if depth >= 15:		
				sam15x += 1
			if depth >= 10:		
				sam10x += 1
	sam50x_cov = sam50x / samlength
	sam40x_cov = sam40x / samlength
	sam20x_cov = sam20x / samlength
	sam15x_cov = sam15x / samlength
	sam10x_cov = sam10x / samlength
	sample_rd_avg = sampledepth / samlength
	fout.write('sample {0} total read depth avg: {1:.2f}\n'.format
				(sample,sample_rd_avg))
	fout.write('sample {0} coverages:\n10x\t15x\t20x\t40x\t50x\n'.format(sample))
	fout.write('{0:.2%}\t{1:.2%}\t{2:.2%}\t{3:.2%}\t{4:.2%}\n'.format
				(sam10x_cov,sam15x_cov,sam20x_cov,sam40x_cov,sam50x_cov))
	fout.close()

if len(sys.argv) != 2: 
	print('Usage: python avg_rd_and_cov.py [sample #]')
	sys.exit(1)

sam = sys.argv[1]
samlength = sam_length(sam)
avg_rd_and_coverage(sam,samlength)
