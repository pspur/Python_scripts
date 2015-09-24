# Create .csv of sample avg read depths, given .out files created
# by WGSavg_rd_and_cov.py

from __future__ import print_function
import glob
import re
import subprocess

fout = open('./read_depths.csv','w')
fout.write('SampleID,Read_Depth\n')
for file in glob.glob('*.out'):
	fin = open(file)
	lines = fin.readlines()
	rd = re.search('(\d+)\n$',lines[115])
	fout.write('{0}_AAA,{1}.0\n'.format(file[0:4],rd.group(1)))
