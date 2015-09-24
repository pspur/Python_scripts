# Create .csv of sample avg read depths, given .out files created
# by WESavg_rd_and_cov.py

from __future__ import print_function
import glob
import re
import subprocess

fout = open('./WESread_depths.csv','w')
fout.write('SampleID,Read_Depth\n')
for file in glob.glob('*.out'):
	fin = open(file)
	lines = fin.readlines()
	rd = re.search('(\d+)\n$',lines[0])
	fout.write('{0}_AAA,{1}.0\n'.format(file[0:5],rd.group(1)))
