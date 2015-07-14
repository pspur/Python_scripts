# Create .csv of sample coverages, given .out created by avg_readdepth.py

from __future__ import print_function
import glob
import re

fout = open('./coverages.csv','w')
fout.write('SampleID,10x,15x,20x,40x,50x\n')
for file in glob.glob('*.out'):
	fin = open(file)
	lines = fin.readlines()
	cov = re.search('(\d+.\d+)%\t(\d+.\d+)%\t(\d+.\d+)%\t(\d+.\d+)%\t(\d+.\d+)',lines[118])
	fout.write('{0}_AAA,{1},{2},{3},{4},{5}\n'.format(file[0:4],
		   cov.group(1),cov.group(2),cov.group(3),cov.group(4),
		   cov.group(5)))
