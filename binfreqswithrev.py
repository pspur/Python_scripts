from __future__ import division
from __future__ import print_function
import re
import sys
import csv

infile = sys.argv[1]
#outfile = sys.argv[2]
bin_size = int(sys.argv[2])
interval = sys.argv[3]
start_loc = int(interval.split('-')[0])
end_loc = int(interval.split('-')[1])
bins = {}
var_total = 0

    
with open(infile,'r') as fin:
#  w = csv.writer(fout)
  for line in fin:
    m = re.search(r'^\d+\s+(\d+)',line)
    if m:
      pos = int(m.group(1))
      if start_loc < end_loc: # forwards interval
        if pos >= start_loc and pos <= end_loc:
          bin = start_loc + bin_size*(int((pos - start_loc) / bin_size))
          bins[bin] = bins.get(bin, 0) + 1
          var_total += 1
      else: # backwards interval
        if pos <= start_loc and pos >= end_loc:
          bin = start_loc + bin_size*int(((pos - start_loc) / bin_size))
          bins[bin] = bins.get(bin, 0) + 1
          var_total += 1
  for count in (bins.values()):
    count = round(count/var_total,3)
#  w.writerow(bins.keys())
#  w.writerow(bins.values())
print('# Source file: {0}, {1} total variants, bin size of {2}, interval: {3}\n'
      '# bin start position\t% of total variants in bin'.format(infile,var_total,bin_size,interval))
'''for bin in (bins.keys()):
  print('{0},'.format(bin), end = '')
print()
for count in (bins.values()):
  print('{0},'.format(round(count/var_total,3)), end = '')
print()'''
for bin,count in sorted(bins.items()):
  print('{0},{1}'.format(bin,round(count/var_total,3)))
