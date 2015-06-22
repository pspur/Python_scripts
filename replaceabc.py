# Create file with 'abc' (already been calculated) scores replaced
# with correct scores contained in earlier lines

import re

with open('/home/paulspur/R_scripts/pwise_scores_complete.txt','w') as fout:
  with open('/home/paulspur/R_scripts/pwise_scores.txt','r') as fin:
    for line in fin:
      m = re.search('(\d{4}) <> (\d{4}) = (\d+\.\d+|abc)',line)
      if m.group(3) != 'abc':
        fout.write('{0} <> {1} = {2}\n'.format(m.group(1),m.group(2),m.group(3)))
      else:
        for line in open('/home/paulspur/R_scripts/pwise_scores.txt','r'):
          mflip = re.search('(\d{4}) <> (\d{4}) = (\d+\.\d+)',line)
          if mflip and mflip.group(1) == m.group(2) and mflip.group(2) == m.group(1):
            fout.write('{0} <> {1} = {2}\n'.format(m.group(1),m.group(2),mflip.group(3)))
	
