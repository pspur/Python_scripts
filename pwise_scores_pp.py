import re

sampleid = []

with open('/home/paulspur/R_scripts/vcfs.txt','r') as vcf_paths:
  for line in vcf_paths:
    m = re.search(r'Sample_(\d+)',line)
    sampleid.append(m.group(1))
sampleid.sort()

with open('/home/paulspur/R_scripts/pwise_scores_pp.txt','w') as fout:
  fout.write('sID  ')
  for sid in sampleid:
    fout.write('{0},'.format(sid))
  with open('/home/paulspur/R_scripts/pwise_scores_complete.txt','r') as fin:
    for line in fin:
      line = line.strip()
      m = re.search('(\d{4}) <> (\d{4}) = (\d?\.\d+)',line)
      if m.group(2) == '1344':
        fout.write('\n{0},{1}'.format(m.group(1),m.group(3)))
      else: fout.write(',{0}'.format(m.group(3)))
