### Final vcf produced by pipeline is placed in svelter working directory
### Slurm logs are placed in a log directory created under svelter working directory

### TODO: 1. Change svelter location in subprocess calls to global one
###          once it's properly set up.

import os
import sys
import time
import subprocess
import argparse

def createNullModel(bam,sam,workdir,logdir,scriptdir):
    with open(r'{0}/sv_nullmodel_{1}.sh'.format(scriptdir,sam),'w') as nm:
        nm.write('#!/bin/sh\n\n')
        nm.write('#SBATCH --partition=general-compute\n')
        nm.write('#SBATCH --time=30:00\n')
        nm.write('#SBATCH --nodes=1\n')
        nm.write('#SBATCH --mem=8000\n')
        nm.write('#SBATCH --ntasks-per-node=1\n')
        nm.write('#SBATCH --job-name=sv_nullmodel_{0}\n'.format(sam))
        nm.write('#SBATCH --output={0}/nullmodel/sv_nullmodel_{1}.log'
                 '\n'.format(logdir,sam))
        nm.write('#SBATCH --qos=supporters\n')
        nm.write('#SBATCH --account=big\n')
        nm.write('echo "SLURM_JOBID="$SLURM_JOBID\n')
        nm.write('echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST\n')
        nm.write('echo "SLURM_NNODES"=$SLURM_NNODES\n')
        nm.write('echo "SLURMTMPDIR="$SLURMTMPDIR\n')
        nm.write('echo "working directory = "$SLURM_SUBMIT_DIR\n\n')
        nm.write('ulimit -s unlimited\n')
        nm.write('module load python\n')
        nm.write('module load samtools\n')
        nm.write('module list\n\n')
        nm.write('echo "Launch job"\n\n')
        nm.write('python svelter.py NullModel --sample {0} '
                                             '--workdir {1} '
                                             '--null-copyneutral-perc 0.01 '
                                             '\n\n'.format(bam,workdir))
        nm.write('echo "All Done!"')

def createBPSearch(bam,sam,workdir,logdir,scriptdir,chroms,job):
    for chrom in chroms:
        with open(r'{0}/sv_bpsearch_{1}chr{2}.sh'.format(scriptdir,sam,chrom),'w') as bps:
                bps.write('#!/bin/sh\n\n')
                bps.write('#SBATCH --partition=general-compute\n')
                if type(chrom) is int and chrom <= 2:
                    bps.write('#SBATCH --time=10:00:00\n')
                    bps.write('#SBATCH --nodes=1\n')
                    bps.write('#SBATCH --mem=48000\n')
                elif type(chrom) is int and chrom <= 10:
                    bps.write('#SBATCH --time=10:00:00\n')
                    bps.write('#SBATCH --nodes=1\n')
                    bps.write('#SBATCH --mem=24000\n')
                else:
                    bps.write('#SBATCH --time=5:00:00\n')
                    bps.write('#SBATCH --nodes=1\n')
                    bps.write('#SBATCH --mem=16000\n')
                bps.write('#SBATCH --ntasks-per-node=1\n')
                bps.write('#SBATCH --job-name=sv_bpsearch_{0}chr{1}\n'.format(sam,chrom))
                bps.write('#SBATCH --output={0}/bpsearch/sv_bpsearch_{1}chr{2}.log'
                          '\n'.format(logdir,sam,chrom))
                bps.write('#SBATCH --qos=supporters\n')
                bps.write('#SBATCH --account=big\n')
                if job:
                    bps.write('#SBATCH --dependency=afterok:{0}\n'.format(job))
                bps.write('echo "SLURM_JOBID="$SLURM_JOBID\n')
                bps.write('echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST\n')
                bps.write('echo "SLURM_NNODES"=$SLURM_NNODES\n')
                bps.write('echo "SLURMTMPDIR="$SLURMTMPDIR\n')
                bps.write('echo "working directory = "$SLURM_SUBMIT_DIR\n\n')
                bps.write('ulimit -s unlimited\n')
                bps.write('module load python\n')
                bps.write('module load samtools\n')
                bps.write('module list\n\n')
                bps.write('echo "Launch job"\n\n')
                bps.write('python svelter.py BPSearch --sample {0} '
                                                     '--workdir {1} '
                                                     '--chromosome {2}\n\n'.format(bam,workdir,chrom))
                bps.write('echo "All Done!"')

def createBPIntegrate(bam,sam,workdir,logdir,scriptdir,jobs):
    with open(r'{0}/sv_bpintegrate_{1}.sh'.format(scriptdir,sam),'w') as bpi:
        bpi.write('#!/bin/sh\n\n')
        bpi.write('#SBATCH --partition=general-compute\n')
        bpi.write('#SBATCH --time=30:00\n')
        bpi.write('#SBATCH --nodes=1\n')
        bpi.write('#SBATCH --mem=8000\n')
        bpi.write('#SBATCH --ntasks-per-node=1\n')
        bpi.write('#SBATCH --job-name=sv_bpintegrate_{0}\n'.format(sam))
        bpi.write('#SBATCH --output={0}/bpintegrate/sv_bpintegrate_{1}.log'
                  '\n'.format(logdir,sam))
        bpi.write('#SBATCH --qos=supporters\n')
        bpi.write('#SBATCH --account=big\n')
        if jobs:
            bpi.write('#SBATCH --dependency=afterok:{0}\n'.format(':'.join(jobs)))
        bpi.write('echo "SLURM_JOBID="$SLURM_JOBID\n')
        bpi.write('echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST\n')
        bpi.write('echo "SLURM_NNODES"=$SLURM_NNODES\n')
        bpi.write('echo "SLURMTMPDIR="$SLURMTMPDIR\n')
        bpi.write('echo "working directory = "$SLURM_SUBMIT_DIR\n\n')
        bpi.write('ulimit -s unlimited\n')
        bpi.write('module load python\n')
        bpi.write('module load samtools\n')
        bpi.write('module list\n\n')
        bpi.write('echo "Launch job"\n\n')
        bpi.write('python svelter.py BPIntegrate --sample {0} '
                                                '--workdir {1} '
                                                '\n\n'.format(bam,workdir))
        bpi.write('echo "All Done!"')
        
def createSVPredict(bam,sam,workdir,logdir,scriptdir,job):
    bamname = bam.split('/')[-1]
    bamtextname = bamname.replace('.bam','.txt')
    with open(r'{0}/sv_svpredict_{1}.sh'.format(scriptdir,sam),'w') as svp:
        svp.write('#!/bin/sh\n\n')
        svp.write('#SBATCH --partition=general-compute\n')
        svp.write('#SBATCH --time=6:00:00\n')
        svp.write('#SBATCH --nodes=1\n')
        svp.write('#SBATCH --mem=16000\n')
        svp.write('#SBATCH --ntasks-per-node=1\n')
        svp.write('#SBATCH --job-name=sv_svpredict_{0}\n'.format(sam))
        svp.write('#SBATCH --output={0}/svpredict/sv_svpredict_{1}.log'
                  '\n'.format(logdir,sam))
        svp.write('#SBATCH --qos=supporters\n')
        svp.write('#SBATCH --account=big\n')
        if job:
            svp.write('#SBATCH --dependency=afterok:{0}\n'.format(job))
        svp.write('echo "SLURM_JOBID="$SLURM_JOBID\n')
        svp.write('echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST\n')
        svp.write('echo "SLURM_NNODES"=$SLURM_NNODES\n')
        svp.write('echo "SLURMTMPDIR="$SLURMTMPDIR\n')
        svp.write('echo "working directory = "$SLURM_SUBMIT_DIR\n\n')
        svp.write('ulimit -s unlimited\n')
        svp.write('module load python\n')
        svp.write('module load samtools\n')
        svp.write('module list\n\n')
        svp.write('echo "Launch job"\n\n')
        svp.write('python svelter.py SVPredict --sample {0} '
                                              '--workdir {1} '
                                              '--bp-file {1}/bp_files.{2}/{3}'
                                              '\n\n'.format(bam,workdir,bamname,bamtextname))
        svp.write('echo "All Done!"')
        
def createSVIntegrate(bam,sam,workdir,logdir,scriptdir,job):
    bamname = bam.split('/')[-1]
    outprefix = bamname.replace('.bam','_svelter')
    with open(r'{0}/sv_svintegrate_{1}.sh'.format(scriptdir,sam),'w') as svp:
        svp.write('#!/bin/sh\n\n')
        svp.write('#SBATCH --partition=general-compute\n')
        svp.write('#SBATCH --time=10:00\n')
        svp.write('#SBATCH --nodes=1\n')
        svp.write('#SBATCH --mem=8000\n')
        svp.write('#SBATCH --ntasks-per-node=1\n')
        svp.write('#SBATCH --job-name=sv_svintegrate_{0}\n'.format(sam))
        svp.write('#SBATCH --output={0}/svintegrate/sv_svintegrate_{1}.log'
                  '\n'.format(logdir,sam))
        svp.write('#SBATCH --qos=supporters\n')
        svp.write('#SBATCH --account=big\n')
        if job:
            svp.write('#SBATCH --dependency=afterok:{0}\n'.format(job))
        svp.write('echo "SLURM_JOBID="$SLURM_JOBID\n')
        svp.write('echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST\n')
        svp.write('echo "SLURM_NNODES"=$SLURM_NNODES\n')
        svp.write('echo "SLURMTMPDIR="$SLURMTMPDIR\n')
        svp.write('echo "working directory = "$SLURM_SUBMIT_DIR\n\n')
        svp.write('ulimit -s unlimited\n')
        svp.write('module load python\n')
        svp.write('module list\n\n')
        svp.write('echo "Launch job"\n\n')
        svp.write('python svelter.py SVIntegrate --prefix {1}/{0} '
                                              '--workdir {1} '
                                              '--input-path {1}/bp_files.{2}/'
                                              '\n\n'.format(outprefix,workdir,bamname))
        svp.write('echo "All Done!"')

# Execute slurm script, return job # for use as dependency
def runSlurmScript(script):
    sh = subprocess.Popen(['sbatch', '-D', 
                           '/projects/academic/big/paulspur/svelter/svelter_sv/',
                            script], stdout=subprocess.PIPE)
    slout, err = sh.communicate()
    job = slout.split(' ')[-1].strip()
    return(job)    

def createlogDirs(bams,workdir):
    logdir = workdir + '/' + bams.split('/')[-1].split('.')[0] + 'logs'
    sublogdirs = ['bpintegrate','bpsearch','nullmodel','svintegrate','svpredict']
    for sublogdir in sublogdirs:
        try:
            os.makedirs(logdir + '/' + sublogdir)
        except OSError:
            if not os.path.isdir(logdir + '/' + sublogdir):
                raise
    return(logdir)

# Run svelter steps in sequence. Only run steps which are in step list, which allows
# rerunning pipeline starting at specific step.
def main(args):
    numjobs = 0
    #chroms = [21,22]
    chroms = [n for n in range(1,23)]
    chroms.append('X')
    chroms.append('Y')
    chroms.append('MT')

    steps = ['nullmodel','bpsearch','bpintegrate','svpredict','svintegrate']
    steps = steps[steps.index(args.step):]
    logdir = createlogDirs(args.bams,args.work)
    
    # get username for use in upcoming squeue subprocess calls
    runid = subprocess.Popen(['id', '-u', '-n'],stdout=subprocess.PIPE)
    runidout, err = runid.communicate()
    uname = runidout.strip()

    # allow for a single bam to be run without txt file
    bamlist = []
    if not args.bams.endswith('.bam'):
        with open(args.bams,'r') as fin:
            for line in fin:
                if not line.startswith('#') and line.strip() != '':
                    bamlist.append(line.strip())
    else:
        bamlist.append(args.bams)

    print('Starting from step: {0}'.format(args.step))
    for i,bam in enumerate(bamlist, start=1):
        sam = bam.split('/')[-1].split('_')[0]
        
        # 29 jobs submitted per bam, keep submitted jobs under 990
        while(numjobs >= 960):
            time.sleep(30)
            squeue = subprocess.Popen(['squeue', '-h', '-u', '{0}'.format(uname)],
                                      stdout=subprocess.PIPE)
            squeueout, err = squeue.communicate()
            numjobs = len(squeueout.split('\n')) - 1

        print('{0}/{1}: Submitting jobs for sample {2}'.format(i,len(bamlist),sam))

        if 'nullmodel' in steps:
            createNullModel(bam,sam,args.work,logdir,args.scripts)
            script = '{0}/sv_nullmodel_{1}.sh'.format(args.scripts,sam)
            job = runSlurmScript(script)
            time.sleep(1)
        
        if 'bpsearch' in steps:
            if steps.index('bpsearch') == 0:
                job = False
            createBPSearch(bam,sam,args.work,logdir,args.scripts,chroms,job)
            jobs = []
            for chrom in chroms:
                script = '{0}/sv_bpsearch_{1}chr{2}.sh'.format(args.scripts,sam,chrom)
                job = runSlurmScript(script)
                jobs.append(job)
                time.sleep(1)
        
        if 'bpintegrate' in steps:
            if steps.index('bpintegrate') == 0:
                jobs = False
            createBPIntegrate(bam,sam,args.work,logdir,args.scripts,jobs)
            script = '{0}/sv_bpintegrate_{1}.sh'.format(args.scripts,sam)
            job = runSlurmScript(script)
            time.sleep(1)
        
        if 'svpredict' in steps:
            if steps.index('svpredict') == 0:
                job = False
            createSVPredict(bam,sam,args.work,logdir,args.scripts,job)
            script = '{0}/sv_svpredict_{1}.sh'.format(args.scripts,sam)
            job = runSlurmScript(script)
            time.sleep(1)
        
        if 'svintegrate' in steps:
            if steps.index('svintegrate') == 0:
                job = False
            createSVIntegrate(bam,sam,args.work,logdir,args.scripts,job)
            script = '{0}/sv_svintegrate_{1}.sh'.format(args.scripts,sam)
            job = runSlurmScript(script)
            time.sleep(1)
        
        squeue = subprocess.Popen(['squeue', '-h', '-u', '{0}'.format(uname)],
                              stdout=subprocess.PIPE)
        squeueout, err = squeue.communicate()
        numjobs = len(squeueout.split('\n')) - 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run svelter on linux using SLURM')
    parser.add_argument('-b', '--bams', 
                        action='store',
                        required=True,
                        help='Either single .bam or file containing bam locations, '
                             '1 bam per line')

    parser.add_argument('--scripts',
        	       	    action='store',
        	       	    required=True,
        	       	    help='Desired path for slurm scrips')

    parser.add_argument('-w', '--work',
        	       	    action='store',
        	       	    required=True,
        	       	    help='Svelter working directory, also where final vcfs reside')

    parser.add_argument('-s', '--step',
        	       	    default='nullmodel',
                        choices=['nullmodel', 
                                 'bpsearch', 
                                 'bpintegrate', 
                                 'svpredict', 
                                 'svintegrate'],
                        action='store',
       	        	    help='Desired starting step, default: nullmodel.')

    args = parser.parse_args()
    main(args)
