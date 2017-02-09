#!/usr/bin/env python
__version__ = '1.0'
__author__  = 'H. Rakotoarivelo'

import os
import sys
import time
import optparse
import subprocess

if __name__ == '__main__':
    
    # extract params
    parser = optparse.OptionParser()
    parser.add_option('-p', dest='cores' , default=1, type='int',help='nb of cores')
    parser.add_option('-c', dest='cxx', default='intel', help='specify compiler [intel|gnu]')
    parser.add_option('-n', dest='no_numa', action='store_true',help='non NUMA-aware')
    parser.add_option('-m', dest='mcdram', help='allocate on MCDRAM (Intel KNL)')
    parser.add_option('-n', dest='N'     , default='80000000', help='stream array size')
    parser.add_option('-t', dest='ntimes', default='100', help='ntimes')

    # retrieve params
    (param,args) = parser.parse_args()
    assert(param.cores > 0) 

    # generate Makefile
    if param.cxx == 'intel':
        with open('Makefile','w') as f:
            f.write('CC=icc\n')
            f.write('CFLAGS= -O3 -xCORE-AVX2 -qopenmp -lnuma -DSTREAM_ARRAY_SIZE='+param.N)
            f.write(' -DNTIMES='+param.ntimes+' -ffreestanding ')
            s = ' -DNON_NUMA\n' if param.no_numa else '\n'
            f.write(s+'\n') 
            f.write('all: numa_stream\n')    
            f.write('numa_stream: stream.c\n')    
            f.write('\t$(CC) $(CFLAGS) stream.c -o numa_stream\n\n')
            f.write('clean:\n')
            f.write('\trm -rf numa_stream *.o')

    # set affinity (1 thread/core)
    os.environ['KMP_AFFINITY']='granularity=fine,scatter'
    subprocess.call(['make clean'])
    subprocess.call(['make'])

    for i in range(param.cores):
        os.environ['OMP_NUM_THREADS']=str(i+1)
        subprocess.call(['./numa_stream'])
