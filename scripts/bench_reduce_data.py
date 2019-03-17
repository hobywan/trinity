#!/usr/bin/env python
__author__    = 'Hoby Rakotoarivelo'
__license__   = 'GPL v3'
__version__   = '1.0'
__copyright__ = 'Copyright 2016, The trinity project'

import os
import sys
import time
import optparse
import subprocess
from .bench_export_gnuplot import *

#
def compact_benchmark_data(param):

    testcase = param.name
    nb_cores = int(param.cores)
    nb_field = int(param.fields)

    # for each kernel
    for k in range(1,5):
        target = 'results/benchs/'+testcase+'_'+str(k)+'.dat'
        # skip if already exists
        if os.path.exists(target):
            continue

        with open(target,'w') as output:
            prefix = 'results/_'+str(k)+'/'+testcase+'_w_'
            for rank in range(nb_cores):
                path = prefix + str(rank+1).rjust(2,'0')+'.dat'
                # skip if file does not exist
                if not os.path.exists(path):
                    continue

                count = 0
                stats = [0] * nb_field
                with open(path,'r') as f:
                    for line in f:
                        values = line.split()
                        for i in range(nb_field):
                            stats[i] += int(values[i])
                        count += 1

                assert(count > 0)
                for i in range(nb_field):
                    stats[i] /= count
                    output.write('{:01d}'.format(int(stats[i]))+'\t')
                output.write('\n')
        #
        print('>> '+target+'')

#
def clean_data(param):

    testcase = param.name
    nb_cores = int(param.cores)
    eps = [''] * 2

    for k in range(1,5):
        #
        prefix = 'results/_'+str(k)+'/'+testcase+"_w_"
        for rank in range(nb_cores):
            path = prefix + str(rank+1).rjust(2,'0')+'.dat'
            # skip if file does not exist
            if os.path.exists(path):
                print('>> rm '+path)
                os.remove(path)
        #
        data = 'results/benchs/'+testcase+'_'+str(k)+'.dat'
        if os.path.exists(data):
            print('>> rm '+data)
            os.remove(data)

    #
    script = 'results/benchs/'+testcase+'.gnuplot'
    eps[0] = 'results/'+testcase+'_runtime.eps'
    eps[1] = 'results/'+testcase+'_speedup.eps'
    if os.path.exists(script):
        print('>> rm '+script)
        os.remove(script)

    for i in range(0,2):
        if os.path.exists(eps[i]):
            print('>> rm '+eps[i])
            os.remove(eps[i])

    sys.exit(0)


# sum over the 3 architectures
def reduce_knl_stats(param):

    testcase = param.name
    nb_cores = int(param.cores)
    nb_field = 7
    nb_runs  = 7
    out = ['']*3

    #
    stats = [[]]
    run = ['']*3
    run[0] = 'knl'
    run[1] = 'hyp'
    run[2] = 'llc'

    for mode in range(3):
        # flush
        stats = [[0 for y in range(nb_field)] for x in range(nb_cores)]
        # set output path
        out[mode] = 'results/benchs/'+testcase+'_'+run[mode]+'_reduc.dat'
        # reduce over kernels
        for k in range(4):
            data = 'results/benchs/'+testcase+'_'+run[mode]+'_'+str(k+1)+'.dat'
            assert(os.path.exists(data))
            with open(data,'r') as f:
                c = 0
                for line in f:
                    values = line.split()
                    for i in range(nb_field):
                        stats[c][i] += int(values[i])
                    c = c+1
                assert(c == nb_runs)

        # write to file
        with open(out[mode],'w') as f:
            for r in range(nb_runs):
                for i in range(nb_field):
                    if i == 0:
                        stats[r][i] /= 4
                    f.write('{:01d}'.format(int(stats[r][i]))+'\t')
                f.write('\n')

# sum over the 3 architectures
def reduce_hasw_stats(param):

    testcase = param.name
    nb_cores = int(param.cores)
    nb_field = 7
    nb_runs  = 6
    out = ['']*3

    #
    stats = [[]]
    run = ['']*3
    run[0] = 'nor'
    run[1] = 'hyp'

    for mode in range(2):
        # flush
        stats = [[0 for y in range(nb_field)] for x in range(nb_cores)]
        # set output path
        out[mode] = 'results/benchs/'+testcase+'_'+run[mode]+'_reduc.dat'
        # reduce over kernels
        for k in range(4):
            data = 'results/benchs/'+testcase+'_'+run[mode]+'_'+str(k+1)+'.dat'
            assert(os.path.exists(data))
            with open(data,'r') as f:
                c = 0
                for line in f:
                    values = line.split()
                    for i in range(nb_field):
                        stats[c][i] += int(values[i])
                    c = c+1
                assert(c == nb_runs)

        # write to file
        with open(out[mode],'w') as f:
            for r in range(nb_runs):
                for i in range(nb_field):
                    if i == 0:
                        stats[r][i] /= 4
                    f.write('{:01d}'.format(int(stats[r][i]))+'\t')
                f.write('\n')

#
if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-n', dest='name'  , default='test',help='run/testcase name')
    parser.add_option('-p', dest='cores' , default=4     ,help='number of cores')
    parser.add_option('-f', dest='fields', default=25    ,help='number of fields')
    parser.add_option('-g', dest='gplot', action='store_true',help='generate gnuplot script')
    parser.add_option('-c', dest='clean', action='store_true',help='clean benchmark data')

    # retrieve and check params
    (param,args) = parser.parse_args()

    # remove data
    #if(param.clean):
    #    clean_data(param)

    # compact data
    #compact_benchmark_data(param)

    # create gnuplot file
    #if(param.gplot):
    #    generate_gnuplot_script(param)


