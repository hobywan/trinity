#!/usr/bin/env python
__version__ = '1.0'
__author__  = 'H. Rakotoarivelo'

import os
import sys
import time
import optparse
import subprocess

#
def compact_data(param):
    
    prefix = ['hasw_nor','hasw_hyp','nehalem','shock_knl','shock_hyp','shock_llc']
    label = ['hsw-nor','hsw-hyp','nehalem','knl-fla','knl-hyp','knl-llc']
    
    for k in range(3):

        data = param.output+'_'+str(k+1)+'.dat'
        if os.path.exists(data):
            os.remove(data)

        for t in range(6):
            path = 'benchs/'+prefix[t]+'_'+str(k+1)+'.dat'
            print(path)
            if not os.path.exists(path):
                continue

            step  = [0] * 4
            ratio = [0] * 4
            total = 0
            with open(path,'r') as f:
                for line in f:
                    values = line.split()
                    for i in range(4):
                        step[i] += int(values[20+i])
            # reduct
            for i in range(4):
                total += step[i]

            # append to file
            with open(data,'a') as f:
                f.write(label[t])
                sum_line = 0
                for i in range(4):
                    # compute ratio
                    ratio[i] = float(step[i]/total) * 100
                    sum_line += step[i]
                    # write to file
                    f.write('\t{0:.3f}'.format(ratio[i]))
                f.write('\t{:06d}'.format(sum_line)+'\n')

#
#def extract_load_imbalance(param):
#
#    prefix = ['hasw_nor','hasw_hyp','nehalem','shock_knl','shock_hyp','shock_llc']
#    label = ['hsw-nor','hsw-hyp','nehalem','knl-fla','knl-hyp','knl-llc']
#
#    for k in range(3):
#        #
#        data = param.imb+'_'+str(k+1)+'.dat'
#        if os.path.exists(data):
#            os.remove(data)
#        #
#        for t in range(6):
#            path = 'benchs/'+prefix[t]+'_'+str(k+1)+'.dat'
#            print(path)
#            if not os.path.exists(path):
#                continue

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-k', dest='kernel', default='1',help='kernel num. [1-4]')
    parser.add_option('-o', dest='output', default='repart', help='output .dat file')

    # retrieve and check params
    (param,args) = parser.parse_args() 

    compact_data(param)
