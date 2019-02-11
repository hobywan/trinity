#!/usr/bin/env python
__author__    = 'Hoby Rakotoarivelo'
__license__   = 'GPL v3'
__version__   = '1.0'
__copyright__ = 'Copyright 2016, The trinity project'

import numpy as np
import scipy.sparse as sps
import scipy.io as io
import matplotlib.pyplot as plt
import sys
import optparse

def parse_matrix_market(matrix, figure):
    
    row = []
    col = []
    data = []
    
    M = io.mmread(matrix)
    
    plt.title('graph sparsity pattern')
    plt.spy(M, precision=0, marker='.', markersize=1)
    plt.savefig(figure)


def parse_scotch_grf(graph, figure):

    row = []
    col = []
    data = []
    
    g = open(graph,'r')
    g.readline()
    line = g.readline()
    graph_size = (line.split())[0]
    print(graph_size)
    g.readline()
    
    for index in range(0, int(graph_size)-1):
        line = g.readline()
        linesplit = line.split()
        degree = linesplit[2]
    
        for icol in range(3, int(degree)+3):
            row.append(index)
            col.append(int(linesplit[icol]))
            data.append(1)
    
    dual = sps.csr_matrix((data,(row,col)))
    #plt.ticklabel_format(axis='x', style='sci', scilimits=(-5,5))
    
    #plt.spy(dual, precision=0, marker='.', markersize=1, color='#4682B4')
    plt.spy(dual, precision=0, marker='.', color='black', markersize=2)
    plt.savefig(figure)
    g.close()

#
if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-i', dest='graph' , help='input grf graph')
    parser.add_option('-o', dest='output', help='output png image')
    (param,args) = parser.parse_args()
    
    parse_scotch_grf(param.graph,param.output)
    

