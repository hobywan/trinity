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
    
    #dual = sps.csr_matrix((data,(row,col)))
    #plt.xlabel('indices')
    #plt.ylabel('Indices')
    plt.title('graph sparsity pattern')
    plt.spy(M, precision=0, marker='.', markersize=1)
    plt.savefig(figure)


def parse_scotch_grf(graph, figure):

    row = []
    col = []
    data = []
    
    g = open(graph,'r')
    #print(f)
    g.readline()
    line = g.readline()
    #print(line)
    numberOfVertices = (line.split())[0]
    print(numberOfVertices)
    g.readline()
    
    for index in range(0, int(numberOfVertices)-1):
        line = g.readline()
        linesplit = line.split()
        numberOfNeighbours = linesplit[2]
    
        for icol in range(3, int(numberOfNeighbours)+3):
            row.append(index)
            col.append(int(linesplit[icol]))
            data.append(1)
    
    dual = sps.csr_matrix((data,(row,col)))
    #plt.xlabel('indices')
    #plt.ylabel('Indices')
    #plt.title('mesh nodes incidence matrix')
    #plt.ticklabel_format(axis='x', style='sci', scilimits=(-5,5))
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(-5,5))
    
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
    

