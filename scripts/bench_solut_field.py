#!/usr/bin/env python
__author__    = 'Hoby Rakotoarivelo'
__license__   = 'GPL v3'
__version__   = '1.0'
__copyright__ = 'Copyright 2016, The trinity project'

import os
import optparse
import subprocess
import math
import bench_parse_mesh

#
def convert_3D_to_2D(param):
#    parser.add_option('-b', dest='beg'   , help='Line range lower bound')
#    parser.add_option('-e', dest='end'   , help='Line range upper bound')
    x = int(param.beg)
    y = int(param.end)

    # verif params
    assert(os.path.exists(param.input))
    assert(x >= 0 and y >= x)

    cmd_eol = "wc -l < "+param.input
    eol = os.popen(cmd_eol).read().split()[0]
    print(">> EOL = "+eol)

#    cmd_beg = "sed -n '1,"+str(x-1)+" p' "+param.input+" > toto.mesh"#+param.output
#    print(cmd_beg)

    # copy header
    with open(param.output,'w') as f:
        f.write('MeshVersionFormatted 1\n')
        f.write('Dimension\n2\n')
        f.write('Vertices\n'+str(y-x+1)+'\n')

    os.system("cat "+param.output)

    # replace pattern
    cmd_cut  = "awk 'NR>="+param.beg+" && NR<="+param.end+"' \""+param.input+"\" | "
    cmd_cut += "awk {'print $1,$2,$4'} >> "+param.output
    print('>> '+cmd_cut)
    os.system(cmd_cut)

    # append the rest
    cmd_app  = "awk 'NR>="+str(y+1)+" && NR<="+eol+"' \""+param.input+"\" >> "+param.output
    print('>> '+cmd_app)
    os.system(cmd_app)



#
def gauss(x,y):
    assert(type(x) is float)
    assert(type(y) is float)
    return math.exp((-10.)*(x*x + y*y))

#
def shock(x,y):
    assert(type(x) is float)
    assert(type(y) is float)
    return 0.1 * math.sin(50.*x) + math.atan(0.1/(math.sin(5*y) - 2*x))

#
def waves(x,y):
    assert(type(x) is float)
    assert(type(y) is float)

    xy   = x * y
    coef = math.pi/50.
    fact = math.sin(xy * 50.)

    if(xy <= (-1.) * coef):
        return 0.1 * fact
    elif(xy <= 2.0 * coef):
        return fact
    else:
        return 0.1 * fact

#
_fields = {
    'shock' : shock,
    'waves' : waves,
    'gauss' : gauss
}

#
def compute_solution_field(param):

    global _fields;

    print(">> Computing solution field ... ")
    #
    mesh = _parse.parse_medit(param.mesh)
    nb_nodes = len(mesh['nodes'])
    assert(nb_nodes > 0)

    # compute and store solution field
    with open(param.solu, 'w') as f:
        f.write("2 1 {:d} 2\n".format(nb_nodes))
        for i in mesh['nodes']:
            x,y = mesh['nodes'][i][:2] # from index 0 to 1
            u = _fields[param.field](x,y)
            f.write("{:.10f}\n".format(u)) # 10 decimales
    #
    print(param.solu)

#
if __name__ == '__main__':

    parser = optparse.OptionParser()
    parser.add_option('-i', dest='mesh' , help='input mesh')
    parser.add_option('-s', dest='solu' , help='solut file')
    parser.add_option('-f', dest='field' , help='field [shock|gauss|waves]')
    # retrieve and check params
    (param,args) = parser.parse_args()

    compute_solution_field(param)

