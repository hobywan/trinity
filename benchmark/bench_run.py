#!/usr/bin/env python
__author__    = 'Hoby Rakotoarivelo'
__license__   = 'GPL v3'
__version__   = '1.0'
__copyright__ = 'Copyright 2016, The trinity project'

import os
import sys
import optparse
import subprocess
import shlex
from .bench_solut_field import *
from .bench_parse_mesh import *

#
def compile_sources(param):
    """Compile source code using CMake
    """

    run_mode = ['benchmark','debug','release']
    opt_prof = ['cache','cycles','tlb','branch']
    compiler = ['intel','gnu','llvm']
    cpu_arch = ['nehalem','haswell','knc','knl']

    # check params
    if param.mode not in run_mode:
        msg = "Invalid mode. Please retry using ["
        for i,item in enumerate(run_mode):
          msg += item
          msg += "|" if i < len(run_mode)-1 else ""
        msg += "]"
        sys.exit(msg)
    #
    if(param.papi not in opt_prof):
        msg = "Invalid profile option. Please retry using ["
        for i,item in enumerate(opt_prof):
          msg += item
          msg += "|" if i < len(opt_prof)-1 else ""
        msg += "]"
        sys.exit(msg)
    #
    if param.cxx not in compiler:
        msg = "Invalid CXX compiler. Please retry using ["
        for i,item in enumerate(compiler):
          msg += item
          msg += "|" if i < len(compiler)-1 else ""
        msg += "]"
        sys.exit(msg)

    #
    if param.arch not in cpu_arch:
        msg = "Invalid CPU architecture. Please retry using ["
        for i,item in enumerate(cpu_arch):
          msg += item
          msg += "|" if i < len(cpu_arch)-1 else ""
        msg += "]"
        sys.exit(msg)

    assert(param.target > 0);

    os.chdir("../build")

    # generate makefile (using cmake)
    subprocess.call(["cmake", "-DCMAKE_BUILD_TYPE="+param.mode.upper(), ".."])

    # process compilation
    proc = subprocess.Popen(shlex.split("make -j4"))
    proc.communicate()
    if(proc.wait() != 0):
        sys.exit(0)

    # stop if compile only mode selected
    if(param.compil):
        sys.exit(0)

#
def adaptive_loop(param,cores):

    # use the same path for the solution field
    mesh   = param.mesh
    solut  = parse.replace_ext(mesh,".bb")
    binary = '../trinity_'+param.cxx
    # adjust num. threads if hyperthreading activated
    num = 4 if param.arch == 'knl' else 2
    threads = cores if not param.hyper else num * cores
    # (!)
    policy = 'compact' if param.hyper else 'scatter'
    os.environ['KMP_AFFINITY'] = 'verbose,granularity=fine,'+policy

    #
    for it_adap in range(0,param.adap):

        # update solution field
        preprocess.compute_solution_field(mesh,solut)

        print(">> "+
              binary+" "+
              str(threads)+" "+
              str(param.maxit)+" "+
              str(param.rmsh)+" "+
              str(it_adap+1)+" "+
              str(param.target)+" "+
              mesh+" "+
              solut+" "+
              param.run+" "+
              param.papi)

        # process remeshing
        if(param.arch == 'knl'):
            assert(param.cxx == 'intel')
            subprocess.call(['numactl',
                '--membind=1',
                binary,
                str(threads),
                str(param.maxit),
                str(param.rmsh),
                str(it_adap+1),
                str(param.target),
                mesh,
                solut,
                param.run,
                param.papi
            ])
        else:
            subprocess.call([binary,
                str(threads),
                str(param.maxit),
                str(param.rmsh),
                str(it_adap+1),
                str(param.target),
                mesh,
                solut,
                param.run,
                param.papi
            ])
        # update couple mesh-solution
        mesh  = "data/"+param.run+"2D.mesh"
        solut = parse.replace_ext(mesh,".bb")

#
def process_benchmark(param):

    # thread binding according to hyperthreading
    policy = 'compact' if param.hyper else 'scatter'

    # set environment variables
    os.environ['OMP_NUM_THREADS'] = str(param.cores)
    os.environ['GOMP_CPU_AFFINITY'] = '0-'+str(param.cores)
    os.environ['KMP_AFFINITY'] = 'verbose,granularity=fine,'+policy

    print('\n')
    print('>> OMP_NUM_THREADS='+os.environ['OMP_NUM_THREADS'])
    print('>> GOMP_CPU_AFFINITY='+os.environ['GOMP_CPU_AFFINITY'])
    print('>> KMP_AFFINITY='+os.environ['KMP_AFFINITY'])
    #
    cur=2
    while cur <= param.cores:
        adaptive_loop(param,cur)
        cur *= 2

#
if __name__ == '__main__':

    # extract params
    parser = optparse.OptionParser()
    parser.add_option('-c', dest="compil", action="store_true", help="compile only sources")
    parser.add_option('-x', dest="cxx"   , default="gnu", help="specify compiler [intel|gnu]")
    parser.add_option('-m', dest="mode"  , default="benchmark", help="select mode [benchmark|debug]")
    parser.add_option('-n', dest="run"   , default="test" , help="run/testcase name", )
    parser.add_option('-p', dest="cores" , default=    4, type="int",help="number of cores")
    parser.add_option('-t', dest="target", default=10000, type="int",help="target num. of vertices")
    parser.add_option('-k', dest="maxit" , default=   15, type="int",help="max iter per kernel")
    parser.add_option('-r', dest="rmsh"  , default=    5, type="int",help="remeshing stages")
    parser.add_option('-a', dest="adap"  , default=    1, type="int",help="adaptive stages")
    parser.add_option('-i', dest="mesh"  , default="data/GRID105.mesh", help="input mesh")
    parser.add_option('-b', dest="bench" , action="store_true", help="benchmark mode")
    parser.add_option('-u', dest="arch"  , default="haswell",
                      help="cpu architecture [nehalem|haswell|knl]")
    parser.add_option('-y', dest="hyper" , action="store_true", default=False,
                      help="activate hyperthreading")
    parser.add_option('-l', dest="papi"  , default="cache",
                      help="hardware profile mode [cache|cycles|tlb|branch]")

    # retrieve params
    (param,args) = parser.parse_args()
    cores = param.cores

    # compile source
    #compile_sources(param)

    # process
    #process_benchmark(param) if param.bench else adaptive_loop(param,cores)

