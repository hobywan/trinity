<img src="docs/figures/logo.png" alt="logo" width="210">

**trinity** is a C++ library and command-line tool for [anisotropic mesh adaptation](https://pdfs.semanticscholar.org/3246/d43a28559273446811da8acc705c37198926.pdf).  
It is targetted to [non-uniform memory access](https://en.wikipedia.org/wiki/Non-uniform_memory_access) multicore and [manycore](https://en.wikipedia.org/wiki/Manycore_processor) processors.  
It was primarly designed for **performance** and hence for [HPC](https://en.wikipedia.org/wiki/Parallel_computing) applications.  
It is intended to be involved within a numerical simulation loop.  

<img src="docs/figures/adaptive_loop.png" alt="adaptive-loop" width="390">

[![Build Status](https://travis-ci.com/hobywan/trinity.svg?branch=master)](https://travis-ci.com/hobywan/trinity)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2ae6bd595ce54105b445e81e2d132eb8)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hobywan/trinity&amp;utm_campaign=Badge_Grade)
[![license](https://img.shields.io/badge/license-apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)

###### Table of contents

- [Build and use](#build-use)
- [Features](#features)
- [Profile and deploy](#benchmarks)
- [How to contribute](#license)

----
### Build and use <a name="build-use"></a>
###### Building the library <a name="build"></a>
[![Build Status](https://travis-ci.com/hobywan/trinity.svg?branch=master)](https://travis-ci.com/hobywan/trinity)

**trinity** is completely standalone.  
It can be built on Linux or macOS using [CMake](https://cmake.org).  
It only requires a [C++14](https://isocpp.org/wiki/faq/cpp14-language) compiler endowed with [OpenMP](https://www.openmp.org).  
It can build [medit](https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf) to render meshes but it is optional though.  
It supports [hwloc](https://www.open-mpi.org/projects/hwloc/) to retrieve and print more information on the host machine.  

``` bash
mkdir build                                          # out-of-source build recommended
cd build                                             #
cmake ..                                             # see build options below
make -j4                                             # use multiple jobs for compilation
make install                                         # optional, can use a prefix
```
| Option           | Description                                                                                      | Default |   
|------------------|--------------------------------------------------------------------------------------------------|---------|   
| `Build_Medit`    | Build [medit](https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf) mesh renderer         | `ON`    |   
| `Build_GTest`    | Build [googletest](https://github.com/google/googletest) for *future* unit tests                 | `OFF`    |   
| `Build_Main`     | Build the command-line tool                                                                      | `ON`    |
| `Build_Examples` | Build provided examples                                                                          | `ON`    |   
| `Use_Deferred`   | Use deferred topology updates scheme in [pragmatic](https://github.com/meshadaptation/pragmatic) | `OFF`   |  

###### Linking to your project 
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2ae6bd595ce54105b445e81e2d132eb8)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=hobywan/trinity&amp;utm_campaign=Badge_Grade)

**trinity** is exported as a package.  
To enable the library, update your CMakeLists.txt with:

``` cmake
find_package(trinity)                                # in build or install trees
target_link_libraries(target PRIVATE trinity)        # replace 'target'
```
And then include `trinity.h` in your application.  
Please take a look at the [examples](examples/) folder for basic usage.

###### Using the tool <a name="using"></a>
The list of command arguments is given by the `-h` option.

``` console
host:~$ bin/trinity -h
Usage: trinity [options]

Options:
  -h, --help            show this help message and exit
  -m CHOICE             select mode [release|benchmark|debug]
  -a CHOICE             cpu architecture [skl|knl|kbl]
  -i STRING             initial mesh file
  -o STRING             result mesh file
  -s STRING             solution field .bb file
  -c INT                number of threads
  -b INT                vertex bucket capacity [64-256]
  -t FLOAT              target resolution factor [0.5-1.0]
  -p INT                metric field L^p norm [0-4]
  -r INT                remeshing rounds [1-5]
  -d INT                max refinement/smoothing depth [1-3]
  -v INT                verbosity level [0-2]
  -P CHOICE             enable papi [cache|cycles|tlb|branch]
```
> For now, only `.mesh` files used in [medit](https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf) are supported.

###### Setting thread-core affinity <a name="cpu-affinity"></a>

For performance reasons, I recommend to explicitly set [thread-core affinity](https://eli.thegreenplace.net/2016/c11-threads-affinity-and-hyperthreading/) before any run.  
Indeed, threads should be statically bound to cores to prevent the OS from migrating them.  
Besides, [simultaneous multithreading](https://en.wikipedia.org/wiki/Simultaneous_multithreading) (or [hyperthreading](https://en.wikipedia.org/wiki/Hyper-threading) on Intel) should be:

- enabled to ease memory latency penalties especially on Intel KNL.
- disabled to reduce shared caches saturation on faster nodes.  

It can be done by setting some environment variables:

```bash
export OMP_PLACES=[cores|threads] OMP_PROC_BIND=close  # with GNU or clang/LLVM
export KMP_AFFINITY=granularity=[core|fine],compact    # with Intel compiler  
```

----
### Features <a name="features"></a>

###### Overview
**trinity**  aims to reduce and equidistribute the interpolation error of a computed physical field **_u_** on a triangulated  
**planar** domain **M** by adapting its discretization with respect to a target number of points **_n_**.  
Basically, it takes (**_u_**, **M**, **_n_**) and outputs a mesh adapted to the variation of the gradient of **_u_** on **M** using **_n_** points.  
It uses [metric tensors](https://en.wikipedia.org/wiki/Metric_tensor) to encode the desired point distribution with respect to the estimated error.  

<img src="docs/figures/principle.png" alt="principle" width="820">

It enables to **resample** and **regularize** a planar triangular mesh **M**.  
It aims to reduce and equidistribute the error of a solution field **_u_** on **M** using **_n_** points.  
For that, it uses five kernels:

-	[metric recover](sources/metric.h): compute a tensor field which encodes desired point density.
-	[refinement](sources/refinement.h): add points on areas where the error of **_u_** is large.
-	[coarsening](sources/coarsening.h): remove points on areas where the error of **_u_** is small.
-   [swapping](sources/swapping.h): flip edges to locally improve cell quality.
-   [smoothing](sources/smoothing.h): relocate points to locally improve cell qualities.

###### Error estimate <a name="error-estimate"></a>
**trinity** uses [metric tensors](https://en.wikipedia.org/wiki/Metric_tensor) to link the error of **_u_** with mesh points distribution.  
A tensor encodes the desired edge length incident to a point, which may be [direction-dependent](https://en.wikipedia.org/wiki/Anisotropy).  
**trinity** enables to tune the sensitivity of the error estimate according to the simulation needs.  
For that, it provides a multi-scale estimate in [**L^_p_** norm](https://en.wikipedia.org/wiki/Lp_space):

-    a small **_p_** will distribute points to capture all scales of the error of **_u_**.
-    a large **_p_** will distribute points mainly on large variation areas (shocks).

It actually implements the [continuous metric](https://link.springer.com/chapter/10.1007/978-3-540-34958-7_12) defined in:
>📄 Fréderic Alauzet, Adrien Loseille, Alain Dervieux and Pascal Frey (2006).  
>["Multi-Dimensional Continuous Metric for Mesh Adaptation"](https://link.springer.com/chapter/10.1007/978-3-540-34958-7_12).  
>In _proceedings of the 15th International Meshing Roundtable_, pp 191-214, Springer Berlin.

To obtain a good mesh, it needs an accurate metric tensor field.  
The latter rely on the computation of the variations of the gradient of **_u_**.  
It is given by its local [hessian matrices](https://en.wikipedia.org/wiki/Hessian_matrix).  
It is computed in **trinity** through a [L^2 projection](https://doi.org/10.1002/nme.2036).  

<img src="docs/figures/multiscale_meshes.png" alt="multiscale_meshes.png" width="790">

###### Fine-grained parallelism <a name="parallelism"></a>
**trinity** enables intra-node parallelism by [multithreading](https://en.wikipedia.org/wiki/Multithreading_(computer_architecture)).  
It relies on a [fork-join](https://en.wikipedia.org/wiki/Fork–join_model) model through OpenMP.  
All kernels are structured into synchronous stages.  
A stage consists of local computation, a [reduction](https://en.wikipedia.org/wiki/Reduce_(parallel_pattern)) in shared-memory, and a [barrier](https://en.wikipedia.org/wiki/Barrier_(computer_science)).  

<table>
  <tr><td><img src="docs/figures/algo_structure.png" alt="algo_structure" width="600"></td></tr>
</table>

>It does **not** rely on [domain partitioning](http://www.cs.cmu.edu/~quake/sc96/node5.html) unlike coarse-grained parallel remeshers.  
>It does **not** rely on [task parallelism](https://en.wikipedia.org/wiki/Task_parallelism) and runtime capabilities such as [Cilk](http://supertech.csail.mit.edu/papers/PPoPP95.pdf), [TBB](https://software.intel.com/en-us/intel-tbb) or [StarPU](http://starpu.gforge.inria.fr) neither.

In fact [manycore](https://en.wikipedia.org/wiki/Manycore_processor) machines have plenty of slow cores with small caches.  
To scale up, one needs plenty of very thin and local tasks to keep them busy.  
In **trinity**, remesh kernels are expressed into a [graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)), except for refinement.  
Runnable tasks are then extracted using multithreaded heuristics:  

-    maximal [stable set](https://en.wikipedia.org/wiki/Maximal_independent_set) for coarsening
-    maximal [matching](https://en.wikipedia.org/wiki/Matching_(graph_theory)) for swapping
-    maximal [coloring](https://en.wikipedia.org/wiki/Graph_coloring) for smoothing

<img src="docs/figures/graph_matching.png" alt="graph_matching.png" width="450">

 **trinity** fixes incidence data only at the end of a round of any kernel.  
It uses an explicit synchronization scheme to fix them.  
It relies on the use of low-level [atomic primitives](https://fr.cppreference.com/w/cpp/atomic).  
It was designed to minimize data movement penalties, especially on [NUMA](https://en.wikipedia.org/wiki/Non-uniform_memory_access) cases.  
For further details, please take a look at:

>📄 Hoby Rakotoarivelo, Franck Ledoux, Franck Pommereau and Nicolas Le-Goff (2017).  
>["Scalable fine-grained metric-based remeshing algorithm for manycore/NUMA architectures"](https://doi.org/10.1007/978-3-319-64203-1_43).  
>In _proceedings of 23rd International European Conference on Parallel and Distributed Computing_, Springer.

----
### Benchmark <a name="benchmarks"></a>
###### Profiling <a name="profiling"></a>
**trinity** is natively instrumented.  
It prints the runtime stats with three verbosity level.  
Here is an output example with the medium level.  

<img src="docs/figures/screenshot.png" alt="screenshot" width="650">

>Stats are exported as tab-separated values and can be easily plotted with [gnuplot](http://www.gnuplot.info) or [matplotlib](https://matplotlib.org).  
>You can use [wrappi](https://github.com/hobywan/wrappi) to profile oncore events such as cycles, caches misses, branch predictions.

###### Deployment on a cluster <a name="deployment"></a>
Preparing a benchmark campaign can be tedious 😩.  
I included some python scripts to help setting it up on a node, enabling to:
<!--(https://blogs.cisco.com/performance/process-and-memory-affinity-why-do-you-care)-->

-    compute a synthetic solution field.
-    rebuild sources and set thread-core affinity.
-    set memory affinity through [`numactl`](https://linux.die.net/man/8/numactl), which is useful on a [Intel KNL](https://colfaxresearch.com/knl-numa/) node.
-    compact profiling data and generate [gnuplot](http://www.gnuplot.info) script for plots.
-    profile memory bandwith of the host machine using [STREAM](https://www.cs.virginia.edu/stream/).
-    plot sparsity pattern of mesh incidence graph.

>They are somewhat outdated, so adapt them to your needs.

----
<a name="license">
  <img src="docs/figures/logo.png" alt="logo" width="210">
</a>

###### Copyright 2016, Hoby Rakotoarivelo.

[![license](https://img.shields.io/badge/license-apache-blue.svg)](https://www.apache.org/licenses/LICENSE-2.0)
  
**trinity** is free and intended for research purposes.  
It was written during my doctorate, so improvements are welcome.  
To get involved, you can:

-    report bugs or request features by submitting an [issue](https://github.com/hobywan/trinity/issues).
-    submit code contributions using feature branches and [pull requests](https://github.com/hobywan/trinity/pulls).

Enjoy! 😊