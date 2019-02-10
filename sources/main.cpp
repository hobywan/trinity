/*
 *                          'main.cpp'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *               Copyright (c) 2016 Hoby Rakotoarivelo.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/* ------------------------------------ */
#include "optparse.h"
#include "tools.h"
#include "mesh.h"
#include "io.h"
#include "metric.h"
#include "refinement.h"
#include "coarsening.h"
#include "smoothing.h"
#include "swapping.h"
#include "matching.h"
#include "random_engine.h"
#include <hwloc.h>
/* ------------------------------------ */
// accessible only in this file
namespace {

int _cores;
int _hw_cores;
int _threads;
int _bucket;
int _rounds;
int _depth;
int _norm;
int _verb;
int _makespan;
int _size[2];

double _target;
double _h_min;
double _h_max;

std::string _arch;
std::string _name;
std::string _input;
std::string _result;
std::string _solut;
}

/* ------------------------------------ */
void parse(int argc, char* argv[]) {

  optparse::Parser parser;
  std::string const mode[] = {"normal", "benchmark", "debug"};
  std::string const arch[] = {"neh", "hsw", "knl", "mac"};
  std::string const papi[] = {"cache", "cycles", "tlb", "branch"};

  // add options
  auto& mode_ = parser.add_option("-m").dest("mode").help("select mode [normal|benchmark|debug]");
  auto& arch_ = parser.add_option("-a").dest("arch").help("cpu architecture [neh|hsw|knl]");
  auto& input = parser.add_option("-i").dest("in").help("initial mesh file (eg. GRID5)");
  auto& rsult = parser.add_option("-o").dest("out").help("result  mesh file (eg. RESUL)");
  auto& solut = parser.add_option("-s").dest("solut").help("solut. field file (eg. SHOCK)");
  auto& cores = parser.add_option("-c").dest("cores").help("number of threads (<= cores)");
  auto& buck_ = parser.add_option("-b").dest("buck").help("vertex bucket capacity [64-256]");
  auto& targ_ = parser.add_option("-t").dest("targ").help("metric field target coef [0.5-1.0]");
  auto& norm_ = parser.add_option("-p").dest("norm").help("metric field L^p norm [0-4]");
  auto& round = parser.add_option("-r").dest("round").help("remeshing rounds [1-5]");
  auto& depth = parser.add_option("-d").dest("depth").help("max refinement/smoothing depth [1-3]");
  auto& verb_ = parser.add_option("-v").dest("verb").help("verbosity level [0-2]");
  auto& papi_ = parser.add_option("-P").dest("papi").help("profile [cache|cycles|tlb|branch]");

  cores.set_default(4).type("int");
  buck_.set_default(64).type("int");
  targ_.set_default(1.0).type("float");
  norm_.set_default(2).type("int");
  depth.set_default(3).type("int");
  round.set_default(2).type("int");
  verb_.set_default(1).type("int");
  input.set_default("GRID4");
  solut.set_default("solut/gauss4");
  rsult.set_default("tests/adap");
  mode_.set_default(mode[0]).choices(mode, mode + 3);
  papi_.set_default(papi[0]).choices(papi, papi + 4);
  arch_.set_default(arch[1]).choices(arch, arch + 4);

  // read stdin
  const optparse::Values& param = parser.parse_args(argc, argv);

  _cores   = std::thread::hardware_concurrency();  // C+11
  _threads = std::max(std::atoi(param["cores"]), 1);
  _name    = trinity::tools::basename(param["solut"]);
  _arch    = param["arch"];
  _bucket  = std::max(std::min(std::atoi(param["buck"]), 256), 64);
  _rounds  = std::max(std::min(std::atoi(param["round"]), 5), 1);
  _depth   = std::max(std::min(std::atoi(param["depth"]), 3), 1);
  _target  = std::max(std::min(std::atof(param["targ"]), 1.), 0.5);   // 0.5*1e4
  _norm    = std::max(std::min(std::atoi(param["norm"]), 4), 0);
  _verb    = std::max(std::min(std::atoi(param["verb"]), 2), 0);
  _input   = "data/" + std::string(param["in"]) + ".mesh";
  _result  = "data/" + std::string(param["out"]) + ".mesh";
  _solut   = "data/" + std::string(param["solut"]) + ".bb";
  _h_min   = EPSILON;
  _h_max   = 0.4;

  if (_threads <= _cores) {
    _hw_cores = 0;
    hwloc_topology_t topology;

    if (hwloc_topology_init(&topology) == 0 and hwloc_topology_load(topology) == 0) {
      _hw_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);
      hwloc_topology_destroy(topology);
      omp_set_num_threads(_threads);
    } else
      trinity::tools::abort('c', "unable to retrieve the number of physical cores", parser);
  } else
    trinity::tools::abort('c', "not enough available logical cores", parser);


  if (trinity::tools::exists(_input)) {
    std::ifstream file(_input, std::ios::in);
    assert(file.good());
    _size[0] = trinity::io::find("Vertices", file);
    _size[1] = trinity::io::find("Triangles", file);
    file.close();
  } else
    trinity::tools::abort('i', "invalid mesh file", parser);

  if (trinity::tools::exists(_solut)) {
    std::ifstream file(_solut, std::ios::in);
    assert(file.good());
    int col[] = {0, 0, 0, 0};
    file >> col[0] >> col[1] >> col[2] >> col[3];
    file.close();
    if (col[2] not_eq _size[0])
      trinity::tools::abort('s', "wrong number of vertices in solut. file", parser);
  } else
    trinity::tools::abort('s', "invalid solut. file", parser);

}

/* ------------------------------------ */
void showDesc() {

  std::string compil;

#ifdef __INTEL_COMPILER
  compil = "Intel compiler";
#else
#ifdef __clang_major__
  compil = "Clang/LLVM";
#else
#ifdef __GNUG__
  compil = "GNU compiler";
#endif
#endif
#endif

  char const symbol = static_cast<char>(1 == _threads ? '\0' : 's');

  std::printf("\n\t= trinity =\n\t");
  std::printf("(c) 2015 H. Rakotoarivelo\n\t");
  std::printf("compiled with %s (%s) on %s, %s\n\t", compil.data(), __VERSION__, __DATE__, __TIME__);
  std::printf("using %d thread%c on %d core%c (%s)\n\n", _threads, symbol, std::min(_hw_cores, _threads),
              symbol, _threads > _hw_cores ? "hyperthreading" : "native");
}

/* ------------------------------------ */
void recap(trinity::Stats* stat) {

  _makespan = 0;
  int rate[] = {0, 0, 0, 0, 0};
  int form[] = {0, 0, 0, 0};

  for (int i = 0; i < 5; ++i) {
    _makespan += stat[i].elap;
    rate[i] = (int) std::floor(stat[i].task / (stat[i].elap * 1e-3));
    form[0] = std::max(form[0], trinity::tools::format(stat[i].elap));  // float
    form[1] = std::max(form[1], trinity::tools::format(rate[i]));
    form[2] = std::max(form[2], trinity::tools::format(stat[i].task));
  }

  std::printf("\n\n= recap: %d rounds, %d threads (%.1f sec)\n", _rounds, _threads, (float) _makespan / 1e3);
  std::printf("=  %2d %% metric (%*.2f s), %*d  calc/sec (%*d done)\n",
              (int) std::floor(stat[0].elap * 100 / _makespan), form[0], (float) stat[0].elap / 1e3,
              form[1], rate[0], form[2], stat[0].task);
  std::printf("=  %2d %% refine (%*.2f s), %*d split/sec (%*d done)\n",
              (int) std::floor(stat[1].elap * 100 / _makespan), form[0], (float) stat[1].elap / 1e3,
              form[1], rate[1], form[2], stat[1].task);
  std::printf("=  %2d %% coarse (%*.2f s), %*d merge/sec (%*d done)\n",
              (int) std::floor(stat[2].elap * 100 / _makespan), form[0], (float) stat[2].elap / 1e3,
              form[1], rate[2], form[2], stat[2].task);
  std::printf("=  %2d %% swap   (%*.2f s), %*d  flip/sec (%*d done)\n",
              (int) std::floor(stat[3].elap * 100 / _makespan), form[0], (float) stat[3].elap / 1e3,
              form[1], rate[3], form[2], stat[3].task);
  std::printf("=  %2d %% smooth (%*.2f s), %*d  move/sec (%*d done)\n\n",
              (int) std::floor(stat[4].elap * 100 / _makespan), form[0], (float) stat[4].elap / 1e3,
              form[1], rate[4], form[2], stat[4].task);

  std::memset(form, 0, sizeof(int) * 4);
  for (int i = 0; i < 4; ++i)
    for (int j : stat[i].step)
      form[i] = std::max(form[i], trinity::tools::format(j));

  std::printf("|%6s%-15s | %5s%-15s | %5s%-15s | %5s%-15s |\n",
              "", "refinement", "", "contraction", "", "swapping", "", "smoothing");

  std::printf("|%3d %% filter \e[0m(%*.2f s)\e[0m |"
              "%3d %% filter \e[0m(%*.2f s)\e[0m |"
              "%3d %% qualit \e[0m(%*.2f s)\e[0m |"
              "%3d %% primal \e[0m(%*.2f s)\e[0m |\n",
              stat[1].step[0] * 100 / stat[1].elap, form[0], (float) stat[1].step[0] / 1e3,
              stat[2].step[0] * 100 / stat[2].elap, form[1], (float) stat[2].step[0] / 1e3,
              stat[3].step[0] * 100 / stat[3].elap, form[2], (float) stat[3].step[0] / 1e3,
              stat[4].step[0] * 100 / stat[4].elap, form[3], (float) stat[4].step[0] / 1e3);
  std::printf("|%3d %% stein  \e[0m(%*.2f s)\e[0m |"
              "%3d %% primal \e[0m(%*.2f s)\e[0m |"
              "%3d %% dual   \e[0m(%*.2f s)\e[0m |"
              "%3d %% color  \e[0m(%*.2f s)\e[0m |\n",
              stat[1].step[1] * 100 / stat[1].elap, form[0], (float) stat[1].step[1] / 1e3,
              stat[2].step[1] * 100 / stat[2].elap, form[1], (float) stat[2].step[1] / 1e3,
              stat[3].step[1] * 100 / stat[3].elap, form[2], (float) stat[3].step[1] / 1e3,
              stat[4].step[1] * 100 / stat[4].elap, form[3], (float) stat[4].step[1] / 1e3);
  std::printf("|%3d %% kernel \e[0m(%*.2f s)\e[0m |"
              "%3d %% indep  \e[0m(%*.2f s)\e[0m |"
              "%3d %% match  \e[0m(%*.2f s)\e[0m |"
              "%3d %% qualit \e[0m(%*.2f s)\e[0m |\n",
              stat[1].step[2] * 100 / stat[1].elap, form[0], (float) stat[1].step[2] / 1e3,
              stat[2].step[2] * 100 / stat[2].elap, form[1], (float) stat[2].step[2] / 1e3,
              stat[3].step[2] * 100 / stat[3].elap, form[2], (float) stat[3].step[2] / 1e3,
              stat[4].step[2] * 100 / stat[4].elap, form[3], (float) stat[4].step[2] / 1e3);
  std::printf("|%3d %% fixes  \e[0m(%*.2f s)\e[0m |"
              "%3d %% kernel \e[0m(%*.2f s)\e[0m |"
              "%3d %% kernel \e[0m(%*.2f s)\e[0m |"
              "%3d %% kernel \e[0m(%*.2f s)\e[0m |\n",
              stat[1].step[3] * 100 / stat[1].elap, form[0], (float) stat[1].step[3] / 1e3,
              stat[2].step[3] * 100 / stat[2].elap, form[1], (float) stat[2].step[3] / 1e3,
              stat[3].step[3] * 100 / stat[3].elap, form[2], (float) stat[3].step[3] / 1e3,
              stat[4].step[3] * 100 / stat[4].elap, form[3], (float) stat[4].step[3] / 1e3);
  std::printf("|%6s%-15s |"
              "%3d %% fixes  \e[0m(%*.2f s)\e[0m |"
              "%3d %% fixes  \e[0m(%*.2f s)\e[0m |"
              "%6s%-15s |\n\n",
              "", "N/A",
              stat[2].step[4] * 100 / stat[2].elap, form[1], (float) stat[2].step[4] / 1e3,
              stat[3].step[4] * 100 / stat[3].elap, form[2], (float) stat[3].step[4] / 1e3,
              "", "N/A");
}

/* ------------------------------------ */
void expor(trinity::Stats* stat) {

  _name.pop_back();
  std::string path;

#ifdef DEFERRED_UPDATES
  std::string const suffix = "def_"+_arch+"_"+_name+".dat";
#else
  std::string const suffix = "perf_" + _arch + "_" + _name + ".dat";
#endif

  for (int i = 1; i < 5; ++i) {
    path = "profile/_" + std::to_string(i) + "/" + suffix;

    std::FILE* file = std::fopen(path.data(), "a");
    std::fprintf(file, "%2d \t%d \t%d \t%d \t%3.2f \t%8d \t%8d \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f\n",
                 _threads,
                 _rounds,
                 _size[0],
                 _size[1],
                 _target,
                 stat[i].task,
                 stat[i].eval,
                 (float) _makespan / 1e3,
                 (float) stat[i].elap / 1e3,
                 (float) stat[i].step[0] / 1e3,
                 (float) stat[i].step[1] / 1e3,
                 (float) stat[i].step[2] / 1e3,
                 (float) stat[i].step[3] / 1e3,
                 (float) stat[i].step[4] / 1e3);
    std::fclose(file);
  }
  std::printf("= '%s' exported\n", suffix.data());
}

/* ------------------------------------ */
int main(int argc, char* argv[]) {

  parse(argc, argv);

  // -----
  showDesc();
  // -----

  trinity::Mesh    mesh(_size, _bucket, _depth, _verb, _rounds);
  trinity::Metrics metric(&mesh, _target, _norm, _h_min, _h_max);
  trinity::Partit  heuris(mesh.getCapaNode(), 8);
  trinity::Refine  refine(&mesh, _depth);
  trinity::Swap    swap(&mesh);
  trinity::Coarse  coarse(&mesh, &heuris);
  trinity::Smooth  smooth(&mesh, &heuris, _depth);
  trinity::Stats   stat[5];

  mesh.loadFrom(_input, _solut);
  metric.computeTensorField(stat);
  metric.clear();

  for (int iter = 0; iter < _rounds; ++iter) {
    refine.run(stat + 1);
    coarse.run(stat + 2);
    swap.run(stat + 3);
    smooth.run(stat + 4);
  }

  recap(stat);
  expor(stat);

  mesh.storeTo(_result);
}
