/*
 *                          'parser.cpp'
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

/* ------------------------------------*/
#include "trinity/parser.h"
#include "trinity/io.h"
#include <hwloc.h>
/* ------------------------------------*/
namespace trinity {
/* ------------------------------------*/
Parser::Parser(int argc, char* argv[]) {

  optparse::Parser parser;
  std::string const mode[] = {"normal", "benchmark", "debug"};
  std::string const papi[] = {"cache", "cycles", "tlb", "branch"};

  // add options
  auto& modes = parser.add_option("-m").dest("mode") .help("select mode [normal|benchmark|debug]");
  auto& arch  = parser.add_option("-a").dest("arch") .help("cpu architecture [skl|knl|kbl]");
  auto& input = parser.add_option("-i").dest("in")   .help("initial mesh file");
  auto& rsult = parser.add_option("-o").dest("out")  .help("result  mesh file");
  auto& solut = parser.add_option("-s").dest("solut").help("solution field .bb file");
  auto& cores = parser.add_option("-c").dest("cores").help("number of threads");
  auto& buck  = parser.add_option("-b").dest("buck") .help("vertex bucket capacity [64-256]");
  auto& targ  = parser.add_option("-t").dest("targ") .help("target resolution factor [0.5-1.0]");
  auto& norm  = parser.add_option("-p").dest("norm") .help("metric field L^p norm [0-4]");
  auto& round = parser.add_option("-r").dest("round").help("remeshing rounds [1-5]");
  auto& depth = parser.add_option("-d").dest("depth").help("max refinement/smoothing depth [1-3]");
  auto& verb  = parser.add_option("-v").dest("verb") .help("verbosity level [0-2]");
  auto& papis = parser.add_option("-P").dest("papi") .help("enable papi [cache|cycles|tlb|branch]");

  cores.set_default(4).type("int");
  buck .set_default(64).type("int");
  targ .set_default(1.0).type("float");
  norm .set_default(2).type("int");
  depth.set_default(3).type("int");
  round.set_default(8).type("int");
  verb .set_default(1).type("int");
  input.set_default(std::string(DEFAULT_INPUT_DIR) + "/mesh/GRID4.mesh");
  solut.set_default(std::string(DEFAULT_INPUT_DIR) + "/solut/shock4.bb");
  rsult.set_default(std::string(DEFAULT_BUILD_DIR) + "/data/adapted.mesh");
  modes.set_default(mode[0]).choices(mode, mode + 3);
  papis.set_default(papi[0]).choices(papi, papi + 4);
  arch .set_default("kbl");

  // read stdin
  const auto& params = parser.parse_args(argc, argv);

  param.cores   = std::thread::hardware_concurrency();  // C+11
  param.threads = std::max(std::atoi(params["cores"]), 1);
  param.name    = tools::basename(params["solut"]);
  param.arch    = params["arch"];
  param.bucket  = std::max(std::min(std::atoi(params["buck"]), 256), 64);
  param.rounds  = std::max(std::min(std::atoi(params["round"]), 5), 1);
  param.depth   = std::max(std::min(std::atoi(params["depth"]), 3), 1);
  param.target  = std::max(std::min(std::atof(params["targ"]), 1.), 0.5);   // 0.5*1e4
  param.norm    = std::max(std::min(std::atoi(params["norm"]), 4), 0);
  param.verb    = std::max(std::min(std::atoi(params["verb"]), 2), 0);
  param.input   = std::string(params["in"]);
  param.result  = std::string(params["out"]);
  param.solut   = std::string(params["solut"]);
  param.h_min   = EPSILON;
  param.h_max   = 0.4;

  if (param.threads <= param.cores) {
    param.hw_cores = 0;
    #if HAVE_HWLOC
      hwloc_topology_t topology;
      int const ok = 0;
      if (hwloc_topology_init(&topology) == ok and hwloc_topology_load(topology) == ok) {
        param.hw_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);
        hwloc_topology_destroy(topology);
        omp_set_num_threads(param.threads);
      } else {
        tools::abort('c', "unable to retrieve the number of physical cores", parser);
      }
    #else
      param.hw_cores = param.cores;
    #endif
  } else
    tools::abort('c', "not enough available logical cores", parser);


  if (tools::exists(param.input)) {
    std::ifstream file(param.input, std::ios::in);
    assert(file.good());
    param.size[0] = io::find("Vertices", file);
    param.size[1] = io::find("Triangles", file);
    file.close();
  } else {
    std::printf("input mesh file: %s\n", param.input.data());
    tools::abort('i', "invalid mesh file", parser);
  }


  if (tools::exists(param.solut)) {
    std::ifstream file(param.solut, std::ios::in);
    assert(file.good());
    int col[] = {0, 0, 0, 0};
    file >> col[0] >> col[1] >> col[2] >> col[3];
    file.close();
    if (col[2] not_eq param.size[0])
      tools::abort('s', "wrong number of vertices in solut. file", parser);
  } else
    tools::abort('s', "invalid solut. file", parser);

  // show description
  showDesc();
}

/* --------------------------------------------------------------------------- */
void Parser::recap(Stats* stat) {

  int rate[] = {0, 0, 0, 0, 0};
  int form[] = {0, 0, 0, 0};

  for (int i = 0; i < 5; ++i) {
    param.makespan += stat[i].elap;
    rate[i] = (int) std::floor(stat[i].task / (stat[i].elap * 1e-3));
    form[0] = std::max(form[0], tools::format(stat[i].elap));  // float
    form[1] = std::max(form[1], tools::format(rate[i]));
    form[2] = std::max(form[2], tools::format(stat[i].task));
  }

  std::printf("\n\n= recap: %d rounds, %d threads (%.1f sec)\n",
              param.rounds, param.threads, (float) param.makespan / 1e3);
  std::printf("=  %2d %% metric (%*.2f s), %*d  calc/sec (%*d done)\n",
              (int) std::floor(stat[0].elap * 100 / param.makespan), form[0], (float) stat[0].elap / 1e3,
              form[1], rate[0], form[2], stat[0].task);
  std::printf("=  %2d %% refine (%*.2f s), %*d split/sec (%*d done)\n",
              (int) std::floor(stat[1].elap * 100 / param.makespan), form[0], (float) stat[1].elap / 1e3,
              form[1], rate[1], form[2], stat[1].task);
  std::printf("=  %2d %% coarse (%*.2f s), %*d merge/sec (%*d done)\n",
              (int) std::floor(stat[2].elap * 100 / param.makespan), form[0], (float) stat[2].elap / 1e3,
              form[1], rate[2], form[2], stat[2].task);
  std::printf("=  %2d %% swap   (%*.2f s), %*d  flip/sec (%*d done)\n",
              (int) std::floor(stat[3].elap * 100 / param.makespan), form[0], (float) stat[3].elap / 1e3,
              form[1], rate[3], form[2], stat[3].task);
  std::printf("=  %2d %% smooth (%*.2f s), %*d  move/sec (%*d done)\n\n",
              (int) std::floor(stat[4].elap * 100 / param.makespan), form[0], (float) stat[4].elap / 1e3,
              form[1], rate[4], form[2], stat[4].task);

  std::memset(form, 0, sizeof(int) * 4);
  for (int i = 0; i < 4; ++i)
    for (int j : stat[i].step)
      form[i] = std::max(form[i], tools::format(j));

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
              "", "",
              stat[2].step[4] * 100 / stat[2].elap, form[1], (float) stat[2].step[4] / 1e3,
              stat[3].step[4] * 100 / stat[3].elap, form[2], (float) stat[3].step[4] / 1e3,
              "", "");
}

/* --------------------------------------------------------------------------- */
void Parser::dump(Stats* stat) {

#ifdef DEFERRED_UPDATES
  std::string const suffix = "def_"+_arch+"_"+_name+".dat";
#else
  std::string const suffix = "perf_" + param.arch + "_" + param.name + ".dat";
#endif

  std::string path;

  for (int i = 1; i < 5; ++i) {
    path = "../profile/_" + std::to_string(i) + "/" + suffix;

    std::FILE* file = std::fopen(path.data(), "a");
    std::fprintf(file, "%2d \t%d \t%d \t%d \t%3.2f \t%8d \t%8d \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f\n",
                 param.threads,
                 param.rounds,
                 param.size[0],
                 param.size[1],
                 param.target,
                 stat[i].task,
                 stat[i].eval,
                 (float) param.makespan / 1e3,
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

/* --------------------------------------------------------------------------- */
void Parser::showDesc() {

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

  char const symbol = static_cast<char>(1 == param.threads ? '\0' : 's');

  std::printf("\n\t= trinity =\n\t");
  std::printf("(c) 2016 H. Rakotoarivelo\n\t");
  std::printf("compiled with %s (%s) on %s, %s\n\t", compil.data(), __VERSION__, __DATE__, __TIME__);
  std::printf("using %d thread%c on %d core%c (%s)\n\n",
              param.threads, symbol, std::min(param.hw_cores, param.threads),
              symbol, param.threads > param.hw_cores ? "hyperthreading" : "native");
}

/* --------------------------------------------------------------------------- */
}