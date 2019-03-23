/*
 *                          'parser.cpp'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *                Copyright 2016, Hoby Rakotoarivelo
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/* -------------------------------------------------------------------------- */
#include "trinity/parser.h"
#include "trinity/io.h"

#if HAVE_HWLOC
  #include <hwloc.h>
#endif
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
Parser::Parser(int argc, char* argv[]) {

  optparse::Parser parser;
  std::string const run_mode[] = {"normal", "benchmark", "debug"};
  std::string const papi_mode[] = {"cache", "cycles", "tlb", "branch"};

  // add options
  auto& modes  = parser.add_option("-m").dest("mode");
  auto& archi  = parser.add_option("-a").dest("arch");
  auto& input  = parser.add_option("-i").dest("in");
  auto& result = parser.add_option("-o").dest("out");
  auto& solut  = parser.add_option("-s").dest("solut");
  auto& cores  = parser.add_option("-c").dest("cores");
  auto& bucket = parser.add_option("-b").dest("buck");
  auto& target = parser.add_option("-t").dest("targ");
  auto& norm   = parser.add_option("-p").dest("norm");
  auto& rounds = parser.add_option("-r").dest("round");
  auto& depth  = parser.add_option("-d").dest("depth");
  auto& verbo  = parser.add_option("-v").dest("verb");
  auto& papis  = parser.add_option("-P").dest("papi");

  modes .help("select mode [normal|benchmark|debug]");
  archi .help("cpu architecture [skl|knl|kbl]");
  input .help("initial mesh file");
  result.help("result  mesh file");
  solut .help("solution field .bb file");
  cores .help("number of threads");
  bucket.help("vertex bucket capacity [64-256]");
  target.help("target resolution factor [0.5-1.0]");
  norm  .help("metric field L^p norm [0-4]");
  rounds.help("remeshing rounds [1-5]");
  depth .help("max refinement/smoothing depth [1-3]");
  verbo .help("verbosity level [0-2]");
  papis .help("enable papi [cache|cycles|tlb|branch]");

  modes .set_default(run_mode[0]).choices(run_mode, run_mode + 3);
  archi .set_default("kbl");
  input .set_default(std::string(TRINITY_EXAMP_DIR) + "/mesh/GRID4.mesh");
  solut .set_default(std::string(TRINITY_EXAMP_DIR) + "/solut/shock4.bb");
  result.set_default(std::string(TRINITY_BUILD_DIR) + "/data/adapted.mesh");
  cores .set_default(4).type("int");
  bucket.set_default(64).type("int");
  target.set_default(1.0).type("float");
  norm  .set_default(2).type("int");
  rounds.set_default(8).type("int");
  depth .set_default(3).type("int");
  verbo .set_default(1).type("int");
  papis .set_default(papi_mode[0]).choices(papi_mode, papi_mode + 4);

  // read stdin
  const auto& params = parser.parse_args(argc, argv);

  param.cores   = std::thread::hardware_concurrency();  // C+11
  param.threads = std::max(std::atoi(params["cores"]), 1);
  param.name    = tools::testcase(params["solut"]);
  param.arch    = "intel_" + std::string(params["arch"]);
  param.bucket  = std::max(std::min(std::atoi(params["buck"]), 256), 64);
  param.rounds  = std::max(std::min(std::atoi(params["round"]), 5), 1);
  param.depth   = std::max(std::min(std::atoi(params["depth"]), 3), 1);
  param.target  = std::max(std::min(std::atof(params["targ"]), 1.), 0.5);
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
      if (hwloc_topology_init(&topology) == ok
      and hwloc_topology_load(topology) == ok) {
        param.hw_cores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);
        hwloc_topology_destroy(topology);
        omp_set_num_threads(param.threads);
      } else {
        tools::abort('c', "unable to retrieve number of physical cores", parser);
      }
    #else
      param.hw_cores = param.cores;
    #endif
  } else {
    tools::abort('c', "not enough available logical cores", parser);
  }


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
    if (col[2] != param.size[0]) {
      tools::abort('s', "wrong number of vertices in solut. file", parser);
    }
  } else {
    tools::abort('s', "invalid solut. file", parser);
  }

  // show description
  showDesc();
}

/* -------------------------------------------------------------------------- */
void Parser::recap(Stats* stat) {

  if (stat != nullptr) {

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
                (int) std::floor(stat[0].elap * 100 / param.makespan), form[0],
                (float) stat[0].elap / 1e3,
                form[1], rate[0], form[2], stat[0].task);
    std::printf("=  %2d %% refine (%*.2f s), %*d split/sec (%*d done)\n",
                (int) std::floor(stat[1].elap * 100 / param.makespan), form[0],
                (float) stat[1].elap / 1e3,
                form[1], rate[1], form[2], stat[1].task);
    std::printf("=  %2d %% coarse (%*.2f s), %*d merge/sec (%*d done)\n",
                (int) std::floor(stat[2].elap * 100 / param.makespan), form[0],
                (float) stat[2].elap / 1e3,
                form[1], rate[2], form[2], stat[2].task);
    std::printf("=  %2d %% swap   (%*.2f s), %*d  flip/sec (%*d done)\n",
                (int) std::floor(stat[3].elap * 100 / param.makespan), form[0],
                (float) stat[3].elap / 1e3,
                form[1], rate[3], form[2], stat[3].task);
    std::printf("=  %2d %% smooth (%*.2f s), %*d  move/sec (%*d done)\n\n",
                (int) std::floor(stat[4].elap * 100 / param.makespan), form[0],
                (float) stat[4].elap / 1e3,
                form[1], rate[4], form[2], stat[4].task);

    std::memset(form, 0, sizeof(int) * 4);
    for (int i = 0; i < 4; ++i)
      for (int j : stat[i].step)
        form[i] = std::max(form[i], tools::format(j));

    std::printf("|%6s%-15s | %5s%-15s | %5s%-15s | %5s%-15s |\n",
                "", "refinement",
                "", "contraction",
                "", "swapping",
                "", "smoothing");

    std::printf("|%3d %% filter \e[0m(%*.2f s)\e[0m |"
                "%3d %% filter \e[0m(%*.2f s)\e[0m |"
                "%3d %% qualit \e[0m(%*.2f s)\e[0m |"
                "%3d %% primal \e[0m(%*.2f s)\e[0m |\n",
                stat[1].step[0] * 100 / stat[1].elap,
                form[0], (float) stat[1].step[0] / 1e3,
                stat[2].step[0] * 100 / stat[2].elap,
                form[1], (float) stat[2].step[0] / 1e3,
                stat[3].step[0] * 100 / stat[3].elap,
                form[2], (float) stat[3].step[0] / 1e3,
                stat[4].step[0] * 100 / stat[4].elap,
                form[3], (float) stat[4].step[0] / 1e3);

    std::printf("|%3d %% stein  \e[0m(%*.2f s)\e[0m |"
                "%3d %% primal \e[0m(%*.2f s)\e[0m |"
                "%3d %% dual   \e[0m(%*.2f s)\e[0m |"
                "%3d %% color  \e[0m(%*.2f s)\e[0m |\n",
                stat[1].step[1] * 100 / stat[1].elap,
                form[0], (float) stat[1].step[1] / 1e3,
                stat[2].step[1] * 100 / stat[2].elap,
                form[1], (float) stat[2].step[1] / 1e3,
                stat[3].step[1] * 100 / stat[3].elap,
                form[2], (float) stat[3].step[1] / 1e3,
                stat[4].step[1] * 100 / stat[4].elap,
                form[3], (float) stat[4].step[1] / 1e3);

    std::printf("|%3d %% kernel \e[0m(%*.2f s)\e[0m |"
                "%3d %% indep  \e[0m(%*.2f s)\e[0m |"
                "%3d %% match  \e[0m(%*.2f s)\e[0m |"
                "%3d %% qualit \e[0m(%*.2f s)\e[0m |\n",
                stat[1].step[2] * 100 / stat[1].elap,
                form[0], (float) stat[1].step[2] / 1e3,
                stat[2].step[2] * 100 / stat[2].elap,
                form[1], (float) stat[2].step[2] / 1e3,
                stat[3].step[2] * 100 / stat[3].elap,
                form[2], (float) stat[3].step[2] / 1e3,
                stat[4].step[2] * 100 / stat[4].elap,
                form[3], (float) stat[4].step[2] / 1e3);

    std::printf("|%3d %% fixes  \e[0m(%*.2f s)\e[0m |"
                "%3d %% kernel \e[0m(%*.2f s)\e[0m |"
                "%3d %% kernel \e[0m(%*.2f s)\e[0m |"
                "%3d %% kernel \e[0m(%*.2f s)\e[0m |\n",
                stat[1].step[3] * 100 / stat[1].elap,
                form[0], (float) stat[1].step[3] / 1e3,
                stat[2].step[3] * 100 / stat[2].elap,
                form[1], (float) stat[2].step[3] / 1e3,
                stat[3].step[3] * 100 / stat[3].elap,
                form[2], (float) stat[3].step[3] / 1e3,
                stat[4].step[3] * 100 / stat[4].elap,
                form[3], (float) stat[4].step[3] / 1e3);

    std::printf("|%6s%-15s |"
                "%3d %% fixes  \e[0m(%*.2f s)\e[0m |"
                "%3d %% fixes  \e[0m(%*.2f s)\e[0m |"
                "%6s%-15s |\n\n",
                "", "",
                stat[2].step[4] * 100 / stat[2].elap,
                form[1], (float) stat[2].step[4] / 1e3,
                stat[3].step[4] * 100 / stat[3].elap,
                form[2], (float) stat[3].step[4] / 1e3,
                "", "");
  }
}

/* -------------------------------------------------------------------------- */
void Parser::dump(Stats* stat) {

  if (stat != nullptr) {
    auto dumpKernel = [=](std::string kernel, int i) {
      if (param.verb) {
        std::printf("Exporting ... ");
        std::fflush(stdout);
      }

      auto start = timer::now();

      auto prefix = std::string(TRINITY_BUILD_DIR);
      auto suffix = param.name + "_" + param.arch + "_" + kernel + ".dat";
      auto path   = prefix + "/data/" + suffix;
      auto file   = std::fopen(path.data(), "a");
      std::fprintf(file,
                   "%2d \t%d \t%d \t%d \t%3.2f \t%8d \t%8d \t%.3f"
                   "\t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f\n",
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

      auto output = path.erase(0, prefix.length() + 1).data();
      if (param.verb) {
        std::printf("%-35s \e[32m(%.2f s)\e[0m\n",
                    output,
                    (float) timer::elapsed_ms(start) / 1e3);
      } else {
        std::printf("= '%-35s' exported\n", output);
      }
    };

    dumpKernel("refine", 1);
    dumpKernel("coarse", 2);
    dumpKernel("swap"  , 3);
    dumpKernel("smooth", 4);
  }
}

/* -------------------------------------------------------------------------- */
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

  char const symbol = char(1 == param.threads ? '\0' : 's');

  std::printf("\n\t= trinity =\n\t");
  std::printf("(c) 2016 H. Rakotoarivelo\n\t");
  std::printf("compiled with %s (%s) on %s, %s\n\t",
              compil.data(), __VERSION__, __DATE__, __TIME__);
  std::printf("using %d thread%c on %d core%c (%s)\n\n", param.threads,
              symbol, std::min(param.hw_cores, param.threads), symbol,
              param.threads > param.hw_cores ? "hyperthreading" : "native");
}

/* -------------------------------------------------------------------------- */
}
