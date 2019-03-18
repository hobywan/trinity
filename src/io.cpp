/*
 *                          'io.cpp'
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

#include "trinity/io.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------- */
int io::find(const std::string key, std::ifstream& file) {

  assert(file.is_open());
  assert(file.good());

  int count = -1;
  std::string line, token;
  std::stringstream buffer;

  //file.seekg(std::ios::beg);
  while (std::getline(file, line)) {
    // skip comment
    tools::ltrim(line);
    if (line[0] == '#' or line[0] == '%')
      continue;

    // tokenize
    buffer.clear();
    buffer.str(line);
    if (buffer >> token and token == key)
      break;
  }
  file >> count;
  return count;
}

/* --------------------------------------------------------------------------- */
void io::recap(Stats* stat, Parser& parser) {

  int rate[] = {0, 0, 0, 0, 0};
  int form[] = {0, 0, 0, 0};

  for (int i = 0; i < 5; ++i) {
    parser.param.makespan += stat[i].elap;
    rate[i] = (int) std::floor(stat[i].task / (stat[i].elap * 1e-3));
    form[0] = std::max(form[0], tools::format(stat[i].elap));  // float
    form[1] = std::max(form[1], tools::format(rate[i]));
    form[2] = std::max(form[2], tools::format(stat[i].task));
  }

  std::printf("\n\n= recap: %d rounds, %d threads (%.1f sec)\n",
              parser.param.rounds, parser.param.threads, (float) parser.param.makespan / 1e3);
  std::printf("=  %2d %% metric (%*.2f s), %*d  calc/sec (%*d done)\n",
              (int) std::floor(stat[0].elap * 100 / parser.param.makespan), form[0], (float) stat[0].elap / 1e3,
              form[1], rate[0], form[2], stat[0].task);
  std::printf("=  %2d %% refine (%*.2f s), %*d split/sec (%*d done)\n",
              (int) std::floor(stat[1].elap * 100 / parser.param.makespan), form[0], (float) stat[1].elap / 1e3,
              form[1], rate[1], form[2], stat[1].task);
  std::printf("=  %2d %% coarse (%*.2f s), %*d merge/sec (%*d done)\n",
              (int) std::floor(stat[2].elap * 100 / parser.param.makespan), form[0], (float) stat[2].elap / 1e3,
              form[1], rate[2], form[2], stat[2].task);
  std::printf("=  %2d %% swap   (%*.2f s), %*d  flip/sec (%*d done)\n",
              (int) std::floor(stat[3].elap * 100 / parser.param.makespan), form[0], (float) stat[3].elap / 1e3,
              form[1], rate[3], form[2], stat[3].task);
  std::printf("=  %2d %% smooth (%*.2f s), %*d  move/sec (%*d done)\n\n",
              (int) std::floor(stat[4].elap * 100 / parser.param.makespan), form[0], (float) stat[4].elap / 1e3,
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
void io::dump(Stats* stat, Parser const& parser) {

#ifdef DEFERRED_UPDATES
  std::string const suffix = "def_"+_arch+"_"+_name+".dat";
#else
  std::string const suffix = "perf_" + parser.param.arch + "_" + parser.param.name + ".dat";
#endif

  std::string path;

  for (int i = 1; i < 5; ++i) {
    path = "../profile/_" + std::to_string(i) + "/" + suffix;

    std::FILE* file = std::fopen(path.data(), "a");
    std::fprintf(file, "%2d \t%d \t%d \t%d \t%3.2f \t%8d \t%8d \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f \t%.3f\n",
                 parser.param.threads,
                 parser.param.rounds,
                 parser.param.size[0],
                 parser.param.size[1],
                 parser.param.target,
                 stat[i].task,
                 stat[i].eval,
                 (float) parser.param.makespan / 1e3,
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
void Mesh::load(const std::string& path, const std::string& solu) {

  if (param.verb) {
    std::printf("Parsing mesh/solu ... ");
    std::fflush(stdout);
  }
  auto start = timer::now();

  assert(nb.nodes);
  assert(nb.elems);

  std::ifstream input(path, std::ios::in);
  assert(input.is_open());
  assert(input.good());

  int k = 0;
  int count = 0;
  int tag = 0;
  int v[] = {-1, -1, -1};

  k = 0;
  count = io::find("Vertices", input);
  assert(count == nb.nodes);
  while (k < count and input >> geom.points[k * 2] >> geom.points[k * 2 + 1] >> tag) k++;

  // todo: RequiredVertices
  k = 0;
  count = io::find("Edges", input);
  while (k < count and input >> v[0] >> v[1] >> tag) {
    sync.tags[v[0] - 1] = mask::bound;
    sync.tags[v[1] - 1] = mask::bound;
    ++k;
  }

  k = 0;
  count = io::find("Triangles", input);
  assert(count == nb.elems);
  while (k < count
    and input >> topo.elems[k * 3]
              >> topo.elems[k * 3 + 1]
              >> topo.elems[k * 3 + 2] >> tag)
  { k++; }

#pragma omp parallel for
  for (int i = 0; i < nb.elems * 3; ++i)
    --(topo.elems[i]);

  k = 0;
  count = io::find("Corners", input);
  while (k++ < count and input >> *v)
    sync.tags[*v - 1] = (mask::bound | mask::corner);

  input.close();

  // -----------
  // parse solution field input
  input.open(solu, std::ios::in);

  k = 0;
  tools::seekToLine(2, input);
  while (k < nb.nodes and input >> geom.solut[k++]);

  input.close();
  if (param.verb) {
    std::printf("done. \e[32m(%.2f s)\e[0m\n", (float) timer::elapsed_ms(start) / 1e3);
    std::printf("Rebuild topology  ... ");
  } else
    std::printf("\r= Preprocess ...  66 %% =");

  std::fflush(stdout);
  start = timer::now();

#pragma omp parallel
  rebuildTopology();

  if (param.verb)
    std::printf("done. \e[32m(%.2f s)\e[0m\n", (float) timer::elapsed_ms(start) / 1e3);
  else
    std::printf("\r= Preprocess ... 100 %% =");

  std::fflush(stdout);

  int nb_bounds = 0;

#pragma omp parallel for reduction(+:nb_bounds)
  for (int i = 0; i < nb.nodes; ++i)
    if (isBoundary(i))
      nb_bounds++;

  const int form[] = {tools::format(nb.elems), tools::format((int) capa.elem)};

  if (param.verb == 1)
    std::printf("\n");

  else if (param.verb == 2) {
    std::printf("= %*d nodes, %*lu max. \e[0m(x%d)\e[0m, %d bound. \e[0m(%.1f %%)\e[0m\n",
                form[0], nb.nodes, form[1], capa.node, capa.scale, nb_bounds,
                (float) nb_bounds * 100 / nb.nodes);
    std::printf("= %*d elems, %*lu max. \e[0m(x%d)\e[0m\n\n",
                form[0], nb.elems, form[0], capa.elem, capa.scale);
  }
}

/* --------------------------------------------------------------------------- */
void Mesh::store(const std::string& path) const {

  if (param.verb) {
    std::printf("Exporting ... ");
    std::fflush(stdout);
  }
  auto start = timer::now();

  std::ofstream file(path, std::ios::out | std::ios::trunc);
  assert(file.is_open());
  assert(file.good());

  // -- temp --
  int real_nb_nodes = 0;
  int real_nb_elems = 0;

#pragma omp parallel
  {
#pragma omp for reduction(+:real_nb_nodes) nowait
    for (int i = 0; i < nb.nodes; ++i)
      if (isActiveNode(i))
        real_nb_nodes++;

#pragma omp for reduction(+:real_nb_elems) nowait
    for (int i = 0; i < nb.elems; ++i)
      if (isActiveElem(i))
        real_nb_elems++;
  }
  // ----

  file << "MeshVersionFormatted 1" << std::endl;
  file << "Dimension 2" << std::endl;
  file << "Vertices" << std::endl;
  file << real_nb_nodes << std::endl;

  int k = 0;
  char buffer[50];

  int* index = new int[nb.nodes];
  std::memset(index, -1, nb.nodes * sizeof(int));

  std::vector<int> required;
  required.reserve((size_t) nb.nodes / 4);

  for (int i = 0; i < nb.nodes; ++i) {
    if (__builtin_expect(isActiveNode(i), 1)) {
      std::sprintf(buffer, "%.8f\t%.8f\t0\n", geom.points[i * 2], geom.points[i * 2 + 1]);
      file << buffer;
      index[i] = ++k;
      if (__builtin_expect(isBoundary(i), 0))
        required.push_back(k);
    }
  }
  file << std::endl;

  file << "Triangles" << std::endl;
  file << real_nb_elems << std::endl;

  for (int i = 0; i < nb.elems; ++i) {
    const int* n = getElem(i);
    if (*n > -1) {
      std::sprintf(buffer, "%d\t%d\t%d\t0\n", index[n[0]], index[n[1]], index[n[2]]);
      file << buffer;
    }
  }
  file << std::endl;

  file << "RequiredVertices" << std::endl;
  file << required.size() << std::endl;

  for (auto v = required.begin(); v < required.end(); ++v)
    file << *v << std::endl;
  file << std::endl;

  file << "End" << std::endl;
  file.close();

  delete[] index;

  if (param.verb)
    std::printf("'%s' \e[32m(%.2f s)\e[0m\n",
                tools::basename(path).data(), (float) timer::elapsed_ms(start) / 1e3);
  else
    std::printf("= '%s' exported\n", tools::basename(path).data());
}

/* --------------------------------------------------------------------------- */
void Mesh::storePrimalGraph(const std::string& path) const {

  if (param.verb) {
    std::printf("Exporting graphs ... ");
    std::fflush(stdout);
  }
  auto start = timer::now();

  std::ofstream file(path, std::ios::out | std::ios::trunc);
  assert(file.is_open());
  assert(file.good());

  int nb_edges = 0;

#pragma omp parallel for reduction(+:nb_edges)
  for (int i = 0; i < nb.nodes; ++i)
    nb_edges += (int) topo.vicin[i].size();

  // 1: version number
  // 2: |V| and deg_max
  // 3: begin index=1, vertex label=yes, dart weight=no, vertex weight=no
  file << 0 << std::endl;
  file << nb.nodes << "\t" << nb_edges << std::endl;
  file << 0 << "\t" << 101 << std::endl;

  for (int i = 0; i < nb.nodes; ++i) {
    file << i << "\t" << 1 << "\t" << topo.vicin[i].size() << "\t";
    for (const int& w : topo.vicin[i])
      file << w << "\t";
    file << std::endl;
  }

  file << std::endl;
  file.close();

  if (param.verb)
    std::printf("done. \e[32m(%d ms)\e[0m\n", timer::elapsed_ms(start));
  else
    std::printf("= '%s' exported\n", tools::basename(path).data());
}
} // namespace trinity
