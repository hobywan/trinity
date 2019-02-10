/* ------------------------------------ */
#include "io.h"
/* ------------------------------------ */
using namespace trinity;

/* -------------------------------- */
int io::find(const std::string key, std::ifstream& file) {

  assert(file.is_open());
  assert(file.good());

  int nb = -1;
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
  file >> nb;
  return nb;
}


/* ------------------------------------ */
void Mesh::loadFrom(const std::string& path, const std::string& solu) {

  if (_verb) {
    std::printf("Parsing mesh/solu ... ");
    std::fflush(stdout);
  }
  auto start = timer::now();

  assert(nb_nodes_);
  assert(nb_elems_);

  std::ifstream input(path, std::ios::in);
  assert(input.is_open());
  assert(input.good());

  int k = 0;
  int nb = 0;
  int tag = 0;
  int v[] = {-1, -1, -1};

  k = 0;
  nb = io::find("Vertices", input);
  assert(nb == nb_nodes_);
  while (k < nb and input >> points_[k * 2] >> points_[k * 2 + 1] >> tag) k++;

  // todo: RequiredVertices
  k = 0;
  nb = io::find("Edges", input);
  while (k < nb and input >> v[0] >> v[1] >> tag) {
    tags_[v[0] - 1] = mask::bound;
    tags_[v[1] - 1] = mask::bound;
    ++k;
  }

  k = 0;
  nb = io::find("Triangles", input);
  assert(nb == nb_elems_);
  while (k < nb and input >> elems_[k * 3] >> elems_[k * 3 + 1] >> elems_[k * 3 + 2] >> tag) k++;

#pragma omp parallel for
  for (int i = 0; i < nb_elems_ * 3; ++i)
    --(elems_[i]);

  k = 0;
  nb = io::find("Corners", input);
  while (k++ < nb and input >> *v)
    tags_[*v - 1] = (mask::bound | mask::corner);

  input.close();

  // -----------
  // parse solution field input
  input.open(solu, std::ios::in);

  k = 0;
  tools::seek_to_line(2, input);
  while (k < nb_nodes_ and input >> solut_[k++]);

  input.close();
  if (_verb) {
    std::printf("done. \e[32m(%.2f s)\e[0m\n", (float) timer::elapsed_ms(start) / 1e3);
    std::printf("Rebuild topology  ... ");
  } else
    std::printf("\r= Preprocess ...  66 %% =");

  std::fflush(stdout);
  start = timer::now();

#pragma omp parallel
  rebuildTopology();

  if (_verb)
    std::printf("done. \e[32m(%.2f s)\e[0m\n", (float) timer::elapsed_ms(start) / 1e3);
  else
    std::printf("\r= Preprocess ... 100 %% =");

  std::fflush(stdout);

  int nb_bounds = 0;

#pragma omp parallel for reduction(+:nb_bounds)
  for (int i = 0; i < nb_nodes_; ++i)
    if (isBoundary(i))
      nb_bounds++;

  int form[] = {tools::format(nb_elems_), tools::format(max_elem)};

  if (_verb == 1)
    std::printf("\n");

  else if (_verb == 2) {
    std::printf("= %*d nodes, %*lu max. \e[0m(x%d)\e[0m, %d isBoundary. \e[0m(%.1f %%)\e[0m\n",
                form[0], nb_nodes_, form[1], max_node, _scale, nb_bounds,
                (float) nb_bounds * 100 / nb_nodes_);
    std::printf("= %*d elems, %*lu max. \e[0m(x%d)\e[0m\n\n",
                form[0], nb_elems_, form[0], max_elem, _scale);
  }
}

/* ------------------------------------ */
void Mesh::storeTo(const std::string& path) const {

  if (_verb) {
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
    for (int i = 0; i < nb_nodes_; ++i)
      if (isActiveNode(i))
        real_nb_nodes++;

#pragma omp for reduction(+:real_nb_elems) nowait
    for (int i = 0; i < nb_elems_; ++i)
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

  int* index = new int[nb_nodes_];
  std::memset(index, -1, nb_nodes_ * sizeof(int));

  std::vector<int> required;
  required.reserve(nb_nodes_ / 4);

  for (int i = 0; i < nb_nodes_; ++i) {
    if (__builtin_expect(isActiveNode(i), 1)) {
      std::sprintf(buffer, "%.8f\t%.8f\t0\n", points_[i * 2], points_[i * 2 + 1]);
      //const double& x = points[i*2];
      //const double& y = points[i*2+1];
      //std::sprintf(buffer, "%.8f\t%.8f\t%.8f\t0\n", x, y, std::exp((-10.)*(x*x + y*y)));
      file << buffer;
      index[i] = ++k;
      if (__builtin_expect(isBoundary(i), 0))
        required.push_back(k);
    }
  }
  file << std::endl;

  file << "Triangles" << std::endl;
  file << real_nb_elems << std::endl;

  for (int i = 0; i < nb_elems_; ++i) {
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

  if (_verb)
    std::printf("'%s' \e[32m(%.2f s)\e[0m\n",
                tools::basename(path).data(), (float) timer::elapsed_ms(start) / 1e3);
  else
    std::printf("= '%s' exported\n", tools::basename(path).data());
}

/* ------------------------------------ */
void Mesh::storePrimalGraph(const std::string& path) const {

  if (_verb) {
    std::printf("Exporting graphs ... ");
    std::fflush(stdout);
  }
  auto start = timer::now();

  std::ofstream file(path, std::ios::out | std::ios::trunc);
  assert(file.is_open());
  assert(file.good());

  int nb_edges = 0;

#pragma omp parallel for reduction(+:nb_edges)
  for (int i = 0; i < nb_nodes_; ++i)
    nb_edges += (int) vicin_[i].size();

  // 1: version number
  // 2: |V| and deg_max
  // 3: begin index=1, vertex label=yes, dart weight=no, vertex weight=no
  file << 0 << std::endl;
  file << nb_nodes_ << "\t" << nb_edges << std::endl;
  file << 0 << "\t" << 101 << std::endl;

  for (int i = 0; i < nb_nodes_; ++i) {
    file << i << "\t" << 1 << "\t" << vicin_[i].size() << "\t";
    for (const int& w : vicin_[i])
      file << w << "\t";
    file << std::endl;
  }

  file << std::endl;
  file.close();

  if (_verb)
    std::printf("done. \e[32m(%d ms)\e[0m\n", timer::elapsed_ms(start));
  else
    std::printf("= '%s' exported\n", tools::basename(path).data());
}

