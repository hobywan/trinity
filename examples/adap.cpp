/*
 *                       'examples/adap.cpp'
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

#include <trinity.h>

int main(int argc, char* argv[]) {

  int const threads = std::thread::hardware_concurrency();
  int const verbose = 1;
  int const bucket  = 64;
  int const depth   = 3;
  int const rounds  = 8;
  int const max_col = 8;

  auto const norm   = 2;
  auto const target = 1.0;
  auto const h_min  = 1.E-9;
  auto const h_max  = 0.4;

  auto const input  = std::string(SOURCE_DIR) + "/examples/mesh/GRID4.mesh";
  auto const solut  = std::string(SOURCE_DIR) + "/examples/solut/shock4.mesh";
  auto const result = std::string( BUILD_DIR) + "/examples/adapted.mesh";

  assert(trinity::tools::exists(input));
  assert(trinity::tools::exists(solut));

  std::ifstream file(input, std::ios::in);
  assert(file.good());
  int size[] = {
    trinity::io::find("Vertices",  file),
    trinity::io::find("Triangles", file)
  };
  file.close();

  trinity::Mesh    mesh   (size, bucket, depth, verbose, rounds);
  trinity::Metrics metric (&mesh, target, norm, h_min, h_max);
  trinity::Partit  heuris (&mesh, max_col);
  trinity::Refine  refine (&mesh, depth);
  trinity::Swap    swap   (&mesh);
  trinity::Coarse  coarse (&mesh, &heuris);
  trinity::Smooth  smooth (&mesh, &heuris, depth);

  mesh.load(input, solut);
  metric.run();

  for (int iter = 0; iter < rounds; ++iter) {
    refine.run();
    coarse.run();
      swap.run();
    smooth.run();
  }

  mesh.store(result);
}