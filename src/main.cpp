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

/* --------------------------------------------------------------------------- */
#include "trinity.h"

/* --------------------------------------------------------------------------- */
int main(int argc, char* argv[]) {

  trinity::Parser parser(argc, argv);

  auto size   = parser.param.size;
  auto bucket = parser.param.bucket;
  auto depth  = parser.param.depth;
  auto verb   = parser.param.verb;
  auto rounds = parser.param.rounds;
  auto target = parser.param.target;
  auto norm   = parser.param.norm;
  auto h_min  = parser.param.h_min;
  auto h_max  = parser.param.h_max;
  auto input  = parser.param.input;
  auto solut  = parser.param.solut;
  auto result = parser.param.result;

  trinity::Mesh    mesh  (size, bucket, depth, verb, rounds);
  trinity::Metrics metric(&mesh, target, norm, h_min, h_max);
  trinity::Partit  heuris(&mesh, 8);
  trinity::Refine  refine(&mesh, depth);
  trinity::Swap    swap  (&mesh);
  trinity::Coarse  coarse(&mesh, &heuris);
  trinity::Smooth  smooth(&mesh, &heuris, depth);
  trinity::Stats   stat[5];

  mesh.load(input, solut);
  metric.run(stat);

  for (int iter = 0; iter < parser.param.rounds; ++iter) {
    refine.run(stat + 1);
    coarse.run(stat + 2);
      swap.run(stat + 3);
    smooth.run(stat + 4);
  }

  parser.recap(stat);
  parser.dump(stat);
  mesh.store(result);
}
