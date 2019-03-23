/*
 *                          'main.cpp'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *                Copyright 2016, Hoby Rakotoarivelo.
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
#include "trinity.h"
/* -------------------------------------------------------------------------- */
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
