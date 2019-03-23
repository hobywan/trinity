/*
 *                         'parser.h'
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

#pragma once
/* -------------------------------------------------------------------------- */
#include "tools.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
class Parser {

public:
  struct {
    int cores    = 0;
    int hw_cores = 0;
    int threads  = 0;
    int bucket   = 0;
    int rounds   = 0;
    int depth    = 0;
    int norm     = 0;
    int verb     = 0;
    int makespan = 0;
    int size[2]  = {0,0};

    double target = 0.;
    double h_min  = 0.;
    double h_max  = 0.;

    std::string arch   = "";
    std::string name   = "";
    std::string input  = "";
    std::string result = "";
    std::string solut  = "";
  } param;

   Parser() = default;
   Parser(int argc, char* argv[]);
  ~Parser() = default;

  void recap(Stats* stat);
  void dump(Stats* stat);

private:
  void showDesc();
};
/* -------------------------------------------------------------------------- */
}