/*
 *                          'io.h'
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

#pragma once
/* ------------------------------------*/
#include "tools.h"
/* ------------------------------------*/
namespace trinity {
/* ------------------------------------*/
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
/* ------------------------------------*/
}