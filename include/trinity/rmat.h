/*
 *                          'rmat.h'
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
#include "header.h"
#include "timer.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
class RMAT {

  friend class Partit;

public:

  // rule of five
  RMAT();
  RMAT(const RMAT& other) = delete;
  RMAT& operator=(RMAT other) = delete;
  RMAT(RMAT&& other) noexcept = delete;
  RMAT& operator=(RMAT&& other) noexcept = delete;
  ~RMAT();

  // utils
  void reset();
  void load(std::string path);
  void info(std::string name);

private:

  // utils
  void saveChrono();
  int elapsed();

  Graph graph;

  struct {
    int nodes;
    int edges;
    int rounds;
    int error;
    int color;
  } nb;

  struct {
    int max;
    int avg;
  } deg;

  struct { double ratio; } stat;
  struct { Time start; }   time;
};
/* -------------------------------------------------------------------------- */
} // namespace trinity
