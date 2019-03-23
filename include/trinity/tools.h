/*
 *                          'tools.h'
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
#include "optparse.h"
/* -------------------------------------------------------------------------- */
namespace trinity { namespace tools {
/* -------------------------------------------------------------------------- */
// hash provided key
uint32_t hash(const uint32_t id);

// get number of digits to print 'num'
int format(int num);

// print a separator line
void separator();

// trim left blank space
void ltrim(std::string& line);

// get end file name of a given path
std::string basename(const std::string& s);

// check if a file exists
bool exists(const std::string& path);

// abort if given command option is not supported
void abort(char option, const char* msg, const optparse::Parser& parser);

// check equality of two c-strings
bool equals(const char* s1, const char* s2);

// retrieve extension of a file
std::string getExt(const char* path);

// seek line pointer to a given line number
std::ifstream& seekToLine(int nb, std::ifstream& file);

// check if provided arg is a digit
bool isDigit(const char* arg);

// get file name without extension
std::string rootOf(const std::string& path);

// get testcase name from path
std::string testcase(const std::string& path);

// replace a file extension by another
std::string replaceExt(const std::string& fname, const std::string& ext);

// display elapsed time
void showElapsed(Time& tic, const char* msg, int step);

/* -------------------------------------------------------------------------- */
template<typename type_t>
void display(const std::vector<type_t>& list) {

  std::stringstream buffer;
  buffer << "[";
  for (auto it = list.begin(); it != list.end(); ++it) {
    buffer << *it;
    if (it + 1 != list.end()) { buffer << ","; }
  }
  buffer << "]";
  std::printf("%s\n", buffer.str().data());
}

/* -------------------------------------------------------------------------- */
template<typename type_t>
void erase(type_t needle, std::vector<type_t>& list) {
  auto found = std::find(list.begin(), list.end(), needle);
  assert(found != list.end());
  std::swap(*found, list.back());
  list.pop_back();
}

/* -------------------------------------------------------------------------- */
}} // namespace trinity::tools
