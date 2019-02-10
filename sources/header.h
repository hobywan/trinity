/* ------------------------------------*/
#pragma once
/* ------------------------------------ */
#include <stack>
#include <set>
#include <map>
#include <algorithm>
/* ------------------------------------ */
#include <omp.h>
#include <thread>
#include <atomic>
#include <memory>
/* ------------------------------------ */
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
/* ------------------------------------ */
#include <cmath>
#include <random>
#include <limits>
/* ------------------------------------ */
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <csignal>
#include <cassert>
#include <unistd.h>
/* ------------------------------------ */
#define EPSILON 1.e-9
/* ------------------------------------ */
// bitmask for multiple tags
namespace trinity { namespace mask {

uint8_t const unset  = (0);
uint8_t const active = (1 << 0);
uint8_t const reeval = (1 << 1);
uint8_t const bound  = (1 << 2);
uint8_t const corner = (1 << 3);
uint8_t const frag   = (1 << 4);   // for mesh
uint8_t const full   = (1 << 5);

}} // namespace trinity::mask

namespace trinity {

double const l_min = 0.707106781186547;
double const l_max = 1.414213562373095;
double const q_norm = 6.928203230275510;
double const q_min = 0.9;
double const max_f = (double) 1 / std::numeric_limits<uint32_t>::max();

struct Stats {
  int eval = 0;
  int task = 0;
  int elap = 0;
  int step[5] = {0, 0, 0, 0, 0};
};

// shortcuts
using Graph = std::vector<std::vector<int> >;
struct Patch { std::vector<int> node; std::vector<int> elem; };

} // namespace trinity

