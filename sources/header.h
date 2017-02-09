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
#include <math.h>
#include <cmath>
#include <random>
#include <limits>
/* ------------------------------------ */
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <csignal>
#include <assert.h>
#include <unistd.h>
/* ------------------------------------ */
#define EPSILON 1.e-9
/* ------------------------------------ */
namespace trigen {

  double const l_min  = 0.707106781186547;
  double const l_max  = 1.414213562373095;
  double const q_norm = 6.928203230275510;
  double const q_min  = 0.9;
  double const max_f  = (double) 1/std::numeric_limits<uint32_t>::max();

  // bitmask for multiple tags
  namespace mask {
    uint8_t const unset  = (0);
    uint8_t const active = (1 << 0);
    uint8_t const reeval = (1 << 1);
    uint8_t const bound  = (1 << 2);
    uint8_t const corner = (1 << 3);
    uint8_t const frag   = (1 << 4);   // for mesh
    uint8_t const full   = (1 << 5);
  }

  struct stats_t {
    int eval = 0;
    int task = 0;
    int elap = 0;
    int step[5] = {0,0,0,0,0};
  };    

  // shortcuts
  using  graph_t = std::vector<std::vector<int> >;  
  struct entry_t { int key; int val; };
  struct patch_t { std::vector<int> node; std::vector<int> elem; };
}
