/*
 *                          'sync.h'
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
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
/* ------------------------------------ */
namespace trinity { namespace sync {
/* ------------------------------------ */
template<typename type_t>
inline bool compareAndSwap(type_t* flag, int expected, int value) {
  return __sync_bool_compare_and_swap(flag, expected, value);
}
/* ------------------------------------ */
template<typename type_t>
inline type_t fetchAndAdd(type_t* shared, int value) {
  return __sync_fetch_and_add(shared, value);
}
/* ------------------------------------ */
template<typename type_t>
inline type_t fetchAndSub(type_t* shared, int value) {
  return __sync_fetch_and_sub(shared, value);
}

/* ------------------------------------ */
inline void reduceTasks(int* array, std::vector<int>* heap, int* count, int stride) {
  size_t nb = heap->size();
  if (nb) {
    const int offset = sync::fetchAndAdd(count, int(nb / stride));
    std::memcpy(array + (offset * stride), heap->data(), nb * sizeof(int));
    heap->clear();
  }
#pragma omp barrier
}

/* ------------------------------------ */
inline int kernelPrefixSum(int* begin, int* end, size_t grain_size) {
  size_t len = end - begin;
  if (len < grain_size) {
    for (size_t i = 1; i < len; ++i)
      begin[i] += begin[i - 1];
  } else {

    size_t mid = len / 2;
    int sum = 0;

#pragma omp task shared(sum)
    sum = kernelPrefixSum(begin, begin + mid, grain_size);
#pragma omp task
    kernelPrefixSum(begin + mid, end, grain_size);
#pragma omp taskwait

#pragma omp parallel for
    for (size_t i = mid; i < len; ++i)
      begin[i] += sum;
  }
  return begin[len - 1];
}

/* ------------------------------------ */
inline void prefixSum(int* values, size_t nb, size_t grain_size) {
#pragma omp single
  kernelPrefixSum(values, values + nb, grain_size);
}

/* ------------------------------------ */
inline void reduceTasks(int* array, std::vector<int>* heap, int* count, int* off) {
  size_t nb = heap->size();
  int tid = omp_get_thread_num();
  int cores = omp_get_num_threads();  // negligible cost

  if (tid < cores - 1) {
    off[tid + 1] = (int) nb;
  } else {
    off[0] = 0;
    *count = (int) nb;
  }
#pragma omp barrier

  prefixSum(off, (size_t) cores, 10);
  std::memcpy(array + off[tid], heap->data(), sizeof(int) * nb);
  heap->clear();

#pragma omp single
  *count += off[cores - 1]; // off[p-1] + N[p-1]
}

/* ------------------------------------ */
inline void reallocBucket(std::vector<int>* bucket, int index, size_t chunk, int verbose) {
  if (__builtin_expect(bucket[index].size() <= chunk, 0))
    // at least one thread will perform the reallocation ('single' is not ok)
#pragma omp critical(resize)
  {
    size_t old = bucket[index].size();
    // (!) recheck
    if (old <= chunk) {
      bucket[index].resize(old * 2);
      std::fill(bucket[index].begin() + old, bucket[index].end(), -1);
      if (verbose)
        std::fprintf(stderr, "warning: stenc[%d] was reallocated\n", index);
    }
  }
}
/* ------------------------------------ */
}} // namespace trinity::sync
