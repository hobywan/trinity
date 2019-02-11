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
/* --------------------------------------------------------------------------- */
#include "header.h"
#include "timer.h"
/* --------------------------------------------------------------------------- */
namespace trinity { namespace sync {
/* --------------------------------------------------------------------------- */

void reduceTasks(int* array, std::vector<int>* heap, int* count, int stride);
void reduceTasks(int* array, std::vector<int>* heap, int* count, int* off);
void prefixSum(int* values, size_t nb, size_t grain_size);
void reallocBucket(std::vector<int>* bucket, int index, size_t chunk, int verbose);

/* --------------------------------------------------------------------------- */
// internal method
namespace { int kernelPrefixSum(int* begin, int* end, size_t grain_size); }
/* --------------------------------------------------------------------------- */
template<typename type_t>
bool compareAndSwap(type_t* flag, int expected, int value) {
  return __sync_bool_compare_and_swap(flag, expected, value);
}

template<typename type_t>
type_t fetchAndAdd(type_t* shared, int value) {
  return __sync_fetch_and_add(shared, value);
}

template<typename type_t>
type_t fetchAndSub(type_t* shared, int value) {
  return __sync_fetch_and_sub(shared, value);
}
/* --------------------------------------------------------------------------- */
}} // namespace trinity::sync