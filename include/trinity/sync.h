/*
 *                          'sync.h'
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
namespace trinity { namespace sync {
/* -------------------------------------------------------------------------- */
void reduceTasks(int* array, std::vector<int>* heap, int* count, int stride);
void reduceTasks(int* array, std::vector<int>* heap, int* count, int* off);
void prefixSum(int* values, size_t nb, size_t grain_size);
void reallocBucket(std::vector<int>* bucket, int index, size_t chunk, int verb);

/* -------------------------------------------------------------------------- */
template <typename type_t,
          typename = std::enable_if_t<std::is_integral<type_t>::value> >
bool compareAndSwap(type_t* flag, int expected, int value) {
  assert(flag != nullptr);
  return __sync_bool_compare_and_swap(flag, expected, value);
}

/* -------------------------------------------------------------------------- */
template <typename type_t,
          typename = std::enable_if_t<std::is_integral<type_t>::value> >
type_t fetchAndAdd(type_t* shared, int value) {
  assert(shared != nullptr);
  return __sync_fetch_and_add(shared, value);
}

/* -------------------------------------------------------------------------- */
template <typename type_t,
          typename = std::enable_if_t<std::is_integral<type_t>::value> >
type_t fetchAndSub(type_t* shared, int value) {
  assert(shared != nullptr);
  return __sync_fetch_and_sub(shared, value);
}
/* -------------------------------------------------------------------------- */
}} // namespace trinity::sync