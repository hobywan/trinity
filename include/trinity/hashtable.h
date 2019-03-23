/*
 *                          'hashtable.h'
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
#include "sync.h"
#include "tools.h"
/* -------------------------------------------------------------------------- */
namespace trinity {
/* -------------------------------------------------------------------------- */
template <typename type_t,
          typename = std::enable_if_t<std::is_integral<type_t>::value> >
class Hashtable {

public:

  // rule of five
  Hashtable() = default;
  Hashtable(const Hashtable& other) = delete;
  Hashtable& operator=(Hashtable other) = delete;
  Hashtable(Hashtable&& other) noexcept = delete;
  Hashtable& operator=(Hashtable&& other) noexcept = delete;
  Hashtable(size_t table_size, size_t bucket_size, size_t bucket_stride);
  ~Hashtable();

  type_t generateKey(type_t i, type_t j, size_t scale = 1) const;
  void push(type_t key, const std::initializer_list<type_t>& val);
  type_t getValue(type_t v1, type_t v2, bool use_hash = false) const;
  size_t getCapacity() const;
  void reset();

private:

  type_t** bucket   = nullptr;
  int*     offset   = nullptr;
  size_t   capacity = 1;
  size_t   size     = 0;
  size_t   stride   = 1;
  int      nb_cores = 1;
};
/* -------------------------------------------------------------------------- */
} // namespace trinity

#include "hashtable.impl.h"