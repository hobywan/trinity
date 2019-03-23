/*
 *                      'hashtable.impl.h'
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
namespace trinity {
/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
Hashtable<type_t, flag_t>::Hashtable(size_t table_size,
                                     size_t bucket_size,
                                     size_t bucket_stride) {
  assert(table_size);
  assert(bucket_size);
  assert(bucket_stride);

  nb_cores = omp_get_max_threads();
  size     = table_size;
  capacity = bucket_size;
  stride   = bucket_stride;
  offset   = new int[table_size];
  bucket   = new type_t* [table_size];

#pragma omp parallel for
  for (int i = 0; i < (int) table_size; ++i) {
    bucket[i] = new type_t[capacity];
  }
}

/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
Hashtable<type_t, flag_t>::~Hashtable() {

  for (int i = 0; i < (int) size; ++i) {
    delete[] bucket[i];
  }
  delete[] bucket;
  delete[] offset;
}

/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
type_t Hashtable<type_t, flag_t>::generateKey(type_t i, type_t j, size_t scale) const {

  assert(scale);
  assert(nb_cores);

  auto min_key = (uint32_t) std::min(i, j);
  return (type_t) tools::hash(min_key) % (scale * nb_cores);
}

/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
size_t Hashtable<type_t, flag_t>::getCapacity() const { return size; }

/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
void Hashtable<type_t, flag_t>::push(type_t key, const std::initializer_list<type_t>& val) {
  assert(val.size() == stride);

  auto j = sync::fetchAndAdd(offset + key, (int) stride);
  assert((j + stride) < capacity);

  for (auto i = 0; i < (int) stride; ++i)
    bucket[key][j + i] = *(val.begin() + i);
}

/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
type_t Hashtable<type_t, flag_t>::getValue(type_t v1, type_t v2, bool use_hash) const {

  assert(stride == 2);
  type_t key  = (use_hash ? generateKey(v1,v2) : std::min(v1, v2));
  type_t hint = std::max(v1, v2);

  for (int k = 0; k < offset[key] - 1; k += 2) {
    if (bucket[key][k] == hint) {
      return bucket[key][k + 1];
    }
  }
  return (type_t) -1; // not found
}

/* -------------------------------------------------------------------------- */
template <typename type_t, typename flag_t>
void Hashtable<type_t, flag_t>::reset() {
#pragma omp for
  for (int i = 0; i < (int) size; ++i) {
    std::memset(bucket[i], -1, capacity * sizeof(int));
  }

#pragma omp for
  for (int i = 0; i < (int) size; ++i) {
    offset[i] = 0;
  }
}
/* -------------------------------------------------------------------------- */
} // namespace trinity