/*
 *                          'table.h'
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
#include "sync.h"
#include "tools.h"
/* --------------------------------------------------------------------------- */
namespace trinity {
/* --------------------------------------------------------------------------- */
template<typename type_t>
class Hashtable {

public:

  // rule of five
  Hashtable() = default;
  Hashtable(const Hashtable& other) = delete;
  Hashtable& operator=(Hashtable other) = delete;
  Hashtable(Hashtable&& other) noexcept = delete;
  Hashtable& operator=(Hashtable&& other) noexcept = delete;

  Hashtable(size_t size, size_t bucket, size_t stride) {

    this->size     = size;
    this->capacity = bucket;
    this->stride   = stride;
    this->offset   = new int[size];
    this->bucket   = new type_t* [size];

    #pragma omp parallel for
    for (int i = 0; i < size; ++i){
      this->bucket[i] = new type_t[capacity];
    }
  }

  /* --------------------------------------------------------------------------- */
  ~Hashtable() {
    for (int i = 0; i < size; ++i){
      delete[] bucket[i];
    }
    delete[] bucket;
    delete[] offset;
  }

  /* --------------------------------------------------------------------------- */
  inline int generateKey(int i, int j, int scale, int nb_cores) const {

    auto min_key = static_cast<uint32_t>(std::min(i, j));
    return tools::hash(min_key) % (scale * nb_cores);
  }

  /* --------------------------------------------------------------------------- */
  inline size_t getCapacity() const { return size; }

  /* --------------------------------------------------------------------------- */
  void push(int key, const std::initializer_list<type_t>& val) {
    assert(val.size() == stride);
    int j = sync::fetchAndAdd(offset + key, (int) stride);
    assert((j + stride) < capacity);

    for (int i = 0; i < stride; ++i)
      bucket[key][j + i] = *(val.begin() + i);
  }

  /* --------------------------------------------------------------------------- */
  inline int getValue(int v1, int v2) const {

    const int index[] = {std::min(v1, v2), std::max(v1, v2)};

    for (int k = 0; k < offset[*index] - 1; k += 2)
      if (bucket[*index][k] == *(index + 1))
        return bucket[*index][k + 1];
    // not found
    return -1;
  }

  /* --------------------------------------------------------------------------- */
  inline void reset() {
#pragma omp for
    for (int i = 0; i < size; ++i)
      std::memset(bucket[i], -1, capacity * sizeof(int));
#pragma omp for
    for (int i = 0; i < size; ++i)
      offset[i] = 0;
  }

private:

  type_t** bucket;
  int*     offset;
  size_t   capacity;
  size_t   size;
  size_t   stride;

};
/* --------------------------------------------------------------------------- */
} // namespace trinity
