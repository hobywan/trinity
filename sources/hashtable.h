/*
 *                          'hashtable.h'
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
template<typename type_t = int>
class Hashtable {

public:

  // rule of five
  Hashtable() = default;
  Hashtable(const Hashtable& other) = delete;
  Hashtable& operator=(Hashtable other) = delete;
  Hashtable(Hashtable&& other) noexcept = delete;
  Hashtable& operator=(Hashtable&& other) noexcept = delete;
  Hashtable(size_t size, size_t bucket, size_t stride);
  ~Hashtable();

  int generateKey(int i, int j, int scale, int nb_cores) const;
  void push(int key, const std::initializer_list<type_t>& val);
  int getValue(int v1, int v2) const;
  size_t getCapacity() const;
  void reset();

private:

  type_t** bucket = nullptr;
  int*     offset = nullptr;
  size_t   capacity = 0;
  size_t   size = 0;
  size_t   stride = 0;

};
/* --------------------------------------------------------------------------- */
} // namespace trinity

#include "hashtable.impl.h"