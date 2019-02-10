/* ------------------------------------*/
#pragma once
/* ------------------------------------*/
#include "sync.h"
#include "tools.h"
/* ------------------------------------*/
namespace trinity {

template<typename type_t>
class Table {

public:
  /* ------------------------------------*/
  Table() = default;
  /* ------------------------------------*/
  Table(size_t size, size_t bucket, size_t stride)
    : size_(size),
      capacity_(bucket),
      stride_(stride)
  {
    offset_ = new int[size_];
    bucket_ = new type_t* [size_];

#pragma omp parallel for
    for (int i = 0; i < size_; ++i)
      bucket_[i] = new type_t[capacity_];
  }

  /* ------------------------------------*/
  ~Table() {
    for (int i = 0; i < size_; ++i)
      delete[] bucket_[i];
    delete[] bucket_;
    delete[] offset_;
  }

  /* ------------------------------------*/
  inline int generateKey(int i, int j, int scale, int nb_cores) const {
    auto min_key = static_cast<uint32_t>(std::min(i, j));
    return tools::hash(min_key) % (scale * nb_cores);
  }

  /* ------------------------------------*/
  inline size_t getCapacity() const { return size_; }

  /* ------------------------------------*/
  void push(int key, const std::initializer_list<type_t>& val) {
    assert(val.size() == stride_);
    int j = sync::fetchAndAdd(offset_ + key, (int) stride_);
    assert((j + stride_) < capacity_);

    for (int i = 0; i < stride_; ++i)
      bucket_[key][j + i] = *(val.begin() + i);
  }

  /* ------------------------------------*/
  inline int getValue(int v1, int v2) const {

    const int index[] = {std::min(v1, v2), std::max(v1, v2)};

    for (int k = 0; k < offset_[*index] - 1; k += 2)
      if (bucket_[*index][k] == *(index + 1))
        return bucket_[*index][k + 1];
    // not found
    return -1;
  }

  /* ------------------------------------*/
  inline void reset() {
#pragma omp for
    for (int i = 0; i < size_; ++i)
      std::memset(bucket_[i], -1, capacity_ * sizeof(int));
#pragma omp for
    for (int i = 0; i < size_; ++i)
      offset_[i] = 0;
  }

private:

  type_t** bucket_;
  int*     offset_;
  size_t capacity_;
  size_t     size_;
  size_t   stride_;

};
}
