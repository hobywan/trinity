/* ------------------------------------ */
#pragma once
/* ------------------------------------ */
#include "header.h"
#include "timer.h"
/* -------------------------------- */
namespace trinity { namespace sync {

template<typename type_t>
inline bool compare_and_swap(type_t *flag, int expected, int value) {
  return __sync_bool_compare_and_swap(flag, expected, value);
}

template<typename type_t>
inline int fetch_and_add(type_t *shared, int value) {
  return __sync_fetch_and_add(shared, value);
}

template<typename type_t>
inline int fetch_and_sub(type_t *shared, int value) {
  return __sync_fetch_and_sub(shared, value);
}
/* -------------------------------- */
inline void task_reduction(int *array, std::vector<int> *heap, int *count, int stride) {
  size_t nb = heap->size();
  if (nb) {
    const int offset = sync::fetch_and_add(count, int(nb / stride));
    std::memcpy(array + (offset * stride), heap->data(), nb * sizeof(int));
    heap->clear();
  }
#pragma omp barrier
}

/* -------------------------------- */
inline void task_reduction(int *array, std::vector<int> *heap, int *count) {
  task_reduction(array, heap, count, 1);
}

/* -------------------------------- */
inline int kernel_prefix_sum(int *begin, int *end, size_t grain_size) {
  size_t len = end - begin;
  if (len < grain_size) {
    for (size_t i = 1; i < len; ++i)
      begin[i] += begin[i - 1];
  } else {

    size_t mid = len / 2;
    int sum = 0;

#pragma omp task shared(sum)
    sum = kernel_prefix_sum(begin, begin + mid, grain_size);
#pragma omp task
    kernel_prefix_sum(begin + mid, end, grain_size);
#pragma omp taskwait

#pragma omp parallel for
    for (size_t i = mid; i < len; ++i)
      begin[i] += sum;
  }
  return begin[len - 1];
}

/* -------------------------------- */
inline void prefix_sum(int *values, size_t nb, size_t grain_size) {
#pragma omp single
  kernel_prefix_sum(values, values + nb, grain_size);
}

/* -------------------------------- */
inline void task_reduction(int *array, std::vector<int> *heap, int *count, int *off) {
  size_t nb = heap->size();
  int tid = omp_get_thread_num();
  int cores = omp_get_num_threads();  // negligible cost

  if (tid < cores - 1) {
    off[tid + 1] = nb;
  } else {
    off[0] = 0;
    *count = nb;
  }
#pragma omp barrier

  prefix_sum(off, cores, 10);
  std::memcpy(array + off[tid], heap->data(), sizeof(int) * nb);
  heap->clear();

#pragma omp single
  *count += off[cores - 1]; // off[p-1] + N[p-1]
}

/* -------------------------------- */
inline void check_realloc(std::vector<int> *stenc, int index, size_t chunk, int verbose) {
  if (__builtin_expect(stenc[index].size() <= chunk, 0))
    // at least one thread will perform the realloc ('single' is not ok)
#pragma omp critical(resize)
  {
    size_t old = stenc[index].size();
    // (!) recheck
    if (old <= chunk) {
      stenc[index].resize(old * 2);
      std::fill(stenc[index].begin() + old, stenc[index].end(), -1);
      if (verbose)
        std::fprintf(stderr, "warning: stenc[%d] was reallocated\n", index);
    }
  }
}


}} // namespace trinity
