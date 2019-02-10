
#ifndef _LOCK_H_
#define _LOCK_H_
/* ------------------------------------
 * memory_order_relaxed : no synch/ordering constraints, only atomicity
 * memory_order_acq_rel : no mem accesses in the current thread
 * can be reorderered before loadFile and after storeFile : all writes in other
 * threads that release the same atomic variable are visible before the modif,
 * and the modif is visible in other threads that acquire the same variable.
 * taken from 'pragmatic' with slight modifs
 */
#include <atomic>
/* ------------------------------------ */
class Lock {

private :
  // integral type for atomic builtins (fetch_or..)
  std::atomic<uint8_t> _state;

public :

  // construct
  Lock() : _state(0) {}
  Lock(const uint8_t value) : _state(value) {}
  Lock(const Lock& other) : _state(other._state.load(std::memory_order_relaxed)) {}

  ~Lock(){}

  // acquire
  inline bool try_lock(){
    uint8_t val = _state.load(std::memory_order_relaxed);
    if(val == 1)
      return false;
    val = _state.fetch_or(1,std::memory_order_acq_rel);
    return (val == 0);
  }

  // spin
  inline void lock(){
    uint8_t val = _state.load(std::memory_order_relaxed);
    if(val == 1){
      while(!try_lock());
      return;
    }
    bool unlocked = _state.compare_exchange_weak(val,1,std::memory_order_acq_rel);
    // spin if not unlocked
    if(!unlocked){
      while(!try_lock());
      return;
    }
  }

  // check
  inline bool locked(){
    return (_state.load(std::memory_order_acquire) == 1);
  }

  inline void unlock(){
    _state.store(0,std::memory_order_release);
  }
};
#endif
