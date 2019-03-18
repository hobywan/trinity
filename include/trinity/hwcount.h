/*
 *                          'hwcount.h'
 *            This file is part of the "trinity" project.
 *               (https://github.com/hobywan/trinity)
 *               Copyright (c) 2016 Hoby Rakotoarivelo.
 *            adapted from 'papi-wrapper' of Sean Chester.
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
/**
 *         Library for updating/accessing standard PAPI counters.
 *                  Copyright (C) 2017 Sean Chester
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * @note The first time you use this library on a new architecture, it may be
 * worth running \c papi_avail and \c papi_native_avail to determine which
 * hardware counters are available on the machine that you are using.
 */
#if HAVE_PAPI

#include <header.h>
#include <papi.h>

/**
 * The base class that defines the general behaviour related to any
 * set of hardware counters.
 * @tparam NUM_COUNTERS The number of hardware counters used in the given subset.
 *
 * A PAPI object defines all the methods related to hardware
 * counters without knowing how many there are (the template parametre
 * NUM_COUNTERS) and without knowing what they are. These are specialised
 * within the specific subclasses.
 */
class PAPI {

public:

  virtual ~PAPI() {}

  /** Starts tracking/reset all the hardware counters in this subset. */
  void start();
  /** Resets all counters to zero. */
  void reset();
  /**
   * Stops tracking all the hardware counters in this subset.
   * @post The values accessible in this object have been updated
   * to reflect the change in the actual, physical hardware counters
   * since start() was invoked.
   * @note Without calling stop(), none of the member variables of
   * this object will have been updated.
   */
  void stop();

  /**
   * Returns a string of tab-separated human-readable names
   * for this set of counters, in the same order that the
   * counters would be printed to an output stream.
   * @return The string of tab-separated headers.
   */
  std::string headers();

  // operator overloads
  PAPI& operator=(const PAPI& other);
  PAPI& operator+=(const PAPI& other);
  PAPI& operator-=(const PAPI& other);
  PAPI& operator/=(const uint32_t scalar);
  PAPI& operator*=(const uint32_t scalar);

  // static functions

  /**
   * Starts up one set of PAPI counters for each thread.
   * @param counter_array An array of sets of PAPI counters, one for
   * each thread.
   * @post All threads will now start tracking PAPI statistics.
   */
  static void start(std::vector<PAPI*>& counters);
  static void stop(std::vector<PAPI*>& counters);
  static void sum(std::vector<PAPI*>& counters);

  // static public constants

  /** Indicates not to use PAPI counters. */
  static const uint32_t papi_mode_off = 0;
  /** Indicates that pre-built L2 and L3 cache performance counters should be used. */
  static const uint32_t papi_mode_cache = 1;
  /** Indicates that pre-built branch prediction counters should be used. */
  static const uint32_t papi_mode_branch = 2;
  /** Indicates that prebuilt throughput instruction counters should be used. */
  static const uint32_t papi_mode_cycle = 3;
  /** Indicates that prebuilt tlb counters should be used. */
  static const uint32_t papi_mode_tlb = 4;
  /** Indicates that custom (i.e., user-selected) PAPI events will be specified. */
  static const uint32_t papi_mode_custom = 5;

  friend std::ostream& operator<<(std::ostream& os, const PAPI& p);

protected:

  /**
   * Registers a new counter with this object based on a PAPI event. Accepts both
   * preset and native event names.
   * @param counter_name The ASCII name of the PAPI event to be registered
   * (decoding to an int takes place within this method).
   * @param header A human-readable string indicating what is being tracked by
   * this hardware counter.
   * @post This object is modified to now track a new PAPI event.
   */
  void register_counter(const std::string& counter_name, const std::string& header);

  /**
   * An array containing the actual hardware counter values (recorded at stop()
   * invocations).
   */
  std::vector<long long> values_; // accessed by friended stream overloading.


private:

  /** Writes the polymorphic object to an output stream */
  virtual void to_stream(std::ostream& os) const = 0;

  /**
   * A flag indicating whether or not this set of papi counters is currently
   * being tracked.
   */
  bool papi_started_ = false;
  /** An array mapping the indexes used in this object to PAPI hardware counter ids. */
  std::vector<int> counters_;
  /** An array that maps PAPI hardware counters onto human-readable strings. */
  std::vector<std::string> headers_;
};


/**
 * A subset of PAPI hardware counters related to the distribution of instructions.
 */
class PAPI_Instructions : public PAPI {
public:

  PAPI_Instructions();

  /** Returns the number of load instructions issued. */
  long long inline load_instructions() const;
  /** Returns the number of store instructions issued. */
  long long inline store_instructions() const;
  /** Returns the number of branch instructions issued. */
  long long inline branch_instructions() const;
  /** Returns the total number of instructions issued. */
  long long inline total_instructions() const;

  /** Returns the number of loadFile instructions issued as a fraction of the total instructions. */
  double load_instructions_ratio() const;
  /** Returns the number of storeFile instructions issued as a fraction of the total instructions. */
  double store_instructions_ratio() const;
  /** Returns the number of branch instructions issued as a fraction of the total instructions. */
  double branch_instructions_ratio() const;

private:
  void to_stream(std::ostream& os) const override;
};


/**
 * A subset of PAPI hardware counters related to the distribution of cycles.
 */
class PAPI_Cycles : public PAPI {

public:

  PAPI_Cycles();

  /** Returns the total number of cycles spent idling. */
  inline long long idle_cycles() const;
  /** Returns the total number of cycles that are utilised. */
  inline long long utilised_cycles() const;
  /** Returns the total number of cycles that are stalled. */
  inline long long stalled_cycles() const;
  /** Returns the total number of cycles. */
  inline long long cycles() const;
  /** Returns the total number of instructions retired. */
  inline long long instructions() const;
  /** Returns the fraction of cycles that are stalled due to any resource. */
  double stalled_cycles_ratio() const;
  /** Returns the fraction of cycles that are idled. */
  double idle_cycles_ratio() const;
  /** Returns the fraction of cycles that are maximally utilised. */
  double utilised_cycles_ratio() const;
  /** Returns the average number of cycles spent retiring an instruction (CPI). */
  double cycles_per_instruction() const;

private:
  void to_stream(std::ostream& os) const override;
};

/**
 * A subset of PAPI hardware counters related to cache hit performance.
 */
class PAPI_Cache : public PAPI {

public:
  PAPI_Cache();

  /** Returns the number of Level 2 total cache misses. */
  long long inline l2_cache_misses() const;
  /** Returns the number of Level 2 total cache accesses. */
  long long inline l2_cache_accesses() const;
  /** Returns the number of Level 3 total cache misses. */
  long long inline l3_cache_misses() const;
  /** Returns the number of Level 3 total cache misses. */
  long long inline l3_cache_accesses() const;

  /**
   * Returns the fraction of L2/L3 cache accesses (both data and instruction) that were
   * missed (because the resource was not available in L2/L3 cache).
   */
  double l2_miss_ratio() const;
  double l3_miss_ratio() const;

//private:
  void to_stream(std::ostream& os) const override;
};

/**
 * A subset of PAPI hardware counters related to branch prediction performance.
 */
class PAPI_Branch : public PAPI {

public:
  PAPI_Branch();

  /** Returns the number of branch instructions issued. */
  long long inline branch_instructions() const;
  /** Returns the number of branch instructions mispredicted. */
  long long inline branch_misses() const;
  /** Returns the number of branch instructions correctly predicted. */
  long long inline branch_hits() const;
  /** Returns the fraction of conditional branches that were incorrectly predicted. */
  double misprediction_ratio() const;
  /** Returns the fraction of conditional branches that were predicted correctly. */
  double prediction_ratio() const;

private:
  void to_stream(std::ostream& os) const override;
};

/**
 * A subset of PAPI hardware counters related to transaction lookaside buffer performance.
 */
class PAPI_TLB : public PAPI {

public:
  PAPI_TLB();

  /** Returns the total number of data tlb misses. */
  long long inline data_tlbs() const;
  /** Returns the total number of instruction tlb misses. */
  long long inline instruction_tlbs() const;

private:
  void to_stream(std::ostream& os) const override;
};

/**
 * An set of PAPI hardware counters that use custom-defined PAPI events.
 */
class PAPI_Custom : public PAPI {

public:
  /**
   * Constructs a new papi_custom set using a specified set of PAPI events.
   * @param event_names Pairs of event names corresponding to the PAPI event
   * names and the human-readable headers for the custom events that should
   * be tracked
   * @post Constructs a new instance of a papi_custom set.
   */
  PAPI_Custom(const std::vector<std::pair<std::string, std::string> >& event_names);

private:
  void to_stream(std::ostream& os) const override;
};

namespace trinity { namespace papi {
//
void init(int num_cores, int mode, std::vector<PAPI*>* counters);
//
void finalize(std::vector<PAPI*> papi_counters[4]);
}} // namespace trinity::papi
#endif
