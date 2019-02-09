/**
 * Library for updating/accessing standard PAPI counters
 *
 * @note The first time you use this library on a new architecture, it may be
 * worth running \c papi_avail and \c papi_native_avail to determine which
 * hardware counters are available on the machine that you are using.
 */
#pragma once

#ifdef HAVE_PAPI
#include <header.h>
#include <papi.h>

/**
 * The base class that defines the general behaviour related to any
 * set of hardware counters.
 * @tparam NUM_COUNTERS The number of hardware counters used in the given subset.
 *
 * A papi_t object defines all the methods related to hardware
 * counters without knowing how many there are (the template parametre
 * NUM_COUNTERS) and without knowing what they are. These are specialised
 * within the specific subclasses.
 */
class papi_t {

public:

	virtual ~papi_t() {}

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
	papi_t & operator= (const papi_t & other);
	papi_t & operator+=(const papi_t & other);
	papi_t & operator-=(const papi_t & other);
	papi_t & operator/=(const uint32_t scalar);
	papi_t & operator*=(const uint32_t scalar);

	// static functions

	/**
	 * Starts up one set of PAPI counters for each thread.
	 * @param counter_array An array of sets of PAPI counters, one for
	 * each thread.
	 * @post All threads will now start tracking PAPI statistics.
	 */
	static void start(std::vector<papi_t*> &counters);
	static void stop (std::vector<papi_t*> &counters);
	static void sum  (std::vector<papi_t*> &counters);

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

	friend std::ostream& operator<<(std::ostream& os, const papi_t & p);

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
	void register_counter( const std::string &counter_name, const std::string &header );

	/**
	 * An array containing the actual hardware counter values (recorded at stop()
	 * invocations).
	 */
	std::vector< long long > values_; // accessed by friended stream overloading.


private:

	/** Writes the polymorphic object to an output stream */
	virtual void to_stream( std::ostream& os ) const = 0;

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
class papi_instructions : public papi_t {
public:

	papi_instructions();

	/** Returns the number of load instructions issued. */
	long long inline load_instructions() const;
	/** Returns the number of store instructions issued. */
	long long inline store_instructions() const;
	/** Returns the number of branch instructions issued. */
	long long inline branch_instructions() const;
	/** Returns the total number of instructions issued. */
	long long inline total_instructions() const;

	/** Returns the number of load instructions issued as a fraction of the total instructions. */
	double load_instructions_ratio() const;
	/** Returns the number of store instructions issued as a fraction of the total instructions. */
	double store_instructions_ratio() const;
	/** Returns the number of branch instructions issued as a fraction of the total instructions. */
	double branch_instructions_ratio() const;

private:
	void to_stream( std::ostream& os ) const;
};


/**
 * A subset of PAPI hardware counters related to the distribution of cycles.
 */
class papi_cycles : public papi_t {

public:

	papi_cycles();

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
	void to_stream( std::ostream& os ) const;
};

/**
 * A subset of PAPI hardware counters related to cache hit performance.
 */
class papi_cache : public papi_t {

public:
	papi_cache();

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
	void to_stream( std::ostream& os ) const;
};

/**
 * A subset of PAPI hardware counters related to branch prediction performance.
 */
class papi_branch : public papi_t {

public:
	papi_branch();

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
	void to_stream( std::ostream& os ) const;
};

/**
 * A subset of PAPI hardware counters related to transaction lookaside buffer performance.
 */
class papi_tlb : public papi_t {

public:
	papi_tlb();

	/** Returns the total number of data tlb misses. */
	long long inline data_tlbs() const;
	/** Returns the total number of instruction tlb misses. */
	long long inline instruction_tlbs() const;

private:
	void to_stream( std::ostream& os ) const;
};

/**
 * An set of PAPI hardware counters that use custom-defined PAPI events.
 */
class papi_custom : public papi_t {

public:
	/**
	 * Constructs a new papi_custom set using a specified set of PAPI events.
	 * @param event_names Pairs of event names corresponding to the PAPI event
	 * names and the human-readable headers for the custom events that should
	 * be tracked
	 * @post Constructs a new instance of a papi_custom set.
	 */
	papi_custom( const std::vector< std::pair< std::string, std::string > > &event_names );

private:
	void to_stream( std::ostream& os ) const;
};

namespace trinity {
  //
  void papi_init(const uint32_t num_cores,
                 const uint32_t papi_mode,
                 std::vector<papi_t*> papi_counters[4]);
  //
  void papi_finalize(std::vector<papi_t*> papi_counters[4]);
}
#endif
