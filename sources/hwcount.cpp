/**
 * Implementation of papi-wrapper library.
 */
#ifdef HAVE_PAPI
#include "hwcount.h"

/* ------------------------------------ */
void papi_t::start(std::vector<papi_t*>& counter_array) {
#pragma omp parallel num_threads( counter_array.size())
  {
    counter_array[omp_get_thread_num()]->start();
  }
}

/* ------------------------------------ */
void papi_t::stop(std::vector<papi_t*>& counter_array) {
#pragma omp parallel num_threads( counter_array.size())
  {
    counter_array[omp_get_thread_num()]->stop();
  }
}

/* ------------------------------------ */
void papi_t::sum(std::vector<papi_t*>& counter_array) {
  for (uint32_t i = 1; i < counter_array.size(); ++i) {
    *(counter_array[0]) += *(counter_array[i]);
  }
}

/* ------------------------------------ */
std::ostream& operator<<(std::ostream& os, const papi_t& b) {
  b.to_stream(os);
  return os;
}

/* -- print --------------------------- */

void papi_branch::to_stream(std::ostream& os) const {
  os << branch_hits() << "\t" << branch_misses() << "\t" << branch_instructions()
     << "\t" << misprediction_ratio();
}

/* ------------------------------------ */
void papi_cycles::to_stream(std::ostream& os) const {
  os << instructions() << "\t" << cycles() << "\t" << cycles_per_instruction();
}

/* ------------------------------------ */
void papi_instructions::to_stream(std::ostream& os) const {
  os << load_instructions_ratio() << "\t" << store_instructions_ratio()
     << "\t" << branch_instructions_ratio();
}

/* ------------------------------------ */
void papi_tlb::to_stream(std::ostream& os) const {
  os << data_tlbs() << "\t" << instruction_tlbs();
}

/* ------------------------------------ */
void papi_cache::to_stream(std::ostream& os) const {
  //first output L2 statistics: misses, accesses, and miss ratio
  //then output L3 statistics: misses, accesses, and miss ratio
  os << l2_cache_misses() << "\t" << l2_cache_accesses() << "\t" << l2_miss_ratio()
     << "\t"
     << l3_cache_misses() << "\t" << l3_cache_accesses() << "\t" << l3_miss_ratio();
}

/* ------------------------------------ */
void papi_custom::to_stream(std::ostream& os) const {
  if (values_.size() > 0) {
    os << values_[0];
    for (uint32_t i = 1; i < values_.size(); ++i) {
      os << "\t" << values_[i];
    }
  }
}

/* ---------- helper functions -------- */
long long absolute_difference(const long long v1, const long long v2) {
  return (v1 > v2 ? v1 - v2 : v2 - v1);
}

/* ------------------------------------ */
void handle_error(int retval) {
  std::cerr << "PAPI error " << retval << ": " << PAPI_strerror(retval) << std::endl;
  exit(1);
}

//--------- Base class implementation -------------//
void papi_t::register_counter(const std::string& counter_name, const std::string& header) {
  int counter;
  char c_style_string[counter_name.length() + 1];
  strcpy(c_style_string, counter_name.c_str());
  c_style_string[counter_name.length()] = '\0';
  int retval = PAPI_event_name_to_code(c_style_string, &counter);
  if (retval not_eq PAPI_OK) {
    std::cerr << "Could not decode PAPI counter: " << counter_name << std::endl;
    return;
  } else {
    values_.push_back(0);
    counters_.push_back(counter);
    headers_.push_back(header);
  }
}

/* ------------------------------------ */
void papi_t::reset() {
  for (auto it = values_.begin(); it not_eq values_.end(); ++it) {
    *it = 0;
  }
}

/* ------------------------------------ */
void papi_t::start() {
  std::vector<int> eventsMutable(counters_.data(), counters_.data() + counters_.size());
  int retval = PAPI_start_counters(&eventsMutable[0], counters_.size());
  if (retval == PAPI_OK) {
    papi_started_ = true;
  } else {
    std::cerr << "PAPI error " << retval << ": " << PAPI_strerror(retval) << std::endl;
    papi_started_ = false;
  }
}

/* ------------------------------------ */
void papi_t::stop() {

  if (values_.size() == 0) { return; }
  if (papi_started_) {
    long long v[counters_.size()];
    int retval = PAPI_stop_counters(&v[0], counters_.size());
    if (retval not_eq PAPI_OK) handle_error(retval);
    for (uint32_t i = 0; i < values_.size(); ++i) {
      values_[i] += v[i];
    }
  } else {
    for (auto it = values_.begin(); it not_eq values_.end(); ++it) {
      *it = -1;
    }
  }
  papi_started_ = false;
}

/* ------------------------------------ */
std::string papi_t::headers() {

  if (headers_.size() == 0) { return ""; }
  std::string output(headers_[0]);
  for (uint32_t i = 1; i < headers_.size(); ++i) {
    output = output + "\t" + headers_[i];
  }
  return output;
}

/* ------------------------------------ */
papi_t& papi_t::operator=(const papi_t& other) {
  for (uint32_t i = 0; i < values_.size(); ++i) {
    values_[i] = other.values_[i];
  }
  return *this;
}

/* ------------------------------------ */
papi_t& papi_t::operator+=(const papi_t& other) {
  for (uint32_t i = 0; i < values_.size(); ++i) {
    values_[i] += other.values_[i];
  }
  return *this;
}

/* ------------------------------------ */
papi_t& papi_t::operator-=(const papi_t& other) {
  for (uint32_t i = 0; i < values_.size(); ++i) {
    values_[i] = absolute_difference(values_[i], other.values_[i]);
  }
  return *this;
}

/* ------------------------------------ */
papi_t& papi_t::operator/=(const uint32_t scalar) {
  for (auto it = values_.begin(); it not_eq values_.end(); ++it) {
    *it /= scalar;
  }
  return *this;
}

/* ------------------------------------ */
papi_t& papi_t::operator*=(const uint32_t scalar) {
  for (auto it = values_.begin(); it not_eq values_.end(); ++it) {
    *it *= scalar;
  }
  return *this;
}

//--------- Predefined classes implementation -------------//

papi_instructions::papi_instructions() {
  this->register_counter(std::string("PAPI_LD_INS"), std::string("loads")); //NOT on AMD
  this->register_counter(std::string("PAPI_SR_INS"), std::string("stores")); //NOT on AMD
  this->register_counter(std::string("PAPI_BR_INS"), std::string("branches"));
  this->register_counter(std::string("PAPI_TOT_INS"), std::string("total"));
}

/* ------------------------------------ */
long long inline papi_instructions::load_instructions() const { return values_[0]; }
long long inline papi_instructions::store_instructions() const { return values_[1]; }
long long inline papi_instructions::branch_instructions() const { return values_[2]; }
long long inline papi_instructions::total_instructions() const { return values_[3]; }

/* ------------------------------------ */
double papi_instructions::load_instructions_ratio() const {
  return load_instructions() / (double) total_instructions();
}

/* ------------------------------------ */
double papi_instructions::store_instructions_ratio() const {
  return store_instructions() / (double) total_instructions();
}

/* ------------------------------------ */
double papi_instructions::branch_instructions_ratio() const {
  return branch_instructions() / (double) total_instructions();
}

/* ------------------------------------ */
papi_cycles::papi_cycles() {
  this->register_counter("PAPI_STL_ICY", "cycles_no_instructions_issue");
  this->register_counter("PAPI_FUL_CCY", "cycles_max_instructions_completed"); //NOT on AMD
  this->register_counter("PAPI_RES_STL", "cycles_stalled_any_resource");
  this->register_counter("PAPI_TOT_CYC", "total_cycles");
  this->register_counter("PAPI_TOT_INS", "total_instructions");
}

/* ------------------------------------ */
inline long long papi_cycles::idle_cycles() const { return values_[0]; }
inline long long papi_cycles::utilised_cycles() const { return values_[1]; }
inline long long papi_cycles::stalled_cycles() const { return values_[2]; }
inline long long papi_cycles::cycles() const { return values_[3]; }
inline long long papi_cycles::instructions() const { return values_[4]; }

/* ------------------------------------ */
double papi_cycles::stalled_cycles_ratio() const {
  return stalled_cycles() / (double) cycles();
}

/* ------------------------------------ */
double papi_cycles::idle_cycles_ratio() const {
  return idle_cycles() / (double) cycles();
}

/* ------------------------------------ */
double papi_cycles::utilised_cycles_ratio() const {
  return utilised_cycles() / (double) cycles();
}

/* ------------------------------------ */
double papi_cycles::cycles_per_instruction() const {
  return cycles() / (double) instructions();
}

/* ------------------------------------ */
papi_cache::papi_cache() {
  this->register_counter("PAPI_L2_TCM", "L2_total_cache_misses");
  this->register_counter("PAPI_L2_TCA", "L2_total_cache_accesses");
  this->register_counter("PAPI_L3_TCM", "L3_total_cache_misses");
  this->register_counter("PAPI_L3_TCA", "L3_total_cache_accesses");
}

/* ------------------------------------ */
long long inline papi_cache::l2_cache_misses() const { return values_[0]; }
long long inline papi_cache::l2_cache_accesses() const { return values_[1]; }
long long inline papi_cache::l3_cache_misses() const { return values_[2]; }
long long inline papi_cache::l3_cache_accesses() const { return values_[3]; }

/* ------------------------------------ */
double papi_cache::l2_miss_ratio() const {
  return l2_cache_misses() / (double) l2_cache_accesses();
}

/* ------------------------------------ */
double papi_cache::l3_miss_ratio() const {
  return l3_cache_misses() / (double) l3_cache_accesses();
}

/* ------------------------------------ */
papi_branch::papi_branch() {
  this->register_counter("PAPI_BR_PRC", "conditional_branches_correctly_predicted"); //NOT on AMD
  this->register_counter("PAPI_BR_MSP", "conditional_branches_mispredicted");
}

/* ------------------------------------ */
long long inline papi_branch::branch_instructions() const { return values_[0] + values_[1]; }
long long inline papi_branch::branch_misses() const { return values_[1]; }
long long inline papi_branch::branch_hits() const { return values_[0]; }

/* ------------------------------------ */
double papi_branch::misprediction_ratio() const {
  return branch_misses() / (double) branch_instructions();
}

/* ------------------------------------ */
double papi_branch::prediction_ratio() const {
  return branch_hits() / (double) branch_instructions();
}

/* ------------------------------------ */
papi_tlb::papi_tlb() {
  this->register_counter("PAPI_TLB_DM", "dtlb_misses");
  this->register_counter("PAPI_TLB_IM", "itlb_misses");
}

/* ------------------------------------ */
long long inline papi_tlb::data_tlbs() const { return values_[0]; }
long long inline papi_tlb::instruction_tlbs() const { return values_[1]; }

/* ------------------------------------ */
papi_custom::papi_custom(const std::vector <std::pair<std::string, std::string>>& events) {
  for (auto it = events.begin(); it not_eq events.end(); ++it) {
    this->register_counter(it->first, it->second);
  }
}

/**
 * Preprocessing stage
 * @param num_cores the number of cores (1 thread/core)
 * @param papi_mode the type of papi counter to track
 * @param papi_counters array of polymorphic papi counter sets per thread and per kernel
 */
void trinity::papi_init(const uint32_t num_cores,
                        const uint32_t papi_mode,
                        std::vector<papi_t*> hw_counters[4]) {
  assert(num_cores);
  assert(hw_counters not_eq nullptr);
  // init random seed
  srand((unsigned) time(0));

  // Initialize PAPI library for each thread.
  if (PAPI_is_initialized() == PAPI_NOT_INITED) {
    PAPI_library_init(PAPI_VER_CURRENT);
#pragma omp parallel num_threads(num_cores)
    if (PAPI_thread_init(pthread_self) not_eq PAPI_OK) {
      exit(0);
    }
  }

  // step 2: initialise and start the papi counter sets, one for each thread
  if (papi_mode not_eq papi_t::papi_mode_off) {
    for (int k = 0; k < 4; ++k) {
      for (uint32_t t = 0; t < num_cores; ++t) {
        switch (papi_mode) {
          case papi_t::papi_mode_cache :
            hw_counters[k].push_back(new papi_cache());
            break;

          case papi_t::papi_mode_branch :
            hw_counters[k].push_back(new papi_branch());
            break;

          case papi_t::papi_mode_cycle :
            hw_counters[k].push_back(new papi_cycles());
            break;

          case papi_t::papi_mode_tlb :
            hw_counters[k].push_back(new papi_tlb());
            break;
        }
      }
    }
  }
}

/**
 * Report the PAPI counters' values, reduce over all threads, and clean up.
 */
void trinity::papi_finalize(std::vector<papi_t*> hw_counters[4]) {
  for (int k = 0; k < 4; ++k) {
    // 1: reduce counters for each kernel
    papi_t::sum(hw_counters[k]);
    std::cout << *(hw_counters[k][0]) << std::endl;
    // 2: clean up
    while (!hw_counters[k].empty())
      hw_counters[k].pop_back();   // delete performed here
  }
}
#endif
