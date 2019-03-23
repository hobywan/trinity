
/* -------------------------------------------------------------------------- *
 * Code style guidelines by example (or see .clang-format for an approximation
 * of how to automate part of this)
 * -------------------------------------------------------------------------- */

/*
 * Use "struct" for plain data without method.
 * Use "class" otherwise.
 */
struct Param {
  int depth;
  int verbose;
  int iter;
  int rounds;
};

class Base {
public:
  virtual void bar(int a) = 0;
};

template <typename T, typename U = void>
class Example : Base {

  /*
   * Do not use "typedef", always use "using".
   * The "using" construct actually gets the order sane.
   * Make type aliases whenever it seems useful. Do not hold back,
   * especially if the type is long or *could ever change*.
   */
   using Graph = std::vector<std::vector<int>>;

public:
  /*
   * Use 'defaults' or 'delete' when possible and needed.
   */
  Example() = default;
  Example(MyClassName&&) = default;
  Example(MyClassName const&) = delete;  // no allowed copy construct
  explicit Example(int in_data);

  /*
   * Often your destructor should be virtual
   */
  virtual ~Example() = default;

  /*
   * If it is above 80 cols, wrap to next line with regular indenting
   */
  static inline typename Abc<T,X>::type foo(
    int a, int b, int c, int d
  );

  /*
   * Use the override keyword, it can only help you
   */
  void bar(int a) override {

    /*
     * Initialize stack variables as much as possible
     * (and static when possible)
     */
    int a = 0, tmp = 0;
    int my_local_count = 0;
    std::vector<int> my_vec = {};

    /*
     * Don't be afraid to use "and" or "or". Often it's more
     * readable. It's fully C++ standard compliant.
     */
    if (my_local_count != 0 and a == 0) {
      a = 10;
    }

    // Braces instead of semicolon when empty
    while (true) { }

    // Prefer spacing here
    for (int i = 0; i < 10; i++) {

    }

    /*
     * Prefer for each loop with && reference whenever possible.
     */
    for (auto&& elm : my_vec) { }
  }

  // if function returns bool, it should start with "is" or "has"
  bool isValid() const;

private:

  /*
   * Semantically pack member data within a anonymous struct whenever possible.
   * Provide few comments if member data names are not enough explicit.
   */
  struct {
    int    scale;                     // capacity scale factor
    size_t bucket;                    // bucket capacity
    size_t node;                      // node capacity
    size_t elem;                      // elem capacity
  } capacity;

  int my_data;                        // no trailing underscore for private data
  bool is_done;        // if data is boolean, it should start with "is" or "has"
};

/*
 * Wrap typename args when they get long, start ">" on next line, indented as
 * so. Use "typename" unless you have to use "class": exception being template
 * templates. Space goes here: "template <". Ellipses on on LHS for parameter
 * packs
 */
template <
  typename CollectionT, typename IndexT, typename TupleT, typename RetT,
  typename... Args
>
class MyFunctor {

  /*
   * Wrap function args when they get long, indented as so.
   */
  RetT operator()(
    CollectionT const& col, IndexT const& idx, TupleT&& tup,
    std::tuple<Args...> tup, std::unique_ptr<int> ptr
  ) {
    return {};
  }
};

/*
 * Use an enum class when possible. Degrading automatically (happens with
 * non-strongly typed enums) may introduce bugs
 */
enum class State : int {
  unset  = 0,
  active = 1,
  reeval = 2
};

// A new line should be after the template, even in cases like this:
template <typename T>
class A;

/*
 * Do not indent code under a namespace. Always combine multiple nested
 * namespaces in one line (makes it easier to identify the full namespace path)
 * unless you need to wrap, then try to semantically group them
 */
namespace X { namespace Y {

}} // end namespace X::Y
