
BIN=./bin/bench_graph RMAT_ER 4 1 16000000 128000000

perf record -g -s -F 99 $BIN
perf script | $FLAME_DIR/stackcollapse-perf.pl --tid | $FLAME_DIR/flamegraph.pl > perf.svg
perf script | $FLAME_DIR/stackcollapse-perf.pl --chrono | $FLAME_DIR/flamegraph.pl --tid > perf.svg --hash --width=1600

perf stat -B -e cycles,idle-cycles-frontend,idle-cycles-backend,instructions,cache-references,cache-misses $binaire


