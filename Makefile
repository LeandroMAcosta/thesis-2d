# 2D event-driven gas simulation — CPU + OpenMP build.
#
# Independent repo (no parent). Apply the same philosophy as the 1D
# event-driven implementation: postcondition-based testing, per-particle
# event integration, no CUDA (event-driven matches CPU model better).

CC      := gcc
CSTD    := -std=c99
CWARN   := -Wall -Wextra -Werror -pedantic -Wno-unused-result
CINC    := -Iinclude -Ithird_party/tomlc99
CLIBS   := -lm
OMP     := -fopenmp

CFLAGS_RELEASE := $(CSTD) $(CWARN) $(CINC) -O3 -march=native -ffast-math
CFLAGS_DEBUG   := $(CSTD) $(CWARN) $(CINC) -O0 -g3 -fno-omit-frame-pointer \
                  -fsanitize=address,undefined

TOML_CFLAGS := $(CSTD) $(CINC) -O2 -w

# --- Sources ----------------------------------------------------------------

COMMON_SRC := src/main.c src/state.c src/io.c
COMMON_OBJ := $(patsubst src/%.c,build/%.o,$(COMMON_SRC))
TOML_OBJ   := build/toml.o

OMP_SRC := src/physics_event.c src/histogram_omp.c src/backend_omp.c
OMP_OBJ := $(patsubst src/%.c,build/%-omp.o,$(OMP_SRC))

# --- Targets ----------------------------------------------------------------

.PHONY: all debug clean distclean test verify periodicity isotropy independence \
        format main main-event

all: main-event main

main-event: $(OMP_OBJ) $(COMMON_OBJ) $(TOML_OBJ)
	$(CC) $(CFLAGS_RELEASE) $(OMP) -o $@ $^ $(CLIBS)

main: main-event
	@ln -sf main-event main

build/%-omp.o: src/%.c | build
	$(CC) $(CFLAGS_RELEASE) $(OMP) -c -o $@ $<

build/%.o: src/%.c | build
	$(CC) $(CFLAGS_RELEASE) -c -o $@ $<

build/toml.o: third_party/tomlc99/toml.c | build
	$(CC) $(TOML_CFLAGS) -c -o $@ $<

build:
	mkdir -p build

debug: CFLAGS_RELEASE := $(CFLAGS_DEBUG)
debug: clean main-event
	@echo "Built debug binary (main-event) with -O0 -g + ASan/UBSan."

# --- Tests ------------------------------------------------------------------

test: main-event
	@bash tests/regression.sh

verify: main-event
	@bash tests/verify.sh

periodicity: main-event
	@bash tests/periodicity.sh

isotropy: main-event
	@bash tests/isotropy.sh

independence: main-event
	@bash tests/independence.sh

# Run all check scripts in sequence
test-all: test verify periodicity isotropy independence
	@echo "test-all: all suites passed."

format:
	clang-format -style=Microsoft -i src/*.c include/*.h

clean:
	rm -rf build main main-event

distclean: clean
	rm -f graba.dmp graba.dmp.* X*.dat
