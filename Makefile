SRC = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.h)

RTM_OBJS = $(SRC:src/%.cc=build/rtm/%.o)
RTM_OPTIM_OBJS = $(SRC:src/%.cc=build/rtm_optim/%.o)
RTM_DEBUG_OBJS = $(SRC:src/%.cc=build/rtm_debug/%.o)
RTM_HANSON_OBJS = $(SRC:src/%.cc=build/rtm_hanson/%.o)
RTM_MORRICONE_OBJS = $(SRC:src/%.cc=build/rtm_morricone/%.o)
RTM_HPC_OBJS = $(SRC:src/%.cc=build/rtm_hpc2/%.o)

LDFLAGS := $(LDFLAGS) -lgsl -lgslcblas
CXXFLAGS := $(CXXFLAGS) -g -DHAVE_INLINE -std=c++11

rtm_debug: $(RTM_DEBUG_OBJS) $(HEADERS)
	$(CXX) $(RTM_DEBUG_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_debug


rtm: $(RTM_OBJS) $(HEADERS)
	$(CXX) $(RTM_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm

rtm_optim: $(RTM_OPTIM_OBJS) $(HEADERS)
	$(CXX) -fopenmp $(RTM_OPTIM_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_optim

rtm_hanson: $(RTM_HANSON_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_hanson

rtm_morricone: $(RTM_MORRICONE_OBJS) $(HEADERS)
	$(CXX) -fopenmp $(RTM_MORRICONE_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_morricone

rtm_hpc2: $(RTM_HPC_OBJS) $(HEADERS)
	icpc $^ -fopenmp $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_hpc2

.PHONY: all
all: rtm rtm_debug rtm_optim rtm_hanson rtm_morricone rtm_hpc2

.PHONY: clean
clean:
	rm -rf build
	rm -f rtm rtm_debug rtm_optim rtm_hanson rtm_morricone rtm_hpc2

build/rtm_debug/%.o: src/%.cc
	@mkdir -p build/rtm_debug
	$(CXX) -c -o $@ $< $(CXXFLAGS)

build/rtm/%.o: src/%.cc
	@mkdir -p build/rtm
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS -O3 -march=native

build/rtm_optim/%.o: src/%.cc
	@mkdir -p build/rtm_optim
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS -O3 -march=native -DNDEBUG

build/rtm_hanson/%.o: src/%.cc
	@mkdir -p build/rtm_hanson
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS -O3 -march=native

build/rtm_morricone/%.o: src/%.cc
	@mkdir -p build/rtm_morricone
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS -O3 -march=native

build/rtm_hpc2/%.o: src/%.cc
	@mkdir -p build/rtm_hpc2
	icpc -g -Ofast -xHOST -fopenmp -DUSE_THREADS -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -DNDEBUG -c -o $@ $<
