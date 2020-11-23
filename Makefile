SRC = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.h)

RTM_OBJS = $(SRC:src/%.cc=build/rtm/%.o)
RTM_OPTIM_OBJS = $(SRC:src/%.cc=build/rtm_optim/%.o)
RTM_DEBUG_OBJS = $(SRC:src/%.cc=build/rtm_debug/%.o)

LDFLAGS := $(LDFLAGS) -lgsl -lgslcblas -lgomp -L/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/gsl-2.4-z4fspadciir53pifcikqs3zbce67d7ru/lib -lm
CXXFLAGS := $(CXXFLAGS) -g -DHAVE_INLINE -I/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/gsl-2.4-z4fspadciir53pifcikqs3zbce67d7ru/include

rtm_debug: $(RTM_DEBUG_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_debug


rtm: $(RTM_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm

rtm_optim: $(RTM_OPTIM_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_optim

.PHONY: all
all: rtm rtm_debug rtm_optim

.PHONY: clean
clean:
	rm -rf build
	rm -f rtm rtm_debug rtm_optim

build/rtm_debug/%.o: src/%.cc
	@mkdir -p build/rtm_debug
	$(CXX) -c -o $@ $< $(CXXFLAGS)

build/rtm/%.o: src/%.cc
	@mkdir -p build/rtm
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS -O3 -march=native

build/rtm_optim/%.o: src/%.cc
	@mkdir -p build/rtm_optim
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS -O3 -march=native
