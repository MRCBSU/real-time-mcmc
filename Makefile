SRC = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.h)

RTM_OBJS = $(SRC:src/%.cc=build/rtm/%.o)
RTM_OPTIM_OBJS = $(SRC:src/%.cc=build/rtm_optim/%.o)
RTM_DEBUG_OBJS = $(SRC:src/%.cc=build/rtm_debug/%.o)
RTM_PROF_DEBUG_OBJS = $(SRC:src/%.cc=build/rtm_prof_debug/%.o)
RTM_INTEL_DEBUG_OBJS = $(SRC:src/%.cc=build/rtm_intel_debug/%.o)
RTM_HANSON_OBJS = $(SRC:src/%.cc=build/rtm_hanson/%.o)
RTM_MORRICONE_OBJS = $(SRC:src/%.cc=build/rtm_morricone/%.o)
RTM_HPC_OBJS = $(SRC:src/%.cc=build/rtm_hpc2/%.o)

# LDFLAGS := $(LDFLAGS) -L/usr/local/packages/OpenBLAS/0.3.17/lib64 -L/usr/local/packages/openmp/12.0.1/lib -L/usr/local/packages/gsl/2.7/lib -lgsl -lgslcblas -lgfortran -lm -lpthread
LDFLAGS := $(LDFLAGS) -lgsl -lgslcblas -lgfortran -lm -lpthread
## LDFLAGS := $(LDFLAGS) -lgsl -lgslcblas
# CXXFLAGS := $(CXXFLAGS) -g -std=c++11 -DHAVE_INLINE -I/usr/local/packages/gsl/2.7/include -I/usr/local/packages/openmp/12.0.1/include -L/usr/local/packages/openmp/12.0.1/lib -I/usr/local/packages/OpenBLAS/0.3.17/include 
CXXFLAGS := $(CXXFLAGS) -g -std=c++11 -DHAVE_INLINE

rtm_debug: $(RTM_DEBUG_OBJS) $(HEADERS)
	g++ -g $(RTM_DEBUG_OBJS)  -I/usr/local/packages/gsl/2.7/include -L/usr/local/packages/gsl/2.7/lib -I/usr/local/packages/OpenBLAS/0.3.17/include -L/usr/local/packages/OpenBLAS/0.3.17/lib64 -fopenmp $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_debug

rtm_prof_debug: $(RTM_PROF_DEBUG_OBJS) $(HEADERS)
	g++ -g -pg $(RTM_PROF_DEBUG_OBJS) -I/usr/local/packages/gsl/2.7/include -L/usr/local/packages/gsl/2.7/lib -I/usr/local/packages/OpenBLAS/0.3.17/include -L/usr/local/packages/OpenBLAS/0.3.17/lib64 -fopenmp $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_prof_debug

rtm_intel_debug: $(RTM_INTEL_DEBUG_OBJS) $(HEADERS)
	icpc -pg -g $(RTM_INTEL_DEBUG_OBJS) -fopenmp $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_intel_debug

rtm: $(RTM_OBJS) $(HEADERS)
	$(CXX) $(RTM_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm

rtm_optim: $(RTM_OPTIM_OBJS) $(HEADERS)
	$(CXX) -fopenmp $(RTM_OPTIM_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_optim

rtm_hanson: $(RTM_HANSON_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_hanson

rtm_morricone: $(RTM_MORRICONE_OBJS) $(HEADERS)
	$(CXX) -fopenmp $(RTM_MORRICONE_OBJS) $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_morricone

rtm_hpc2: $(RTM_HPC_OBJS) $(HEADERS)
	icpc $(RTM_HPC_OBJS) -fopenmp $(LDFLAGS) $(LOADLIBES) $(LDLIBS)  -o rtm_hpc2

.PHONY: all
all: rtm rtm_debug rtm_optim rtm_hanson rtm_morricone rtm_hpc2 rtm_intel_debug rtm_prof_debug

.PHONY: clean
clean:
	rm -rf build
	rm -f rtm rtm_debug rtm_optim rtm_hanson rtm_morricone rtm_hpc2 rtm_intel_debug rtm_prof_debug

build/rtm_debug/%.o: src/%.cc
	@mkdir -p build/rtm_debug
	g++ -g -std=c++11 -L/usr/local/packages/gsl/2.7/lib -I/usr/local/packages/gsl/2.7/include -I/usr/local/packages/OpenBLAS/0.3.17/include -O0 -fopenmp -c -o $@ $< $(CXXFLAGS)

build/rtm_prof_debug/%.o: src/%.cc
	@mkdir -p build/rtm_prof_debug
	g++ -g -std=c++11 -L/usr/local/packages/gsl/2.7/lib -I/usr/local/packages/gsl/2.7/include -I/usr/local/packages/OpenBLAS/0.3.17/include -O0 -fopenmp -c -o $@ $< $(CXXFLAGS)


build/rtm_intel_debug/%.o: src/%.cc
	@mkdir -p build/rtm_intel_debug
	icpc -pg -g -std=c++11 -I/usr/local/packages/gsl/2.7/include -I/usr/local/packages/OpenBLAS/0.3.17/include -O0 -fopenmp -c -o $@ $< $(CXXFLAGS)

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
	icpc -g -std=c++11 -Ofast -fopenmp -march=native -DGSL_RANGE_CHECK_OFF -DNDEBUG -DUSE_THREADS -DHAVE_INLINE  -c -o $@ $<
	
## g++ -g -std=c++11 -Ofast -xHOST -fopenmp -DUSE_THREADS -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -DNDEBUG -c -o $@ $<
## icpc -g -std=c++11 -Ofast -I/usr/local/packages/gsl/2.7/include -I/usr/local/packages/OpenBLAS/0.3.17/include -fopenmp -DGSL_RANGE_CHECK_OFF -DNDEBUG -DUSE_THREADS -DHAVE_INLINE  -c -o $@ $<
	
