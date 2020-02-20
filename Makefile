SRC = $(wildcard src/*.cc)
HEADERS = $(wildcard src/*.h)

RTM_OBJS = $(SRC:src/%.cc=build/rtm/%.o)
RTM_DEBUG_OBJS = $(SRC:src/%.cc=build/rtm_debug/%.o)

LDFLAGS := $(LDFLAGS) -lgsl -lgslcblas -lgomp -lstdc++fs
CXXFLAGS := $(CXXFLAGS) -g

rtm_debug: $(RTM_DEBUG_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_debug


rtm: $(RTM_OBJS) $(HEADERS)
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm

.PHONY: all
all: rtm rtm_debug rtm_optim

.PHONY: clean
clean:
	rm -rf build
	rm -f rtm rtm_debug

build/rtm_debug/%.o: src/%.cc
	@mkdir -p build/rtm_debug
	$(CXX) -c -o $@ $< $(CXXFLAGS)

build/rtm/%.o: src/%.cc
	@mkdir -p build/rtm
	$(CXX) -c -o $@ $< $(CXXFLAGS) -fopenmp -DUSE_THREADS