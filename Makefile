SRC = $(wildcard *.cc)
OBJ = $(SRC:%.cc:build/%.o)
rtm_gnu: CXXFLAGS := $(CXXFLAGS) -fopenmp


LDFLAGS := $(LDFLAGS) -lgsl -lgslcblas -lgomp -lstdc++fs



rtm_gnu_debug	: build/RTM_RegRTM_Beta3.0DEBUG.o build/RTM_Inputs.o build/RTM_StructAllocFree.o build/RTM_MixingMatrices.o build/RTM_WithinRegion.o build/RTM_LikelihoodsDEBUG.o build/RTM_MetropHast.o build/RTM_flagclass.o build/RTM_MvNorm.o build/gp_param_patch.o build/R_like_fns.o build/file_string_fns.o build/string_fns.o build/gsl_vector_exts.o build/gsl_matrix_exts.o src/distributions.h src/gsl_mat_ext.h src/gsl_vec_ext.h src/R_like_fns.h src/RTM_flagclass.h src/RTM_FunctDefs.h src/RTM_Header.h src/RTM_StructAssign.h src/RTM_StructDefs.h src/string_fns.h src/RTM_modelstate.h build/RTM_modelstate.o build/rtm_data_class.o src/rtm_data_class.h src/stl_vec_ext.h build/stl_vector_exts.o
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_gnu_debug


rtm_gnu	: build/RTM_RegRTM_Beta3.0.o build/RTM_Inputs.o build/RTM_StructAllocFree.o build/RTM_MixingMatrices.o build/RTM_WithinRegion.o build/RTM_Likelihoods.o build/RTM_MetropHast.o build/RTM_flagclass.o build/RTM_MvNorm.o build/gp_param_patch.o build/R_like_fns.o build/file_string_fns.o build/string_fns.o build/gsl_vector_exts.o build/gsl_matrix_exts.o src/distributions.h src/gsl_mat_ext.h src/gsl_vec_ext.h src/R_like_fns.h src/RTM_flagclass.h src/RTM_FunctDefs.h src/RTM_Header.h src/RTM_StructAssign.h src/RTM_StructDefs.h src/string_fns.h src/RTM_modelstate.h build/RTM_modelstate.o build/rtm_data_class.o src/rtm_data_class.h src/stl_vec_ext.h build/stl_vector_exts.o
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o rtm_gnu


.PHONY: all
all: rtm_gnu_debug rtm_gnu rtm_gnu_empirical rtm_gnu_empirical_debug

.PHONY: clean
clean:
	rm -rf build
	rm -f rtm_gnu rtm_gnu_debug rtm_gnu_empirical rtm_gnu_empirical_debug

build/%.o: src/%.cc
	@mkdir -p build
	$(CXX) -c -o $@ $< $(CXXFLAGS)

