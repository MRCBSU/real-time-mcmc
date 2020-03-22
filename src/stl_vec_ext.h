// FUNCTIONS IN stl_vector_exts.cc

#ifndef HEADER_stlvec_
#define HEADER_stlvec_

#include <string>
#include <vector>

using namespace std;
using std::string;

template <typename D>
void stl_vector_sscanf(const string &, vector<D> &,
                       const char *chr_delim = ",;");

#endif
