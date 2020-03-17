// LOADING IN USER-DEFINED HEADER FILE, CONTAINING ALL STANDARD HEADER FILES REQUIRED
#include "stl_vec_ext.h"
#include "string_fns.h"

template <typename D>
void stl_vector_sscanf(const string& in_string, vector<D>& out_vec, const char* chr_delim)
{

  // How long should the vector be?
  int int_vec_length = count_delims_in_string(in_string, ",") + 1;

  // Read in the new vector, element by element and add to the end of out_vec
  int indx = 0;
  for(int i = 0; i < int_vec_length; i++)
    out_vec.push_back(read_from_delim_string<D>(in_string, chr_delim, indx));
}

template void stl_vector_sscanf<int>(const string&, vector<int>&, const char*);
template void stl_vector_sscanf<double>(const string&, vector<double>&, const char*);
