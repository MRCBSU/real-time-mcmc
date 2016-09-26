#ifndef HEADER_string_fn_
#define HEADER_string_fn_

// WHEN I CAN FIGURE OUT THE CORRECT FILES TO INCLUDE, THE BELOW INCLUDE STATEMENT SHOULD BE REMOVED AND REPLACES WITH STANDARD C++ LIBRARIES, AND THIS SHOULD THEN BE PORTABLE CODE
#include "RTM_Header.h"

using namespace std;
using std::string;

#define STRING_DELIMS " ,\t\n\r"

// DECLARATION OF FUNCTIONS IN string_fns.cc
template <typename D>D read_from_delim_string(const string&, const char*, int&);
int count_delims_in_string(const string, const char* delims = ",;");
int count_instances_in_string(const string, const char*);
int nth_instance_in_string(const string, const char*, const int);
string TrimStr(const string& Src, const string& c = STRING_DELIMS);
string int_to_string(const int n);
string dbl_to_string(const double);

// DECLARATION OF FUNCTIONS IN file_string_fns.cc
void fn_load_file(string*, const char*);
int count_instances_in_file(const char *, const char *);
void nth_instance(char *&, const char *, const char *, const int);
void read_double_input(double&, const char*, const char*);
void read_int_input(int&, const char*, const char*);
void read_char_input(char*, const char*, const char*, const int);

#endif
