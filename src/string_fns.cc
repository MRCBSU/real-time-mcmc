#include <iomanip>

#include "string_fns.h"

using namespace std;
using std::string;

// Read in a delimited character from a string.
template <typename D>
D read_from_delim_string(const string& str_in, const char* str_delim, int& indx)
{
  int int_temp = str_in.find_first_of(str_delim, indx);

  string temp_string = str_in.substr(indx, int_temp - indx);

  indx = int_temp + 1;

  D out_D;
  stringstream outStream(TrimStr(temp_string));
  outStream >> out_D;
  return out_D;

}

template string read_from_delim_string<string>(const string&, const char*, int&);
template int read_from_delim_string<int>(const string&, const char*, int&);
template double read_from_delim_string<double>(const string&, const char*, int&);

int count_delims_in_string(const string in_string, const char *delims)
{
  int int_marker = 0, counter = 0;
  for(counter = 0; int_marker != string::npos; )
      int_marker = in_string.find_first_of(delims, int_marker + 1 - ((int) counter++ == 0));
  return --counter;
}


int count_instances_in_string(const string in_string, const char *sep_string)
{

  int cat_appearances = 0;
  int i_start, i_end;

  for(int i = in_string.find(sep_string, 0); i != string::npos; i = in_string.find(sep_string, i))
    {
      // IS THIS INSTANCE A WORD? OR A SUB-WORD?
      i_start = (i == 0) ? -1 : in_string.find_first_of(STRING_DELIMS, i - 1);
      i_end = in_string.find_first_of(STRING_DELIMS, i + strlen(sep_string));
      if((i == (i_start + 1)) && ((i + strlen(sep_string)) == i_end)) // THEN WE HAVE A MATCH FOR THE WORD
	cat_appearances++;
      i++;
    }

  return cat_appearances;

}

int nth_instance_in_string(const string in_string, const char *sep_string, const int n_instance){

  // RETURNS THE POSITION OF THE N^TH INSTANCE OF SEP_STRING AS A WORD IN THE INPUT STRING IN_STRING
  if(n_instance <= 0){

    // INVALID NUMBER OF INSTANCES - RETURN A NULL VALUE
    return -1;

  } else {

    // COUNT THE TOTAL NUMBER OF INSTANCES
    int tot_appearances = count_instances_in_string(in_string, sep_string);

    // CHECK THERE ARE SUFFICIENT INSTANCES IN THE FILE
    if(n_instance > tot_appearances)
      { // ERROR
	return -1;
      }

    // ELSE RETURN A VALUE
    int num_instance = 0;
    int i = 0, i_start, i_end;
    for(; num_instance < n_instance; ++i)
      {
	i = in_string.find(sep_string, i);
	i_start = (i == 0) ? -1 : in_string.find_first_of(STRING_DELIMS, i - 1);
	i_end = in_string.find_first_of(STRING_DELIMS, i + strlen(sep_string));
	if((i == (i_start + 1)) && ((i + strlen(sep_string)) == i_end)) // THEN WE HAVE A MATCH FOR THE WORD
	  num_instance++;
      }

    return --i;

  }

}

// TRIM WHITE SPACE AT BEGINNING AND END OF A STRING
string TrimStr(const string& Src, const string& c)
{
  int p2 = Src.find_last_not_of(c);
  if(p2 == string::npos) return string();
  int p1 = Src.find_first_not_of(c);
  if(p1 == string::npos) p1 = 0;
  return Src.substr(p1, (p2 - p1) + 1);
}

string int_to_string(const int n)
{
  string out_string;
  stringstream i_to_str;
  i_to_str << n;
  i_to_str >> out_string;
  return out_string;
}

string dbl_to_string(const double x)
{
  string out_string;
  stringstream i_to_str;
  i_to_str << std::fixed << x;
  i_to_str >> out_string;
  return out_string;
}
