#ifndef HEADER_flagclass_
#define HEADER_flagclass_

#include "RTM_Header.h"

using namespace std;
using std::string;

// DEFINE A CLASS FOR HANDLING THE VARIOUS REGIONAL UPDATE FLAGS

class flagclass
{
 private:
  bool *regional_update_flags;
  static unsigned int instances;
  static unsigned int iSize;
  static string* rmp_member_names;
  void flagclass_static_new();
  void new_and_set_all_flags(const bool&);
 public:
  flagclass();
  string UPItoRMP(const int&);
  bool getFlag(const string&);
  void switchFlag(const string&);
  flagclass(const int&);
  ~flagclass();
};

#endif
