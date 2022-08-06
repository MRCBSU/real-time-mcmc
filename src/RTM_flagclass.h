#ifndef HEADER_flagclass_
#define HEADER_flagclass_

#include "RTM_Header.h"

using namespace std;
using std::string;

// DEFINE A CLASS FOR HANDLING THE VARIOUS REGIONAL UPDATE FLAGS

class flagclass
{
 private:
  static unsigned int instances;
  static unsigned int iSize;
  static string* rmp_member_names;
  void flagclass_static_new();
  void new_and_set_all_flags(const bool&);
public:
  // CCS: Temporarily public for testing block updates
  bool *regional_update_flags;

  flagclass();
  string UPItoRMP(const int&);
  bool getFlag(const string&);
  unsigned int getSize();
  void switchFlag(const string&);
  flagclass(const int&);
  ~flagclass();
};

#endif
