#ifndef HEADER_BlockUpdate_
#define HEADER_BlockUpdate_

#include "RTM_updParams.h"



//extern paramBlock globalParams;
//extern std::vector<paramBlock> localParams;


void updateGlobalParams(const gslVector &globalParamValues, const std::vector<updParam> &globalParams);
void updateRegionParams(int region, const gslVector &localParamValues, const std::vector<updParam> &localParams);

#endif // HEADER_BlockUpdate_
