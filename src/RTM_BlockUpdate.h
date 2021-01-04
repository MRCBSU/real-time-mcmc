#ifndef HEADER_BlockUpdate_
#define HEADER_BlockUpdate_

#include "params.h"

void updateGlobalParams(const gslVector &globalParamValues, const std::vector<updateableParam> &globalParams);
void updateRegionParams(int region, const gslVector &localParamValues, const std::vector<updateableParam> &localParams);

#endif // HEADER_BlockUpdate_
