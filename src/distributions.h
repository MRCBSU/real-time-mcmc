#ifndef HEADER_distributions_
#define HEADER_distributions_

// GLOBAL ENUMERATIONS
typedef enum { cNO_TRUNC, cTRUNC, cTRUNC_LO, cTRUNC_UP} trunc_flags;
typedef enum { cCONSTANT = 1, cGAMMA, cBETA, cNORMAL, cMVNORMAL, cHALFNORMAL, cUNIFORM } distribution_type;
typedef enum { cIDENTITY, cLOG, cLOGIT } link_function;
typedef enum { cPOISSON, cNEGBIN, cBINOMIAL } data_type;

#endif
