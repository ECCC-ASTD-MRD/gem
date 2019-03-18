#include <rpnmacros.h>

#define   COMPRESS    1
#define UNCOMPRESS    2

#define NIL               0
#define SAMPLE            1
#define MOYENNE           2
#define MINIMUM           3
#define PARALLELOGRAM     4
#define PARALLELOGRAM32   5

#define BEST    1
#define FAST    0

#define FTN2C(i,j,ni)  (unsigned int)((ni) * (j-1) + i-1)
#define  ROUND(x)               (int)(x + 0.5)

#define min(x,y) (x < y ? x : y)
#define max(x,y) (x > y ? x : y)

typedef struct
{
  int code_methode;
  int ni,nj,nk;
  int ni_coarse, nj_coarse, nk_coarse;
  int nbits_asked, nbits_got;
  int bitstreams[4];
  int nstreams;
  int interp_degree;
  int step;
  float min,max,delta;
  int total_record_size;
  int *streams[4];
  float *coarse;
  float *fld;
  char *bzfld;
  unsigned int bzsize;
} _zfst;

typedef struct
{
#if defined(Little_Endian)
  word predictor_type:4, degree:3, step:3, nbits:5, levels: 3, version: 6, reserved3:8;
#else
  word reserved3:8, version:6, levels:3, nbits:5, step:3, degree:3, predictor_type:4;
#endif
} _fstzip;

typedef union 
  {
  INT_32 i;
  float f;
  } _floatint;


#define stuffmore(token, availableWordPtr, wordSize, bitSizeOfPackedToken, lastWordShifted, spaceInLastWord)  \
{                                                                                                             \
stuff(token, availableWordPtr, wordSize, bitSizeOfPackedToken, lastWordShifted, spaceInLastWord);             \
if (spaceInLastWord < wordSize )                                                                              \
  {                                                                                                           \
  *availableWordPtr = ( lastWordShifted << spaceInLastWord) | ( *availableWordPtr & ~(-1 << spaceInLastWord));                          \
  }                                                                                                           \
}                                                                                                             \


