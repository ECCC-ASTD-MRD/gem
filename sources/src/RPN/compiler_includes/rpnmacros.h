#if !defined (_RPN_MACROS_)
#define _RPN_MACROS_ _rpn_macros_

#include <rpn_macros_arch.h>

#if !defined(F2Cl)
#define F2Cl int
#endif

#if defined (f77_name)
#include <rpnmacros_global.h>

#else


#define UNIX unix

#if defined(_AIX)
#define __AIX__ 
#define AIX IBM_AIX
#endif

#define FLOAT_PACK 1
#define FLOAT_UNPACK 2

#if defined (mips) || defined (__mips)
#    define C910
#if defined (sgi) &&  (defined (__STDC__) || defined (__EXTENSIONS__))
#    define f77name(x) x##_
#else
#    define f77name(x) x/**/_
#endif
#endif
#if defined (C4680)
#    define f77name(x) x/**/_
#endif

/*  PASSE CROCHE SEEK_CUR remplace par 1 sur linux pour tell */
#if defined (i386)
#    define tell(fdesc) lseek(fdesc,0,1)
#    define f77name(x) x##_
#    define Little_Endian
#endif

#if defined (SUN)
#    define f77name(x) x##_
#endif

#if defined (__hp9000s800) || defined (__hp9000s700)
#    define HP
#    define f77name(x) x
#endif

#if defined (__AIX__)
#    define f77name(x) x/**/_
#endif

#if defined(__amd64__) || defined(__amd64) || defined(__x86_64) || defined(__x86_64__)
#    define Amd64 Amd64
#    define f77name(x) x##_
#    define tell(fdesc) lseek(fdesc,0,1)
#    define Little_Endian
#endif

#if defined(__alpha__) || defined(__alpha)
#    define Alpha Alpha
#    define f77name(x) x##_
#    define tell(fdesc) lseek(fdesc,0,1)
#    define Little_Endian
#endif

#if defined(_FLOAT0) ||  defined(_FLOAT1) ||  defined(_FLOAT2) || defined (__uxpv__)
#define NEC
#define f77name(x) x##_
#if defined (_FLOAT2)
    error " wrong float mode : must be FLOAT1 or FLOAT0"
#endif
#if defined (WORD_64)
# define NEC64 nec64
#else
# if defined (NEC64)
#  define WORD_64 word_64
# else
#  define NEC32 nec32
# endif
#endif
#endif

#if defined(_FLOAT0) || defined(Alpha) || defined (__uxpv__) || defined (__AIX__) || defined(__amd64__) || defined(__amd64) || defined(__x86_64) || defined(__x86_64__)
#    define INT_32 int
#    define INT_64 long
#else
#    define INT_32 long
#    define INT_64 long long
#endif

/* en fait toute machine dans un mode ou les pointeurs ont 64 bits */

#if defined (NEC) || defined(Alpha) || defined (__AIX__) || defined(__amd64__) || defined(__amd64) || defined(__x86_64) || defined(__x86_64__)
#define PTR_AS_INT INT_64
#else
#define PTR_AS_INT INT_32
#endif

#    define word unsigned INT_32

#if defined (WORD_64)
#    define wordint INT_64
#    define ftnword INT_64
#    define ftnfloat double
#    define bytesperword 8
#    define wordfloat double
#else
#    define wordfloat float
#    define wordint INT_32
#    define ftnword INT_32
#    define ftnfloat float
#    define bytesperword 4
#endif

#if defined(HP) || defined (i386) || defined(Alpha) || (defined (SGI) && !defined (F77_O32)) || defined (__AIX__) || defined(__amd64__) || defined(__amd64) || defined(__x86_64) || defined(__x86_64__)
#    define D77MULT bytesperword
#else
#    define D77MULT 1
#endif

#if defined(Little_Endian)
    typedef struct {
      unsigned INT_32 mantis3:29, mantis2:3, mantis1:20, expo:11, sign:1;
    }IEEE_FLOAT_DOUBLE;
    typedef struct {
      unsigned INT_32 mantis:23, expo:8 , sign:1;
    }IEEE_FLOAT;
#else
    typedef struct {
      unsigned INT_32 sign:1, expo:11 , mantis1:20,mantis2:3,mantis3:29;
    }IEEE_FLOAT_DOUBLE;
    typedef struct {
      unsigned INT_32 sign:1, expo:8 , mantis:23;
    }IEEE_FLOAT;
#endif
typedef struct {
   unsigned INT_32 sign:1, expo:7 , mantis:24;
}IBM_FLOAT;
typedef struct {
   unsigned INT_32 sign:1, expo:7 , mantis1:24 , mantis2:32;
}IBM_FLOAT_DOUBLE;
typedef union {
  float X;
  double XD;
  unsigned INT_32 U;
  IEEE_FLOAT_DOUBLE MD;
  IEEE_FLOAT M;
  IBM_FLOAT I;
  IBM_FLOAT_DOUBLE ID;
}ALL_FLOAT;





/*****************************************************************************
 *                                                                           *
 *  Objective : extract a token from a stream of packed token                *
 *  Arguments :                                                              *
 *      OUT     packedToken           token extracted                        *
 *   IN/OUT     packedWordPtr         pointer to the token to be extracted   *
 *   IN         wordSize              size of a word                         *
 *   IN         bitSizeOfPackedToken  size of a token in bit                 *
 *   IN/OUT     packedWord            word holding the desired token         *
 *   IN/OUT     bitPackInWord         no. of bits remained packed            *
 *                                    in the packedWord                      * 
 *                                                                           *
 ****************************************************************************/
#define extract( packedToken, packedWordPtr, wordSize, bitSizeOfPackedToken,                 \
                 packedWord, bitPackInWord)                                                  \
  {                                                                                          \
        /**  obtain the integer representation  **/                                          \
        if (  bitPackInWord >= bitSizeOfPackedToken )                                        \
          {                                                                                  \
                                                                                             \
            packedToken = (packedWord >> ( wordSize - bitSizeOfPackedToken ) );              \
            packedWord <<= bitSizeOfPackedToken;                                             \
            bitPackInWord -= bitSizeOfPackedToken;                                           \
          }                                                                                  \
        else                                                                                 \
          {                                                                                  \
            packedToken = (packedWord >> ( wordSize - bitSizeOfPackedToken ));               \
            packedWordPtr++;                                                                 \
            packedWord = *packedWordPtr;                                                     \
            packedToken |= ( packedWord >> ( wordSize - (bitSizeOfPackedToken-bitPackInWord)));\
            packedWord <<= ( bitSizeOfPackedToken - bitPackInWord );                         \
            bitPackInWord = wordSize - (bitSizeOfPackedToken - bitPackInWord);               \
          }; /* if */                                                                        \
          /*    handle special condition    */                                               \
        if ( bitPackInWord == 0 )                                                            \
          {                                                                                  \
            packedWordPtr++;                                                                 \
            packedWord = *packedWordPtr;                                                     \
            bitPackInWord = wordSize;                                                        \
          };/* if */                                                                         \
  }                                                                                          \

/*****************************************************************************
 *                                                                           *
 *  Objective : discard unwanted bits                                        *
 *  Arguments :                                                              *
 *   IN/OUT     packedWordPtr         pointer to the token to be extracted   *
 *   IN         wordSize              size of a word                         *
 *   IN         discardBit            bit to be discarded                    *
 *   IN/OUT     packedWord            word holding the desired token         *
 *   IN/OUT     bitPackInWord         no. of bits remained packed            *
 *                                    in the packedWord                      * 
 *                                                                           *
 ****************************************************************************/
#define discard(packedWordPtr, wordSize, discardBit, packedWord, bitPackInWord)       \
  {                                                                                   \
            if ( bitPackInWord > discardBit )                                         \
              {                                                                       \
                packedWord <<= discardBit;                                            \
                bitPackInWord -= discardBit;                                          \
              }                                                                       \
            else                                                                      \
              {                                                                       \
                packedWordPtr++;                                                      \
                packedWord = *packedWordPtr;                                          \
                packedWord <<= ( discardBit - bitPackInWord);                         \
                bitPackInWord = wordSize - ( discardBit - bitPackInWord );            \
              };                                                                      \
            /*******  special condition ********/                                     \
            if ( bitPackInWord == 0 )                                                 \
              {                                                                       \
                packedWordPtr++;                                                      \
                packedWord = *packedWordPtr;                                          \
                bitPackInWord = wordSize;                                             \
              };                                                                      \
  }                                                                                   \

/*************************************************************************************
 *                                                                                   *
 *  Objective : pack a token into a stream of packed token                           *
 *  Arguments :                                                                      *
 *   IN         token                 token to be packed                             *
 *   IN/OUT     availableWordPtr      pointer to the word position available for     *
 *                                    packing                                        *
 *   IN         wordSize              size of a word                                 *
 *   IN         bitSizeOfPackedToken  bit size of the token                          *
 *   IN/OUT     lastWordShifted       word designated to hold the token              *
 *   IN/OUT     spaceInLastWord       no. of bits available for packing in last word *
 *                                                                                   *
 ************************************************************************************/
#define stuff(token, availableWordPtr, wordSize, bitSizeOfPackedToken, lastWordShifted,    \
              spaceInLastWord)                                                             \
  {                                                                                        \
              if ( spaceInLastWord >= bitSizeOfPackedToken )                               \
                /* integer token fits into space left */                                   \
                {                                                                          \
                  lastWordShifted = ( lastWordShifted << bitSizeOfPackedToken ) | token;   \
                  spaceInLastWord = spaceInLastWord - bitSizeOfPackedToken;                \
                }                                                                          \
              else                                                                         \
                /* integer token can not fit into space left */                            \
                {                                                                          \
                  *availableWordPtr = ((lastWordShifted << spaceInLastWord) |              \
                                             ( token >>                                    \
                                               ( bitSizeOfPackedToken - spaceInLastWord)));\
                  lastWordShifted = token &                                                \
                    ( -1 >> ( wordSize - ( bitSizeOfPackedToken - spaceInLastWord)));      \
                  spaceInLastWord = wordSize - ( bitSizeOfPackedToken - spaceInLastWord ); \
                  availableWordPtr++;                                                      \
                  lastSlot++;                                                              \
                }; /* if */                                                                \
  }
#define swap_word_endianness(mot) { register unsigned INT_32 tmp =(unsigned INT_32)mot; \
   mot = (tmp>>24) | (tmp<<24) | ((tmp>>8)&0xFF00) | ((tmp&0xFF00)<<8); }
#define swap_buffer_endianness(buff,nwds) {\
    unsigned INT_32 *buf=(unsigned INT_32 *)buff ;\
   register INT_32 nwords=nwds ;\
   while(nwords--) { swap_word_endianness(*buf);buf++; };\
   }

#endif

#endif
