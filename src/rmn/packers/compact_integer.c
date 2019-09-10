/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <rmnlib.h>




/***********************************************************
 *                                                         *
 *  Objecttive : find min and max of an integer array      *
 *                                                         *
 *  Argument   :                                           *
 *   IN    arrayOfUnpacked    integer array                *
 *   OUT   min                minimun integer              *
 *   OUT   max                maximun integer              *
 *                                                         *
 **********************************************************/
#define FindMinMax(arrayOfUnpacked, min, max)             \
{                                                         \
                                                          \
  min = arrayOfUnpacked[0];                               \
  max = arrayOfUnpacked[0];                               \
                                                          \
  for(i=stride; i < intCount*stride ; i+=stride)          \
    {                                                     \
      if ( arrayOfUnpacked[i] < min )                     \
        {                                                 \
          min = arrayOfUnpacked[i];                       \
        }                                                 \
      else if ( arrayOfUnpacked[i] > max )                \
        {                                                 \
          max = arrayOfUnpacked[i];                       \
        };                                                \
    };                                                    \
}                                                         \


/*********************************************************************
 *                                                                   *
 *  Objective : obtain and stuff a chain of 32 bit word into an array*
 *              ( macro tobe used in another macro "pack" )          *
 *                                                                   *
 *  Argument  :                                                      *
 *   IN   arrayOfUnpacked       integer array                        *
 *   IN   min                   minimun integer                      *
 *                                                                   *
 ********************************************************************/
#define pack32Bit(arrayOfUnpacked, min)                                  \
{                                                                        \
          i = intCount;                                                  \
          if ( packHeader != NULL )                                      \
            {                                                            \
             while ( i-- )                                               \
              {                                                          \
               *arrayPtr = ((*arrayOfUnpacked - min ) >> shiftRequired );\
               arrayPtr++;                                               \
               arrayOfUnpacked+=stride;                                  \
              };                                                         \
            }                                                            \
          else                                                           \
            {                                                            \
             while ( i-- )                                               \
              {                                                          \
               *arrayPtr = *arrayOfUnpacked + positiveMask;              \
               arrayPtr++;                                               \
               arrayOfUnpacked+=stride;                                  \
              };                                                         \
            };                                                           \
};                                                                       \


/*************************************************************************
 *                                                                       *
 *  Objective : obtain and stuff a chain of word(<32 bit) into an array  *
 *              ( tobe used in another macro "pack" )                    *
 *                                                                       *
 *  Argument  :                                                          *
 *   IN   arrayOfUnpacked       integer array                            *
 *   IN   min                   minimun integer                          *
 *                                                                       *
 ************************************************************************/
#define packBit(arrayOfUnpacked, min)                                    \
        {                                                                \
          if (packHeader !=NULL)                                         \
            {                                                            \
             for ( i = 0; i < intCount*stride; i+=stride)                \
              {                                                          \
               tempInt = ((arrayOfUnpacked[i]- min ) >> shiftRequired ); \
               tempInt &= cleanupMask;                                   \
               stuff(tempInt, arrayPtr, wordSize, bitSizeOfPackedToken,  \
                     lastWordShifted, spaceInLastWord);                  \
              }; /* for */                                               \
            }                                                            \
          else                                                           \
            {                                                            \
              for ( i = 0; i < intCount*stride; i+=stride)               \
              {                                                          \
              tempInt = arrayOfUnpacked[i] + positiveMask;               \
              tempInt &= cleanupMask;                                    \
              stuff(tempInt, arrayPtr, wordSize, bitSizeOfPackedToken,   \
                    lastWordShifted, spaceInLastWord);                   \
              };                                                         \
            };                                                           \
        };                                                               \

/*********************************************************************
 *                                                                   *
 *  Objective : pack a chain of word into an array                   *
 *                                                                   *
 *  Argument  :                                                      *
 *   IN   arrayOfUnpacked       integer array                        *
 *   IN   min                   minimun integer                      *
 *                                                                   *
 ********************************************************************/
#define Pack(arrayOfUnpacked, min)                                       \
{                                                                        \
     /* initialize local variables*/                                     \
                                                                         \
      lastPackBit = off_set;                                             \
      spaceInLastWord =  wordSize - ( lastPackBit % wordSize );          \
      lastSlot = ( lastPackBit / wordSize );                             \
                                                                         \
      if ( spaceInLastWord == wordSize )                                 \
        {                                                                \
          lastWordShifted = 0;                                           \
        }                                                                \
      else                                                               \
        {                                                                \
          lastWordShifted = arrayOfPacked[lastSlot] >> spaceInLastWord ; \
        };                                                               \
                                                                         \
      arrayPtr = &arrayOfPacked[lastSlot];                               \
      if (( spaceInLastWord == wordSize ) &&                             \
          ( bitSizeOfPackedToken == wordSize ))                          \
        {                                                                \
          pack32Bit(arrayOfUnpacked, min);                               \
        }                                                                \
      else                                                               \
        {                                                                \
          packBit(arrayOfUnpacked, min);                                 \
        };                                                               \
                                                                         \
       /* squeezes hole left in the integer array */                     \
                                                                         \
      if ( spaceInLastWord < wordSize )                                  \
        {                                                                \
         *arrayPtr = ( lastWordShifted << spaceInLastWord) |             \
                     ( *arrayPtr & ~(-1 << spaceInLastWord));            \
         };                                                              \
                                                                         \
                                                                         \
      return bitSizeOfPackedToken;                                       \
}                                                                        \

/*******************************************************************
 *                                                                 *
 *  Objective : construct pack header                              *
 *                                                                 *
 *  Argument  :                                                    *
 *    IN     arrayOfUnpacked        integer array                  *
 *    IN     min                    mininum integer                *
 *    IN     max                    maximum integer                *
 *                                                                 *
 ******************************************************************/
#define constructHeader(arrayOfUnpacked, min, max)                       \
 {                                                                       \
  word tempUnsignedMax, tempUnsignedMin;                                 \
                                                                         \
  /*      obtain min, max & range         */                             \
                                                                         \
   FindMinMax(arrayOfUnpacked, min, max);                                \
   if ( (max > 0) && (min < 0) )                                         \
     /* prevent signed overflow */                                       \
     {                                                                   \
       tempUnsignedMax = max;                                            \
       tempUnsignedMin = -min;                                           \
       maxRange = tempUnsignedMax + tempUnsignedMin;                     \
     }                                                                   \
   else                                                                  \
     maxRange = max - min;                                               \
                                                                         \
   /*  compute shift required */                                         \
   bitRequiredForRange = 0;                                              \
   while ( maxRange != 0 )                                               \
     {                                                                   \
       maxRange = maxRange >> 1;                                         \
       bitRequiredForRange++;                                            \
     };                                                                  \
                                                                         \
   shiftRequired = 0;                                                    \
   if ( bitSizeOfPackedToken == -1 )                                     \
     {                                                                   \
       bitSizeOfPackedToken = bitRequiredForRange;                       \
     }                                                                   \
   else                                                                  \
     {                                                                   \
       while ( (bitRequiredForRange-bitSizeOfPackedToken) > 0 )          \
         {                                                               \
           shiftRequired++;                                              \
           bitRequiredForRange--;                                        \
         };                                                              \
     };                                                                  \
                                                                         \
   /*   construct pack header */                                         \
   packHeader[0] = 0xFD000000 | (shiftRequired <<6)|bitSizeOfPackedToken;\
   packHeader[1] = intCount;                                             \
   packHeader[2] = min;                                                  \
   packHeader[3] = max;                                                  \
                                                                         \
 };                                                                      \
        
/*******************************************************************
 *                                                                 *
 *  Objective : unpack an array                                    *
 *                                                                 *
 *  Argument  :                                                    *
 *    IN     arrayOfUnpacked        integer array                  *
 *    IN     requiredShift          shift required                 *
 *    IN     bitSizeOfPackedToken   pack token size                *
 *    IN     min                    mininum integer                *
 *    IN     intCount               element in the packed array    *
 *                                                                 *
 ******************************************************************/
#define Unpack(arrayOfUnpacked, requiredShift, tokenSize, min, intCount)                  \
 {                                                                                        \
                                                                                          \
    /*    initialize variables         */                                                 \
    firstPackBit = off_set;                                                               \
    bitPackInFirstWord =  wordSize - ( firstPackBit % wordSize );                         \
    currentSlot = ( firstPackBit / wordSize );                                            \
    currentWord = arrayOfPacked[currentSlot] << ( wordSize - bitPackInFirstWord );        \
    positiveMask = -positiveMask;                                                         \
                                                                                          \
    if ( tokenSize > wordSize )                                                           \
      {                                                                                   \
        significantBit = wordSize;                                                        \
        inSignificantBit = tokenSize - wordSize;                                          \
      }                                                                                   \
    else                                                                                  \
      {                                                                                   \
        significantBit = tokenSize;                                                       \
        inSignificantBit = 0;                                                             \
      };                                                                                  \
                                                                                          \
    /**    unpack integer numbers             */                                          \
    arrayPtr = &arrayOfPacked[currentSlot];                                               \
    if ( packHeader != NULL )                                                             \
      {                                                                                   \
        for ( i = 0; i < intCount*stride; i+=stride)                                      \
          {                                                                               \
            extract(packInt, arrayPtr, wordSize, significantBit,                          \
                    currentWord, bitPackInFirstWord);                                     \
                                                                                          \
            /*      truncate extra bit        */                                          \
            if ( inSignificantBit > 0 )                                                   \
              {                                                                           \
                discard(arrayPtr, wordSize, inSignificantBit,                             \
                        currentWord, bitPackInFirstWord);                                 \
              }                                                                           \
            arrayOfUnpacked[i] =  (packInt << requiredShift ) + min;                      \
          }                                                                               \
      }                                                                                   \
    else                                                                                  \
      {                                                                                   \
        for ( i = 0; i < intCount*stride; i+=stride)                                      \
          {                                                                               \
            extract(packInt, arrayPtr, wordSize, significantBit,                          \
                    currentWord, bitPackInFirstWord);                                     \
                                                                                          \
            /*      truncate extra bit if necessary           */                          \
            if ( inSignificantBit > 0 )                                                   \
              {                                                                           \
                discard(arrayPtr, wordSize, inSignificantBit, currentWord,                \
                        bitPackInFirstWord);                                              \
              }                                                                           \
                                                                                          \
            arrayOfUnpacked[i] =  packInt + positiveMask;                                 \
          }                                                                               \
      };                                                                                  \
                                                                                          \
    return bitSizeOfPackedToken;                                                          \
 };                                                                                       \

/***********************************************************************************************
 *                                                                                             *
 * Author   : Jianhui He, 1997                                                                 *
 *                                                                                             *
 * Objective: transform integer array between pack and unpacked format                         *
 *                                                                                             *
 * Arguments:                                                                                  *
 *    IN/OUT  unpackedArrayOfInt     unpacked integer array                                    *
 *    IN/OUT  packedArrayOfInt       packed integer array                                      *
 *    IN      elementCount           total count of element in unpacked integer array          *
 *    IN      bitSizeOfPackedToken   packed integer size in bit                                *
 *                                   when equal "-1", autodetection of bit length is on        *
 *    IN      offset                 in packing   : the last bit of integer packed inside array*
 *                                   in unpacking : the first bit of integer packed inside array*
 *    IN      stride                 unpacked integer spacing indicator                        *
 *    IN      opCode                 1: unsigned pack                                          *
 *                                   2: unsigned unpack                                        *
 *                                   3: signed pack                                            *
 *                                   4: signed unpack                                          *
 **********************************************************************************************/
int  compact_integer( void *unpackedArrayOfInt, void *packedHeader, void *packedArrayOfInt, 
                       int elementCount, int bitSizeOfPackedToken, int off_set, 
                       int stride, int opCode)
{

    typedef struct
    {
#if defined(Little_Endian)
      word numOfBitsPerToken : 6, SHIFT : 6, unused : 12, ID : 8; 
#else
      word ID : 8, unused : 12, SHIFT : 6, numOfBitsPerToken : 6;
#endif
      word numOfPackedToken : 32;
      word minValue         : 32;
      word maxValue         : 32;
    }integer_header;

    int wordSize;
    word cleanupMask;
    int minSignedInteger=0, maxSignedInteger=0;
    word minUnsignedInteger=0, maxUnsignedInteger=0;
    word maxRange;
    word maxSpan;
   

    int positiveMask;
    word *arrayOfUnsignedUnpacked;
    int  *arrayOfSignedUnpacked;
    word *arrayOfPacked;
    int i, k;
    int intCount;
    int bitRequiredForRange, shiftRequired = 0;

    

    /****************************************
     *                                      *
     *     variables used by the packer     *
     *                                      *
     ***************************************/
    int lastPackBit, spaceInLastWord, lastSlot;
    word lastWordShifted, tempInt;
    word *packHeader;
    
    /***************************************
     *                                     *
     *    variables used by the unpacker   *
     *                                     *
     **************************************/
    int firstPackBit, bitPackInFirstWord, currentSlot;
    word currentWord, packInt;
    int significantBit, inSignificantBit;
    word *arrayPtr;
    int tokenSize, ShiftIntended, elementCountFromHeader;
    int  minSigned;
    word minUnsigned;
    integer_header *theHeader;


/*  printf("minSignedInteger=%d minUnsignedInteger=%d \n",minSignedInteger,minUnsignedInteger); */

    /********************************
     *                              *
     *   handle abnormal condition  *
     *                              *
     ********************************/
    /* token size is 0 */
    if ( bitSizeOfPackedToken == 0 )
      { 
        return 0;
      };  





   
    /********************************************************
     *                                                      *
     *    determine wordsize and others                     * 
     *                                                      *
     ********************************************************/
    wordSize                 = 8 * sizeof(word);
    arrayOfUnsignedUnpacked  = (word *)unpackedArrayOfInt;
    arrayOfSignedUnpacked    = (int *)unpackedArrayOfInt;
    theHeader                = (integer_header *)packedHeader;
    packHeader               = (word *)packedHeader;
    arrayOfPacked            = (word  *)packedArrayOfInt;
    intCount                 = elementCount;
    cleanupMask              = ((word)(~0)>>(wordSize-bitSizeOfPackedToken));
      


   if ( (opCode==1) || (opCode==3) )
    /************************************************
     *                                              *
     *         collect info for the packing         *
     *                                              *
     ***********************************************/   
   {
     if ( packedHeader != NULL )
       /*******************************************************************
        *                                                                  *
        *  pack header is required, (X - Xmin) is used as packInt          *
        *                                                                  *
        *******************************************************************/
       {
         if ( opCode == 1 )
           { 
             constructHeader(arrayOfUnsignedUnpacked, minUnsignedInteger, maxUnsignedInteger);
           }
         else
           {
             constructHeader(arrayOfSignedUnpacked, minSignedInteger, maxSignedInteger); 
           }
      }
    else
      /*************************************************************
       *                                                           *
       *   pack header not required, X itself is used as packInt,  *
       *   determines bitSizeOfPackedToken, if not available       *
       *                                                           *
       ************************************************************/
      {
        if ( bitSizeOfPackedToken == -1 )
          {
            /**********************************************************
             *                                                        *
             *   obtain minimum, maximun, span                        *
             *                                                        *
             *********************************************************/
            if ( opCode == 3 )
              {
                /* signed integer number */
               
                FindMinMax(arrayOfSignedUnpacked, minSignedInteger, maxSignedInteger);

                maxSpan    = ( abs(minSignedInteger) > maxSignedInteger ) ? abs(minSignedInteger) : 
                             maxSignedInteger;            
              }
            else if (  opCode == 1 )
              {
                /* unsigned integer number */
                maxSpan = arrayOfUnsignedUnpacked[0];

                for(i=stride; i < intCount*stride ; i+=stride)
                  {
                    maxSpan |= arrayOfUnsignedUnpacked[i];
                  }; 
              };
          
            /************************************************************
             *                                                          *
             *           derive bitSizeOfPackedToken                    *
             *                                                          *
             ***********************************************************/
            bitSizeOfPackedToken = 0;
            while ( maxSpan != 0 )
              {
                maxSpan = maxSpan >> 1;
                bitSizeOfPackedToken++;
              };
            if ( opCode == 3 )
              {/* accomodate the signed bit */
                bitSizeOfPackedToken++;
              };
            cleanupMask = ((word)(~0)>>(wordSize-bitSizeOfPackedToken));
          };/* if */



      };/* else */
   } 
   else/* opCode == 2 or 4 */
    /************************************************
     *                                              *
     *         collect info for the unpacking       *
     *                                              *
     ***********************************************/
     {
       if ( packHeader != NULL )
         {
           theHeader     = (integer_header *)packedHeader;
           tokenSize     = theHeader->numOfBitsPerToken; 
           ShiftIntended = theHeader->SHIFT;
           intCount      = theHeader->numOfPackedToken;
           minSigned     = theHeader->minValue;
           minUnsigned   = theHeader->minValue;

         }
       else
         {
           tokenSize     = bitSizeOfPackedToken;
           ShiftIntended = 0;
           intCount      = intCount;
           minSigned     = minSignedInteger;
           minUnsigned   = minUnsignedInteger;
         };
     };
 
   

   /**********************************************
    *                                            *
    *  compute signed int adjustment,            *
    * since 1,2,3,4 without header all needs it  *
    *                                            *
    *********************************************/
   positiveMask = ( opCode < 3 ) ? 0 : ( 1 << ( bitSizeOfPackedToken - 1 ));


   /***********************************************
    *                                             *
    *   pack                                      *
    *                                             *
    **********************************************/
   if ( opCode == 1 )
     {
       Pack(arrayOfUnsignedUnpacked, minUnsignedInteger); 
     }
   else if ( opCode == 3 )
     {
       Pack(arrayOfSignedUnpacked, minSignedInteger);     
     }
   /***********************************************
    *                                             *
    *   unpack                                    *
    *                                             *
    **********************************************/
   else if ( opCode == 2 )
     {
       Unpack(arrayOfUnsignedUnpacked, ShiftIntended, tokenSize, 
              minUnsigned, intCount);
     }
   else if ( opCode == 4 )
     {
       Unpack(arrayOfSignedUnpacked, ShiftIntended, tokenSize, 
              minSigned, intCount);
     }
   else
     {
       printf("\n opCode:%d is not defined \n", opCode);
       return 0;
     };/* if */
 
  return intCount;  /* unused, function must return something */
 

} /* end compact_integer */

/***********************************************************************************************
 *                                                                                             *
 * Author   : M. Lepine, sept 2005                                                             *
 *                                                                                             *
 * Objective: transform array of short integer between pack and unpacked format                *
 *                                                                                             *
 * Arguments:                                                                                  *
 *    IN/OUT  unpackedArrayOfShort   unpacked array of short integer                           *
 *    IN/OUT  packedArrayOfInt       packed integer array                                      *
 *    IN      elementCount           total count of element in unpacked integer array          *
 *    IN      bitSizeOfPackedToken   packed integer size in bit                                *
 *                                   when equal "-1", autodetection of bit length is on        *
 *    IN      offset                 in packing   : the last bit of integer packed inside array*
 *                                   in unpacking : the first bit of integer packed inside array*
 *    IN      stride                 unpacked integer spacing indicator                        *
 *    IN      opCode                 5: unsigned short pack                                    *
 *                                   6: unsigned short unpack                                  *
 *                                   7: signed short pack    (not implemented)                 *
 *                                   8: signed short unpack  (not implemented)                 *
 **********************************************************************************************/
int  compact_short( void *unpackedArrayOfShort, void *packedHeader, void *packedArrayOfInt, 
                       int elementCount, int bitSizeOfPackedToken, int off_set, 
                       int stride, int opCode)
{

    typedef struct
    {
#if defined(Little_Endian)
      word numOfBitsPerToken : 6, SHIFT : 6, unused : 12, ID : 8; 
#else
      word ID : 8, unused : 12, SHIFT : 6, numOfBitsPerToken : 6;
#endif
      word numOfPackedToken : 32;
      word minValue         : 32;
      word maxValue         : 32;
    }integer_header;

    int wordSize;
    word cleanupMask;
    int minSignedInteger=0, maxSignedInteger=0;
    word minUnsignedInteger=0, maxUnsignedInteger=0;
    word maxRange;
    word maxSpan;
   

    int positiveMask;
    unsigned short *arrayOfUnsignedShort;
    word *arrayOfPacked;
    int i, k;
    int intCount;
    int bitRequiredForRange, shiftRequired = 0;

    

    /****************************************
     *                                      *
     *     variables used by the packer     *
     *                                      *
     ***************************************/
    int lastPackBit, spaceInLastWord, lastSlot;
    word lastWordShifted, tempInt;
    word *packHeader;
    
    /***************************************
     *                                     *
     *    variables used by the unpacker   *
     *                                     *
     **************************************/
    int firstPackBit, bitPackInFirstWord, currentSlot;
    word currentWord, packInt;
    int significantBit, inSignificantBit;
    word *arrayPtr;
    int tokenSize, ShiftIntended, elementCountFromHeader;
    int  minSigned;
    word minUnsigned;
    integer_header *theHeader;


/*  printf("minSignedInteger=%d minUnsignedInteger=%d \n",minSignedInteger,minUnsignedInteger); */

    /********************************
     *                              *
     *   handle abnormal condition  *
     *                              *
     ********************************/
    /* token size is 0 */
    if ( bitSizeOfPackedToken == 0 )
      { 
        return 0;
      };  





   
    /********************************************************
     *                                                      *
     *    determine wordsize and others                     * 
     *                                                      *
     ********************************************************/
    wordSize                 = 8 * sizeof(word);
    arrayOfUnsignedShort     = (short *)unpackedArrayOfShort;
    theHeader                = (integer_header *)packedHeader;
    packHeader               = (word *)packedHeader;
    arrayOfPacked            = (word  *)packedArrayOfInt;
    intCount                 = elementCount;
    cleanupMask              = ((word)(~0)>>(wordSize-bitSizeOfPackedToken));
      


   if (opCode==5)
    /************************************************
     *                                              *
     *         collect info for the packing         *
     *                                              *
     ***********************************************/   
   {
     if ( packedHeader != NULL )
       /*******************************************************************
        *                                                                  *
        *  pack header is required, (X - Xmin) is used as packInt          *
        *                                                                  *
        *******************************************************************/
       {
             constructHeader(arrayOfUnsignedShort, minUnsignedInteger, maxUnsignedInteger);
       }
    else
      /*************************************************************
       *                                                           *
       *   pack header not required, X itself is used as packInt,  *
       *   determines bitSizeOfPackedToken, if not available       *
       *                                                           *
       ************************************************************/
      {
        if ( bitSizeOfPackedToken == -1 )
          {
            /**********************************************************
             *                                                        *
             *   obtain minimum, maximun, span                        *
             *                                                        *
             *********************************************************/
            if ( opCode == 5 )
              {
                /* unsigned integer number */
                maxSpan = arrayOfUnsignedShort[0];

                for(i=stride; i < intCount*stride ; i+=stride)
                  {
                    maxSpan |= arrayOfUnsignedShort[i];
                  }; 
              };
          
            /************************************************************
             *                                                          *
             *           derive bitSizeOfPackedToken                    *
             *                                                          *
             ***********************************************************/
            bitSizeOfPackedToken = 0;
            while ( maxSpan != 0 )
              {
                maxSpan = maxSpan >> 1;
                bitSizeOfPackedToken++;
              };
            if ( opCode == 3 )
              {/* accomodate the signed bit */
                bitSizeOfPackedToken++;
              };
            cleanupMask = ((word)(~0)>>(wordSize-bitSizeOfPackedToken));
          };/* if */



      };/* else */
   } 
   else/* opCode == 6 or 8 */
    /************************************************
     *                                              *
     *         collect info for the unpacking       *
     *                                              *
     ***********************************************/
     {
       if ( packHeader != NULL )
         {
           theHeader     = (integer_header *)packedHeader;
           tokenSize     = theHeader->numOfBitsPerToken; 
           ShiftIntended = theHeader->SHIFT;
           intCount      = theHeader->numOfPackedToken;
           minSigned     = theHeader->minValue;
           minUnsigned   = theHeader->minValue;

         }
       else
         {
           tokenSize     = bitSizeOfPackedToken;
           ShiftIntended = 0;
           intCount      = intCount;
           minSigned     = minSignedInteger;
           minUnsigned   = minUnsignedInteger;
         };
     };
 
   

   /**********************************************
    *                                            *
    *  compute signed int adjustment,            *
    * since 1,2,3,4 without header all needs it  *
    *                                            *
    *********************************************/
   positiveMask = ( opCode < 7 ) ? 0 : ( 1 << ( bitSizeOfPackedToken - 1 ));


   /***********************************************
    *                                             *
    *   pack                                      *
    *                                             *
    **********************************************/
   if ( opCode == 5 )
     {
       Pack(arrayOfUnsignedShort, 0); 
     }
   /***********************************************
    *                                             *
    *   unpack                                    *
    *                                             *
    **********************************************/
   else if ( opCode == 6 )
     {
       Unpack(arrayOfUnsignedShort, ShiftIntended, tokenSize, 
              0, intCount);
     }
   else
     {
       printf("\n opCode:%d is not defined \n", opCode);
       return 0;
     };/* if */
 
  return intCount;  /* unused, function must return something */
 

} /* end compact_short */

/***********************************************************************************************
 *                                                                                             *
 * Author   : M. Lepine, sept 2005                                                             *
 *                                                                                             *
 * Objective: transform array of short integer between pack and unpacked format                *
 *                                                                                             *
 * Arguments:                                                                                  *
 *    IN/OUT  unpackedArrayOfBytes   unpacked array of bytes                                   *
 *    IN/OUT  packedArrayOfInt       packed integer array                                      *
 *    IN      elementCount           total count of element in unpacked integer array          *
 *    IN      bitSizeOfPackedToken   packed integer size in bit                                *
 *                                   when equal "-1", autodetection of bit length is on        *
 *    IN      offset                 in packing   : the last bit of integer packed inside array*
 *                                   in unpacking : the first bit of integer packed inside array*
 *    IN      stride                 unpacked integer spacing indicator                        *
 *    IN      opCode                 9: unsigned char pack                                     *
 *                                  10: unsigned char unpack                                   *
 **********************************************************************************************/
int  compact_char( void *unpackedArrayOfBytes, void *packedHeader, void *packedArrayOfInt, 
                       int elementCount, int bitSizeOfPackedToken, int off_set, 
                       int stride, int opCode)
{

    typedef struct
    {
#if defined(Little_Endian)
      word numOfBitsPerToken : 6, SHIFT : 6, unused : 12, ID : 8; 
#else
      word ID : 8, unused : 12, SHIFT : 6, numOfBitsPerToken : 6;
#endif
      word numOfPackedToken : 32;
      word minValue         : 32;
      word maxValue         : 32;
    }integer_header;

    int wordSize;
    word cleanupMask;
    int minSignedInteger=0, maxSignedInteger=0;
    word minUnsignedInteger=0, maxUnsignedInteger=0;
    word maxRange;
    word maxSpan;
   

    int positiveMask;
    unsigned char *arrayOfUnsignedChar;
    word *arrayOfPacked;
    int i, k;
    int intCount;
    int bitRequiredForRange, shiftRequired = 0;

    

    /****************************************
     *                                      *
     *     variables used by the packer     *
     *                                      *
     ***************************************/
    int lastPackBit, spaceInLastWord, lastSlot;
    word lastWordShifted, tempInt;
    word *packHeader;
    
    /***************************************
     *                                     *
     *    variables used by the unpacker   *
     *                                     *
     **************************************/
    int firstPackBit, bitPackInFirstWord, currentSlot;
    word currentWord, packInt;
    int significantBit, inSignificantBit;
    word *arrayPtr;
    int tokenSize, ShiftIntended, elementCountFromHeader;
    int  minSigned;
    word minUnsigned;
    integer_header *theHeader;


/*  printf("minSignedInteger=%d minUnsignedInteger=%d \n",minSignedInteger,minUnsignedInteger); */

    /********************************
     *                              *
     *   handle abnormal condition  *
     *                              *
     ********************************/
    /* token size is 0 */
    if ( bitSizeOfPackedToken == 0 )
      { 
        return 0;
      };  





   
    /********************************************************
     *                                                      *
     *    determine wordsize and others                     * 
     *                                                      *
     ********************************************************/
    wordSize                 = 8 * sizeof(word);
    arrayOfUnsignedChar      = (unsigned char *)unpackedArrayOfBytes;
    theHeader                = (integer_header *)packedHeader;
    packHeader               = (word *)packedHeader;
    arrayOfPacked            = (word  *)packedArrayOfInt;
    intCount                 = elementCount;
    cleanupMask              = ((word)(~0)>>(wordSize-bitSizeOfPackedToken));
      


   if (opCode==9)
    /************************************************
     *                                              *
     *         collect info for the packing         *
     *                                              *
     ***********************************************/   
   {
     if ( packedHeader != NULL )
       /*******************************************************************
        *                                                                  *
        *  pack header is required, (X - Xmin) is used as packInt          *
        *                                                                  *
        *******************************************************************/
       {
             constructHeader(arrayOfUnsignedChar, minUnsignedInteger, maxUnsignedInteger);
       }
    else
      /*************************************************************
       *                                                           *
       *   pack header not required, X itself is used as packInt,  *
       *   determines bitSizeOfPackedToken, if not available       *
       *                                                           *
       ************************************************************/
      {
        if ( bitSizeOfPackedToken == -1 )
          {
            /**********************************************************
             *                                                        *
             *   obtain minimum, maximun, span                        *
             *                                                        *
             *********************************************************/
            if ( opCode == 9 )
              {
                /* unsigned integer number */
                maxSpan = arrayOfUnsignedChar[0];

                for(i=stride; i < intCount*stride ; i+=stride)
                  {
                    maxSpan |= arrayOfUnsignedChar[i];
                  }; 
              };
          
            /************************************************************
             *                                                          *
             *           derive bitSizeOfPackedToken                    *
             *                                                          *
             ***********************************************************/
            bitSizeOfPackedToken = 0;
            while ( maxSpan != 0 )
              {
                maxSpan = maxSpan >> 1;
                bitSizeOfPackedToken++;
              };
            if ( opCode == 3 )
              {/* accomodate the signed bit */
                bitSizeOfPackedToken++;
              };
            cleanupMask = ((word)(~0)>>(wordSize-bitSizeOfPackedToken));
          };/* if */



      };/* else */
   } 
   else/* opCode == 10 */
    /************************************************
     *                                              *
     *         collect info for the unpacking       *
     *                                              *
     ***********************************************/
     {
       if ( packHeader != NULL )
         {
           theHeader     = (integer_header *)packedHeader;
           tokenSize     = theHeader->numOfBitsPerToken; 
           ShiftIntended = theHeader->SHIFT;
           intCount      = theHeader->numOfPackedToken;
           minSigned     = theHeader->minValue;
           minUnsigned   = theHeader->minValue;

         }
       else
         {
           tokenSize     = bitSizeOfPackedToken;
           ShiftIntended = 0;
           intCount      = intCount;
           minSigned     = minSignedInteger;
           minUnsigned   = minUnsignedInteger;
         };
     };
 
   

   /**********************************************
    *                                            *
    *  compute signed int adjustment,            *
    * since 1,2,3,4 without header all needs it  *
    *                                            *
    *********************************************/
   positiveMask = ( opCode < 11 ) ? 0 : ( 1 << ( bitSizeOfPackedToken - 1 ));


   /***********************************************
    *                                             *
    *   pack                                      *
    *                                             *
    **********************************************/
   if ( opCode == 9 )
     {
       Pack(arrayOfUnsignedChar, 0); 
     }
   /***********************************************
    *                                             *
    *   unpack                                    *
    *                                             *
    **********************************************/
   else if ( opCode == 10 )
     {
       Unpack(arrayOfUnsignedChar, ShiftIntended, tokenSize, 
              0, intCount);
     }
   else
     {
       printf("\n opCode:%d is not defined \n", opCode);
       return 0;
     };/* if */
 
  return intCount;  /* unused, function must return something */
 

} /* end compact_char */

