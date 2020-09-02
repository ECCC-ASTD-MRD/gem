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
#include "rmnlib.h"






/***********************************************************************************************
 *                                                                                             *
 * Author   : Jianhui He, 1997                                                                 *
 *                                                                                             *
 * Objective:                                                                                  *
 *            pack IEEE floating point numbers( float and double ),                            *
 *            into a reduced bit size "IEEE" floating point numbers                            *
 *            ( a routine to be included by compact_IEEEblock.c )                              *
 *                                                                                             *
 * Arguments:                                                                                  *
 *    IN/OUT  unpackedArrayOfFloat   array of floating point numbers                           *
 *    IN/OUT  packedHeader           format information of packed integer numbers              *
 *    IN/OUT  packedArrayOfInt       array Of Integers                                         *
 *    IN      elementCount           total count of element in floating point array            *
 *    IN      bitSizeOfPackedToken   token size in bit                                         *
 *    IN      bitSizeOfPackedExpo    packed exponent size in bit                               *
 *    IN      offset                 in packing   : the last bit of integer packed inside array*
 *                                   in unpacking : the first bit of integer packed into array *
 *    IN      stride                 floating point number spacing indicator                   *
 *    IN      opCode                 operator(FLOAT_PACK, FLOAT_UNPACK)                        *
 *    IN      hasMissing             1: inicate existence of missing value in                  *
 *                                      floating point array                                   *
 *                                   0: otherwise                                              *
 *    IN      missingTag             missing value identifier                                  *
 **********************************************************************************************/
void *compact_IEEEblock_FLOAT_4_8(void *unpackedArrayOfFloat, void *packedHeader, 
                                  void *packedArrayOfInt, int elementCount, 
                                  int bitSizeOfPackedToken, int bitSizeOfPackedExpo, int off_set,
                                  int stride, int opCode, int hasMissing, void *missingTag )
{


  /*******************************************
   *                                         *
   *    declare header type                  *
   *                                         *
   ******************************************/
    typedef struct
    { 
#if defined(Little_Endian)
      word bitSizeOfExpo : 5, bitSizeOfToken : 7, maxExpo : 8, marker : 12, count : 32;
#else
      word marker : 12, maxExpo : 8, bitSizeOfToken : 7, bitSizeOfExpo : 5, count : 32;
#endif
    }IEEEblock_struct_data;


    /****************************************
     *                                      *
     *     variables used by the packer     *
     *                                      *
     ***************************************/
    int   wordSize;
    FLOAT_4_8 *arrayOfFloat;
    word *packHeader, *arrayOfInt;
    int   i, k;
    word  floatCount;
    float maxFloat;
    ALL_FLOAT maxFloatTemplate, floatTemplate;
    word  maxFloatExpo, maxPackedExpo, expoAlignmentFactor;
    word  packedSign, packedExpo, packedMantisa;
    int   expoDifference;
    int   lastPackBit, spaceInLastWord, lastSlot;
    word  lastWordShifted, tempInt;
    word *arrayPtr, *arrayOfUnpacked;
    FLOAT_4_8 *arrayOfTempFloat;
    float tempFloat;




    /*****************************************
     *                                       *
     *    variables used by the unpacker     *
     *                                       *
     ****************************************/
    IEEEblock_struct_data *theHeader;
    word currentWord;
    word intCount;
    int  firstPackBit;
    word bitPackInFirstWord;
    int  currentSlot;
    word packInt;
    int  tempSignedInt;
    int  significantBit, inSignificantBit;
    int  packedTokenSize, packedMantisaSize, packedExpoSize;

    
    /****************************************
     *                                      *
     *   obtain an array of power of 2      *    
     *                                      *
     ***************************************/
    if ( ! powerOf2sInitialized )
      {
        powerOf2s[0] = 1.0;
        for ( i = 1; i < powerSpan; i++)
          {
            powerOf2s[i] = 2.0 *powerOf2s[i-1];      
          };
        powerOf2sInitialized = 1;
      };


    /****************************************
     *                                      *
     *   handle abnormal condition          *
     *                                      *
     ***************************************/
    /* token size is 0 */
    if ( bitSizeOfPackedToken == 0 )
      { 
        return NULL;
      };  
    /* missing value handling routine no available yet */
    if ( hasMissing ) 
      {
        return NULL;
      };
    /* element count exceeds allowable limits in header field */
    if ( elementCount > powerOf2s[32] -1 )
      {
        return NULL;
      };

      
    /****************************************
     *                                      *
     *    determine wordsize                * 
     *                                      *
     ***************************************/
    wordSize         = 8 * sizeof(word);
    
 


    


























    if ( opCode == FLOAT_PACK )   
      /********************************************************************
       *                                                                  *
       *       compact a floating point array into an word array          *
       *                                                                  *
       ********************************************************************/
    {
      /*********************************************************
       *                                                       *
       *         initialize variables                          *
       *                                                       *
       *********************************************************/
      arrayOfFloat = (FLOAT_4_8 *)unpackedArrayOfFloat;
      packHeader   = (word  *)packedHeader;
      arrayOfInt   = (word  *)packedArrayOfInt;
      floatCount   = elementCount;


      /*********************************************************
       *                                                       *
       *  determine the bitSizeOfMantisa                       *
       *                                                       *
       ********************************************************/
      packedMantisaSize = bitSizeOfPackedToken - 1 - bitSizeOfPackedExpo; 
      maxPackedExpo     = powerOf2s[bitSizeOfPackedExpo] - 1;

     


      /*********************************************************
       *                                                       *
       *      obtain the maximun float                         *
       *                                                       *
       ********************************************************/
      maxFloat = arrayOfFloat[0];
        
      for(i=stride; i < floatCount*stride ; i+=stride)
        {
          if ( arrayOfFloat[i] > maxFloat )
            {
              maxFloat = arrayOfFloat[i];
            }; 
        }; 

      maxFloatTemplate.X = maxFloat;



  
    /*************************************************************
     *                                                           *
     *      obtain the alignmentFactor                           *
     *                                                           *
     *************************************************************/
      maxFloatExpo        = maxFloatTemplate.M.expo;
      expoAlignmentFactor = maxFloatExpo - maxPackedExpo;

    
    /*************************************************************************
     *                                                                       *
     *              initialize the header of the integer array               *
     *                                                                       *
     ************************************************************************/
    packHeader[0] = 0xfb0 << 20 | maxFloatExpo << 12 | bitSizeOfPackedToken << 5 
                   | bitSizeOfPackedExpo;

    packHeader[1] = floatCount;


    /***********************************************************************
     *                                                                     *
     *    transform the floating point into                                *
     *    the desired representation                                       *
     *                                                                     *
     **********************************************************************/ 

  
    /******************************
     *                            *
     * initialize other variables *
     *                            *
     *****************************/
    lastPackBit     = off_set;
    spaceInLastWord =  wordSize - ( lastPackBit % wordSize );
    lastSlot        = ( lastPackBit / wordSize );

    if ( spaceInLastWord == wordSize )
      { 
        lastWordShifted = 0;
      }
    else
      { 
        lastWordShifted = arrayOfInt[lastSlot] >> spaceInLastWord ;
      };

    arrayPtr = &arrayOfInt[lastSlot];
    if ( bitSizeOfPackedToken == 64 )
      /**************************
       *                        *
       *   64 bits direct copy  *
       *                        *
       *************************/
      {
        
/* CHC/NRC       arrayOfTempFloat = (double *)arrayOfInt; */
        arrayOfTempFloat = (FLOAT_4_8 *)arrayOfInt; 
        for ( i = 0; i < floatCount*stride; i+=stride)
         {
           arrayOfTempFloat[i] = arrayOfFloat[i];
         };
        return arrayOfTempFloat;

      }
    else if (( spaceInLastWord == wordSize ) && ( bitSizeOfPackedToken == wordSize ))
      /**************************
       *                        *
       *   32 bits direct copy  *
       *                        *
       *************************/
      {

        for ( i = 0; i < floatCount*stride; i+=stride)
         {
           floatTemplate.X = arrayOfFloat[i];
           arrayOfInt[i]   = floatTemplate.U;

         };
        return arrayOfTempFloat;

      }
    else
      /**************************
       *                        *
       *  bit by bit shuffle    *
       *                        *
       *************************/
      {
        for ( i = 0; i < floatCount*stride; i+=stride)
         {

            /*****************************
             *                           *
             *     determine tempInt     *
             *                           *
             ****************************/
            floatTemplate.X = arrayOfFloat[i];
            expoDifference = floatTemplate.M.expo - expoAlignmentFactor;
            if ( expoDifference >= 0 )
              {
                packedSign = floatTemplate.M.sign;
                packedExpo = expoDifference;

                
                if ( packedMantisaSize >=24 )
                  {
                    packedMantisa = floatTemplate.M.mantis;
                  }
                else
                  {
                    packedMantisa=floatTemplate.M.mantis>>(24-packedMantisaSize);
                  };
                  
                   
                tempInt = packedSign << ( bitSizeOfPackedToken - 1)  |
                          packedExpo << packedMantisaSize      |
                          packedMantisa; 

              }
            else
              {
                tempInt = 0;
              };

            /*****************************
             *                           *
             *  pack tempInt             *
             *                           *
             ****************************/
            stuff(tempInt, arrayPtr, wordSize, bitSizeOfPackedToken, lastWordShifted, 
                  spaceInLastWord); 
           
         }; /* for */
      };

        /********************************
         *                              *
         *  squeezes hole left in       *
         *   the integer array          *
         *                              *
         *******************************/
    if ( spaceInLastWord < wordSize )
      {
        *arrayPtr = ( lastWordShifted << spaceInLastWord) |
                    ( *arrayPtr & ~(-1 << spaceInLastWord));   
      };
    return (word *)arrayOfInt;
    }



    
 































  else if ( opCode == FLOAT_UNPACK ) 
    /**********************************************************************************
     *                                                                                *
     *       u n p a c k e d                                                          *
     *                                                                                *
     *********************************************************************************/
  {
   
    /********************************************************
     *                                                      *
     *   initialize variables                               *
     *   retrieve information from the header               *
     *                                                      *
     *******************************************************/

    arrayOfFloat        = (FLOAT_4_8 *)unpackedArrayOfFloat;
    theHeader           = ( IEEEblock_struct_data *) packedHeader;
    arrayOfInt          = (word  *)packedArrayOfInt;
    packedTokenSize     = theHeader->bitSizeOfToken;
    packedExpoSize      = theHeader->bitSizeOfExpo;
    maxPackedExpo       = powerOf2s[packedExpoSize] - 1;
    expoAlignmentFactor = theHeader->maxExpo  - maxPackedExpo;
    intCount            = theHeader->count;
    packedMantisaSize   = packedTokenSize - 1 - packedExpoSize;
    firstPackBit        = off_set;
    bitPackInFirstWord  =  wordSize - ( firstPackBit % wordSize );
    currentSlot         = ( firstPackBit / wordSize );
    currentWord         = arrayOfInt[currentSlot] << ( wordSize - bitPackInFirstWord );

    if ( packedTokenSize > wordSize )
      {
        significantBit   = wordSize;
        inSignificantBit = packedTokenSize - wordSize;
      }
    else
      {
        significantBit   = packedTokenSize;
        inSignificantBit = 0;
      };

        

    /******************************************************
     *                                                    *
     *    unpack floating point numbers                   *
     *                                                    *
     *****************************************************/
    arrayPtr = &arrayOfInt[currentSlot];

    if ( packedTokenSize == 64 )
      /******************************************
       *                                        *
       *   64 bit target IEEE size              *
       *                                        *
       *****************************************/
      {
           
/* CHC/NRC        arrayOfTempFloat = (double *)arrayOfInt; */
        arrayOfTempFloat = (FLOAT_4_8 *)arrayOfInt;
        for ( i = 0; i < intCount*stride; i+=stride)
          {
            arrayOfFloat[i] = arrayOfTempFloat[i];
          };
        return arrayOfFloat;
      }
    else if ((significantBit == wordSize) && ( bitPackInFirstWord == wordSize))
      /******************************************
       *                                        *
       *   32 bit target IEEE size              *
       *                                        *
       *****************************************/
      {
        for ( i = 0; i < intCount*stride; i+=stride)
          {
            floatTemplate.U = arrayOfInt[i];
            arrayOfFloat[i] = floatTemplate.X;

          };
        return arrayOfFloat;
      }
    else
      /******************************************
       *                                        *
       *   >32 bit target IEEE size             *
       *                                        *
       *****************************************/
      {
        for ( i = 0; i < intCount*stride; i+=stride)
          {     
            extract(packInt, arrayPtr, wordSize, significantBit,                          
                    currentWord, bitPackInFirstWord); 

            /******************************************
             *                                        *
             *   truncate extra bit if necessary      *
             *                                        *
             *****************************************/
            if ( inSignificantBit > 0 )
              {
                discard(arrayPtr, wordSize, inSignificantBit, currentWord, 
                        bitPackInFirstWord);     
              };   
            
            
            /******************************************
             *                                        *
             *   recover the IEEE numbers             *
             *                                        *
             *****************************************/
            if ( packInt == 0 )
              {
                arrayOfFloat[i] = 0;
              }
            else
              {
                floatTemplate.M.sign    = (packInt << (wordSize - packedTokenSize)) 
                                                   >> ( wordSize -1);
                tempSignedInt           = (packInt << (wordSize - packedTokenSize + 1)) 
                                                   >> (wordSize - packedExpoSize);
                floatTemplate.M.expo    = tempSignedInt + expoAlignmentFactor ;
                floatTemplate.M.mantis  = (packInt << 
                                           (wordSize - packedTokenSize + 1 + packedExpoSize)) 
                                                   >> (wordSize - packedMantisaSize );

                arrayOfFloat[i] = floatTemplate.X;
              };
      


          };/* for */
      };/* else */




    return ((word *)arrayOfFloat);

  }
  else
  {
    printf("\n opCode is not defined \n");
    return NULL;
  };/* if */
 
 

} /* end compact_float */
#undef compact_IEEEblock_FLOAT_4_8
#undef isDouble
#undef FLOAT_4_8
