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



#define powerSpan  33
static double powerOf2s[powerSpan];
static int powerOf2sInitialized = 0;





/***********************************************************************************************
 *                                                                                             *
 * Author   : Jianhui He, 1997                                                                 *
 *                                                                                             *
 * Objective: transform integer array between pack and unpacked format in RLE fashion          *
 *                                                                                             *
 * Arguments:                                                                                  *
 *    IN/OUT  unpackedArrayOfInt     unpacked integer array                                    *
 *    IN/OUT  packedHeader           information about the packed integer array                *
 *    IN/OUT  packedArrayOfInt       packed integer array                                      *
 *    IN      max, min               maximun integer and minimun integer, they are ignored     *
 *                                   if max is equal min                                       *
 *    IN      elementCount           total count of element in unpacked integer array          *
 *    IN      bitSizeOfPackedToken   packed integer size in bit                                *
 *                                   when equal "-1", autodetection of bit length is on        *
 *    IN      offset                 in packing   : the last bit of integer packed inside array*
 *                                   in unpacking : the first bit of integer packed inside array*
 *    IN      stride                 unpacked integer spacing indicator                        *
 *    IN      opCode                 1: pack                                                   *
 *                                   2: unpack                                                 *
 **********************************************************************************************/
int  compact_rle( void *unpackedArrayOfInt, void *packedHeader, void *packedArrayOfInt, 
                  int max, int min,
                  int elementCount, int bitSizeOfPackedToken, int off_set, 
                  int stride, int opCode)
{

    typedef struct
    {
#if defined(Little_Endian)
      word RESERVED :16,RLE_TYP :4,LNG :4,ID :8,numOfPackedToken :26,numOfBitsPerToken :6;
#else
      word ID :8,LNG :4,RLE_TYP :4,RESERVED :16,numOfBitsPerToken :6,numOfPackedToken :26;
#endif
      word minInteger :32,maxRange :32;
    }rle_header;


    int wordSize;
    word *arrayOfUnpacked;
    word *arrayOfPacked;
    int i, k;
    int intCount;

    /****************************************
     *                                      *
     *     variables used by the packer     *
     *                                      *
     ***************************************/
    int lastPackBit, spaceInLastWord, lastSlot;
    word lastWordShifted, tempInt;

    word minInteger, maxInteger;
    word maxRange;
    word currentToken, previousToken;
    
    int previousTokenCount, currentTokenCount;
   
    int currentPosition, totalPackedElementCount, maxRepeatPerToken;
    word *packHeader;
    int lastTokenPosition;


    word *arrayPosition;

    
    /***************************************
     *                                     *
     *    variables used by the unpacker   *
     *                                     *
     **************************************/
    int firstPackBit, bitPackInFirstWord, currentSlot;
    word currentWord, packInt;
    int significantBit, inSignificantBit;
    word *tempPackedArray;
    rle_header *theHeader;
    int repeatCount;

   
    /********************************
     *                              *
     *   handle abnormal condition  *
     *                              *
     ********************************/
    /* token size is 0 */
    if ( bitSizeOfPackedToken < 3 )
      { 
        return 0;
      };  
   

    /*************************************************
     *                                               *
     *      initialize the power ladder of base 2    *
     *                                               *
     ************************************************/
    if ( ! powerOf2sInitialized )
      {
        powerOf2s[0] = 1.0;
        for ( i = 1; i < powerSpan; i++)
          {
            powerOf2s[i] = 2.0 *powerOf2s[i-1];      
          };
        powerOf2sInitialized = 1;
      };

      




    /********************************************************
     *                                                      *
     *    determine wordsize and others                     * 
     *                                                      *
     ********************************************************/
    wordSize        = 8 * sizeof(word);
    arrayOfUnpacked = (word *)unpackedArrayOfInt;
    packHeader      = (word *) packedHeader;
    arrayOfPacked   = (word *)packedArrayOfInt;
    intCount        = elementCount;


    



    if ( opCode == 1 )
    {


    /**********************************************************
     *                                                        *
     *   obtain minimum, maximun, range to derive             *
     *      maxRepeatPerToken                                 *
     *                                                        *
     *********************************************************/
    if ( max == min )
      {
        maxInteger = arrayOfUnpacked[0];
        minInteger = arrayOfUnpacked[0];
        tempInt = intCount*stride;
        i = stride;
        while ( i < tempInt )
          {
            
            if ( arrayOfUnpacked[i] < minInteger )
              {      
                minInteger = arrayOfUnpacked[i];
              }
            else if ( arrayOfUnpacked[i] > maxInteger )
              {
                maxInteger = arrayOfUnpacked[i];
              }; 
                          
            i+=stride;
            
          };
       
      }
    else
      {
        maxInteger = max;
        minInteger = min;
      };
    maxRange       = maxInteger - minInteger;
    maxRepeatPerToken = powerOf2s[bitSizeOfPackedToken] - 1 - maxRange;


     /************************************************
      *                                              *
      *          initialize local variables          *
      *                                              *
      ************************************************/
      lastPackBit = off_set;
      spaceInLastWord =  wordSize - ( lastPackBit % wordSize );
      lastSlot = ( lastPackBit / wordSize );
      if ( spaceInLastWord == wordSize )
        { 
          lastWordShifted = 0;
        }
      else
        { 
          lastWordShifted = arrayOfPacked[lastSlot] >> spaceInLastWord ;
        };

      /*************************************************
       *                                               *
       *       packing                                 *
       *                                               *
       ************************************************/
      totalPackedElementCount = 0;      
      previousToken           = -1;
      previousTokenCount      = 1;
      lastTokenPosition       = (intCount-1)*stride; 
      currentPosition         = 0;
      
        
      arrayPosition = &arrayOfPacked[lastSlot];
      for ( i = 0; i < lastTokenPosition; i += stride )
        {
        
          
          currentToken = arrayOfUnpacked[i] - minInteger;
          


          /*******************************************************
           *                                                     *
           *   construct packing info                            *
           *                                                     *
           ******************************************************/      
          if ( currentToken != previousToken ) 
            /** new token encountered **/
            {
              if (previousTokenCount == 1)
                {
                  stuff(currentToken, arrayPosition, wordSize, bitSizeOfPackedToken, 
                        lastWordShifted, spaceInLastWord);
                  totalPackedElementCount++;
                }
              else if ( previousTokenCount == 2 )
                {
                  stuff(previousToken, arrayPosition, wordSize, bitSizeOfPackedToken, 
                        lastWordShifted, spaceInLastWord);
                  stuff(currentToken, arrayPosition, wordSize, bitSizeOfPackedToken, 
                        lastWordShifted, spaceInLastWord);
                  totalPackedElementCount+=2;
                }
              else
                {
                  tempInt   = maxRange+previousTokenCount;
                  stuff(tempInt, arrayPosition, wordSize, bitSizeOfPackedToken, 
                        lastWordShifted, spaceInLastWord);
                  stuff(currentToken, arrayPosition, wordSize, bitSizeOfPackedToken, 
                        lastWordShifted, spaceInLastWord);
                  totalPackedElementCount+=2;
                }
              previousToken      = currentToken;
              previousTokenCount = 1;
            }
          else
            /** a repeated token **/
            {
              previousTokenCount++;
              if ( previousTokenCount == maxRepeatPerToken )
                {
                  tempInt = maxRange + maxRepeatPerToken;
                  stuff(tempInt, arrayPosition, wordSize, bitSizeOfPackedToken, 
                        lastWordShifted, spaceInLastWord);
                  totalPackedElementCount++;
                  /* treat it same as starting condition */
                  previousToken       = -1; 
                  previousTokenCount = 1;
                }
            }
         
        };/* for */
     

      /**************************************************
       *                                                *
       * handle the last token                          *
       *                                                *
       *************************************************/
      currentToken = arrayOfUnpacked[lastTokenPosition] - minInteger;
      if ( currentToken != previousToken )
        /* new token */
        {
          stuff(currentToken, arrayPosition, wordSize, bitSizeOfPackedToken, 
                lastWordShifted, spaceInLastWord);
        }
      else
        /* repeated token */
        {
          previousTokenCount++;
          if ( previousTokenCount == 2 )
            {
              stuff(currentToken, arrayPosition, wordSize, bitSizeOfPackedToken, 
                    lastWordShifted, spaceInLastWord);  
            }
          else
            {
                      
              tempInt = maxRange + previousTokenCount;
              stuff(tempInt, arrayPosition, wordSize, bitSizeOfPackedToken, 
                    lastWordShifted, spaceInLastWord);
            };    
        };
      totalPackedElementCount++;

      /**************************************************
       *                                                *
       *     squeezes hole left in the integer array    *
       *                                                *
       **************************************************/
      if ( spaceInLastWord < wordSize )
        {
          *arrayPosition = ( lastWordShifted << spaceInLastWord) |
            ( *arrayPosition & ~(-1 << spaceInLastWord)); 
        };

      /*****************************
       *                           *
       * construct pack header     *
       *                           *
       ****************************/
      packHeader[0] = 0xF0000000 | (1 << (32-8-4));
      packHeader[1] = ( bitSizeOfPackedToken << 26 ) | totalPackedElementCount;
      packHeader[2] = minInteger;
      packHeader[3] = maxRange;


      return totalPackedElementCount;
    }


  else if ( opCode == 2 ) 
    /***********************************
     *       u n p a c k e d           *
     *                                 *
     ***********************************/
  {

   
   
    /***********************************
     *                                 *
     *    initialize variables         *
     *                                 *
     **********************************/
    firstPackBit       = off_set;
    bitPackInFirstWord =  wordSize - ( firstPackBit % wordSize );
    currentSlot        = ( firstPackBit / wordSize );
    currentWord        = arrayOfPacked[currentSlot] << ( wordSize - bitPackInFirstWord );
    intCount           = elementCount;


    /*************************************
     *                                   *
     *  recover info from the header     *
     *                                   *
     ************************************/
    theHeader = (rle_header *)packedHeader;
    minInteger = theHeader->minInteger;
    maxInteger = minInteger + theHeader->maxRange;
    totalPackedElementCount = theHeader->numOfPackedToken;
    maxRange = maxInteger - minInteger;

    /*
    significantBit = bitSizeOfPackedToken;
    */
    significantBit = theHeader->numOfBitsPerToken;


    
    /*****************************************
     *                                       *
     *    unpack integer numbers             *
     *                                       *
     *****************************************/
    currentPosition = 0;
    arrayPosition   = &arrayOfPacked[currentSlot];
    for ( i = 0; i < totalPackedElementCount; i++)
      {     
        extract( tempInt, arrayPosition, wordSize, significantBit,
                 currentWord, bitPackInFirstWord);






        currentToken = tempInt + minInteger;
        if (tempInt <= maxRange )
          {
            arrayOfUnpacked[currentPosition] = currentToken;
            currentPosition += stride;

          }
        else
          {
            
            repeatCount = tempInt-maxRange-1;
            
            /******************************************
             *                                        *
             * while loop is more efficient           *
             * when used inside a for loop            *
             *                                        *
             *****************************************/ 
            
            while ( repeatCount--)
              {
                arrayOfUnpacked[currentPosition] = previousToken;
                currentPosition += stride;
              };
                      
          };
        previousToken = currentToken;

      };/* for */

    return totalPackedElementCount;
  }
  else
  {
    printf("\n opCode:%d is not defined \n", opCode);
    return NULL;
  };/* if */
 
  

} /* end compact_rle */









