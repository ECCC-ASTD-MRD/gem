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

#include <rmnlib.h>



#define notSpecified -1





/***********************************************************************************************
 *                                                                                             *
 * Author   : Jianhui He, 1997                                                                 *
 *                                                                                             *
 * Objective: unpack routine wrapper                                                           *
 *                                                                                             *
 * Arguments:                                                                                  *
 *    OUT     unpackedArray          array of unpacked  numbers                                *
 *    IN      packedHeader           format information of packed integer numbers              *
 *    IN      packedArrayOfInt       array Of packed Integers                                  *
 *    IN      stride                 spacing indicator                                         *
 *    IN      missingTag             missing value identifier                                  *
 **********************************************************************************************/

void unpackWrapper(void *unpackedArray, void *packedHeader, void *packedArray,
                   int stride, void *missingValueTag)
{

  int   elementCount, bitSizeOfPackedToken, bitSizeOfPackedExpo, off_set, opCode, hasMissing;
  int   min, max;
  unsigned INT_32  headerType;
  unsigned INT_32 *packHeader;
  


  elementCount         = notSpecified;
  bitSizeOfPackedToken = notSpecified;
/*  bitSizeOfPackedToken = 0; */
  off_set              = 128;
  min                  = notSpecified;
  max                  = notSpecified;
  packHeader           = (unsigned INT_32 *)packedHeader;
  headerType           = packHeader[0] >> 24;



  if (( headerType == 0x0000007f ) || ( headerType == 0x000000ff ))
    /***********************************************
     *                                             *
     *    floating point, without missing value    *
     *                                             *
     **********************************************/
    {
      opCode     = 2;
      hasMissing = 0;
      compact_float(unpackedArray, packedHeader, packedArray, 
                    elementCount, bitSizeOfPackedToken, off_set, stride, opCode, 
                    hasMissing, missingValueTag);
    }
  else if (( headerType == 0x0000007e ) || ( headerType == 0x000000fe ))
    /***********************************************
     *                                             *
     *    floating point, with missing value       *
     *                                             *
     **********************************************/
    {
      opCode     = 2;
      hasMissing = 1;
      compact_float(unpackedArray, packedHeader, packedArray, 
                    elementCount, bitSizeOfPackedToken, off_set, stride, opCode, 
                    hasMissing, missingValueTag);
    }
  else if ( headerType == 0x000000fb )
    /***********************************************
     *                                             *
     *    IEEE block                               *
     *                                             *
     **********************************************/
    {
      opCode              = 2;
      hasMissing          = 0;
      bitSizeOfPackedExpo = notSpecified;
      off_set             = 0;
      compact_IEEEblock_float(unpackedArray, packedHeader, packedArray, 
                              elementCount, bitSizeOfPackedToken, bitSizeOfPackedExpo,
                              off_set, stride, opCode, hasMissing, missingValueTag);
    }
  else if ( headerType == 0x000000fd )
    /***********************************************
     *                                             *
     *    integer                                  *
     *                                             *
     **********************************************/
    {
      opCode = 2;
      
      compact_integer(unpackedArray, packedHeader, packedArray,
                      elementCount, bitSizeOfPackedToken, off_set, stride, 
                      opCode); 
    }
  else if ( headerType == 0x000000f0 )
    /***********************************************
     *                                             *
     *    run length encoding                      *
     *                                             *
     **********************************************/
    { 
      opCode = 2;
      compact_rle(unpackedArray, packedHeader, packedArray,
                  max, min,
                  elementCount, bitSizeOfPackedToken, off_set, 
                  stride, opCode); 
    }
  else 
    {
      printf("\n %8.8x not a valid header \n", headerType);
    };


}

