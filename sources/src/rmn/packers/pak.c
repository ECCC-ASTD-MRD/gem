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

#include<stdio.h>
#include<stdlib.h>
#include "rmnlib.h"
static PackFunctionPointer pfp;






/******************************************************************************************
 *                                                                                        *
 * Author    : Jianhui He, 1997                                                           *
 *                                                                                        *
 * Objective : FORTRAN interface to integer packer/unpacker written in C                  *
 *                                                                                        *
 * Arguments :                                                                            *
 *  IN/OUT     xunpacked     unpacked integer array                                       *
 *  IN/OUT     xpacked       packed integer array                                         *
 *  IN         ni            total count of element in unpacked integer array             *
 *  IN         nj            unpacked integer spacing indicator                           *
 *  IN         nBits         packed integer size in bit                                   *
 *  IN         NB            in packing   : the last bit of integer packed inside array   *
 *                           in unpacking : the first bit of integer packed inside array  *
 *  IN         OP            1   : unsigned pack                                          *
 *                           2   : unsigned unpack                                        *
 *                           3   : signed pack                                            *
 *                           4   : signed unpack                                          *
 *                           +10 : use offset                                             *
 *                           +100: NINJ = NJ and stride = NJ                              *
 *                           +1000: pack header is used                                   *
 *                                                                                        *
 *****************************************************************************************/      
void f77name (iipak) (void *xunpacked, void *xpacked, ftnword *ni, ftnword *nj, ftnword *nBits, 
                 ftnword *NB, ftnword *OP)
{

  /***************************************
   *                                     *
   *   local variables                   *
   *                                     *
   **************************************/ 
  int i;
  int stride;
  int offset;
  int ninj;
  int lnBits;
  int OPER = *OP;
  int lengthOfPackedArray;
  word *tempPackedArray,  *dummyArray, *tempUnpackedArray;
  word *ptrXpacked, *ptrXUnpacked, *transitPtr;
  ftnword *ftnPtrXUnpacked;
  int returnBitSizeOfToken;
  int packHeaderUsed;


  /***************************************
   *                                     *
   * fix lnBits                          *
   *                                     *
   **************************************/
  if(*nBits > 1)
     {
      if(1 > (32 / *nBits))
	 lnBits = 1;
      else
	 lnBits = (32 / *nBits);
      }
  else if (*nBits < 0) 
     lnBits = *nBits * (-1);
  else 
     lnBits = 32;

  /**************************************
   *                                    *
   *   check if pack header             *
   *        is used                     *
   *                                    *
   *************************************/
  if ( (OPER % 10000) >=1000 )
    {
      packHeaderUsed = 1;
    }
  else
    {
      packHeaderUsed = 0;
    };
  OPER = OPER % 1000;

  /***************************************
   *                                     *
   *  handle the stride and offset       *
   *  specification from FORTRAN program *
   *                                     *
   **************************************/
  if ((OPER % 1000) >= 100)
     {
     ninj = *ni;
     stride = *nj;
     OPER = OPER % 100;
     }
  else
    {
      ninj = (*ni) * (*nj);
      stride = 1;
    }

  if (OPER >= 10)
    {
      offset = *NB;
      
      OPER = (OPER % 10);
    }
  else
    offset = 0;
   




#if defined(NEC64)
  /*************************************************************************
   *                                                                        *
   *  64 bit mode on NEC                                                    *
   *                                                                        *
   *************************************************************************/ 

      /**********************************************************
      *                                                        *
      *  determine                                             *
      *     stride                                             *
      *     ptrXpacked, ptrXUnpacked                           *
      *     ftnPtrXUnpacked                                    *
      *     tempPackedArray                                    *
      *                                                        *
      *********************************************************/
      stride =stride*2;

      ptrXpacked   = NULL;
      ptrXpacked   = (word *)xpacked;
      ptrXUnpacked = NULL;
      ptrXUnpacked = (word *)xunpacked; 

      transitPtr = NULL;
      transitPtr = &ptrXUnpacked[1];
      ptrXUnpacked = NULL;
      ptrXUnpacked = transitPtr;
      transitPtr   = NULL;

      ftnPtrXUnpacked = NULL;
      ftnPtrXUnpacked = (ftnword *)xunpacked;

      /*******************************************************
       *                                                     *
       *   handle situation when offset exceeds 32 bit       *
       *                                                     *
       ******************************************************/
      if ( offset >= 32 )
        {
          ptrXpacked = &ptrXpacked[offset/32];
          offset = offset % 32;
        }

      lengthOfPackedArray = ((ninj*lnBits)+31+offset) / 32;
      if (packHeaderUsed) lengthOfPackedArray += 4; 
 
      tempPackedArray = NULL;    
      tempPackedArray = ( word *) malloc(lengthOfPackedArray*sizeof(word));


      
      
      if ( (OPER == 1) || (OPER == 3) )
        /*******************************************************
         *                                                     *
         *   pack                                              *
         *                                                     *
         ******************************************************/
        { 
          
          if ( packHeaderUsed )
            {
              returnBitSizeOfToken = compact_integer(ptrXUnpacked, tempPackedArray, 
                                                     &tempPackedArray[4], 
                                                     ninj, lnBits, 0, stride, OPER);
            }
          else
            {
              returnBitSizeOfToken = compact_integer(ptrXUnpacked, NULL, tempPackedArray, 
                                                     ninj, lnBits, offset, stride, OPER);
            };     

          /******************************************
           *                                        *
           *   avoid overriding the info            *
           *   in the leading cell of packed array  *
           *                                        *
           *****************************************/
          ptrXpacked[0] = ( ptrXpacked[0] & ( (-1) << (32-offset) ) ) | 
                                   tempPackedArray[0];

          /******************************************
           *                                        *
           *  shuffle back the packed integer       *
           *                                        *
           *****************************************/
          for ( i = 1; i < lengthOfPackedArray; i++)
            { 
              ptrXpacked[i] =  tempPackedArray[i];
            };
        }
      else if ( (OPER == 2) || (OPER == 4) )
        /*******************************************************
         *                                                     *
         *   unpack                                            *
         *                                                     *
         ******************************************************/
        {
          
          for ( i = 0; i < lengthOfPackedArray; i++)
            { 
              tempPackedArray[i] = ptrXpacked[i];
            };
          
          if ( packHeaderUsed )
            {        
              returnBitSizeOfToken = compact_integer(ptrXUnpacked, tempPackedArray, 
                                                     &tempPackedArray[4], 
                                                     ninj, lnBits, 0, stride, OPER);
            }
          else
            {
              returnBitSizeOfToken = compact_integer(ptrXUnpacked, NULL, tempPackedArray, 
                                                     ninj, lnBits, offset, stride, OPER);
            };
         
          /******************************************
           *                                        *
           *  shuffle back the unpacked integer     *
           *                                        *
           *****************************************/
          if ( lnBits == 32 )
            /**********************
             *                    *
             *   32 bit           *
             *                    *
             *********************/
            {
              for ( i = 0; i < ninj; i++)
                {
                  
                  ftnPtrXUnpacked[i] =  ptrXUnpacked[i*stride];
                };
            }
          else
            /**********************
             *                    *
             *  less than 32 bit  *
             *                    *
             *********************/
            {
              for ( i = 0; i < ninj; i++)
                {
              
                  if ( (ptrXUnpacked[i*stride] & 0x80000000 ) == 0 )
                    {
                  
                      ftnPtrXUnpacked[i] =  ptrXUnpacked[i*stride];
                    }
                  else
                    {
                      ftnPtrXUnpacked[i] = (-1) * abs(ptrXUnpacked[i*stride]);
                    };
                };

            };
            
        };
       
      free(tempPackedArray);  

#else

  /*******************************************************************
  *                                                                  *
  *   HP or 32 bit on NEC                                            *
  *                                                                  *
  *******************************************************************/

  ptrXpacked = NULL;
  ptrXpacked = xpacked;


  if ( xunpacked == xpacked )
    /*******************************************************
     *                                                     *
     *   in place                                          *
     *                                                     *
     ******************************************************/
    {  
     
      lengthOfPackedArray = ((ninj*lnBits)+31) / 32;
      if (packHeaderUsed) lengthOfPackedArray += 4;                  
      tempPackedArray = NULL;    
      
      tempPackedArray = ( word *) malloc(lengthOfPackedArray*sizeof(word));  
      

      if ( (OPER == 1) || (OPER == 3) )
        /****************************************
         *                                      *
         *   pack                               *
         *                                      *
         ***************************************/
        { 

          if ( packHeaderUsed )
            {
             
              returnBitSizeOfToken = compact_integer(xunpacked, tempPackedArray, tempPackedArray, 
                                                     ninj, lnBits, 128, stride, OPER);
            }
          else
            {
              returnBitSizeOfToken = compact_integer(xunpacked, NULL, tempPackedArray, 
                                                     ninj, lnBits, offset, stride, OPER);
            };
                                        
          for ( i = 0; i < lengthOfPackedArray; i++)
            { 
              ptrXpacked[i] = tempPackedArray[i];  
              
            };
            

        }
      else if ( (OPER == 2) || (OPER == 4) )
        /****************************************
         *                                      *
         *   unpack                             *
         *                                      *
         ***************************************/
        {
          for ( i = 0; i < lengthOfPackedArray; i++)
            { 
              tempPackedArray[i] = ptrXpacked[i] ;  
            };

          if ( packHeaderUsed )
            {            
              returnBitSizeOfToken = compact_integer(xunpacked, tempPackedArray, tempPackedArray, 
                                                     ninj, lnBits, 128, stride, OPER);
            }
          else
            {
              returnBitSizeOfToken = compact_integer(xunpacked, NULL, tempPackedArray, 
                                                     ninj, lnBits, offset, stride, OPER);
            };
          
          
        };
      
      free(tempPackedArray);  
      tempPackedArray = NULL;
      
    }
  else
    /*******************************************************
     *                                                     *
     *   non in place                                      *
     *                                                     *
     ******************************************************/
    {      
      if ( packHeaderUsed )
        {
          returnBitSizeOfToken = compact_integer(xunpacked, xpacked, xpacked, 
                                                 ninj, lnBits, 128, stride, OPER);
        }
      else
        {
          returnBitSizeOfToken = compact_integer(xunpacked, NULL, xpacked, 
                                                 ninj, lnBits, offset, stride, OPER);
        }; 
    };

   ptrXpacked = NULL;
   tempUnpackedArray = NULL;
#endif   
   
}

   

















/******************************************************************************************
 *                                                                                        *
 * Author    : Jianhui He, 1997                                                           *
 *                                                                                        *
 * Objective : FORTRAN interface to floating point packer/unpacker written in C           *
 *                                                                                        *
 * Arguments :                                                                            *
 *  IN/OUT     xunpacked     unpacked floating point array                                *
 *  IN/OUT     xpacked       packed integer array                                         *
 *  IN         ni            total count of element in unpacked array                     *
 *  IN         nj            unpacked array element spacing indicator                     *
 *  IN         nBits         packed integer size in bit                                   *
 *  IN         NB            in packing   : the last bit of integer packed inside array   *
 *                           in unpacking : the first bit of integer packed inside array  *
 *  IN         OP            1   : pack float                                             *
 *                           2   : unpack float                                           *
 *                           5   : pack double                                            *
 *                           6   : unpack double                                          *
 *                                                                                        *
 *****************************************************************************************/
void f77name (xxpak) (void *xunpacked, void *xpacked,ftnword *ni, ftnword *nj, ftnword *nBits,
                      ftnword *NB, ftnword *OP)
{
  /**************************************************************************
   *                                                                        *
   *   variable declaration                                                 *
   *                                                                        *
   *************************************************************************/
  int i;
  int stride = 1;
  int offset = 24;
  int ninj = ((*ni) * (*nj));
  int lnBits = (*nBits);
  int OPER = (*OP);
  word *tempIntArray;
  int lengthOfIntArray;
  word *tempFloatArray;
  double tempFloat=9999.0000;
  word *ptrXpacked;
 


  /***************************************************************************
   *                                                                         *
   *    determine lnBits                                                     *
   *                                                                         *
   **************************************************************************/
  if(*nBits > 1)
     {
      if(1 > (32 / *nBits))
	 lnBits = 1;
      else
	 lnBits = (32 / (int)*nBits);
      }
  else if (*nBits < 0) 
     lnBits = (int)*nBits * (-1);
  else 
     lnBits = 32;



  /***************************************************************************
   *                                                                         *
   *   allows 1 to represent pack doubles under NEC64 mode                   *
   *          2 to represent unpack doubles under NEC64 mode                 *
   *                                                                         *
   **************************************************************************/
#if defined (NEC64)
  if ( OPER == 1 )
    {
      OPER = 5;
    }
  else if (OPER == 2)
    {
      OPER = 6;
    };

#endif


  /*************************************************************************
   *                                                                       *
   *  determine function pointer and make OPER uniform for pack and unpack *
   *                                                                       *
   ************************************************************************/
  pfp = OPER > 3 ? &compact_double : &compact_float;
  OPER = OPER > 3 ? OPER - 4 : OPER;
  ptrXpacked = NULL;
  ptrXpacked = xpacked;


  /************************************************************************
   *                                                                      *
   *  call the appropriate pack/unpack routine written in C               *
   *                                                                      *
   ***********************************************************************/
  if ( xunpacked == xpacked )
    /*****************************************************
     *                                                   *
     *   inplace stuff                                   *
     ****************************************************/
    {  
    
      lengthOfIntArray = (ninj*lnBits) / 32 + 6;
      tempIntArray = NULL;    
      tempIntArray = ( word *) malloc(lengthOfIntArray*sizeof(word));  
      

      if (OPER == 1)
        /**************
         *  pack      *
         *************/
        {
          tempFloatArray = (*pfp)(xunpacked, &tempIntArray[0], &tempIntArray[3], ninj, lnBits, 
                                         offset, stride, OPER, 0, &tempFloat);                               
          for ( i = 0; i < lengthOfIntArray; i++)
            { 
              ptrXpacked[i] = tempIntArray[i];  
            };
        }
      else if ( OPER == 2 )
        /**************
         *  unpack    *
         *************/
        {
          for ( i = 0; i < lengthOfIntArray; i++)
            { 
              tempIntArray[i] = ptrXpacked[i] ;  
            };
          tempFloatArray = (*pfp)(xunpacked, &tempIntArray[0], &tempIntArray[3], ninj, lnBits, 
                                         offset, stride, OPER, 0, &tempFloat);
        };

      free(tempIntArray);  
      tempIntArray = NULL;     
    }
  else
    /***************************************************
     *                                                 *
     *  regular stuff                                  *
     *                                                 *
     **************************************************/
    {
      tempFloatArray = (*pfp)(xunpacked, &ptrXpacked[0], &ptrXpacked[3], ninj, lnBits, 
                                              offset, stride, OPER, 0, &tempFloat);
    };
    
   ptrXpacked = NULL;

}








