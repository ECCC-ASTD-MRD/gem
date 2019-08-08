/* useful routines for C and FORTRAN programming
 * Copyright (C) 2013  centre ESCER, Universite du Quebec A Montreal
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
*/
/* 
 * This is a set of routines used to manage the "missing data" feature of the standard file package.
 * These routines can also be used by other packages, provided that the datatype indicator
 * uses values that are the same as the standard file package, namely:
 * 
 *   data type codes
 *   0: binary, transparent
 *   1: floating point
 *   2: unsigned integer
 *   3: character (R4A in an integer)
 *   4: signed integer
 *   5: IEEE floating point
 *   6: floating point 
 *   7: character string
 *   8: complex IEEE
 * +128 = with second stage packer
 * +64  = missing data flagging is used
 * 
 * to activate the package:
 * - set environment variable MISSING_VALUE_FLAGS
 *   (the order of the values is important)
 *   (real, integer, unsigned integer, 8 byte real, 2 byte integer, unsigned 2 byte integer)
 *   example:
 *   export MISSING_VALUE_FLAGS="1.234E37 123456 9999 4.567E+300 -5432 65432"
 *   NOTE: if less than 6 values are given, the missing ones are set to internal default values
 *            export MISSING_VALUE_FLAGS="   " 
 *         would activate the package with the internal default values
 * 
 * - call subroutine set_missing_value_flags 
 * 
 * There are C and FORTRAN callable routines/functions
 * 
 * =================================================================================
 *   missing_value_used, ForceMissingValueUsage
 * =================================================================================
 * 
 *   get the state of the package
 *   this integer function will return 0 if feature is not active, 1 if it is
 * 
 * C interface:
 *   int missing_value_used()
 * FORTRAN interface:
 *   integer, external :: missing_value_used
 * 
 *   forcibly set the state of the package
 * 
 * C interface:
 *   oldmode = ForceMissingValueUsage(flag)
 * 
 * FORTRAN interface:
 *   oldmode = Force_Missing_ValueUsage(flag)
 * =================================================================================
 *   get_missing_value_flags
 * =================================================================================
 * 
 *   this integer function is used to get the special values used to flag "missing" data for
 *   real, integer, unsigned integer, 8 byte real, 2 byte integer, unsigned 2 byte integer
 *   this integer function will return 0 if feature is not active, 1 if it is
 * 
 * C interface:
 *   int get_missing_value_flags(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
 *                               signed char *b, unsigned char *ub)
 * FORTRAN interface:
 *   integer function get_missing_value_flags(f,i,ui,d,s,us,b,ub)
 *   real, intent(OUT) :: f
 *   integer, intent(OUT) :: i, ui
 *   real *8, intent(OUT) :: d
 *   integer *2, intent(OUT) :: s, us
 *   integer *1, intent(OUT) :: b, ub
 *   NOTE: us, ui, us are unsigned numbers, which can be confusing in FORTRAN
 * 
 * =================================================================================
 *   set_missing_value_flags
 * =================================================================================
 * 
 *   this subroutine is used to set the special values used to flag "missing" data for
 *   real, integer, unsigned integer, 8 byte real, 2 byte integer, unsigned 2 byte integer
 *   calling this subroutine will also activate the package
 * 
 * C interface:
 *   void set_missing_value_flags(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
 *                               signed char *b, unsigned char *ub)
 * FORTRAN interface:
 *   subroutine set_missing_value_flags(f,i,ui,d,s,us,b,ub)
 *   real, intent(IN) :: f
 *   integer, intent(IN) :: i, ui
 *   real *8, intent(IN) :: d
 *   integer *2, intent(IN) :: s, us
 *   integer *1, intent(IN) :: b, ub
 *   NOTE: us, ui, us are unsigned numbers, which can be confusing in FORTRAN
 * 
 * =================================================================================
 * EncodeMissingValue / encode_missing_value
 * =================================================================================
 * 
 * this function copies data from field2 into field replacing values that have the "missing data" value 
 * with values that are non detrimental to further packing by the "standard file" packers. field 2 is left
 * as is, field will contain the output.
 * if proper encoding is not possible, the function returns zero immediately, leaving field untouched.
 * 
 * the return value of the function is the number of missing values found ( 0 if unsuccessful for any reason)
 * 
 * C interface:
 *   int EncodeMissingValue(void *field,void *field2,int nvalues,int datatype,int nbits,int is_byte,int is_short,int is_double);
 * FORTRAN interface:
 *   integer function encode_missing_value(field,field2,nvalues,datatype,nbits,is_byte,is_short,is_double)
 *   "any_type", intent(OUT) :: field
 *   "any_type", intent(IN) :: field2
 *   integer, intent(IN) :: nvalues,datatype,nbits,is_byte,is_short,is_double
 * 
 * field2   : input array (unspecified type)
 * field    : output array (assumed to have the same type as field2)
 * nvalues  : number of elements in arrays field and field2
 * datatype : what kind of data is found in array (see above)
 * nbits    : number of bits per element to be used for packing (used to compute appropriate flag values)
 * is_byte  : non zero if array element length is one byte
 * is_short : non zero if integer (signed or not) data element length is two bytes
 * is_double: non zero if real data element length is eight bytes (double / real *8)
 * 
 * the missing value will be replaced by a value greater than the largest value in field2 and 
 * large enough not to create an after packing collision with the value of the actual largest
 * useful value in field2 (nbits is needed for this calculation)
 * =================================================================================
 * SetMissingValueMapping (set_missing_value_mapping from Fortran)
 * set the functions used for encoding missing values and/or decoding encoded fields
 * 
 * C interface:
 *   SetMissingValueMapping(int what, int datatype, void *processor, int is_byte, int is_short, int is_double)
 * 
 * FORTRAN interface:
 *   call set_missing_value_mapping(what, datatype, processor, is_byte, is_short, is_double)
 *   integer, intent(IN) :: what, datatype, is_byte, is_short, is_double
 * 
 * processor is the address of the processing function called as 
 *   decode_function(void *data, int npoints)
 *   n_missing_points = encode_function(void *dest, void *src, int npoints)
 * NOTE:
 *   the npoints argument is passed BY VALUE (Fortran beware, iso_c_binding necessary)
 * =================================================================================
 * =================================================================================
 * DecodeMissingValue / decode_missing_value
 * =================================================================================
 * 
 * this subroutine will replace all occurrences the largest value in array field by the "missing data"
 * value for this kind (datatype) of data. it assumes that EncodeMissingValue / encode_missing_value
 * has been called before to set the largest value in an adequate manner.
 * 
 * C interface:
 *   void DecodeMissingValue(void *field,int nvalues,int datatype,int is_byte,int is_short,int is_double);
 * 
 * FORTRAN interface:
 *   subroutine decode_missing_value(field,nvalues,datatype,is_byte,is_short,is_double)
 *   "any_type", intent(INOUT) :: field
 *   integer, intent(IN) :: nvalues,datatype,is_byte,is_short,is_double
 * 
 * see EncodeMissingValue / encode_missing_value for the meaning of the arguments
 * =================================================================================
 * 
 * NOTE: this set of routines does not implement full defensive coding. SAFE calls from the
 *       standard file package are expected. USE WITH CARE.
*/
#include <stdlib.h>
#include <stdio.h>

#include <dlfcn.h>
void *DlOpen(const char *filename, int flag);
void *DlSym(void *handle, const char *symbol);

int msg_level;

#include "qstdir.h"

static int plugmode=-1;                             /* -1 not initialized, 1 active, 0 not active */
/* default values for missing values */
static float float_missing_val=-1.0E+38;          /* very large float */
static double double_missing_val=-1.0E+38;       /* very large double */
static int int_missing_val=0x80000000;              /* largest negative 32 bit integer */
static short short_missing_val=0x8000;              /* largest negative 16 bit integer */
static signed char byte_missing_val=0x80;           /* largest negative  8 bit integer */
static unsigned int uint_missing_val=0xFFFFFFFF;    /* largest 32 bit unsigned integer */
static unsigned short ushort_missing_val=0xFFFF;    /* largest 16 bit unsigned integer */
static unsigned char ubyte_missing_val=0xFF;       /* largest  8 bit unsigned integer */

void SetMissingValueMapping(int what, int datatype, void *processor_, int is_byte, int is_short, int is_double);

static void (*set_plugin_missing_value_flags)() = NULL ;

/* Fortran and C callable versions get get "magic values" used flag */
/* set magic values used to flag missing data points from environment variable FST_MISSING_VALUE if it exists */
/* returns 1 if mode active, 0 otherwise */
#pragma weak missing_value_used_ = missing_value_used
#pragma weak missing_value_used__ = missing_value_used
int missing_value_used_();
int missing_value_used__();
int missing_value_used()  /* return 1 if missing value detected, 0 otherwise */
{
  char *text;
  void *handle;
  if(plugmode==-1){
    text=getenv("MISSING_VALUE_FLAGS");
    if(text==NULL){
      plugmode = 0 ;
    }else{
      plugmode = 1 ;
      sscanf(text,"%g %d %u %lg %hd %hu %hhd %hhu",
	     &float_missing_val,
	     &int_missing_val,
	     &uint_missing_val,
	     &double_missing_val,
	     &short_missing_val,
	     &ushort_missing_val,
	     &byte_missing_val,
	     &ubyte_missing_val
	    );
    }
    /* process dynamic plugins here to get encoding/decoding functions from dynamic library */
	/* environment variable MISSING_VALUE_PLUGINS gives the name of the library */
	/* 
	   functions must have the following names in the shared library :
	   
       float_encode
       double_encode
       int_encode
       short_encode
       byte_encode
       uint_encode
       ushort_encode
       ubyte_encode
  
       float_decode
       double_decode
       int_decode
       short_decode
       byte_decode
       uint_decode
       ushort_decode
       ubyte_decode
       
       set_plugin_missing_value_flags
	 */
	text=getenv("MISSING_VALUE_PLUGINS");
	if(text != NULL ){
	  fprintf(stderr,"INFO: opening plugin library '%s'\n",text);
	  handle = DlOpen(text,RTLD_NOW);
	  if(handle != NULL){
		SetMissingValueMapping(1,1,DlSym(handle,"float_decode"),0,0,0);
		SetMissingValueMapping(1,1,DlSym(handle,"double_decode"),0,0,1);
		SetMissingValueMapping(1,2,DlSym(handle,"uint_decode"),0,0,0);
		SetMissingValueMapping(1,2,DlSym(handle,"ubyte_decode"),1,0,0);
		SetMissingValueMapping(1,2,DlSym(handle,"ushort_decode"),0,1,0);
		SetMissingValueMapping(1,4,DlSym(handle,"int_decode"),0,0,0);
		SetMissingValueMapping(1,4,DlSym(handle,"byte_decode"),1,0,0);
		SetMissingValueMapping(1,4,DlSym(handle,"short_decode"),0,1,0);

		SetMissingValueMapping(2,1,DlSym(handle,"float_encode"),0,0,0);
		SetMissingValueMapping(2,1,DlSym(handle,"double_encode"),0,0,1);
		SetMissingValueMapping(2,2,DlSym(handle,"uint_encode"),0,0,0);
		SetMissingValueMapping(2,2,DlSym(handle,"ubyte_encode"),1,0,0);
		SetMissingValueMapping(2,2,DlSym(handle,"ushort_encode"),0,1,0);
		SetMissingValueMapping(2,4,DlSym(handle,"int_encode"),0,0,0);
		SetMissingValueMapping(2,4,DlSym(handle,"byte_encode"),1,0,0);
		SetMissingValueMapping(2,4,DlSym(handle,"short_encode"),0,1,0);
                set_plugin_missing_value_flags = (void(*)()) DlSym(handle,"set_plugin_missing_value_flags");
	  }
	  else {
		fprintf(stderr,"WARNING: plugin library '%s' not found\n",text);
	  }
	}
    if(set_plugin_missing_value_flags != NULL ) (*set_plugin_missing_value_flags)(
	     &float_missing_val,
	     &int_missing_val,
	     &uint_missing_val,
	     &double_missing_val,
	     &short_missing_val,
	     &ushort_missing_val,
	     &byte_missing_val,
	     &ubyte_missing_val
												  );
  }
  return(plugmode);
}

/* C entry point */
/* forcibly set(activate) or reset(deactivate) the missing value mode flag */
int ForceMissingValueUsage(int flag) {
  plugmode = (flag!=0) ? 1 : 0;
  return(plugmode);
}

/* Fortran entry points below (calling C entry point ForceMissingValueUsage) */
/* forcibly set(activate) (flag !=0) or reset(deactivate) (flag=0) the "magic value" mode */
#pragma weak force_missing_value_used_ = force_missing_value_used
#pragma weak force_missing_value_used__ = force_missing_value_used
int force_missing_value_used_(int *flag);
int force_missing_value_used__(int *flag);
int force_missing_value_used(int *flag) { return(ForceMissingValueUsage(*flag)) ; }

/* Fortran and C callable versions */
/* get magic values used to flag missing data points from environment variable FST_MISSING_VALUE */
/* format is "real int unsigne_dint double */
/* return value of function tells whether the the "magic value" mode is active (1) or not(0) */
#pragma weak get_missing_value_flags_  = get_missing_value_flags
#pragma weak get_missing_value_flags__ = get_missing_value_flags
int get_missing_value_flags_(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub);
int get_missing_value_flags__(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub);
int get_missing_value_flags(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub)
{
  int status=missing_value_used();
  *f  = float_missing_val;
  *i  = int_missing_val;
  *ui = uint_missing_val;
  *d  = double_missing_val;
  *s  = short_missing_val;
  *us = ushort_missing_val;
  *b  = byte_missing_val;
  *ub = ubyte_missing_val;
  return(status);
}

/* Fortran and C callable versions, set magic values used to flag missing data points */
#pragma weak set_missing_value_flags_  = set_missing_value_flags
#pragma weak set_missing_value_flags__ = set_missing_value_flags
void set_missing_value_flags_(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub);
void set_missing_value_flags__(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub);
void set_missing_value_flags(float *f, int *i, unsigned int *ui, double *d, short *s, unsigned short *us,
                                      signed char *b, unsigned char *ub)
{
  float_missing_val  = *f ;
  int_missing_val    = *i ;
  uint_missing_val   = *ui ;
  double_missing_val = *d ;
  short_missing_val  = *s ;
  ushort_missing_val = *us ;
  byte_missing_val   = *b ;
  ubyte_missing_val  = *ub ;
  plugmode = 1 ;  /* values have been set, activate plug mode */
  if(set_plugin_missing_value_flags != NULL ) (*set_plugin_missing_value_flags)(f, i, ui, d, s, us,b, ub);
}

#ifdef USE_MACROS
/*
 * functions to find minimum and maximum values of a field while ignoring the missing value code
 * these functions return the number of missing values detected
*/

#define FLD_TYPE_ANAL(name,kind,flag)                          \
static int name(kind *z, int n , kind *zmax, kind *zmin)  {    \
  kind zma, zmi;                                               \
  int i, missing;                                              \
  i = 0; missing = 0;                                          \
  while( i<n && z[i]==flag ) {i++ ; missing++ ; }              \
  zma = zmi = z[i];                                            \
  while(++i < n) {                                             \
    if( z[i] == flag ) {                                       \
      missing++;                                               \
    }else{                                                     \
      zmi = z[i] < zmi ? z[i] : zmi;                           \
      zma = z[i] > zma ? z[i] : zma;                           \
    }                                                          \
  }                                                            \
  *zmax = zma ; *zmin = zmi ;                                  \
  return(missing);                                             \
}

FLD_TYPE_ANAL(fld_float_anal,  float,          float_missing_val)
FLD_TYPE_ANAL(fld_double_anal, double,         double_missing_val)
FLD_TYPE_ANAL(fld_int_anal,    int,            int_missing_val)
FLD_TYPE_ANAL(fld_uint_anal,   unsigned int,   uint_missing_val)
FLD_TYPE_ANAL(fld_short_anal,  short,          short_missing_val)
FLD_TYPE_ANAL(fld_ushort_anal, unsigned short, ushort_missing_val)
FLD_TYPE_ANAL(fld_byte_anal,   char,           byte_missing_val)
FLD_TYPE_ANAL(fld_ubyte_anal,  unsigned char,  ubyte_missing_val)

/*
 * routines to replace maximum value of field with missing value code
 */

#define FLD_TYPE_DECODE(name,kind,kind_anal,flag)               \
static void name(kind *z, int n) {                              \
  kind zma, zmi;                                                \
  int notused;                                                  \
  if(missing_value_used()==0) return;                           \
  notused = kind_anal(z,n,&zma,&zmi);                           \
  while(n--) { if(*z==zma) *z=flag ; z++ ; }                    \
  return;                                                       \
}

FLD_TYPE_DECODE(fst_double_decode_missing, double,         fld_double_anal, double_missing_val)
FLD_TYPE_DECODE(fst_float_decode_missing,  float,          fld_float_anal,  float_missing_val)
FLD_TYPE_DECODE(fst_int_decode_missing,    int,            fld_int_anal,    int_missing_val)
FLD_TYPE_DECODE(fst_uint_decode_missing,   unsigned int,   fld_uint_anal,   uint_missing_val)
FLD_TYPE_DECODE(fst_short_decode_missing,  short,          fld_short_anal,  short_missing_val)
FLD_TYPE_DECODE(fst_ushort_decode_missing, unsigned short, fld_ushort_anal, ushort_missing_val)
FLD_TYPE_DECODE(fst_byte_decode_missing,   char,           fld_byte_anal,   byte_missing_val)
FLD_TYPE_DECODE(fst_ubyte_decode_missing,  unsigned char,  fld_ubyte_anal,  ubyte_missing_val)

#else

static int fld_float_anal(float *z, int n , float *zmax, float *zmin)  /* float values */
{
  float zma, zmi;
  int missing=0;
  
  if(*z == float_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==float_missing_val) zma=*z;
    if(zmi==float_missing_val) zmi=*z;
    if(*z == float_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_double_anal(double *z, int n , double *zmax, double *zmin)  /* double values */
{
  double zma, zmi;
  int missing=0;
  
  if(*z == double_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==double_missing_val) zma=*z;
    if(zmi==double_missing_val) zmi=*z;
    if(*z == double_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_int_anal(int *z, int n , int *zmax, int *zmin) /* integer values */
{
  int zma, zmi;
  int missing=0;
  
  if(*z == int_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==int_missing_val) zma=*z;
    if(zmi==int_missing_val) zmi=*z;
    if(*z == int_missing_val) { missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_short_anal(short *z, int n , short *zmax, short *zmin) /* short signed integer values */
{
  short zma, zmi;
  int missing=0;
  
  if(*z == short_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==short_missing_val) zma=*z;
    if(zmi==short_missing_val) zmi=*z;
    if(*z == short_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_byte_anal(signed char *z, int n , signed char *zmax, signed char *zmin) /* byte signed integer values */
{
  signed char zma, zmi;
  int missing=0;
  
  if(*z == byte_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==byte_missing_val) zma=*z;
    if(zmi==byte_missing_val) zmi=*z;
    if(*z == byte_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_uint_anal(unsigned int *z, int n , unsigned int *zmax, unsigned int *zmin) /* unsigned integers */
{
  unsigned int zma, zmi;
  int missing=0;
  
  if(*z == uint_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==uint_missing_val) zma=*z;
    if(zmi==uint_missing_val) zmi=*z;
    if(*z == uint_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_ushort_anal(unsigned short *z, int n , unsigned short *zmax, unsigned short *zmin) /* short unsigned integer values */
{
  unsigned short zma, zmi;
  int missing=0;
  
  if(*z == ushort_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==ushort_missing_val) zma=*z;
    if(zmi==ushort_missing_val) zmi=*z;
    if(*z == ushort_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}

static int fld_ubyte_anal(unsigned char *z, int n , unsigned char *zmax, unsigned char *zmin) /* byte unsigned integer values */
{
  unsigned char zma, zmi;
  int missing=0;
  
  if(*z == ubyte_missing_val)missing++;
  zma=*z ;
  zmi = zma ;
  n-- ;
  while(n-- > 0) /* loop from second to last element of z */
  {
    z++ ;
    if(zma==ubyte_missing_val) zma=*z;
    if(zmi==ubyte_missing_val) zmi=*z;
    if(*z == ubyte_missing_val) {missing ++ ; continue; }
    if(*z < zmi) zmi = *z ;
    if(*z > zma) zma = *z ;
  }
//  if(zma == 0 && zmi == 0 ) zma = 1 ;
  *zmax = zma ;
  *zmin = zmi ;
  return(missing);
}
/* routines to replace maximum value of field with missing value code */
static void fst_double_decode_missing(double *z, int n)  /* float values */
{
  double zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_double_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=double_missing_val ; z++ ; } /* float values */
  return;
}

static void fst_float_decode_missing(float *z, int n)  /* float values */
{
  float zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_float_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=float_missing_val ; z++ ; } /* float values */
  return;
}

static void fst_int_decode_missing(int *z, int n)  /* signed ints */
{
  int zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_int_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=int_missing_val ; z++ ; }
  return;
}

static void fst_short_decode_missing(short *z, int n) /* signed shorts */
{
  short zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_short_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=short_missing_val ; z++ ; } /* signed shorts */
  return;
}

static void fst_byte_decode_missing(signed char *z, int n) /* signed bytes */
{
  signed char zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_byte_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=byte_missing_val ; z++ ; } /* signed bytes */
  return;
}

static void fst_uint_decode_missing(unsigned int *z, int n) /* unsigned integers */
{
  unsigned int zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_uint_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=uint_missing_val ; z++ ; } /* unsigned integers */
  return;
}

static void fst_ushort_decode_missing(unsigned short *z, int n) /* unsigned shorts */
{
  unsigned short zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_ushort_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=ushort_missing_val ; z++ ;  } /* unsigned shorts */
  return;
}

static void fst_ubyte_decode_missing(unsigned char *z, int n) /* unsigned bytes */
{
  unsigned char zma, zmi;
  int notused;

  if(missing_value_used()==0) return;
  notused = fld_ubyte_anal(z,n,&zma,&zmi);
  while(n--) { if(*z==zma) *z=ubyte_missing_val ; z++ ;  } /* unsigned bytes */
  return;
}
#endif
/* functions te replace missing value code with value suitable for packers */
static int fst_double_encode_missing(double *z, double *z2, int n, int nbits)  /* double values */
{
  double zma, zmi;
  double plug;
  float factor[9] = { 2.0, 1.1, .55, .28, .14, .07, .035, .017, .01 };
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_double_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);
  
  if(nbits >8) /* >8 bits, add 1% of max-min to max */
    plug = zma + (zma-zmi)*.01;
  else         /* <=8 bits, add appropriate fraction of max-min to max */
    plug = zma + (zma-zmi)*factor[nbits] ;

  if(plug == zma){   /* field is constant, zma == zmi */
    plug = (zma == 0.0) ? 1.0 : zma*2.0 ;
  }
  while(n--) { if(*z2==double_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_float_encode_missing(float *z, float *z2, int n, int nbits)  /* float values */
{
  float zma, zmi;
  float plug;
  float factor[9] = { 2.0, 1.1, .55, .28, .14, .07, .035, .017, .01 };
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_float_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  if(nbits >8) /* >8 bits, add 1% of max-min to max */
    plug = zma + (zma-zmi)*.01;
  else         /* <=8 bits, add appropriate fraction of max-min to max */
    plug = zma + (zma-zmi)*factor[nbits] ;

  if(plug == zma){   /* field is constant, zma == zmi */
    plug = (zma == 0.0) ? 1.0 : zma*2.0 ;
  }
  while(n--) { if(*z2==float_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_int_encode_missing(int *z, int *z2, int n, int nbits) /* signed integers */
{
  int zma, zmi;
  int plug;
  unsigned int plug2;
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_int_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  plug2 = 0xFFFFFFFF;
  if(nbits<32) plug2 = plug2 >> (32-nbits);
  plug2 = plug2 >> 1 ;  /* largest positive number using nbits bits */
  
  plug = plug2 ;  /* this code assumes packing of z with no bias removal */
  if(plug > zma) plug = zma + 1;
  else WARNPRINT fprintf(stderr,"fst_int_encode_missing: << WARNING >> MAXIMUM VALUE >= ENCODED MISSING VALUE FLAG\n");
  while(n--) { if(*z2==int_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
  
  if(nbits<32){  /* this code assumes packing of z-zmi (bias removal) */
    plug = plug2  ;
    if(plug > (zma-zmi)) plug = zma + 1; else plug = zma + (zma-zmi);  /* double the range if zma-zmi+1 will not fit in nbits */
  }else{
    plug = zma + 1 ;  /* best effort !! */
  }
  while(n--) { if(*z2==int_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_short_encode_missing(short *z, short *z2, int n, int nbits) /* signed integers */
{
  short zma, zmi;
  short plug;
  unsigned short plug2;
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_short_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  plug2 = 0xFFFF;
  if(nbits<16) plug2 = plug2 >> (16-nbits);
  plug2 = plug2 >> 1 ;  /* largest positive number using nbits bits */

  plug = plug2 ;  /* this code assumes packing of z with no bias removal */
  if(plug > zma) plug = zma + 1;
  else WARNPRINT fprintf(stderr,"fst_short_encode_missing: << WARNING >> MAXIMUM VALUE >= ENCODED MISSING VALUE FLAG\n");
  while(n--) { if(*z2==short_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);

  if(nbits<16){  /* this code assumes packing of z-zmi (bias removal) */
    plug = plug2 ;
    if(plug > (zma-zmi)) plug = zma + 1; else plug = zma + (zma-zmi);  /* double the range if zma-zmi+1 will not fit in nbits */
  }else{
    plug = zma + 1 ;  /* best effort !! */
  }
  while(n--) { if(*z2==short_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_byte_encode_missing(signed char *z, signed char *z2, int n, int nbits) /* signed bytes */
{
  signed char zma, zmi;
  signed char plug;
  unsigned char plug2;
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_byte_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  plug2 = 0xFF;
  if(nbits<8) plug2 = plug2 >> (8-nbits);
  plug2 = plug2 >> 1 ;  /* largest positive number using nbits bits */

  plug = plug2 ;  /* this code assumes packing of z with no bias removal */
  if(plug > zma) plug = zma + 1;
  else WARNPRINT fprintf(stderr,"fst_byte_encode_missing: << WARNING >> MAXIMUM VALUE >= ENCODED MISSING VALUE FLAG\n");
  while(n--) { if(*z2==byte_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_uint_encode_missing(unsigned int *z, unsigned int *z2, int n, int nbits) /* unsigned integers */
{
  unsigned int zma, zmi;
  unsigned int plug;
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_uint_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  plug = 0xFFFFFFFF;
  if(nbits<32) plug = plug >> (32-nbits);

  if(plug > zma) plug = zma + 1;  /* this code assumes packing of z with no bias removal */
  else WARNPRINT fprintf(stderr,"fst_uint_encode_missing: << WARNING >> MAXIMUM VALUE %u >= ENCODED MISSING VALUE FLAG %u\n",zma,plug);
  while(n--) { if(*z2==uint_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);

  if(nbits<32){  /* this code assumes packing of z-zmi (bias removal) */
    if(plug > (zma-zmi)) plug = zma + 1; else plug = zma + (zma-zmi);  /* double the range if zma-zmi+1 will not fit in nbits */
  }else{
    plug = zma + 1 ;  /* best effort !! */
  }
  while(n--) { if(*z2==uint_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_ushort_encode_missing(unsigned short *z, unsigned short *z2, int n, int nbits) /* unsigned integers */
{
  unsigned short zma, zmi;
  unsigned short plug;
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_ushort_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  plug = 0xFFFF;
  if(nbits<16) plug = plug >> (16-nbits);

  if(plug > zma) plug = zma + 1;  /* this code assumes packing of z with no bias removal */
  else WARNPRINT fprintf(stderr,"fst_ushort_encode_missing: << WARNING >> MAXIMUM VALUE %hu >= ENCODED MISSING VALUE FLAG %hu\n",zma,plug);
  while(n--) { if(*z2==ushort_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);

  if(nbits<16){  /* this code assumes packing of z-zmi (bias removal) */
    if(plug > (zma-zmi)) plug = zma + 1; else plug = zma + (zma-zmi);  /* double the range if zma-zmi+1 will not fit in nbits */
  }else{
    plug = zma + 1 ;
  }
  while(n--) { if(*z2==ushort_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

static int fst_ubyte_encode_missing(unsigned char *z, unsigned char *z2, int n, int nbits) /* unsigned bytes */
{
  unsigned char zma, zmi;
  unsigned char plug;
  int missing;
  
  if(missing_value_used()==0) return(0);
  missing = fld_ubyte_anal(z2,n,&zma,&zmi);
  if(missing==0) return(0);

  plug = 0xFF;
  if(nbits<8) plug = plug >> (8-nbits);

  if(plug > zma) plug = zma + 1;  /* this code assumes packing of z with no bias removal */
  else WARNPRINT fprintf(stderr,"fst_ubyte_encode_missing: << WARNING >> MAXIMUM VALUE %hu >= ENCODED MISSING VALUE FLAG %hu\n",zma,plug);
  while(n--) { if(*z2==ubyte_missing_val) *z=plug ; else *z=*z2 ; z++ ; z2++ ; }
  return(missing);
}

#ifdef USE_MACROS
#else
#endif
static void fst_null_decode_missing(void *z, int n) /* null decoder */
{
/*  fprintf(stderr,"using null decoder\n"); */
  return;
}

static int (*__fst_float_encode_missing)() = fst_float_encode_missing;
static int (*__fst_double_encode_missing)() = fst_double_encode_missing;
static int (*__fst_int_encode_missing)() = fst_int_encode_missing;
static int (*__fst_short_encode_missing)() = fst_short_encode_missing;
static int (*__fst_byte_encode_missing)() = fst_byte_encode_missing;
static int (*__fst_uint_encode_missing)() = fst_uint_encode_missing;
static int (*__fst_ushort_encode_missing)() = fst_ushort_encode_missing;
static int (*__fst_ubyte_encode_missing)() = fst_ubyte_encode_missing;

static void (*__fst_float_decode_missing)() = fst_float_decode_missing;
static void (*__fst_double_decode_missing)() = fst_double_decode_missing;
static void (*__fst_int_decode_missing)() = fst_int_decode_missing;
static void (*__fst_short_decode_missing)() = fst_short_decode_missing;
static void (*__fst_byte_decode_missing)() = fst_byte_decode_missing;
static void (*__fst_uint_decode_missing)() = fst_uint_decode_missing;
static void (*__fst_ushort_decode_missing)() = fst_ushort_decode_missing;
static void (*__fst_ubyte_decode_missing)() = fst_ubyte_decode_missing;

/* C entry point , restore mapping functions to their default value */
void RestoreMissingValueMapping(void)
{
__fst_float_encode_missing = fst_float_encode_missing;
__fst_double_encode_missing = fst_double_encode_missing;
__fst_int_encode_missing = fst_int_encode_missing;
__fst_short_encode_missing = fst_short_encode_missing;
__fst_byte_encode_missing = fst_byte_encode_missing;
__fst_uint_encode_missing = fst_uint_encode_missing;
__fst_ushort_encode_missing = fst_ushort_encode_missing;
__fst_ubyte_encode_missing = fst_ubyte_encode_missing;
  
__fst_float_decode_missing = fst_float_decode_missing;
__fst_double_decode_missing = fst_double_decode_missing;
__fst_int_decode_missing = fst_int_decode_missing;
__fst_short_decode_missing = fst_short_decode_missing;
__fst_byte_decode_missing = fst_byte_decode_missing;
__fst_uint_decode_missing = fst_uint_decode_missing;
__fst_ushort_decode_missing = fst_ushort_decode_missing;
__fst_ubyte_decode_missing = fst_ubyte_decode_missing;
}

/* Fortran entry points (calling C entry point RestoreMissingValueMapping) */
void restore_missing_value_mapping(void) { RestoreMissingValueMapping() ; }
void restore_missing_value_mapping_(void) { RestoreMissingValueMapping() ; }
void restore_missing_value_mapping__(void) { RestoreMissingValueMapping() ; }

/* C entry point */
/* change the missing value mapping function */
/* mode >0, set function to processor */
/* mode <0, set function to original value, ignore processor */
/* abs(mode)==11 selectively deactivate the selected decoding function (uses fst_null_decode_missing), ignore processor */
/* abs(mode)==1 decoding functions */
/* abs(mode)==2 encoding functions */
/* datatype, is_byte, is_short, is_double , same as elsewhere in this code */
void SetMissingValueMapping(int what, int datatype, void *processor_, int is_byte, int is_short, int is_double)
{
    void *processor = processor_;
    int mode;
    if(what > 0 && processor == NULL) return ; /* null pointer to function, return */
    mode = (what>0) ? what : -what;
    if(mode == 11) {  /* mode 11 : deactivate decoder for specified type */
      mode = 1 ;  /* reduce mode to 1  */
      processor = (void *)fst_null_decode_missing ;
    }
    if(mode==1){  /* replace a decoding routine */
      if(datatype==1 || datatype==5 || datatype==6) { /* float or IEEE */
        if(is_double) {
          __fst_double_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_double_decode_missing ; /* real or IEEE */
        }else{
          __fst_float_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_float_decode_missing ; /* real or IEEE */
        }
      }
      if(datatype==4) {  /* signed */
        if(is_short) {
          __fst_short_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_short_decode_missing ;   /* signed shorts */
        }else if(is_byte){
          __fst_byte_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_byte_decode_missing ;   /* signed bytes */
        }else{
          __fst_int_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_int_decode_missing ;   /* signed integers */
        }
      }
      if(datatype==2) {  /* unsigned */
        if(is_short) {
          __fst_ushort_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_ushort_decode_missing ;   /* uns  SetMissingValueMapping(1,4,NULL,0,0,0);
          igned shorts */
        }else if(is_byte){
          __fst_ubyte_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_ubyte_decode_missing ;   /* unsigned bytes */
        }else{
          __fst_uint_decode_missing = (what>0) ? (void(*)()) processor : (void(*)()) fst_uint_decode_missing ;   /* unsigned integers */
        }
      }
    }
    if(mode==2){           /* replace an encoding routine */
      if(datatype==1 || datatype==5 || datatype==6) { /* float or IEEE */
         if(is_double) {
           __fst_double_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_double_encode_missing ; /* real or IEEE */
         }else{
           __fst_float_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_float_encode_missing ; /* real or IEEE */
         }
      }
      if(datatype==4) {  /* signed */
        if(is_short) {
          __fst_short_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_short_encode_missing ;   /* signed shorts */
        }else if(is_byte){
          __fst_byte_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_byte_encode_missing ;   /* signed bytes */
        }else{
          __fst_int_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_int_encode_missing ;   /* signed integers */
        }
      }
      if(datatype==2) {  /* unsigned */
        if(is_short) {
          __fst_ushort_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_ushort_encode_missing ;   /* unsigned shorts */
        }else if(is_byte){
          __fst_ubyte_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_ubyte_encode_missing ;   /* unsigned bytes */
        }else{
          __fst_uint_encode_missing = (what>0) ? (int(*)()) processor : (int(*)()) fst_uint_encode_missing ;   /* unsigned integers */
        }
      }
    }
}

/* Fortran entry points (calling C entry point SetMissingValueMapping) */
#pragma weak set_missing_value_mapping_ = set_missing_value_mapping
#pragma weak set_missing_value_mapping__ = set_missing_value_mapping
void set_missing_value_mapping_(int *what, int *datatype, void *processor_, int *is_byte, int *is_short, int *is_double);
void set_missing_value_mapping__(int *what, int *datatype, void *processor_, int *is_byte, int *is_short, int *is_double);
void set_missing_value_mapping(int *what, int *datatype, void *processor_, int *is_byte, int *is_short, int *is_double) {
  SetMissingValueMapping(*what, *datatype, processor_, *is_byte, *is_short, *is_double);
}

/*
int EncodeMissingValue(void *field,void *field2,int nvalues,int datatype,int nbits,int is_byte,int is_short,int is_double);
void DecodeMissingValue(void *field,int nvalues,int datatype,int is_byte,int is_short,int is_double);
*/

/* C entry point */
int EncodeMissingValue(void *field,void *field2,int nvalues,int datatype,int nbits,int is_byte,int is_short,int is_double){
  int missing = 0;

  if(missing_value_used()==0) return(0);  /* "magic value" mode off, return zero missing value found */
  datatype &= 0xF ;
  if(datatype==0 || datatype==3 || datatype==7 || datatype==8) return(0) ; /* not valid for transparent or character types */
  if(datatype==1 || datatype==5 || datatype==6) { /* float or IEEE */
    if(datatype==5 && nbits==64) is_double=1;
    if(is_double) {
      missing = (*__fst_double_encode_missing)(field,field2,nvalues,nbits); /* real or IEEE */
    }else{
      if(nbits>32)  return(0) ; /* datalength > 32 not supported */
        missing = (*__fst_float_encode_missing)(field,field2,nvalues,nbits); /* real or IEEE */
    }
  }
  if(datatype==4) {  /* signed */
    if(is_short) {
      missing = (*__fst_short_encode_missing)(field,field2,nvalues,nbits);   /* signed shorts */
    }else if(is_byte){
      missing = (*__fst_byte_encode_missing)(field,field2,nvalues,nbits);   /* signed bytes */
    }else{
      missing = (*__fst_int_encode_missing)(field,field2,nvalues,nbits);   /* signed integers */
    }
  }
  if(datatype==2) {  /* unsigned */
    if(is_short) {
      missing = (*__fst_ushort_encode_missing)(field,field2,nvalues,nbits);  /* unsigned shorts */
    }else if(is_byte){
      missing = (*__fst_ubyte_encode_missing)(field,field2,nvalues,nbits);  /* unsigned shorts */
    }else{
      missing = (*__fst_uint_encode_missing)(field,field2,nvalues,nbits);  /* unsigned integers */
    }
  }
  if (0 >= msg_level) {
    fprintf(stderr,"DEBUG: %d missing values in %d data values replaced, base datatype=%d\n",missing,nvalues,datatype);
  }
  return missing; /* number of missing values found */
}

/* Fortran entry points (calling C entry point EncodeMissingValue) */
#pragma weak encode_missing_value__ = encode_missing_value
#pragma weak encode_missing_value_  = encode_missing_value
int encode_missing_value_(void *field,void *field2,int *nvalues,int *datatype,int *nbits,int *is_byte,int *is_short,int *is_double);
int encode_missing_value__(void *field,void *field2,int *nvalues,int *datatype,int *nbits,int *is_byte,int *is_short,int *is_double);
int encode_missing_value(void *field,void *field2,int *nvalues,int *datatype,int *nbits,int *is_byte,int *is_short,int *is_double){
  return( EncodeMissingValue(field,field2,*nvalues,*datatype,*nbits,*is_byte,*is_short,*is_double) );
}

/* C entry point */
void DecodeMissingValue(void *field,int nvalues,int datatype,int is_byte,int is_short,int is_double){
  if(missing_value_used()==0) return ;  /* "magic value" mode off, do nothing  */
  datatype &= 0xF ;
  if(datatype==0 || datatype==3 || datatype==7 || datatype==8) return ; /* not valid for complex, transparent or character types */
  if(datatype==1 || datatype==5 || datatype==6) { /* real or IEEE */
    if(is_double) {
      (*__fst_double_decode_missing)(field,nvalues);
    }else{
      (*__fst_float_decode_missing)(field,nvalues);
    }
  }
  if(datatype==4) {
    if(is_short) {
      (*__fst_short_decode_missing)(field,nvalues);   /* signed shorts */
    }else if(is_byte){
      (*__fst_byte_decode_missing)(field,nvalues);   /* signed shorts */
    }else{
      (*__fst_int_decode_missing)(field,nvalues);     /* signed integers */
    }
  }
  if(datatype==2) {
    if(is_short) {
      (*__fst_ushort_decode_missing)(field,nvalues);  /* unsigned shorts */
    }else if(is_byte){
      (*__fst_ubyte_decode_missing)(field,nvalues);  /* unsigned shorts */
    }else{
      (*__fst_uint_decode_missing)(field,nvalues);    /* unsigned integers */
    }
  }
}

/* Fortran entry points (calling C entry point DecodeMissingValue) */
#pragma weak decode_missing_value__ = decode_missing_value
#pragma weak decode_missing_value_  = decode_missing_value
void decode_missing_value_(void *field,int *nvalues,int *datatype,int *is_byte,int *is_short,int *is_double);
void decode_missing_value__(void *field,int *nvalues,int *datatype,int *is_byte,int *is_short,int *is_double);
void decode_missing_value(void *field,int *nvalues,int *datatype,int *is_byte,int *is_short,int *is_double){
  DecodeMissingValue(field,*nvalues,*datatype,*is_byte,*is_short,*is_double) ;
}

#ifdef SELFTEST
#define ASIZE 16
/* test program , the 8 basic types are encoded/decoded
 * float, int, unsigned int, double, short, unsigned short, char, unsigned char
 * both the C and the FORTRAN interface are used in the test
*/
int main()
{
  float f, fa[ASIZE], fp[ASIZE];
  double d, da[ASIZE], dp[ASIZE];
  int i, ii, ia[ASIZE], ip[ASIZE];
  unsigned int ui, uia[ASIZE], uip[ASIZE];
  int junk;
  short s, sa[ASIZE], sp[ASIZE];
  unsigned short us, usa[ASIZE], usp[ASIZE];
  signed char b, ba[ASIZE], bp[ASIZE];
  unsigned char ub, uba[ASIZE], ubp[ASIZE];
  int asize=ASIZE;
  int zero=0;
  int one=1;
  int two=2;
  int four=4;
  int sixteen=16;
  int fstat,istat,uistat,dstat,sstat,usstat,bstat,ubstat;

  RestoreMissingValueMapping();
#ifdef NODECODE
  SetMissingValueMapping(1,1,fst_null_decode_missing,0,0,0);
  SetMissingValueMapping(1,1,fst_null_decode_missing,0,0,1);
  SetMissingValueMapping(1,2,fst_null_decode_missing,0,0,0);
  SetMissingValueMapping(1,2,fst_null_decode_missing,1,0,0);
  SetMissingValueMapping(1,2,fst_null_decode_missing,0,1,0);
  SetMissingValueMapping(1,4,fst_null_decode_missing,0,0,0);
  SetMissingValueMapping(1,4,fst_null_decode_missing,1,0,0);
  SetMissingValueMapping(1,4,fst_null_decode_missing,0,1,0);
#endif
  SetMissingValueMapping(-1,2,NULL,0,0,0);  /* restore uint decoder */
  SetMissingValueMapping(-1,2,NULL,1,0,0);  /* restore ubyte decoder */
  SetMissingValueMapping(-1,2,NULL,0,1,0);  /* restore ushort decoder */
  for (i=0 ; i<ASIZE ; i++) {
    fa[i]= i*1.0+.5;
    da[i]=fa[i]+1.2345;
    ia[i]=i-13;
    uia[i]=i+22;
    sa[i]=i-13;
    usa[i]=i+22;
    ba[i]=i-13;
    uba[i]=i+22;
    bp[i]=i-13;
    ubp[i]=i+22;
  }
  junk=get_missing_value_flags(&f,&i,&ui,&d,&s,&us,&b,&ub);
  f = -1.0E+37; d = -1.0E+37;
  set_missing_value_flags(&f,&i,&ui,&d,&s,&us,&b,&ub);
  fprintf(stderr,"%d float=%g, int=%d, uint=%u, double=%lg, short=%hd, ushort=%hu, byte=%hhd, ubyte=%hhu\n",junk, f,i,ui,d,s,us,b,ub);
  junk = missing_value_used();
  junk=get_missing_value_flags(&f,&i,&ui,&d,&s,&us,&b,&ub);
  fprintf(stderr,"%d float=%g, int=%d, uint=%u, double=%lg, short=%hd, ushort=%hu, byte=%hhd, ubyte=%hhu\n",junk,f,i,ui,d,s,us,b,ub);
  for (ii=0 ; ii<ASIZE ; ii++) {
    fprintf(stderr,"%2d float=%12g, int=%12d, short=%8hd, byte=%8hhd, uint=%12u, ushort=%8hu, ubyte=%8hhu,  double=%12lg\n",ii,fa[ii],ia[ii],sa[ii],ba[ii],uia[ii],usa[ii],uba[ii],da[ii]);
  }
  fprintf(stderr,"================================================================\n");
  fa[0]=f ; fa[ASIZE-1]=f;
  fstat=encode_missing_value(fp,fa,&asize,&one,&sixteen,&zero,&zero,&zero);
  da[0]=d ; da[ASIZE-1]=d;
  dstat=encode_missing_value(dp,da,&asize,&one,&sixteen,&zero,&zero,&one);
  ia[1]=i ; ia[ASIZE-2]=i;
  istat=EncodeMissingValue(ip,ia,ASIZE,4,8,0,0,0);
  sa[2]=s ; sa[ASIZE-3]=s;
  sstat=EncodeMissingValue(sp,sa,ASIZE,4,8,0,1,0);
  ba[3]=b ; ba[ASIZE-4]=b;
  bstat=EncodeMissingValue(bp,ba,ASIZE,4,8,1,0,0);
  uia[2]=ui ; uia[ASIZE-3]=ui;
  uistat=EncodeMissingValue(uip,uia,ASIZE,2,8,0,0,0);
  usa[3]=us ; usa[ASIZE-4]=us;
  usstat=EncodeMissingValue(usp,usa,ASIZE,2,8,0,1,0);
  uba[4]=ub ; uba[ASIZE-5]=ub;
  ubstat=EncodeMissingValue(ubp,uba,ASIZE,2,8,1,0,0);
  fprintf(stderr,"status values: %d %d %d %d %d %d %d %d\n",fstat,dstat,istat,sstat,bstat,uistat,usstat,ubstat);
  for (ii=-1 ; ++ii<ASIZE ; ) {
    fprintf(stderr,"%2d float=%12g, int=%12d, short=%8hd, byte=%8hhd, uint=%12u, ushort=%8hu, ubyte=%8hhu,  double=%12lg\n",ii,fp[ii],ip[ii],sp[ii],bp[ii],uip[ii],usp[ii],ubp[ii],dp[ii]);
  }
  fprintf(stderr,"================================================================\n");
  DecodeMissingValue(fp,ASIZE,1,0,0,0);
  DecodeMissingValue(dp,ASIZE,1,0,0,1);
  DecodeMissingValue(ip,ASIZE,4,0,0,0);
  decode_missing_value(sp,&asize,&four,&zero,&one,&zero);
  DecodeMissingValue(bp,ASIZE,4,1,0,0);
  DecodeMissingValue(uip,ASIZE,2,0,0,0);
  decode_missing_value(usp,&asize,&two,&zero,&one,&zero);
  DecodeMissingValue(ubp,ASIZE,2,1,0,0);
  for (ii=0 ; ii<ASIZE ; ii++) {
    fprintf(stderr,"%2d float=%12g, int=%12d, short=%8hd, byte=%8hhd, uint=%12u, ushort=%8hu, ubyte=%8hhu,  double=%12lg\n",ii,fp[ii],ip[ii],sp[ii],bp[ii],uip[ii],usp[ii],ubp[ii],dp[ii]);
  }
  return(0);
}
#endif
