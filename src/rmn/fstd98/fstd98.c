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

#define use_old_signed_pack_unpack_code YES

#include <stdio.h>
#ifdef WIN32    /*CHC/NRC*/
#include <string.h>
#include <malloc.h>
#else
#include <unistd.h>
#include <alloca.h>
#endif
#include <stdlib.h>
#include "qstdir.h"
#include "convert_ip.h"
#include "proto.h"
#include <rmnlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <regex.h>

#define Max_Ipvals 50

static int ip_nb[3] = {0,0,0};
static int ip1s_flag = 0;
static int ip2s_flag = 0;
static int ip3s_flag = 0;
static int dejafait_1=0, dejafait_2=0;

static int ips_tab[3][Max_Ipvals];

static word link_list[1024];
static int link_n;
static int stdf_version = 200001;
static int downgrade_32=0;    /* downgrade 64 bit field to 32 bit field when reading */
static char *ARMNLIB=NULL;    /* ARMNLIB environment variable */
static char *debug_filename=NULL;   /* print on debug file when available */
static char *requetes_filename=NULL;  /* filter file, desire/exclure */
static int remap_table[2][10];    /* used for datatype remapping */
static int nb_remap=0;            /* number of datatype remapping,  0 = no remapping */
static char prnt_options[128]="NINJNK+DATESTAMPO+IP1+IG1234";   /* what is printed with fstecr */
static void crack_std_parms(stdf_dir_keys *stdf_entry,
        stdf_special_parms *cracked_parms);
static void print_std_parms(stdf_dir_keys *stdf_entry,
                            char *pre, char *option,int header);

static int kinds_table_init = 1;
static char kind_chars[96];

static char *kinds_table[] =
{
  "??", "??", "??", "??", "??", "??", "??", "??",
  "??", "??", "??", "??", "??", "??", "??", "??",
  "??", "??", "??", "??", "??", "??", "??", "??",
  "??", "??", "??", "??", "??", "??", "??", "??"
} ;


void KindToString(int kind, char *s1, char *s2);  /* fortran routine from comvertip_123 */
int EncodeMissingValue(void *field,void *field2,int nvalues,int datatype,int nbits,int is_byte,int is_short,int is_double);
void DecodeMissingValue(void *field,int nvalues,int datatype,int is_byte,int is_short,int is_double);

/*splitpoint aaaa_comment */
/*****************************************************************************
 *                                                                           *
 *                              F S T D 9 8                                  *
 *                                                                           *
 *             RPN standard file system software version 1998.               *
 *                                                                           *
 *Authors                                                                    *
 *  M. Lepine - M. Valin                                                     *
 *  Revised version based on XDF internal file structure.                    *
 *  Revised feb 2014, revisited ip translation and IP printing.              *
 *        print_std_parms, ConvipPlus, kinds                                 *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************
 *                        C O N V I P _ P L U S                              *
 *                                                                           *
 *Object                                                                     *
 *   Interface to IP conversion package                   .                  *
 *   ignore s and flag from old convip                    .                  *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN/OUT ip_new, kind, level, mode   (see ConvertIp)                       *
 *                                                                           *
 *****************************************************************************/

static void ConvipPlus(int *ip_new,float *level, int *kind, int *mode)
{
  ConvertIp(ip_new,level,kind,*mode);  /* ignore s and flag that are not used anyway */
}

/*****************************************************************************
 *                            K I N D S                                      *
 *                                                                           *
 *Object                                                                     *
 *   translate kind code to 2 char string                   .                *
 *   private function                                       .                *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  kind                                                                 *
 *                                                                           *
 *****************************************************************************/

static char *kinds(int kind){
static char *dummy="??";
  int i;
  if(kinds_table_init){       /* initialize table if first call */
    for(i=0 ; i<=31 ; i++) {
      KindToString(i,&kind_chars[3*i],&kind_chars[3*i+1]);
      kind_chars[3*i+2] = '\0' ;
    }
    kinds_table_init=0;
  }
  return &kind_chars[3*kind];
}



/*splitpoint c_fstapp */
/*****************************************************************************
 *                            C _ F S T A P P                                *
 *                                                                           *
 *Object                                                                     *
 *   Position at the end of a sequential file for an append.                 *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  IN  option  kept for backward compatibility                              *
 *                                                                           *
 *****************************************************************************/

int c_fstapp(int iun, char *option)
{
  int index_fnom, index, width, nw, end_of_file;
  file_table_entry *f;
  xdf_record_header *header;
  seq_dir_keys * seq_entry;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstapp",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstapp",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (!f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not sequential",iun);
    return(error_msg("c_fstapp",ERR_BAD_FTYPE,WARNING));
  }

  end_of_file = 0;
  width = W64TOWD(f->primary_len);

  if (f->fstd_vintage_89) {
    while (!end_of_file) {
      nw = c_waread2(iun,f->head_keys,f->cur_addr,width);
      header = (xdf_record_header *) f->head_keys;
      if (nw < width) {
        end_of_file = 1;
        header->idtyp = 127;
        header->lng = 1;
        break;
      }
      seq_entry = (seq_dir_keys *) f->head_keys;
      if (seq_entry->eof > 0) {
        header->idtyp = 112 + seq_entry->eof;
        header->lng = 1;
        end_of_file = 1;
        break;
      }
      header->lng = ((seq_entry->lng +3) >> 2) + 15;
      f->cur_addr += W64TOWD(header->lng);
    }
  }
  else {
    while (!end_of_file) {
      nw = c_waread2(iun,f->head_keys,f->cur_addr,width);
      header = (xdf_record_header *) f->head_keys;
      if (nw < W64TOWD(1)) {
        end_of_file = 1;
        header->idtyp = 127;
        header->lng = 1;
        break;
      }
      if ((header->idtyp >= 112) && (header->idtyp <= 127)) {
        end_of_file = 1;
        break;
      }
      f->cur_addr += W64TOWD(header->lng);
    }
  }
  f->nxtadr = f->cur_addr;
  return(0);
}

/*splitpoint c_fstckp */
/*****************************************************************************
 *                            C _ F S T C K P                                *
 *                                                                           *
 *Object                                                                     *
 *   Checkpoint. Clear buffers, rewrite headers.                             *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/
int c_fstckp(int iun)
{
  int index_fnom, index, ier;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstckp",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstckp",ERR_NO_FILE,ERROR));
  }

  xdf_checkpoint = 1;           /* checkpoint mode, not a complete close */
  ier = c_xdfcls(iun);
  return(ier);
}

/*splitpoint c_fst_data_length */
/*****************************************************************************
 *                  C _ F S T _ D A T A _ L E N G T H                        *
 *                                                                           *
 *Object                                                                     *
 *   Gives information on data lenght of the elements passed to fstecr       *
 *   and fstlir (double, short integer, byte ...)                            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  length_type     data length kind                                     *
 *                      1: byte                                              *
 *                      2: short (16 bits)                                   *
 *                      4: regular 32 bits                                   *
 *                      8: double (64 bits)                                  *
 *                                                                           *
 *****************************************************************************/
int c_fst_data_length(int length_type)
{

  switch (length_type) {

    case 1:
      xdf_byte = 1;
      xdf_short = 0;
      xdf_double = 0;
      break;

    case 2:
      xdf_byte = 0;
      xdf_short = 1;
      xdf_double = 0;

      break;

    case 4:
      xdf_byte = 0;
      xdf_short = 0;
      xdf_double = 0;
      break;

    case 8:
      xdf_byte = 0;
      xdf_short = 0;
      xdf_double = 1;
      break;

    default:
      fprintf(stderr,"c_fst_data_length invalid length type=%d",length_type);
      xdf_byte = 0;
      xdf_short = 0;
      xdf_double = 0;
      break;
  }
return(0);
}


/*
void *fst_encode_missing_value(void *field,void *field2,int nvalues,int datatype,int nbits,int is_byte,int is_short,int is_double);
void fst_decode_missing_value(void *field,int nvalues,int datatype,int is_byte,int is_short,int is_double);
*/
/*splitpoint c_fstecr */
/*****************************************************************************
 *                           C _ F S T E C R                                 *
 *                                                                           *
 *Object                                                                     *
 *   Writes record to file.                                                  *
 *                                                                           *
 *Revision                                                                   *
 *   Sept 2011 - Deltat, becomes a long long (deet*npas > 32bit)             *
 *   Mar  2013 - M.Valin  missing values, turbo compression, byte/short mods *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  field   field to write to the file                                   *
 *  IN  work    work field (kept for backward compatibility)                 *
 *  IN  npak    number of bits kept for the elements of the field (-npak)    *
 *  IN  iun     unit number associated to the file                           *
 *  IN  date    date time stamp                                              *
 *  IN  deet    length of a time step in seconds                             *
 *  IN  npas    time step number                                             *
 *  IN  ni      first dimension of the data field                            *
 *  IN  nj      second dimension of the data field                           *
 *  IN  nk      third dimension of the data field                            *
 *  IN  ip1     vertical level                                               *
 *  IN  ip2     forecast hour                                                *
 *  IN  ip3     user defined identifier                                      *
 *  IN  typvar  type of field (forecast, analysis, climatology)              *
 *  IN  nomvar  variable name                                                *
 *  IN  etiket  label                                                        *
 *  IN  grtyp   type of geographical projection                              *
 *  IN  ig1     first grid descriptor                                        *
 *  IN  ig2     second grid descriptor                                       *
 *  IN  ig3     third grid descriptor                                        *
 *  IN  ig4     fourth grid descriptor                                       *
 *  IN  datyp   data type of the elements                                    *
 *          0: binary, transparent                                           *
 *          1: floating point                                                *
 *          2: unsigned integer                                              *
 *          3: character (R4A in an integer)                                 *
 *          4: signed integer                                                *
 *          5: IEEE floating point                                           *
 *          6: floating point (16 bit, made for compressor)                  *
 *          7: character string                                              *
 *          8: complex IEEE                                                  *
 *        130: compressed short integer  (128+2)                             *
 *        133: compressed IEEE           (128+5)                             *
 *        134: compressed floating point (128+6)                             *
 *      +128 : second stage packer active                                    *
 *      +64  : missing value convention used                                 *
 *  IN  rewrit  rewrite flag (true=rewrite existing record, false=append)    *
 *                                                                           *
 *****************************************************************************/
int c_fstecr(word *field_in, void * work, int npak,
                        int iun, int date,
                        int deet, int npas,
                        int ni, int nj, int nk,
                        int ip1, int ip2, int ip3,
                        char *in_typvar, char *in_nomvar, char *in_etiket,
                        char *in_grtyp, int ig1, int ig2,
                        int ig3, int ig4,
                        int in_datyp_ori, int rewrit)
{
  word *field = field_in ; /* use field internally in case we have to allocate new array because of missing values */
  word *field3;
  short *s_field;
  signed char *b_field;
  int ier,l1,l2,l3,l4;
  int index, index_fnom, nbits, handle;
  long long deltat;
  unsigned int datev;
  int p1out,p2out,header_size,stream_size;
  int nw, keys_len, one=1, zero=0, njnk;
  int bitmot=32;
  int minus_nbits, i, lngw, compressed_lng;
  int datyp, nbytes;
  int niout, njout, nkout;
  int is_missing ; /*  missing value feature used flag */
  int in_datyp = in_datyp_ori & 0xFFBF ; /* suppress missing value flag (64) */
  int sizefactor ; /* number of bytes per data item */
  int IEEE_64=0  ; /* flag 64 bit IEEE (type 5 or 8) */

  file_table_entry *f;
  stdf_dir_keys *stdf_entry;
  buffer_interface_ptr buffer;
  PackFunctionPointer packfunc;
  ftnword f_datev;
  double nhours, tempfloat=99999.0;
  char string[12];

  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
                   ,'\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};


  is_missing = in_datyp_ori & 64 ; /* will be cancelled later if not supported or no missing values detected */
  if ( (in_datyp&0xF) == 8) {
    if(in_datyp_ori!=8) {
      WARNPRINT fprintf(stderr,"WARNING: compression and/or missing values not supported, data type %d reset to %d (complex)\n",in_datyp_ori,8);
    }
    is_missing = 0;   /* missing values not supported for complex type */
    in_datyp = 8; /* extra compression not supported for complex type */
  }

  l1 = strlen(in_typvar);
  l2 = strlen(in_nomvar);
  l3 = strlen(in_etiket);
  l4 = strlen(in_grtyp);

  string_copy(typvar,in_typvar,l1);
  string_copy(nomvar,in_nomvar,l2);
  string_copy(etiket,in_etiket,l3);
  string_copy(grtyp,in_grtyp,l4);

  if (in_datyp == 801)  /* 512+256+32+1 no interference with turbo pack (128) and missing value (64) flags */
  { datyp = 1; /* xdf_double=1 ; */ }
  else
  { datyp = in_datyp; }

  if ((xdf_double) || (in_datyp == 801))
    packfunc = &compact_double;
  else
    packfunc = &compact_float;


  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstecr",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstecr",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (! f->cur_info->attr.std) {
    sprintf(errmsg,"file (unit=%d) is not a RPN standard file",iun);
    return(error_msg("c_fstecr",ERR_NO_FILE,ERROR));
  }

  if (f->fstd_vintage_89) {
    sprintf(errmsg,"can not write (unit=%d) on an old (version 89) RPN standard file",iun);
    return(error_msg("c_fstecr",ERR_NO_WRITE,ERRFATAL));
  }

  if (f->cur_info->attr.read_only) {
    sprintf(errmsg,"file (unit=%d) not open with write permission",iun);
    return(error_msg("c_fstecr",ERR_NO_WRITE,ERROR));
  }

  if (npak == 0)
    nbits = FTN_Bitmot;
  else
    nbits = (npak < 0) ? -npak : Max(1,FTN_Bitmot/Max(1,npak));
  nk = Max(1,nk);
  minus_nbits = -nbits;
  njnk = nj * nk;

  if ( (in_datyp_ori == 133) && (nbits > 32) ) {
      WARNPRINT fprintf(stderr,"WARNING: extra compression not supported for IEEE when nbits > 32, data type 133 reset to 5 (IEEE)\n");
      in_datyp = 5; /* extra compression not supported */
      datyp = 5;
  }

  if ((in_datyp == 1) && ((nbits == 31) || (nbits == 32)) && !image_mode_copy) {
    datyp = 5;              /* R32 to E32 automatic conversion */
    nbits = 32;
    minus_nbits = -32;
    }

    if ( ((in_datyp & 0xF) == 5) && (nbits == 64) ) IEEE_64=1;  /* 64 bits real IEEE */
    if ( ((in_datyp & 0xF) == 8) && (nbits == 64) ) IEEE_64=1;  /* 64 bits complex IEEE */

  /* validate range of arguments */
  VALID(ni,1,NI_MAX,"ni","c_fstecr")
  VALID(nj,1,NJ_MAX,"nj","c_fstecr")
  VALID(nk,1,NK_MAX,"nk","c_fstecr")
  VALID(deet,0,DEET_MAX,"deet","c_fstecr")
  VALID(npas,0,NPAS_MAX,"npas","c_fstecr")
  VALID(nbits,1,NBITS_MAX,"nbits","c_fstecr")
  VALID(ig1,0,IG1_MAX,"ig1","c_fstecr")
  VALID(ig2,0,IG2_MAX,"ig2","c_fstecr")
  VALID(ig3,0,IG3_MAX,"ig3","c_fstecr")
  VALID(ig4,0,IG4_MAX,"ig4","c_fstecr")
  VALID(ip1,0,IP1_MAX,"ip1","c_fstecr")
  VALID(ip2,0,IP2_MAX,"ip2","c_fstecr")
  VALID(ip3,0,IP3_MAX,"ip3","c_fstecr")
  VALID(ni*nj*nk*nbits/FTN_Bitmot,0,MAX_RECORD_LENGTH,"record length > 128MB","c_fstecr")
  ;

  datev = date;
  f_datev = (ftnword) datev;
  if (( (long long) deet * npas) > 0) {
    deltat = (long long) deet * npas;
    nhours = (double) deltat;
    nhours = nhours / 3600.;
    f77name(incdatr)(&f_datev,&f_datev,&nhours);
    datev = (unsigned int) f_datev;
  }

  if ((npak == 0) || (npak == 1))
    datyp = 0;                     /* no compaction */

  /* allocate and initialize a buffer interface for xdfput */
  /* an extra 512 bytes are allocated for cluster alignment purpose (seq) */

  if (! image_mode_copy)
    for (i=0; i < nb_remap; i++)
      if (datyp == remap_table[0][i]) {
        datyp = remap_table[1][i];
/*        printf("Debug+ remapping %d to %d\n",remap_table[0][i],datyp); */
        }

  if ((nbits > 16) && (datyp != 133)) datyp &= 0x7F;    /* no extra compression if nbits > 16 */
/*  if ((datyp < 128) && (extra_compression > 0) && (nbits <= 16)) datyp += extra_compression; */
  if ((datyp == 6) && (nbits > 24)) {
    if (! dejafait_1) {
      WARNPRINT fprintf(stderr,"c_fstecr WARNING: nbits > 16, writing E32 instead of F%2d\n",nbits);
      dejafait_1 = 1;
      }
    datyp = 5;
    nbits = 32;
    minus_nbits = -32;
    }
  if ((datyp == 6) && (nbits > 16)) {
    if (! dejafait_2) {
/*      printf("msg_level=%d dejafait_2=%d !dejafait_2=%d\n",msg_level,dejafait_2,!dejafait_2); */
      WARNPRINT fprintf(stderr,"c_fstecr WARNING: nbits > 16, writing R%2d instead of F%2d\n",nbits,nbits);
      dejafait_2 = 1;
      }
    datyp = 1;
    }
  switch (datyp) {

    case 6:
      c_float_packer_params(&header_size,&stream_size,&p1out,&p2out,ni*nj*nk);
      nw = ((header_size+stream_size) * 8 + 63) / 64;
      header_size /= sizeof(INT_32);
      stream_size /= sizeof(INT_32);
      break;

    case 8:
      nw = 2 * ((ni * nj * nk *nbits + 63) / 64);
      break;

    case 129:
      nw = (ni*nj*nk*Max(nbits,16) + 128 + 32 + 63) / 64;         /* 120 bits (floatpack header)+8, 32 bits (extra header) */
      break;

    case 130:
      nw = (ni*nj*nk*Max(nbits,16) +32 + 63) / 64;                /* 32 bits (extra header) */
      break;

    case 134:
      c_float_packer_params(&header_size,&stream_size,&p1out,&p2out,ni*nj*nk);
      nw = ((header_size+stream_size) * 8 + 32 + 63) / 64;
      stream_size /= sizeof(INT_32);
      header_size /= sizeof(INT_32);
      break;

    default:
      nw = (ni*nj*nk*nbits + 120 + 63) / 64;
      break;
    } /* end switch */

  nw = W64TOWD(nw);

  keys_len = W64TOWD((f->primary_len + f->info_len));
  buffer = (buffer_interface_ptr) alloca((10+keys_len+nw+128)*sizeof(word));
  if (buffer)
    memset(buffer,0,(10+keys_len+nw+128)*sizeof(word));
  else {
    sprintf(errmsg,"memory is full, was trying to allocate %ld bytes",
            (10+keys_len+nw+128)*sizeof(word));
    return(error_msg("c_fstecr",ERR_MEM_FULL,ERRFATAL));
  }
  buffer->nwords = 10 + keys_len + nw;
  buffer->nbits = (keys_len + nw) * bitmot;
  buffer->record_index = RECADDR;
  buffer->data_index = buffer->record_index +
                       W64TOWD((f->primary_len + f->info_len));
  buffer->iun = iun;
  buffer->aux_index = buffer->record_index + W64TOWD(f->primary_len);
  buffer->data[buffer->aux_index] = 0;
  buffer->data[buffer->aux_index+1] = 0;

 /* set stdf_entry to address of buffer->data for building keys */
  stdf_entry = (stdf_dir_keys *) &(buffer->data);
  stdf_entry->deleted = 0;
  stdf_entry->select = 1;
  stdf_entry->lng = -1;
  stdf_entry->addr = -1;
  stdf_entry->deet = deet;
  stdf_entry->nbits = nbits;
  stdf_entry->ni = ni;
  stdf_entry->gtyp = grtyp[0];
  stdf_entry->nj = nj;
  stdf_entry->datyp = datyp | is_missing;  /* propagate missing values flag        */
         /* this value may be changed later in the code to eliminate improper flags */
  stdf_entry->nk = nk;
  stdf_entry->ubc = 0;
  stdf_entry->npas = npas;
  stdf_entry->pad7 = 0;
  stdf_entry->ig4 = ig4;
  stdf_entry->ig2a = ig2 >> 16;
  stdf_entry->ig1 = ig1;
  stdf_entry->ig2b = ig2 >> 8;
  stdf_entry->ig3 = ig3;
  stdf_entry->ig2c = ig2 & 0xff;
  stdf_entry->etik15 =
    (ascii6(etiket[0]) << 24) |
    (ascii6(etiket[1]) << 18) |
    (ascii6(etiket[2]) << 12) |
    (ascii6(etiket[3]) <<  6) |
    (ascii6(etiket[4]));
  stdf_entry->pad1 = 0;
  stdf_entry->etik6a =
    (ascii6(etiket[5]) << 24) |
    (ascii6(etiket[6]) << 18) |
    (ascii6(etiket[7]) << 12) |
    (ascii6(etiket[8])<<  6) |
    (ascii6(etiket[9]));
  stdf_entry->pad2 = 0;
  stdf_entry->etikbc =
    (ascii6(etiket[10]) <<  6) |
    (ascii6(etiket[11]));
  stdf_entry->typvar =
    (ascii6(typvar[0]) <<  6) |
    (ascii6(typvar[1]));
  stdf_entry->pad3 = 0;
  stdf_entry->nomvar =
    (ascii6(nomvar[0]) << 18) |
    (ascii6(nomvar[1]) << 12) |
    (ascii6(nomvar[2]) <<  6) |
    (ascii6(nomvar[3]));
  stdf_entry->pad4 = 0;
  stdf_entry->ip1 = ip1;
  stdf_entry->levtyp = 0;
  stdf_entry->ip2 = ip2;
  stdf_entry->pad5 = 0;
  stdf_entry->ip3 = ip3;
  stdf_entry->pad6 = 0;
  stdf_entry->date_stamp = 8 * (datev/10) + (datev % 10);

  handle = 0;
  if ((rewrit) && (!f->xdf_seq)) {       /* find handle for rewrite operation */
    handle = c_fstinf(iun,&niout,&njout,&nkout,-1,etiket,ip1,ip2,ip3,typvar,nomvar);
    if (handle < 0)
      handle = 0;      /* append mode for xdfput */
  }

  if (image_mode_copy) {                /* no pack/unpack, used by editfst */
    if (datyp > 128) {
      lngw = field[0];                  /* first element is length */
/*      fprintf(stderr,"Debug+ datyp=%d ecriture mode image lngw=%d\n",datyp,lngw); */
      buffer->nbits = (keys_len + lngw) * bitmot;
      for (i=0; i < lngw+1; i++)
        buffer->data[keys_len+i] = field[i];
     }
    else {
      if (datyp == 6) {
        c_float_packer_params(&header_size,&stream_size,&p1out,&p2out,ni*nj*nk);
        lngw = (header_size+stream_size) * 8 ;
        }
      else
        lngw = ni*nj*nk*nbits;
      if (datyp == 1) lngw += 120;
      if (datyp == 3) lngw = ni*nj * 8;
      lngw = (lngw + bitmot - 1) / bitmot;
      for (i=0; i < lngw; i++)
        buffer->data[keys_len+i] = field[i];
    }
  }
  else {   /* not image mode copy */
    /* time to fudge field if missing value feature is used */
    sizefactor=4 ;
    if(xdf_byte)  sizefactor=1 ;
    if(xdf_short) sizefactor=2 ;
    if(xdf_double | IEEE_64)sizefactor=8 ;
    if(is_missing){    /* put appropriate values into field after allocating it */
      field= (word *) alloca(ni*nj*nk*sizefactor); /* allocate self deallocating scratch field */
      if( 0 == EncodeMissingValue(field,field_in,ni*nj*nk,in_datyp,nbits,xdf_byte,xdf_short,xdf_double) ) {
        field=field_in ;
        INFOPRINT fprintf(stderr,"NO missing value, data type %d reset to %d\n",stdf_entry->datyp,datyp);
        stdf_entry->datyp = datyp;  /* cancel missing data flag in data type */
        is_missing = 0;
      }
    }
    switch (datyp) {

    case 0: case 128:              /* transparent mode */
      if(datyp==128) {
        WARNPRINT fprintf(stderr,"WARNING: extra compression not available, data type %d reset to %d\n",stdf_entry->datyp,0);
        datyp = 0;
        stdf_entry->datyp = 0;
      }
      lngw = ((ni*nj*nk*nbits) + bitmot-1) / bitmot;
      for (i=0; i < lngw; i++)
        buffer->data[keys_len+i] = field[i];
      break;

    case 1: case 129:             /* floating point */
      if ((datyp > 128) && (nbits <= 16)) {      /* use an additional compression scheme */
        packfunc(field,&(buffer->data[keys_len+1]),&(buffer->data[keys_len+5]),
                ni*nj*nk,nbits+64*Max(16,nbits),0,xdf_stride,1,0,&tempfloat);     /* nbits>64 flags a different packing */
        compressed_lng = armn_compress(&(buffer->data[keys_len+5]),ni,nj,nk,nbits,1);
        if (compressed_lng < 0)
           {
           stdf_entry->datyp = 1;
           packfunc(field,&(buffer->data[keys_len]),&(buffer->data[keys_len+3]),
                  ni*nj*nk,nbits,24,xdf_stride,1,0,&tempfloat);
           }
        else
           {
           nbytes = 16+ compressed_lng;
    /*        fprintf(stderr,"Debug+ apres armn_compress nbytes=%d\n",nbytes); */
           nw = (nbytes * 8 + 63) / 64;
           nw = W64TOWD(nw);
           buffer->data[keys_len] = nw;
    /*        fprintf(stderr,"Debug+ pack buffer->data[keys_len]=%d\n",buffer->data[keys_len]); */
           buffer->nbits = (keys_len + nw) * bitmot;
           }
         }
      else
        {
          packfunc(field,&(buffer->data[keys_len]),&(buffer->data[keys_len+3]),
                   ni*nj*nk,nbits,24,xdf_stride,1,0,&tempfloat);
        }
      break;

    case 2: case 130:              /* integer, short integer or byte stream */
      {
      int offset;
//      short *p16;
      offset = (datyp > 128) ? 1 :0;
      if (datyp > 128) {
        if (xdf_short) {
          stdf_entry->nbits = Min(16,nbits);
          nbits = stdf_entry->nbits;
          memcpy(&(buffer->data[keys_len+offset]),field,ni*nj*nk*2);
        }
        else if (xdf_byte) {
          stdf_entry->nbits = Min(8,nbits);
          nbits = stdf_entry->nbits;
          memcpy_8_16(&(buffer->data[keys_len+offset]),field,ni*nj*nk);
        }
        else {
          memcpy_32_16(&(buffer->data[keys_len+offset]),field,nbits,ni*nj*nk);
/*
          ier = compact_integer(field,(void *) NULL,&(buffer->data[keys_len+offset]),
                                ni*nj*nk,nbits,0,xdf_stride,1);
*/
        }
        c_armn_compress_setswap(0);
        compressed_lng = armn_compress(&(buffer->data[keys_len+offset]),ni,nj,nk,nbits,1);
        c_armn_compress_setswap(1);
        if (compressed_lng < 0) {
            stdf_entry->datyp = 2;
            ier = compact_integer(field,(void *) NULL,&(buffer->data[keys_len+offset]),
                                  ni*nj*nk,nbits,0,xdf_stride,1);
          }
        else {
          nbytes = 4 + compressed_lng;
       /*   fprintf(stderr,"Debug+ fstecr armn_compress compressed_lng=%d\n",compressed_lng); */
          nw = (nbytes * 8 + 63) / 64;
          nw = W64TOWD(nw);
          buffer->data[keys_len] = nw;
          buffer->nbits = (keys_len + nw) * bitmot;
        }
      }
      else {
        if (xdf_short) {
          stdf_entry->nbits = Min(16,nbits);
          nbits = stdf_entry->nbits;
          ier = compact_short(field,(void *) NULL,&(buffer->data[keys_len+offset]),
                              ni*nj*nk,nbits,0,xdf_stride,5);
        }
        else if (xdf_byte) {
          ier = compact_char(field,(void *) NULL,&(buffer->data[keys_len]),
                             ni*nj*nk,Min(8,nbits),0,xdf_stride,9);
          stdf_entry->nbits = Min(8,nbits);
          nbits = stdf_entry->nbits;
        }
        else {
          ier = compact_integer(field,(void *) NULL,&(buffer->data[keys_len+offset]),
                                ni*nj*nk,nbits,0,xdf_stride,1);
        }
      }
      break;
      }

    case 3: case 131:              /* character */
      {
        int nc = (ni*nj+3)/4;
        if(datyp==131) {
          WARNPRINT fprintf(stderr,"WARNING: extra compression not available, data type %d reset to %d\n",stdf_entry->datyp,3);
          datyp = 3;
          stdf_entry->datyp = 3;
        }
        ier = compact_integer(field,(void *) NULL,&(buffer->data[keys_len]),nc,
                              32,0,xdf_stride,1);
        stdf_entry->nbits = 8;
        break;
      }

    case 4: case 132:              /* signed integer */
      if(datyp==132) {
        WARNPRINT fprintf(stderr,"WARNING: extra compression not supported, data type %d reset to %d\n",stdf_entry->datyp,is_missing | 4);
        datyp = 4;
      }
      stdf_entry->datyp = is_missing | 4;  /* turbo compression not supported for this type, revert to normal mode */
#ifdef use_old_signed_pack_unpack_code
// fprintf(stderr,"OLD PACK CODE======================================\n");
      field3 = field;
      if(xdf_short || xdf_byte){
        field3=(word *)alloca(ni*nj*nk*sizeof(word));
        s_field = (short *)field; b_field = (signed char *)field;
        if(xdf_short) for (i=0;i<ni*nj*nk;i++) { field3[i]=s_field[i]; } ;
        if(xdf_byte)  for (i=0;i<ni*nj*nk;i++) { field3[i]=b_field[i]; } ;
      }
      ier = compact_integer(field3,(void *) NULL,&(buffer->data[keys_len]),ni*nj*nk,
                            nbits,0,xdf_stride,3);
#else
// fprintf(stderr,"NEW PACK CODE======================================\n");
      if(xdf_short){
        ier = compact_short(field,(void *) NULL,&(buffer->data[keys_len]),ni*nj*nk,
                            nbits,0,xdf_stride,7);
      } else if(xdf_byte){
        ier = compact_char(field,(void *) NULL,&(buffer->data[keys_len]),ni*nj*nk,
                            nbits,0,xdf_stride,11);
      }else{
        ier = compact_integer(field,(void *) NULL,&(buffer->data[keys_len]),ni*nj*nk,
                            nbits,0,xdf_stride,3);
      }
#endif
      break;

    case 5: case 8: case 133:  case 136:            /* IEEE and IEEE complex representation */
      {
        ftnword f_ni = (ftnword) ni;
        ftnword f_njnk = (ftnword) njnk;
        ftnword f_zero = (ftnword) zero;
        ftnword f_one = (ftnword) one;
        ftnword f_minus_nbits = (ftnword) minus_nbits;
        if(datyp==136) {
          WARNPRINT fprintf(stderr,"WARNING: extra compression not available, data type %d reset to %d\n",stdf_entry->datyp,8);
          datyp = 8;
          stdf_entry->datyp = 8;
        }
        if (datyp == 133) {
          /* use an additional compression scheme */
          compressed_lng = c_armn_compress32(&(buffer->data[keys_len+1]), field, ni,nj,nk,nbits);
          if (compressed_lng < 0)
           {
             stdf_entry->datyp = 5;
             f77name(ieeepak)(field,&(buffer->data[keys_len]),&f_ni,&f_njnk,&f_minus_nbits,
                              &f_zero,&f_one);
           }
          else {
            nbytes = 16 + compressed_lng;
            nw = (nbytes * 8 + 63) / 64;
            nw = W64TOWD(nw);
            buffer->data[keys_len] = nw;
            buffer->nbits = (keys_len + nw) * bitmot;
            }
          }
        else {
          if (datyp == 8) f_ni = f_ni * 2;
          f77name(ieeepak)(field,&(buffer->data[keys_len]),&f_ni,&f_njnk,&f_minus_nbits,
                           &f_zero,&f_one);
          }
        break;
      }

    case 6: case 134:             /* floating point, new packers */
      {
        if ((datyp > 128) && (nbits <= 16)) {      /* use an additional compression scheme */
          c_float_packer(field,nbits,&(buffer->data[keys_len+1]),&(buffer->data[keys_len+1+header_size]),ni*nj*nk);
          compressed_lng = armn_compress(&(buffer->data[keys_len+1+header_size]),ni,nj,nk,nbits,1);
          if (compressed_lng < 0)
           {
             stdf_entry->datyp = 6;
             c_float_packer(field,nbits,&(buffer->data[keys_len]),&(buffer->data[keys_len+header_size]),ni*nj*nk);
           }
          else {
            nbytes = 16 + (header_size*4) + compressed_lng;
        /*   fprintf(stderr,"Debug+ apres armn_compress nbytes=%d\n",nbytes); */
            nw = (nbytes * 8 + 63) / 64;
            nw = W64TOWD(nw);
            buffer->data[keys_len] = nw;
        /*   fprintf(stderr,"Debug+ pack buffer->data[keys_len]=%d\n",buffer->data[keys_len]); */
            buffer->nbits = (keys_len + nw) * bitmot;
          }

        }
        else {
          c_float_packer(field,nbits,&(buffer->data[keys_len]),&(buffer->data[keys_len+header_size]),ni*nj*nk);
/*          fprintf(stderr,"Debug+ fstecr apres float_packer buffer->data=%8X\n",buffer->data[keys_len]); */
        }
        break;
      }

    case 7: case 135:              /* character string */
      if(datyp==135) {
        WARNPRINT fprintf(stderr,"WARNING: extra compression not available, data type %d reset to %d\n",stdf_entry->datyp,7);
        datyp = 7;
        stdf_entry->datyp = 7;
      }
      ier = compact_char(field,(void *) NULL,&(buffer->data[keys_len]),
                            ni*nj*nk,8,0,xdf_stride,9);
      break;

    default: sprintf(errmsg,"(unit=%d) invalid datyp=%d",iun,datyp);
      return(error_msg("c_fstecr",ERR_BAD_DATYP,ERROR));

    } /* end switch */
  }   /* end if image mode copy */


  /* write record to file and add entry to directory */
  ier = c_xdfput(iun,handle,buffer);
  if (msg_level <= INFORM) {
    sprintf(string,"Write(%d)",iun);
    print_std_parms(stdf_entry,string,prnt_options,0);
  }
/*  free(buffer);        replaced by alloca */
  xdf_double = 0;
  xdf_short = 0;
  xdf_byte = 0;
  return(ier);
}

/*splitpoint c_fst_edit_dir */
/*****************************************************************************
 *                     C _ F S T _ E D I T _ D I R                           *
 *                                                                           *
 *Object                                                                     *
 *   Edits the directory content of a RPN standard file.                     *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle     handle to the directory entry to edit                     *
 *                                                                           *
 *****************************************************************************/

int c_fst_edit_dir_plus(int handle,
                   unsigned int date, int deet, int npas,
                   int ni, int nj, int nk,
                   int ip1, int ip2, int ip3,
                   char *in_typvar, char *in_nomvar, char *in_etiket,
                   char *in_grtyp, int ig1, int ig2,
                   int ig3, int ig4, int datyp)
{
  int index,width,pageno,recno;
  int l1,l2,l3,l4;
  file_table_entry *f;
  word *entry;
  stdf_dir_keys *stdf_entry;
  stdf_special_parms cracked;
  char string[20];
  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
                   ,'\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};


  index = INDEX_FROM_HANDLE(handle);

  if ((index < 0) || (index >= MAX_XDF_FILES)) {
    sprintf(errmsg,"invalid handle=%d",handle);
    return(error_msg("c_fst_edit_dir",ERR_BAD_HNDL,ERROR));
  }

  f = file_table[index];

  if (! f->cur_info->attr.std) {
    sprintf(errmsg,"file (unit=%d) is not a RPN standard file",f->iun);
    return(error_msg("c_fst_edit_dir",ERR_NO_FILE,ERROR));
  }

  if (f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not a RPN standard file",f->iun);
    return(error_msg("c_fst_edit_dir",ERR_NO_FILE,ERROR));
  }

  l1 = strlen(in_typvar);
  l2 = strlen(in_nomvar);
  l3 = strlen(in_etiket);
  l4 = strlen(in_grtyp);

  string_copy(typvar,in_typvar,l1);
  string_copy(nomvar,in_nomvar,l2);
  string_copy(etiket,in_etiket,l3);
  string_copy(grtyp,in_grtyp,l4);

  pageno = PAGENO_FROM_HANDLE(handle);
  if (pageno > f->npages) {   /* page is not in current file */
    sprintf(errmsg,"invalid handle, invalid page number\n");
    return(error_msg("c_fst_edit_dir",ERR_BAD_PAGENO,ERROR));
  }

  recno = RECORD_FROM_HANDLE(handle);
/*  printf("Debug+ c_fst_edit_dir pageno=%d recno=%d\n",pageno,recno); */
  width = W64TOWD(f->primary_len);
  entry = (f->dir_page[pageno])->dir.entry + recno * width;
  stdf_entry = (stdf_dir_keys *) entry;

/*  stdf_entry->ni = ni; */
/*  stdf_entry->nj = nj; */
/*  stdf_entry->nk = nk; */
  if(grtyp[0] != ' ') stdf_entry->gtyp = grtyp[0];
/*  stdf_entry->datyp = datyp; */

  if (deet != -1) stdf_entry->deet = deet;
  if (npas != -1) stdf_entry->npas = npas;
  if (ig1 != -1) stdf_entry->ig1 = ig1;
  if (ig2 != -1) {
    stdf_entry->ig2a = ig2 >> 16;
    stdf_entry->ig2b = ig2 >> 8;
    stdf_entry->ig2c = ig2 & 0xff;
  }
  if (ig3 != -1) stdf_entry->ig3 = ig3;
  if (ig4 != -1) stdf_entry->ig4 = ig4;
  if (strcmp(etiket,"            ") !=0 ) {
    stdf_entry->etik15 =
      (ascii6(etiket[0]) << 24) |
      (ascii6(etiket[1]) << 18) |
      (ascii6(etiket[2]) << 12) |
      (ascii6(etiket[3]) <<  6) |
      (ascii6(etiket[4]));
    stdf_entry->etik6a =
      (ascii6(etiket[5]) << 24) |
      (ascii6(etiket[6]) << 18) |
      (ascii6(etiket[7]) << 12) |
      (ascii6(etiket[8])<<  6) |
      (ascii6(etiket[9]));
    stdf_entry->etikbc =
      (ascii6(etiket[10]) <<  6) |
      (ascii6(etiket[11]));
  }
  if (strcmp(typvar,"  ") != 0) {
    stdf_entry->typvar =
      (ascii6(typvar[0]) <<  6) |
      (ascii6(typvar[1]));
  }
  if (strcmp(nomvar,"    ") != 0) {
    stdf_entry->nomvar =
      (ascii6(nomvar[0]) << 18) |
      (ascii6(nomvar[1]) << 12) |
      (ascii6(nomvar[2]) <<  6) |
      (ascii6(nomvar[3]));
  }
/*  stdf_entry->levtyp = 0; */
  if (ip1 != -1) stdf_entry->ip1 = ip1;
  if (ip2 != -1) stdf_entry->ip2 = ip2;
  if (ip3 != -1) stdf_entry->ip3 = ip3;
  if (date != -1) stdf_entry->date_stamp = 8 * (date/10) + (date % 10);

  crack_std_parms(stdf_entry,&cracked);
  sprintf(string,"%5d-",recno);
  print_std_parms(stdf_entry,string,"NINJNK+DATEO+IP1+IG1234",1);
  f->dir_page[pageno]->modified = 1;
  f->modified = 1;
  return(0);
}
int c_fst_edit_dir(int handle,
                   unsigned int date, int deet, int npas,
                   int ni, int nj, int nk,
                   int ip1, int ip2, int ip3,
                   char *in_typvar, char *in_nomvar, char *in_etiket,
                   char *in_grtyp, int ig1, int ig2,
                   int ig3, int ig4, int datyp)
{
  return c_fst_edit_dir_plus(handle,date,deet,npas,-1,-1,-1,ip1,ip2,ip3,in_typvar,in_nomvar,in_etiket," ",ig1,ig2,ig3,ig4,-1);
}

/*splitpoint c_fsteff */
/*****************************************************************************
 *                         C _ F S T E F F                                   *
 *                                                                           *
 *Object                                                                     *
 *   Deletes the record associated to handle.                                *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle  handle to the record to delete                               *
 *                                                                           *
 *****************************************************************************/

int c_fsteff(int handle)
{
  int ier, index;
  file_table_entry *f;

  index = INDEX_FROM_HANDLE(handle);
  if ((index < 0) || (index >= MAX_XDF_FILES)) {
    sprintf(errmsg,"invalid handle=%d",handle);
    return(error_msg("c_fsteff",ERR_BAD_HNDL,ERROR));
  }

  f = file_table[index];
  if (f == NULL) {
    sprintf(errmsg,"invalid handle=%d",handle);
    return(error_msg("c_fsteff",ERR_BAD_HNDL,ERROR));
  }

  if (! f->cur_info->attr.std) {
    sprintf(errmsg,"file (unit=%d) is not a RPN standard file",f->iun);
    return(error_msg("c_fsteff",ERR_NO_FILE,ERROR));
  }

  if (f->fstd_vintage_89) {
    sprintf(errmsg,"can not write (unit=%d) on an old (version 89) RPN standard file",f->iun);
    return(error_msg("c_fsteff",ERR_NO_WRITE,ERRFATAL));
  }

  ier = c_xdfdel(handle);
  return(ier);
}

/*splitpoint c_fsteof */
/*****************************************************************************
 *                            C _ F S T E O F                                *
 *                                                                           *
 *Object                                                                     *
 *   Return the level of end of file for the sequential file.                *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/
int c_fsteof(int iun)
{
  int index_fnom, index, eof;
  file_table_entry *f;
  xdf_record_header *header;
  seq_dir_keys *seq_entry;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fsteof",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fsteof",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (!f->xdf_seq)
    return(0);

  eof = 0;
  if (!f->fstd_vintage_89) {
    header = (xdf_record_header *) f->head_keys;
    if ((header->idtyp >= 112) && (header->idtyp <= 127))
      eof = header->idtyp - 112;
  }
  else {
    seq_entry = (seq_dir_keys *) f->head_keys;
    if (seq_entry->eof > 0)
      eof = (seq_entry->eof == 31) ? 15 : seq_entry->eof;
  }

  return(eof);
}

/*splitpoint c_fstfrm */
/*****************************************************************************
 *                         C _ F S T F R M                                   *
 *                                                                           *
 *Object                                                                     *
 *   Closes a RPN standard file.                                             *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/

int c_fstfrm(int iun)
{
  int ier;

  ier = c_xdfcls(iun);
  return(ier);
}

/*splitpoint c_fstinf */
/*****************************************************************************
 *                        C _ F S T I N F                                    *
 *                                                                           *
 *Object                                                                     *
 *   Locate the next record that matches the research keys.                  *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *  IN  datev   valid date                                                   *
 *  IN  etiket  label                                                        *
 *  IN  ip1     vertical level                                               *
 *  IN  ip2     forecast hour                                                *
 *  IN  ip3     user defined identifier                                      *
 *  IN  typvar  type of field                                                *
 *  IN  nomvar  variable name                                                *
 *                                                                           *
 *****************************************************************************/

int c_fstinf(int iun, int *ni, int *nj, int *nk, int datev,char *in_etiket,
                 int ip1, int ip2, int ip3, char *in_typvar, char *in_nomvar)
{
  int err, handle;

  handle = -2;  /* means handle will be discarded */
  err = c_fstinfx(handle,iun,ni,nj,nk,datev,in_etiket,ip1,ip2,ip3,
                      in_typvar,in_nomvar);
  return(err);
}

/*splitpoint c_fstinfx */
/*****************************************************************************
 *                       C _ F S T I N F X                                   *
 *                                                                           *
 *Object                                                                     *
 *   Locate the next record that matches the research keys.                  *
 *   The search begins at the position given by handle.                      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle  handle from which the search begins                          *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *  IN  datev   valid date stamp                                             *
 *              last 3 bits ignored, used to be run number                   *
 *              search resolution drops to 40s                               *
 *  IN  etiket  label                                                        *
 *  IN  ip1     vertical level                                               *
 *  IN  ip2     forecast hour                                                *
 *  IN  ip3     user defined identifier                                      *
 *  IN  typvar  type of field                                                *
 *  IN  nomvar  variable name                                                *
 *                                                                           *
 *****************************************************************************/

int c_fstinfx(int handle, int iun, int *ni, int *nj, int *nk,
                 int datev,char *in_etiket,
                 int ip1, int ip2, int ip3, char *in_typvar, char *in_nomvar)
{
  stdf_dir_keys *stdf_entry, *search_mask;
  word * pkeys, *pmask;
  file_table_entry *f;
  int i, index_fnom, index, index_h;
  int addr, lng, idtyp, ier, l1, l2, l3;
  unsigned int u_datev=datev;

  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '
                   ,'\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};

  l1 = strlen(in_etiket);
  l2 = strlen(in_typvar);
  l3 = strlen(in_nomvar);

  string_copy(etiket,in_etiket,l1);
  string_copy(typvar,in_typvar,l2);
  string_copy(nomvar,in_nomvar,l3);
  if (msg_level < INFORM)
    fprintf(stdout,"Debug fstinf iun %d recherche: datev=%d etiket=[%s] ip1=%d ip2=%d ip3=%d typvar=[%s] nomvar=[%s]\n",iun,datev,etiket,ip1,ip2,ip3,typvar,nomvar);

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstinf",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstinf",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  stdf_entry = (stdf_dir_keys *) calloc(1,sizeof(stdf_dir_keys));
  search_mask = (stdf_dir_keys *) calloc(1,sizeof(stdf_dir_keys));

  pkeys = (word *) stdf_entry;
  pmask = (word *) search_mask;

  for (i = 0; i < (sizeof(stdf_dir_keys) / sizeof(word)); i++)
    pmask[i] = -1;

  search_mask->pad1 = 0;
  search_mask->pad2 = 0;
  search_mask->pad3 = 0;
  search_mask->pad4 = 0;
  search_mask->pad5 = 0;
  search_mask->pad6 = 0;
  search_mask->pad7 = 0;
  search_mask->deleted = 0;
  search_mask->select = 0;
  search_mask->lng = 0;
  search_mask->addr = 0;
  search_mask->deet = 0;
  search_mask->nbits = 0;
  search_mask->ni = 0;
  search_mask->gtyp = 0;
  search_mask->nj = 0;
  search_mask->datyp = 0;
  search_mask->nk = 0;
  search_mask->ubc = 0;
  search_mask->npas = 0;
  search_mask->ig4 = 0;
  search_mask->ig2a = 0;
  search_mask->ig1 = 0;
  search_mask->ig2b = 0;
  search_mask->ig3 = 0;
  search_mask->ig2c = 0;
  search_mask->levtyp = 0;

  stdf_entry->date_stamp = 8 * (u_datev/10) + (u_datev % 10);
  search_mask->date_stamp &= ~(0x7);
  if (datev == -1) search_mask->date_stamp = 0;

  stdf_entry->ip1 = ip1;
  if ((ip1 == -1) || (ip1s_flag)) search_mask->ip1 = 0;

  stdf_entry->ip2 = ip2;
  if ((ip2 == -1) || (ip2s_flag)) search_mask->ip2 = 0;

  stdf_entry->ip3 = ip3;
  if ((ip3 == -1) || (ip3s_flag)) search_mask->ip3 = 0;

  stdf_entry->nomvar = (ascii6(nomvar[0]) << 18) |
                       (ascii6(nomvar[1]) << 12) |
                       (ascii6(nomvar[2]) <<  6) |
                       (ascii6(nomvar[3]));
  if (stdf_entry->nomvar == 0) search_mask->nomvar = 0;

  stdf_entry->typvar = (ascii6(typvar[0]) << 6) |
                       (ascii6(typvar[1]));
  if (stdf_entry->typvar == 0) search_mask->typvar = 0;

  stdf_entry->etik15 = (ascii6(etiket[0]) << 24) |
                       (ascii6(etiket[1]) << 18) |
                       (ascii6(etiket[2]) << 12) |
                       (ascii6(etiket[3]) <<  6) |
                       (ascii6(etiket[4]));

  stdf_entry->etik6a = (ascii6(etiket[5]) << 24) |
                       (ascii6(etiket[6]) << 18) |
                       (ascii6(etiket[7]) << 12) |
                       (ascii6(etiket[8]) <<  6) |
                       (ascii6(etiket[9]));

  stdf_entry->etikbc = (ascii6(etiket[10]) <<  6) |
                       (ascii6(etiket[11]));

  if ((stdf_entry->etik15 == 0) && (stdf_entry->etik6a == 0)) {
    search_mask->etik15 = 0;
    search_mask->etik6a = 0;
    search_mask->etikbc = 0;
  }
  pkeys += W64TOWD(1);
  pmask += W64TOWD(1);
  if (handle == -2)  /* means handle not specified */
    if (f->xdf_seq)
      handle = c_xdfloc2(iun,-1,pkeys,16,pmask);
    else
      handle = c_xdfloc2(iun,0,pkeys,16,pmask);
  else {
    if (handle > 0) {
      index_h = INDEX_FROM_HANDLE(handle);
      if (index_h != index) {
        sprintf(errmsg,"invalid handle=%d, or iun=%d",handle,iun);
        free(stdf_entry);
        free(search_mask);
        return(error_msg("c_fstinfx",ERR_BAD_HNDL,ERROR));
      }
    }
    handle = c_xdfloc2(iun,handle,pkeys,16,pmask);
  }

  if (handle < 0) {
    if (msg_level == TRIVIAL)
      fprintf(stdout,"c_fstinf: (unit=%d) record not found, errcode=%d\n",
              iun,handle);
    if (ip1s_flag || ip2s_flag || ip3s_flag) ier = init_ip_vals();
    free(stdf_entry);
    free(search_mask);
    return(handle);
  }
  ier = c_xdfprm(handle,&addr,&lng,&idtyp,pkeys,16);

  if (ip1s_flag || ip2s_flag || ip3s_flag) {
    int nomatch = 1;
    while ((handle >=  0) && (nomatch)) {
      nomatch = 0;
      if ((ip1s_flag) && (ip1 >= 0))
        if (ip_is_equal(ip1,stdf_entry->ip1,1) == 0)
          nomatch = 1;
        else if ((ip2s_flag) && (ip2 >= 0))
           if (ip_is_equal(ip2,stdf_entry->ip2,2) == 0)
             nomatch = 1;
           else if ((ip3s_flag) && (ip3 >= 0))
             if (ip_is_equal(ip3,stdf_entry->ip3,3) == 0)
              nomatch = 1;
      if (nomatch) {
        handle = c_xdfloc2(iun,-1,pkeys,16,pmask);
        if (handle >= 0) ier = c_xdfprm(handle,&addr,&lng,&idtyp,pkeys,16);
      }
    } /* end while */
  }
  if (ip1s_flag || ip2s_flag || ip3s_flag) {   /* arranger les masques de recherches pour un fstsui */
    if (ip1s_flag) search_mask->ip1 = 0xFFFFFFF;
    if (ip2s_flag) search_mask->ip2 = 0xFFFFFFF;
    if (ip3s_flag) search_mask->ip3 = 0xFFFFFFF;
    f->build_primary(f->target,pkeys,f->cur_mask,pmask,index,1);
    ier = init_ip_vals();
  }
  *ni = stdf_entry->ni;
  *nj = stdf_entry->nj;
  *nk = stdf_entry->nk;
  free(stdf_entry);
  free(search_mask);
  return(handle);
}

/*splitpoint c_fstinl */
/*****************************************************************************
 *                        C _ F S T I N L                                    *
 *                                                                           *
 *Object                                                                     *
 *   Locates all the records that matches the research keys.                 *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *  IN  datev   valid date                                                   *
 *  IN  etiket  label                                                        *
 *  IN  ip1     vertical level                                               *
 *  IN  ip2     forecast hour                                                *
 *  IN  ip3     user defined identifier                                      *
 *  IN  typvar  type of field                                                *
 *  IN  nomvar  variable name                                                *
 *  OUT liste   list of handle to the records                                *
 *  OUT infon   number of elements for the list (number of records found)    *
 *  OUT nmax    dimension of list as given by caller                         *
 *                                                                           *
 *****************************************************************************/

int c_fstinl(int iun, int *ni, int *nj, int *nk, int datev, char *etiket,
                 int ip1, int ip2, int ip3, char *typvar, char *nomvar,
                 word *liste, int *infon, int nmax)
{
  int handle, nfound=0, nimax, njmax, nkmax, nijkmax;

  if (msg_level < INFORM)
    fprintf(stdout,"Debug fstinl iun %d recherche: datev=%d etiket=[%s] ip1=%d ip2=%d ip3=%d typvar=[%s] nomvar=[%s]\n",iun,datev,etiket,ip1,ip2,ip3,typvar,nomvar);

  handle = c_fstinf(iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar);
  nijkmax = (*ni) * (*nj) * (*nk);
  nimax = *ni;
  njmax = *nj;
  nkmax = *nk;

  while ((handle >= 0) && (nfound < nmax)) {
    liste[nfound] = handle;
    nfound++;
    if (nfound >= nmax) break;
    handle = c_fstsui(iun,ni,nj,nk);
    if ( ((*ni) * (*nj) * (*nk)) > nijkmax ) {
      nimax = *ni;
      njmax = *nj;
      nkmax = *nk;
      nijkmax = (*ni) * (*nj) * (*nk);
    }
  }
  *ni = nimax;
  *nj = njmax;
  *nk = nkmax;
  *infon = nfound;
  if (msg_level < INFORM)
    fprintf(stdout,"Debug fstinl nombre trouve=%d nmax=%d\n",nfound,nmax);

  while ( (handle = c_fstsui(iun,ni,nj,nk)) >= 0) nfound++;
  if (nfound > nmax) {
     sprintf(errmsg,"number of records found (%d) > nmax specified (%d)",nfound,nmax);
     return(error_msg("FSTINL",-nfound,ERROR));
  }
  else
    return(0);
}

/*splitpoint c_fstlic */
/*****************************************************************************
 *                        C _ F S T L I C                                    *
 *                                                                           *
 *Object                                                                     *
 *   Search for a record that matches the research keys and check that the   *
 *   remaining parameters match the record descriptors                       *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  OUT field    data field to be read                                       *
 *  IN  iun      unit number associated to the file                          *
 *  IN  niin     dimension 1 of the data field                               *
 *  IN  njin     dimension 2 of the data field                               *
 *  IN  nkin     dimension 3 of the data field                               *
 *  IN  datein   valid date                                                  *
 *  IN  etiketin label                                                       *
 *  IN  ip1in    vertical level                                              *
 *  IN  ip2in    forecast hour                                               *
 *  IN  ip3in    user defined identifier                                     *
 *  IN  typvarin type of field                                               *
 *  IN  nomvarin variable name                                               *
 *  IN  ig1      first grid descriptor                                       *
 *  IN  ig2      second grid descriptor                                      *
 *  IN  ig3      third grid descriptor                                       *
 *  IN  ig4      fourth grid descriptor                                      *
 *  IN  grtypin  type of geographical projection                             *
 *                                                                           *
 *****************************************************************************/

int c_fstlic(word *field, int iun, int niin, int njin, int nkin,
                 int datein, char *etiketin, int ip1in, int ip2in, int ip3in,
                 char *typvarin, char *nomvarin,
                 int ig1in, int ig2in, int ig3in, int ig4in, char *grtypin)
{
  int handle, ier;
  int ni, nj, nk, date, ip1, ip2, ip3, ig1, ig2, ig3, ig4;
  int nbits, swa, ubc, lng, dltf, xtra1, xtra2, xtra3, deet, npas, datyp;
  char etiket[13]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};
  char typvar[3]={' ',' ','\0'};
  char nomvar[5]={' ',' ',' ',' ','\0'};
  char grtyp[2]={' ','\0'};

  handle = c_fstinf(iun,&ni,&nj,&nk,datein,etiketin,
                        ip1in,ip2in,ip3in,typvarin,nomvarin);

  if (handle < 0) {
    sprintf(errmsg,"(unit=%d) record not found",iun);
    return(error_msg("c_fstlic",ERR_NOT_FOUND,ERROR));
  }

  if ((niin != ni) || (njin != nj) || (nkin != nk)) {
    fprintf(stderr,"c_fstlic unit=%d, ni or nj or nk error:",iun);
    fprintf(stderr," user ni=%d, file ni=%d\n",niin,ni);
    fprintf(stderr," user nj=%d, file nj=%d\n",njin,nj);
    fprintf(stderr," user nk=%d, file nk=%d\n",nkin,nk);
    return(-1);
  }

  ier = c_fstprm(handle,&date,&deet,&npas,&ni,&nj,&nk,&nbits,&datyp,&ip1,
                     &ip2,&ip3,typvar,nomvar,etiket,grtyp,&ig1,&ig2,
                     &ig3,&ig4,&swa,&lng,&dltf,&ubc,&xtra1,&xtra2,&xtra3);

  if ((strcmp(grtypin,grtyp) != 0) || (ig1in != ig1) ||
      (ig2in != ig2) || (ig3in != ig3) || (ig4in != ig4)) {
    fprintf(stderr,"c_fstlic unit=%d, grtyp ig1 ig2 ig3 ig4 error:",iun);
    fprintf(stderr," user grtyp=%s, file grtyp=%s\n",grtypin,grtyp);
    fprintf(stderr," user ig1=%d, file ig1=%d\n",ig1in,ig1);
    fprintf(stderr," user ig2=%d, file ig2=%d\n",ig2in,ig2);
    fprintf(stderr," user ig3=%d, file ig3=%d\n",ig3in,ig3);
    fprintf(stderr," user ig4=%d, file ig4=%d\n",ig4in,ig4);
    return(-1);
  }

  ier = c_fstlir(field,iun,&ni,&nj,&nk,datein,etiketin,ip1in,ip2in,
                 ip3in,typvarin,nomvarin);
  return(ier);
}
/*splitpoint c_fstlir */
/*****************************************************************************
 *                        C _ F S T L I R                                    *
 *                                                                           *
 *Object                                                                     *
 *   Reads the next record that matches the research keys.                   *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  OUT field   data field to be read                                        *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *  IN  datev   valid date                                                   *
 *  IN  etiket  label                                                        *
 *  IN  ip1     vertical level                                               *
 *  IN  ip2     forecast hour                                                *
 *  IN  ip3     user defined identifier                                      *
 *  IN  typvar  type of field                                                *
 *  IN  nomvar  variable name                                                *
 *                                                                           *
 *****************************************************************************/

int c_fstlir(word *field, int iun, int *ni, int *nj, int *nk,
                 int datev, char *etiket,
                 int ip1, int ip2, int ip3, char *typvar, char *nomvar)
{
  int key,handle;

  handle = -2;  /* means handle will be discarded */
  key = c_fstlirx(field,handle,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,
                      typvar,nomvar);
  return(key);
}

/*splitpoint c_fstlirx */
/*****************************************************************************
 *                      C _ F S T L I R X                                    *
 *                                                                           *
 *Object                                                                     *
 *   Reads the next record that matches the research keys.                   *
 *   The search begins at the position given by handle.                      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  OUT field   data field to be read                                        *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *  IN  datev   valid date                                                   *
 *  IN  etiket  label                                                        *
 *  IN  ip1     vertical level                                               *
 *  IN  ip2     forecast hour                                                *
 *  IN  ip3     user defined identifier                                      *
 *  IN  typvar  type of field                                                *
 *  IN  nomvar  variable name                                                *
 *                                                                           *
 *****************************************************************************/

int c_fstlirx(word *field, int handle, int iun,
                  int *ni, int *nj, int *nk, int datev, char *etiket,
                  int ip1, int ip2, int ip3, char *typvar, char *nomvar)
{
  int ier;

  handle = c_fstinfx(handle,iun,ni,nj,nk,datev,etiket,
                         ip1,ip2,ip3,typvar,nomvar);
  if (handle < 0) {
    if (msg_level <= WARNING)
      fprintf(stdout,"c_fstlirx: (unit=%d) record not found, errcode=%d\n",
              iun,handle);
    return(handle);
  }

  ier = c_fstluk(field,handle,ni,nj,nk);
  if (ier < 0)
    return(ier);
  else
    return(handle);
}

/*splitpoint c_fstlis */
/*****************************************************************************
 *                         C _ F S T L I S                                   *
 *                                                                           *
 *Object                                                                     *
 *   Reads the next record that matches the last search criterias            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  OUT field   data field to be read                                        *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *                                                                           *
 *****************************************************************************/

int c_fstlis(word *field, int iun, int *ni, int *nj, int *nk)
{
  int index_fnom, index, handle, ier;
  word *primk=NULL;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstlis",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstlis",ERR_NO_FILE,ERROR));
  }

  /* position to the next record that matches the last search criterias */
  handle = c_xdfloc(iun,-1,primk,0); /* find next with handle=-1 and nprim=0 */
  if (handle < 0) {
    if (msg_level <= WARNING)
      fprintf(stdout,"c_fstlis: (unit=%d) record not found, errcode=%d\n",
              iun,handle);
    return(handle);
  }

  ier = c_fstluk(field,handle,ni,nj,nk);
  return(ier);
}


/*splitpoint c_fstluk */
/*****************************************************************************
 *                         C _ F S T L U K                                   *
 *                                                                           *
 *Object                                                                     *
 *   Read the record at position given by handle.                            *
 *Revision                                                                   *
 *   Mar  2013 - M.Valin  missing values, byte/short mods                    *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  OUT field   data field to be read                                        *
 *  IN  handle  positioning information to the record                        *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *                                                                           *
 *****************************************************************************/

int c_fstluk(word *field, int handle, int *ni, int *nj, int *nk)
{
  stdf_dir_keys *stdf_entry;
  word *pkeys, *work_field;
  float *ptr_real;
  double *ptr_double;
  buffer_interface_ptr buf;
  int ier, addr, lng, lng2, idtyp, mode, npak, nelm, i, zero=0, one=1;
  int nbytes, header_size, stream_size, p1out, p2out;
  int bitmot=32;
  int lngw, stdf_aux_keys[2];
  char string[11];
  PackFunctionPointer packfunc;
  double tempfloat=99999.0;
  int has_missing = 0;    /* missing data flag (bit with value 64 in datatype) */
  int *field_out;
  short *s_field_out;
  signed char *b_field_out;

  stdf_entry = (stdf_dir_keys *) calloc(1,sizeof(stdf_dir_keys));
  pkeys = (word *) stdf_entry;
  pkeys += W64TOWD(1);

  ier = c_xdfprm(handle,&addr,&lng,&idtyp,pkeys,16);
  if (ier < 0) return(ier);
  *ni = stdf_entry->ni;
  *nj = stdf_entry->nj;
  *nk = stdf_entry->nk;
  has_missing = stdf_entry->datyp & 64 ;          /* get missing data flag */
  stdf_entry->datyp = stdf_entry->datyp & 0xBF ;  /* suppress missing data flag */
  xdf_datatyp = stdf_entry->datyp ;

  if (xdf_double)
    packfunc = &compact_double;
  else
    packfunc = &compact_float;

  lng = W64TOWD(lng);
  if ((xdf_datatyp == 1) || (xdf_datatyp == 5))
    lng = (xdf_double) ? 2*lng : lng;


  if ((xdf_datatyp == 6) || (xdf_datatyp == 134)) {       /* new packer */
    c_float_packer_params(&header_size,&stream_size,&p1out,&p2out,(*ni)*(*nj)*(*nk));
    header_size /= sizeof(INT_32);
    lng2 = 1 + ((*ni * *nj * *nk* 16 + 32 + 31) / 32) + header_size + 20;
    }
  else if (xdf_datatyp == 133) {  /* compressed ieee */
    lng2 = 1 + lng;
    }
  else if (xdf_datatyp > 128) {
    lng2 = 4 + ((*ni * *nj * *nk * 16 + 32 + 31) / 32) + 20;          /* 16 for 16 bits armn_compress */
/*   fprintf(stderr,"Debug+ fstluk ni=%d nj=%d nk=%d lng=%d lng2=%d\n",*ni,*nj,*nk,lng,lng2); */
    }
  else
    lng2 = lng;

/*  printf("Debug+ fstluk lng2 = %d\n",lng2); */
  /* allocate 8 more bytes in case of realingment for 64 bit data */
  if ((work_field = alloca(8+(lng2+10)*sizeof(word))) == NULL) {
    sprintf(errmsg,"memory is full, was trying to allocate %ld bytes",
            lng*sizeof(word));
    return(error_msg("c_fstluk",ERR_MEM_FULL,ERRFATAL));
  }
  else
    memset(work_field,0,8+(lng2+10)*sizeof(word));

  buf = (buffer_interface_ptr) work_field;
  if ( (((&(buf->data[0]) - &(buf->nwords)) * sizeof(word)) & 0x7) != 0 ) {
    /* realign buf to make sure that buf->data is 64bit align */
    buf = (buffer_interface_ptr) (work_field + 1);
  }
  buf->nwords = -(lng+10);   /* negative value means get data only */
  buf->nbits = -1;
  ier = c_xdfget2(handle,buf,stdf_aux_keys);
  if (ier < 0) return(ier);
  if ((stdf_aux_keys[0] != 0) && (stdf_aux_keys[1] != 0)) {
    printf("c_fstluk aux_keys[0]=%d, aux_keys[1]=%d\n",stdf_aux_keys[0],stdf_aux_keys[1]);
    sprintf(errmsg,"wrong version of fstd98 (%d), recompile with a more recent version",
            stdf_version);
    return(error_msg("c_fstluk",ERR_STDF_VERSION,ERRFATAL));
  }


  nelm = stdf_entry->ni * stdf_entry->nj * stdf_entry->nk;
  if (stdf_entry->datyp == 8) nelm *= 2;

  npak = -(stdf_entry->nbits);
  if (image_mode_copy) {            /* no pack/unpack, used by editfst */
    if (stdf_entry->datyp > 128) {
      lngw = buf->data[0];
/*      fprintf(stderr,"Debug+ lecture mode image lngw=%d\n",lngw); */
      for (i=0; i < lngw+1; i++) {
        field[i] =  buf->data[i];
      }
    }
    else {
      lngw = nelm * stdf_entry->nbits;
      if (stdf_entry->datyp == 1) lngw += 120;
      if (stdf_entry->datyp == 3) lngw = *ni * *nj * 8;
      if (stdf_entry->datyp == 6) {
        c_float_packer_params(&header_size,&stream_size,&p1out,&p2out,nelm);
        lngw = (header_size+stream_size) * 8 ;
        }
      lngw = (lngw + bitmot - 1) / bitmot;
      for (i=0; i < lngw; i++) {
        field[i] =  buf->data[i];
      }
    }
  }
  else {
    switch (stdf_entry->datyp)
      {
      case 0:  /* raw binary */
        lngw = ((nelm * stdf_entry->nbits) + bitmot - 1) / bitmot;
        for (i=0; i < lngw; i++)
          field[i] =  buf->data[i];
        break;

      case 1: case 129:  /* floating point */
        if (stdf_entry->datyp > 128) {
/*          fprintf(stderr,"Debug+ unpack buf->data=%d\n",*(buf->data)); */
          nbytes = armn_compress(buf->data+5,*ni,*nj,*nk,stdf_entry->nbits,2);
/*          fprintf(stderr,"Debug+ buf->data+4+(nbytes/4)-1=%X buf->data+4+(nbytes/4)=%X \n",
              *(buf->data+4+(nbytes/4)-1),*(buf->data+4+(nbytes/4))); */
          packfunc(field,buf->data+1,buf->data+5,
                   nelm,stdf_entry->nbits+64*Max(16,stdf_entry->nbits),0,xdf_stride,FLOAT_UNPACK,0,
                   &tempfloat);
        }
        else
          packfunc(field,buf->data,buf->data+3,
                   nelm,stdf_entry->nbits,24,xdf_stride,FLOAT_UNPACK,0,
                   &tempfloat);
        break;

      case 2: case 130:          /* integer, short integer or byte stream */
        {
        int offset;
        offset = (stdf_entry->datyp > 128) ? 1 : 0;
        if (xdf_short) {
          if (stdf_entry->datyp > 128) {
            c_armn_compress_setswap(0);
            nbytes = armn_compress(buf->data+offset,*ni,*nj,*nk,stdf_entry->nbits,2);
          /*  printf("Debug+ fstluk mode short compress nbytes=%d\n",nbytes); */
            c_armn_compress_setswap(1);
            memcpy(field,buf->data+offset,nbytes);
            }
          else {
            mode = 6;
            ier = compact_short(field,(void *) NULL,buf->data+offset,nelm,
                                stdf_entry->nbits,0,xdf_stride,mode);
            }
          }
        else if (xdf_byte) {
          if (stdf_entry->datyp > 128) {
            c_armn_compress_setswap(0);
            nbytes = armn_compress(buf->data+offset,*ni,*nj,*nk,stdf_entry->nbits,2);
            c_armn_compress_setswap(1);
        /*    printf("Debug+ fstluk xdf_byte armn_compress nbytes=%d nelm=%d\n",nbytes,nelm); */
            memcpy_16_8(field,buf->data+offset,nelm);
          }
          else {
            mode = 10;
            ier = compact_char(field,(void *) NULL,buf->data,nelm,
                               8,0,xdf_stride,mode);
          }
        }
        else {
          if (stdf_entry->datyp > 128) {
            c_armn_compress_setswap(0);
            nbytes = armn_compress(buf->data+offset,*ni,*nj,*nk,stdf_entry->nbits,2);
            c_armn_compress_setswap(1);
          /*  printf("Debug+ fstluk mode int compress nbytes=%d\n",nbytes); */
            memcpy_16_32(field,buf->data+offset,stdf_entry->nbits,nelm);
          }
          else {
            mode = 2;
            ier = compact_integer(field,(void *) NULL,buf->data+offset,nelm,
                                  stdf_entry->nbits,0,xdf_stride,mode);
          }
        }
        break;
        }

      case 3: mode=2;  /* character */
        {
        int nc = (nelm+3)/4;
        ier = compact_integer(field,(void *) NULL,buf->data,nc,
                              32,0,xdf_stride,mode);
        break;
        }

      case 4: mode=4;  /* signed integer */
#ifdef use_old_signed_pack_unpack_code
// fprintf(stderr,"OLD UNPACK CODE ======================================\n");
        if(xdf_short || xdf_byte){
          field_out=alloca(nelm*sizeof(int));
          s_field_out=(short *)field;
          b_field_out=(signed char *)field;
        }else{
          field_out=(int *)field;
        }
        ier = compact_integer(field_out,(void *) NULL,buf->data,nelm,
                              stdf_entry->nbits,0,xdf_stride,mode);
        if(xdf_short){ for (i=0;i<nelm;i++) s_field_out[i]=field_out[i]; } ;
        if(xdf_byte) { for (i=0;i<nelm;i++) b_field_out[i]=field_out[i]; } ;
#else
// fprintf(stderr,"NEW UNPACK CODE ======================================\n");
        if(xdf_short){
          ier = compact_short(field,(void *) NULL,buf->data,nelm,
                              stdf_entry->nbits,0,xdf_stride,8);
        }else if(xdf_byte){
          ier = compact_char(field,(void *) NULL,buf->data,nelm,
                              stdf_entry->nbits,0,xdf_stride,12);
        }else{
          ier = compact_integer(field,(void *) NULL,buf->data,nelm,
                              stdf_entry->nbits,0,xdf_stride,mode);
        }
#endif
        break;

      case 5: case 8: mode=2;  /* IEEE representation */
        {
          register INT_32 temp32,*src,*dest;
          if ((downgrade_32) && (stdf_entry->nbits == 64)) {       /* downgrade 64 bit to 32 bit */
            ptr_real = (float *) field;
            ptr_double = (double *) buf->data;
#if defined(Little_Endian)
            src = (INT_32 *) buf->data;
            dest = (INT_32 *) buf->data;
            for (i=0; i < nelm; i++) {
              temp32 = *src++;
              *dest++ = *src++;
              *dest++ = temp32;
            }
#endif
            for (i=0; i < nelm; i++) {
              *ptr_real++ = *ptr_double++ ;
            }
          }
          else {
            ftnword f_nelm = (ftnword) nelm;
            ftnword f_one = (ftnword) one;
            ftnword f_npak = (ftnword) npak;
            ftnword f_zero = (ftnword) zero;
            ftnword f_mode = (ftnword) mode;
            f77name(ieeepak)(field,buf->data,&f_nelm,&f_one,&f_npak,&f_zero,&f_mode);
          }
        }
        break;

      case 6: case 134:  /* floating point, new packers */
        {
          int nbits;

          if (stdf_entry->datyp > 128) {
            nbytes = armn_compress(buf->data+1+header_size,*ni,*nj,*nk,stdf_entry->nbits,2);
  /*          fprintf(stderr,"Debug+ buf->data+4+(nbytes/4)-1=%X buf->data+4+(nbytes/4)=%X \n",
                *(buf->data+4+(nbytes/4)-1),*(buf->data+4+(nbytes/4))); */
            c_float_unpacker(field,buf->data+1,buf->data+1+header_size,nelm,&nbits);
          }
          else {
            c_float_unpacker(field,buf->data,buf->data+header_size,nelm,&nbits);
            }
          break;
        }

      case 133:  /* floating point, new packers */
          nbytes = c_armn_uncompress32(field, buf->data+1, *ni,*nj,*nk,stdf_entry->nbits);
          break;

      case 7: mode=10;  /* character string */
/*        printf("Debug fstluk compact_char xdf_stride=%d nelm =%d\n",xdf_stride,nelm); */
        ier = compact_char(field,(void *) NULL,buf->data,nelm,
                              8,0,xdf_stride,mode);
        break;

      default: sprintf(errmsg,"invalid datyp=%d",stdf_entry->datyp);
        return(error_msg("c_fstluk",ERR_BAD_DATYP,ERROR));

      } /* end switch */
  }

  if (msg_level <= INFORM) {
    sprintf(string,"Read(%d)",buf->iun);
    stdf_entry->datyp = stdf_entry->datyp | has_missing;
    print_std_parms(stdf_entry,string,prnt_options,0);
  }
  if(has_missing) {
    /* replace "missing" data points with the appropriate values given the type of data (int/float) */
    /* if nbits = 64 and IEEE , set xdf_double */
    if((stdf_entry->datyp & 0xF) == 5 && stdf_entry->nbits == 64 ) xdf_double=1;
    DecodeMissingValue( field , (*ni)*(*nj)*(*nk) , xdf_datatyp&0x3F,xdf_byte,xdf_short,xdf_double ); /* */
  }
  free(stdf_entry);
/* free(work_field);        replaced by alloca */
  xdf_double = 0;
  xdf_short = 0;
  xdf_byte = 0;
  return(handle);
}

/*splitpoint c_fstmsq */
/*****************************************************************************
 *                        C _ F S T M S Q                                    *
 *                                                                           *
 *Object                                                                     *
 *   Mask a portion of the research keys.                                    *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *   IN    iun     unit number associated to the file                        *
 * IN/OUT  mip1    mask for vertical level                                   *
 * IN/OUT  mip2    mask for forecast hour                                    *
 * IN/OUT  mip3    mask for ip3 (user defined identifier)                    *
 * IN/OUT  metiket mask for label                                            *
 *   IN    getmode logical (1: getmode 0:set mode)                           *
 *                                                                           *
 *****************************************************************************/

int c_fstmsq(int iun, int *mip1, int *mip2, int *mip3, char *metiket,
                     int getmode)
{
  stdf_dir_keys *search_mask;
  int index, i;
  file_table_entry *f;

#define isignore(c) ((c == '*') ? 0 : 0x3f)
#define inv_isignore(c) ((c == 0x3f) ? ' ' : '*')

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstmsq",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (! f->cur_info->attr.std) {
    sprintf(errmsg,"file (unit=%d) is not a RPN standard file",iun);
    return(error_msg("c_fstmsq",ERR_NO_FILE,ERROR));
  }

  search_mask = (stdf_dir_keys *) f->srch_mask;
  if (getmode) {
    *mip1 = ~(search_mask->ip1) & 0xfffffff;
    *mip2 = ~(search_mask->ip2) & 0xfffffff;
    *mip3 = ~(search_mask->ip3) & 0xfffffff;
    for (i=0; i <= 4; i++)
      metiket[i] = inv_isignore( ((search_mask->etik15 >> ((4-i)*6)) & 0x3f) );

    for (i=5; i <=9; i++)
      metiket[i] = inv_isignore( ((search_mask->etik6a >> ((9-i)*6)) & 0x3f) );

    metiket[10] = inv_isignore( ((search_mask->etikbc >> 6) & 0x3f) );
    metiket[11] = inv_isignore((search_mask->etikbc  & 0x3f));
    metiket[12] = '\0';
  }
  else {
    search_mask->ip1 = ~(*mip1) & 0xfffffff;
    search_mask->ip2 = ~(*mip2) & 0xfffffff;
    search_mask->ip3 = ~(*mip3) & 0xfffffff;
    search_mask->etik15 =
      (isignore(metiket[0]) << 24) |
      (isignore(metiket[1]) << 18) |
      (isignore(metiket[2]) << 12) |
      (isignore(metiket[3]) <<  6) |
      (isignore(metiket[4]));
    search_mask->etik6a =
      (isignore(metiket[5]) << 24) |
      (isignore(metiket[6]) << 18) |
      (isignore(metiket[7]) << 12) |
      (isignore(metiket[8])<<  6) |
      (isignore(metiket[9]));
    search_mask->etikbc =
      (isignore(metiket[10]) <<  6) |
      (isignore(metiket[11]));
  }
  return(0);
}

/*splitpoint c_fstnbr */
/*****************************************************************************
 *                          C _ F S T N B R                                  *
 *                                                                           *
 *Object                                                                     *
 *   Returns the number of records of the file associated with unit number.  *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/

int c_fstnbr(int iun)
{
  int index, index_fnom, ier, nrec;
  file_table_entry *f;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstnbr",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    /*    sprintf(errmsg,"file (unit=%d) is not open",iun); */
    /*    return(error_msg("c_fstnbr",ERR_NO_FILE,ERROR)); */
    ier = c_fstouv(iun,"RND");
    index = file_index(iun);
    f = file_table[index];
    nrec = f->nrecords;
    ier = c_fstfrm(iun);
    return(nrec);
  }

  f = file_table[index];
  return(f->nrecords);

}

/*splitpoint c_fstnbrv */
/*****************************************************************************
 *                          C _ F S T N B R V                                *
 *                                                                           *
 *Object                                                                     *
 *   Returns the number of valid records (excluding deleted records) of the  *
 *   file associated with unit number.                                       *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/

int c_fstnbrv(int iun)
{
  int index, index_fnom, ier, nrec;
  file_table_entry *f;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstnbrv",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    ier = c_fstouv(iun,"RND");
    index = file_index(iun);
    f = file_table[index];
    nrec = f->header->nrec;
    ier = c_fstfrm(iun);
    return(nrec);
  }

  f = file_table[index];
  return(f->header->nrec);

}

/*splitpoint c_fstopc */
/*****************************************************************************
 *                          C _ F S T O P C                                  *
 *                                                                           *
 *Object                                                                     *
 *   Print out, get, or set a fstd or xdf global variable option.            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *   IN     option   option name to be set/printed                           *
 *   IN     value    option value                                            *
 *   IN     getmode  (1: print option, 0: set option, 2: get option)         *
 *                                                                           *
 *****************************************************************************/
static int turbocomp_mode=0;
static char *comptab[2] = {"FAST","BEST"};
static char *msgtab[7] = {"DEBUG","INFORM","WARNIN","ERRORS","FATALE","SYSTEM","CATAST"};
static int nivmsg[7] = {0,2,4,6,8,10,10};

int c_fstopc(char *option, char *value, int getmode)
{
  int i;
  int val = 0;

  if (strcmp(option,"MSGLVL") == 0) {
    if (getmode){
      if (getmode == 2) val = msg_level;
    }else{
      for (i = 0; i < 7; i++) {
        if (strcmp(msgtab[i],value) == 0) {
          msg_level = i;
          break;
        }
      }
      c_armn_compress_option(option,value);  /* pass msglvl option to compressor */
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopc MSGLVL option = '%s' , %d\n",msgtab[msg_level],nivmsg[msg_level]);
    return val;
  }

  if (strcmp(option,"TOLRNC") == 0) {
    if (getmode){
      if (getmode == 2) val = xdf_toler;
    }else{
      for (i = 0; i < 7; i++){
        if (strcmp(msgtab[i],value) == 0) {
          xdf_toler = i;
          break;
        }
      }
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopc option %s, TOLRNC=%d\n",msgtab[xdf_toler],nivmsg[xdf_toler]);
    return val;
  }

  if (strcmp(option,"PRINTOPT") == 0) {
    if (getmode){
      if (getmode == 2) val = 0;
    }else{
      sprintf(prnt_options,"%s",value);
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopc option PRINTOPT='%s'\n",prnt_options);
    return val;
  }

  if (strcmp(option,"TURBOCOMP") == 0) {
    if (getmode){
      if (getmode == 2) val = turbocomp_mode;
    }else{
      for (i = 0; i < 2; i++) {
        if (strcmp(comptab[i],value) == 0) {
          turbocomp_mode = i;
          break;
        }
      }
      c_armn_compress_setlevel(turbocomp_mode);
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopc option TURBOCOMP=%s\n",comptab[turbocomp_mode]);
    return val;
  }

  fprintf(stderr,"c_fstopc: unknown option %s\n",option);
  return(val);
}

/*splitpoint c_fstopi */
/*****************************************************************************
 *                          C _ F S T O P I                                  *
 *                                                                           *
 *Object                                                                     *
 *   Print out, get, or set a fstd or xdf global variable option.            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *   IN     option   option name to be set/printed                           *
 *   IN     value    option value                                            *
 *   IN     getmode  (1: print option, 0: set option, 2: get option)         *
 *                                                                           *
 *****************************************************************************/
int c_fstopi(char *option, int value, int getmode)
{
  int i;
  int val = 0;

  if (strcmp(option,"MSGLVL") == 0) {
    if (getmode){
      if (getmode == 2) val = msg_level;
    }else{
      for (i = 0; i < 7; i++) {
        if (nivmsg[i] == value) {
          msg_level = i;
          break;
        }
      }
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopi  MSGLVL option = '%s' , %d\n",msgtab[msg_level],nivmsg[msg_level]);
    return val;
  }

  if (strcmp(option,"TOLRNC") == 0) {
    if (getmode){
      if (getmode == 2) val = xdf_toler;
    }else{
      for (i = 0; i < 7; i++){
        if (nivmsg[i] == value) {
          xdf_toler = i;
          break;
        }
      }
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopi option TOLRNC=%d\n",nivmsg[xdf_toler]);
    return val;
  }

  if (strcmp(option,"TURBOCOMP") == 0) {
    if (getmode){
      if (getmode == 2) val = turbocomp_mode;
    }else{
      if(value == 0 || value == 1){
        turbocomp_mode = value;
        c_armn_compress_setlevel(turbocomp_mode);
      }
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopi option TURBOCOMP= '%s' , %d\n",comptab[turbocomp_mode],turbocomp_mode);
    return val;
  }

  fprintf(stderr,"c_fstopi: unknown option %s\n",option);
  return(val);
}

/*splitpoint c_fstopl */
/*****************************************************************************
 *                          C _ F S T O P L                                  *
 *                                                                           *
 *Object                                                                     *
 *   Print out, get, or set a fstd or xdf global variable option.            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *   IN     option   option name to be set/printed                           *
 *   IN     value    option value                                            *
 *   IN     getmode  (1: print option, 0: set option, 2: get option)         *
 *                                                                           *
 *****************************************************************************/
int c_fstopl(char *option, int value, int getmode)
{
  int val = 0;

  if (strcmp(option,"FASTIO") == 0) {
    if (getmode){
      if (getmode == 2) val = 0;
    }else{
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopl FASTIO mode NOT ACTIVE\n");
    return val;
  }

  if (strcmp(option,"IMAGE") == 0) {
    if (getmode){
      if (getmode == 2) val = image_mode_copy;
    }else{
      image_mode_copy = value;
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopl option IMAGE_MODE_COPY=%d\n",image_mode_copy);
    return val;
  }

  if (strcmp(option,"REDUCTION32") == 0) {
    if (getmode){
      if (getmode == 2) val = downgrade_32;
    }else {
      downgrade_32 = value;
    }
    if(getmode == 1 || msg_level <= INFORM) fprintf(stdout,"c_fstopl option REDUCTION32=%d\n",downgrade_32);
    return val;
  }

  fprintf(stderr,"c_fstopi: unknown option %s\n",option);
  return(val);
}


/*splitpoint c_fstopr */
/*****************************************************************************
 *                          C _ F S T O P R                                  *
 *                                                                           *
 *Object                                                                     *
 *   Print out, get, or set a fstd or xdf global variable option.            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *   IN     option   option name to be set/printed                           *
 *   IN     value    option value                                            *
 *   IN     getmode  (1: print option, 0: set option, 2: get option)         *
 *                                                                           *
 *****************************************************************************/
int c_fstopr(char *option, float value, int getmode)
{
  /* no current float variable to be set for now */
  return(0);
}

/*splitpoint c_fstouv */
/*****************************************************************************
 *                         C _ F S T O U V                                   *
 *                                                                           *
 *Object                                                                     *
 *   Opens a RPN standard file.                                              *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  IN  options random or sequential access                                  *
 *                                                                           *
 *****************************************************************************/

int c_fstouv(int iun, char *options)
{
  int ier,nrec,i,iwko;
  static int premiere_fois=1;
//  char *turbo_compression;
  char appl[5];

  if (premiere_fois) {
    premiere_fois = 0;
/*    printf("DEBUG++ fstouv appel a c_env_var_cracker\n"); */
    c_env_var_cracker("FST_OPTIONS", c_fst_env_var, "C");    /* obtain options from environment variable */
    C_requetes_init(requetes_filename,debug_filename);
    ier = init_ip_vals();
  }
  i = fnom_index(iun);
  if (i == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstouv",ERR_NO_FNOM,ERROR));
  }

  if ((strstr(options,"RND"))  || (strstr(options,"rnd")))
    sprintf(appl,"%s","STDR");      /* standard random */
  else
    sprintf(appl,"%s","STDS");      /* standard sequential */

  FGFDT[i].attr.std = 1;  /* force attribute to standard file */
  if (FGFDT[i].attr.remote)
    if ((FGFDT[i].eff_file_size == 0) && (! FGFDT[i].attr.old))
      ier = c_xdfopn(iun,"CREATE",(word_2 *) &stdfkeys,16,(word_2 *) &stdf_info_keys,2,appl);
    else
      ier = c_xdfopn(iun,"R-W",(word_2 *) &stdfkeys,16,(word_2 *) &stdf_info_keys,2,appl);
  else
    if (((iwko=c_wkoffit(FGFDT[i].file_name,strlen(FGFDT[i].file_name))) == -2) &&
       (! FGFDT[i].attr.old)) {
      ier = c_xdfopn(iun,"CREATE",(word_2 *) &stdfkeys,16,(word_2 *) &stdf_info_keys,2,appl);
      }
    else {
      ier = c_xdfopn(iun,"R-W",(word_2 *) &stdfkeys,16,(word_2 *) &stdf_info_keys,2,appl);
      }

  if (ier < 0) return(ier);
  nrec = c_fstnbr(iun);
  return(nrec);
}

/*splitpoint c_fstprm */
/*****************************************************************************
 *                       C _ F S T P R M                                     *
 *                                                                           *
 *Object                                                                     *
 *   Get all the description informations of the record.                     *
 *                                                                           *
 *Rev 001 - M. Lepine - Oct 2002, returns extra1 as the validity datestamp   *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  handle  positioning information to the record                        *
 *  OUT date    date time stamp                                              *
 *  OUT deet    length of a time step in seconds                             *
 *  OUT npas    time step number                                             *
 *  OUT ni      first dimension of the data field                            *
 *  OUT nj      second dimension of the data field                           *
 *  OUT nk      third dimension of the data field                            *
 *  OUT nbits   number of bits kept for the elements of the field            *
 *  OUT datyp   data type of the elements                                    *
 *  OUT ip1     vertical level                                               *
 *  OUT ip2     forecast hour                                                *
 *  OUT ip3     user defined identifier                                      *
 *  OUT typvar  type of field (forecast, analysis, climatology)              *
 *  OUT nomvar  variable name                                                *
 *  OUT etiket  label                                                        *
 *  OUT grtyp   type of geographical projection                              *
 *  OUT ig1     first grid descriptor                                        *
 *  OUT ig2     second grid descriptor                                       *
 *  OUT ig3     third grid descriptor                                        *
 *  OUT ig4     fourth grid descriptor                                       *
 *  OUT swa     starting word address                                        *
 *  OUT lng     record length                                                *
 *  OUT dltf    delete flag                                                  *
 *  OUT ubc     unused bit count                                             *
 *  OUT extra1  extra parameter                                              *
 *  OUT extra2  extra parameter                                              *
 *  OUT extra3  extra parameter                                              *
 *                                                                           *
 *****************************************************************************/

int c_fstprm(int handle,
                     int *dateo, int *deet, int *npas,
                     int *ni, int *nj, int *nk,
                     int *nbits, int *datyp, int *ip1,
                     int *ip2, int *ip3, char *typvar,
                     char *nomvar, char *etiket, char *grtyp,
                     int *ig1, int *ig2, int *ig3,
                     int *ig4, int *swa, int *lng,
                     int *dltf, int *ubc, int *extra1,
                     int *extra2, int *extra3)
{
  stdf_dir_keys *stdf_entry;
  stdf_special_parms cracked;
  word *pkeys;
  int ier, addr, idtyp, xdflng, l1, l2, l3, l4;

  stdf_entry = (stdf_dir_keys *) calloc(1,sizeof(stdf_dir_keys));
  pkeys = (word *) stdf_entry;
  pkeys += W64TOWD(1);

  ier = c_xdfprm(handle,&addr,&xdflng,&idtyp,pkeys,16);

  crack_std_parms(stdf_entry,&cracked);

  *ni = stdf_entry->ni;
  *nj = stdf_entry->nj;
  *nk = stdf_entry->nk;
  *dateo = cracked.date_stamp;
  *deet = stdf_entry->deet;
  *npas = stdf_entry->npas;
  *nbits = stdf_entry->nbits;
  *datyp = stdf_entry->datyp;
  *ip1 = stdf_entry->ip1;
  *ip2 = stdf_entry->ip2;
  *ip3 = stdf_entry->ip3;
  *ig1 = stdf_entry->ig1;
  *ig2 = cracked.ig2;
  *ig3 = stdf_entry->ig3;
  *ig4 = stdf_entry->ig4;
  *swa = addr;
  *lng = W64TOWD(xdflng);
  *dltf = stdf_entry->deleted;
  *ubc = stdf_entry->ubc;
  *extra1 = cracked.date_valid;              /* new, use to be undefined */
  *extra2 = 0;
  *extra3 = 0;
  /*  strcpy(typvar,cracked.typvar);
  strcpy(nomvar,cracked.nomvar);
  strcpy(etiket,cracked.etiket);
  strcpy(grtyp,cracked.gtyp);
  strncpy(typvar,cracked.typvar,1);
  strncpy(nomvar,cracked.nomvar,2);
  strncpy(etiket,cracked.etiket,8);
  strncpy(grtyp,cracked.gtyp,1); */
  for (l1=0; (typvar[l1] != '\0') && (l1 <2); l1++);
  for (l2=0; (nomvar[l2] != '\0') && (l2 <4); l2++);
  for (l3=0; (etiket[l3] != '\0') && (l3 <12); l3++);
  l4 = 1;
  string_copy(typvar,cracked.typvar,l1);
  string_copy(nomvar,cracked.nomvar,l2);
  string_copy(etiket,cracked.etiket,l3);
  string_copy(grtyp,cracked.gtyp,l4);
  free(stdf_entry);
  return(ier);
}

/*splitpoint c_fstreset_ip_flags */
/*****************************************************************************
 *                 C _ F S T R E S E T _ I P _ F L A G S                     *
 *                                                                           *
 *Object                                                                     *
 *   Reset all the flags previously set by ip(1-3)_val                       *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
void c_fstreset_ip_flags()
{
  int ier;
  ier = init_ip_vals();
}

/*splitpoint c_fstrwd */
/*****************************************************************************
 *                           C _ F S T R W D                                 *
 *                                                                           *
 *Object                                                                     *
 *   Rewinds a RPN standard sequential file.                                 *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/
int c_fstrwd(int iun)
{
  int index, index_fnom;
  file_table_entry *f;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstrwd",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstrwd",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (! f->cur_info->attr.std) {
    sprintf(errmsg,"file (unit=%d) is not a RPN standard file",iun);
    return(error_msg("c_fstrwd",ERR_NO_FILE,ERROR));
  }

  if (! f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not sequential",iun);
    return(error_msg("c_fstrwd",ERR_BAD_FTYPE,WARNING));
  }

  f->cur_addr = f->seq_bof;
  f->valid_pos = 0;
  return(0);
}

/*splitpoint c_fstskp */
/*****************************************************************************
 *                            C _ F S T S K P                                *
 *                                                                           *
 *Object                                                                     *
 *   Skip nrec records forward or backward in the sequential file.           *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  IN  nrec    number of records to skip (negative nrec means backward)     *
 *                                                                           *
 *****************************************************************************/
int c_fstskp(int iun, int nrec)
{
  int index_fnom, index, i, nw, cur_pos, dim;
  file_table_entry *f;
  xdf_record_header header64;
  postfix_seq postfix;
  seq_dir_keys seq_entry;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstskp",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstskp",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (!f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not sequential",iun);
    return(error_msg("c_fstskp",ERR_BAD_FTYPE,WARNING));
  }

  if (! f->fstd_vintage_89) {
    if (nrec > 0) {           /* skip forward */
      for (i=0; i < nrec; i++) {
        nw = c_waread2(iun,&header64,f->cur_addr,W64TOWD(1));
        if ((nw != W64TOWD(1)) ||
            ((header64.idtyp >= 112) && (header64.idtyp <=127))) {
          if (msg_level <= INFORM)
            fprintf(stdout,"c_fstskp: (unit %d) skip to end of file\n",iun);
          break;
        }
        cur_pos = f->cur_addr;
        f->cur_addr += W64TOWD(header64.lng);
        c_waread(iun,&postfix,f->cur_addr,W64TOWD(2));
        if ((postfix.idtyp == 0) && (postfix.lng == 2) &&
            (postfix.addr == -1))
          f->cur_addr += W64TOWD(2);  /* skip postfix also */
        else {
          sprintf(errmsg,
                  "file (unit=%d) has invalid or no record postfix",iun);
          return(error_msg("c_fstskp",ERR_NO_POSTFIX,ERRFATAL));
        }
      }
    }
    else {
      nrec = -nrec;
      for (i=0; i < nrec; i++) {
        if ((f->cur_addr - W64TOWD(2)) > f->seq_bof) {
          c_waread(iun,&postfix,(f->cur_addr - W64TOWD(2)),W64TOWD(2));
          if ((postfix.idtyp == 0) && (postfix.lng == 2) &&
              (postfix.addr == -1))
            f->cur_addr = W64TOWD( (postfix.prev_addr-1) )+1;
          else {
            sprintf(errmsg,"file (unit=%d) has no record postfix",iun);
            return(error_msg("c_fstskp",ERR_NO_POSTFIX,ERRFATAL));
          }
          c_waread(iun,&header64,f->cur_addr,W64TOWD(1));
          if (header64.addr != (WDTO64( (f->cur_addr -1) )+1)) {
            sprintf(errmsg,
                    "file (unit=%d), postfix address (%d) not equal to record address (%d) ",iun,(WDTO64( (f->cur_addr -1) )+1),header64.addr);
            return(error_msg("c_fstskp",ERR_NO_POSTFIX,ERRFATAL));
          }
        }
      }
    }
  }
  else {        /* old sequential standard file */
    dim = sizeof(seq_entry) / sizeof(word);
    if (nrec < 0) {    /* skip forward */
      for (i=0; i < nrec; i++) {
        nw = c_waread2(iun,&seq_entry,f->cur_addr,dim);
        if ((nw != dim) || (seq_entry.eof > 0)) {
          if (msg_level <= INFORM)
            fprintf(stderr,"c_fstskp: (unit %d) skip to end of file\n",iun);
          break;
        }
        f->cur_addr += W64TOWD( (((seq_entry.lng + 3) >> 2) + 15) );
      }
    }
    else {             /* skip backward */
      nrec = -nrec;
      for (i=0; i < nrec; i++) {
        if (f->cur_addr <= 1) {
          f->cur_addr = 1;
          break;
        }
        c_waread(iun,&seq_entry,f->cur_addr,dim);
        cur_pos = f->cur_addr;
        f->cur_addr = seq_entry.swa_last * W64TOWD(15) + 1;
        if (f->cur_addr >= cur_pos) {    /* not moving backward anymore */
          f->cur_addr = 1;
          break;
        }
      }
    }
  }
  return(0);
}

/*splitpoint c_fstsui */
/*****************************************************************************
 *                         C _ F S T S U I                                   *
 *                                                                           *
 *Object                                                                     *
 *   Finds the next record that matches the last search criterias            *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  OUT ni      dimension 1 of the data field                                *
 *  OUT nj      dimension 2 of the data field                                *
 *  OUT nk      dimension 3 of the data field                                *
 *                                                                           *
 *****************************************************************************/

int c_fstsui(int iun, int *ni, int *nj, int *nk)
{
  int addr,lng,idtyp,ier,handle;
  stdf_dir_keys *stdf_entry;
  word *pkeys;
  word *primk=NULL;

  /* position to the next record that matches the last search criterias */
  handle = c_xdfloc(iun,-1,primk,0); /* find next with handle=-1 and nprim=0 */
  if (handle < 0) {
    if (msg_level == TRIVIAL)
      fprintf(stdout,"c_fstsui: record not found, errcode=%d\n",handle);
    return(handle);
  }

  stdf_entry = (stdf_dir_keys *) calloc(1,sizeof(stdf_dir_keys));
  pkeys = (word *) stdf_entry;
  pkeys += W64TOWD(1);

  ier = c_xdfprm(handle,&addr,&lng,&idtyp,pkeys,16);
  *ni = stdf_entry->ni;
  *nj = stdf_entry->nj;
  *nk = stdf_entry->nk;
  free(stdf_entry);
  return(handle);

}

/*splitpoint c_fst_version */
/*****************************************************************************
 *                       C _ F S T  _ V E R S I O N                          *
 *                                                                           *
 *Object                                                                     *
 *   Returns package version number.                                          *
 *                                                                           *
 *****************************************************************************/

int c_fst_version()
{
  return(stdf_version);
}

/*splitpoint c_fstvoi */
/*****************************************************************************
 *                       C _ F S T V O I                                     *
 *                                                                           *
 *Object                                                                     *
 *   Prints out the directory content of a RPN standard file.                *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *                                                                           *
 *****************************************************************************/

int c_fstvoi(int iun, char* options) {
   int index, index_fnom, i, j, width, nrec, nw, end_of_file;
   file_table_entry* file_entry;
   xdf_dir_page* curpage;
   word* entry;
   stdf_dir_keys* stdf_entry;
   seq_dir_keys* seq_entry;
   stdf_special_parms cracked;
   xdf_record_header* header;
   char string[20];
   ftnword f_datev;
   double nhours;
   int deet, npas, run;
   unsigned int datexx;
   long long deetnpas, i_nhours;

   index_fnom = fnom_index(iun);
   if (index_fnom == -1) {
      sprintf(errmsg, "file (unit=%d) is not connected with fnom", iun);
      return(error_msg("c_fstvoi", ERR_NO_FNOM, ERROR));
   }

   if ((index = file_index(iun)) == ERR_NO_FILE) {
      sprintf(errmsg, "file (unit=%d) is not open", iun);
      return(error_msg("c_fstvoi", ERR_NO_FILE, ERROR));
   }

   file_entry = file_table[index];

   if (! file_entry->cur_info->attr.std) {
      sprintf(errmsg, "file (unit=%d) is not a RPN standard file", iun);
      return(error_msg("c_fstvoi", ERR_NO_FILE, ERROR));
   }

   nrec = 0;
   width = W64TOWD(file_entry->primary_len);
   if (! file_entry->xdf_seq) {
     for (i=0; i < file_entry->npages; i++) {
       entry = (file_entry->dir_page[i])->dir.entry;
       for (j=0; j < (file_entry->dir_page[i])->dir.nent; j++) {
         header = (xdf_record_header *) entry;
         if (header->idtyp < 112) {
           stdf_entry = (stdf_dir_keys *) entry;
           crack_std_parms(stdf_entry,&cracked);
           sprintf(string,"%5d-",nrec);
           print_std_parms(stdf_entry,string,options,((nrec % 70) == 0));
           nrec++;
         }
         entry += width;
       }
       curpage = &((file_entry->dir_page[i])->dir);
     } /* end for i */
   }
   else { /* xdf sequential */
     end_of_file = 0;
     while (! end_of_file) {
       nw = c_waread2(iun,file_entry->head_keys,file_entry->cur_addr,width);
       header = (xdf_record_header *) file_entry->head_keys;
       if ((header->idtyp >= 112) || (nw < W64TOWD(1))) {
         if ((header->idtyp >= 112) && (header->idtyp < 127)) {
            file_entry->cur_addr += W64TOWD(1);
         }
         end_of_file = 1;
         break;
       }
       if (file_entry->fstd_vintage_89) {   /* old sequential standard */
         if ((stdf_entry = calloc(1,sizeof(stdf_dir_keys))) == NULL) {
           sprintf(errmsg,"memory is full");
           return(error_msg("c_fstvoi",ERR_MEM_FULL,ERRFATAL));
         }
         seq_entry = (seq_dir_keys *) file_entry->head_keys;
         if (seq_entry->dltf) {
           file_entry->cur_addr += W64TOWD( (((seq_entry->lng + 3) >> 2)+15) );
           continue;
         }
         if (seq_entry->eof > 0) {
           if (seq_entry->eof < 15)
           file_entry->cur_addr += W64TOWD(1);
           end_of_file = 1;
           break;
         }
         stdf_entry->deleted = 0;
         stdf_entry->select = 1;
         stdf_entry->lng = ((seq_entry->lng + 3) >> 2) + 15;
         stdf_entry->addr = (seq_entry->swa >> 2) +1;
         stdf_entry->deet = seq_entry->deet;
         stdf_entry->nbits = seq_entry->nbits;
         stdf_entry->ni = seq_entry->ni;
         stdf_entry->gtyp = seq_entry->grtyp;
         stdf_entry->nj = seq_entry->nj;
         stdf_entry->datyp = seq_entry->datyp;
         stdf_entry->nk = seq_entry->nk;
         stdf_entry->ubc = 0;
         stdf_entry->npas = (seq_entry->npas2 << 16) |
           seq_entry->npas1;
         stdf_entry->pad7 = 0;
         stdf_entry->ig4 = seq_entry->ig4;
         stdf_entry->ig2a = 0;
         stdf_entry->ig1 = seq_entry->ig1;
         stdf_entry->ig2b = seq_entry->ig2 >> 8;
         stdf_entry->ig3 = seq_entry->ig3;
         stdf_entry->ig2c = seq_entry->ig2 & 0xff;
         stdf_entry->etik15 =
           (ascii6(seq_entry->etiq14 >> 24) << 24) |
           (ascii6((seq_entry->etiq14 >> 16) & 0xff) << 18) |
           (ascii6((seq_entry->etiq14 >>  8) & 0xff) << 12) |
           (ascii6((seq_entry->etiq14      ) & 0xff) <<  6) |
           (ascii6((seq_entry->etiq56 >>  8) & 0xff));
         stdf_entry->pad1 = 0;
         stdf_entry->etik6a =
           (ascii6((seq_entry->etiq56      ) & 0xff) << 24) |
           (ascii6((seq_entry->etiq78 >>  8) & 0xff) << 18) |
           (ascii6((seq_entry->etiq78      ) & 0xff) << 12);
         stdf_entry->pad2 = 0;
         stdf_entry->etikbc = 0;
         stdf_entry->typvar = ascii6(seq_entry->typvar) << 6;
         stdf_entry->pad3 = 0;
         stdf_entry->nomvar =
           (ascii6((seq_entry->nomvar >>  8) & 0xff) << 18) |
           (ascii6((seq_entry->nomvar      ) & 0xff) << 12);
         stdf_entry->pad4 = 0;
         stdf_entry->ip1 = seq_entry->ip1;
         stdf_entry->levtyp = 0;
         stdf_entry->ip2 = seq_entry->ip2;
         stdf_entry->pad5 = 0;
         stdf_entry->ip3 = seq_entry->ip3;
         stdf_entry->pad6 = 0;
         stdf_entry->date_stamp = seq_entry->date;
         deet = stdf_entry->deet;
         npas = stdf_entry->npas;
         deetnpas = npas ; deetnpas = deetnpas * deet ;
         if ((deetnpas % 3600) != 0) {
            /* recompute datev to take care of rounding used with 1989 version
             * de-octalise the date_stamp */
            run = stdf_entry->date_stamp & 0x7;
            datexx = (stdf_entry->date_stamp >> 3) * 10 + run;

            f_datev = (ftnword) datexx;
            i_nhours = (deetnpas - ((deetnpas+1800)/3600)*3600);
            nhours = i_nhours;
            nhours = (nhours / 3600.0);
            f77name(incdatr)(&f_datev,&f_datev,&nhours);
            datexx = (unsigned int) f_datev;
            /* re-octalise the date_stamp */
            stdf_entry->date_stamp = 8 * (datexx/10) + (datexx % 10);
         }
         sprintf(string,"%5d-",nrec);
         print_std_parms(stdf_entry,string,options,((nrec % 70) == 0));
         nrec++;
         file_entry->cur_addr += W64TOWD( (((seq_entry->lng + 3) >> 2)+15) );
         free(stdf_entry);
       } /* end if fstd_vintage_89 */
       else {
         if ((header->idtyp < 1) || (header->idtyp > 127)) {
           file_entry->cur_addr += W64TOWD(header->lng);
           continue;
         }
         stdf_entry = (stdf_dir_keys *) file_entry->head_keys;
         sprintf(string,"%5d-",nrec);
         print_std_parms(stdf_entry,string,options,((nrec % 70) == 0));
         nrec++;
         file_entry->cur_addr += W64TOWD(header->lng);
       }
     } /* end while */
   }
   fprintf(stdout,"\nSTATISTICS for file %s, unit=%d\n\n",
           FGFDT[index_fnom].file_name,iun);
   if (file_entry->fstd_vintage_89)
     sprintf(string,"Version 1989");
   else
     sprintf(string,"Version 1998");
   if (file_entry->xdf_seq)
     fprintf(stdout,"%d records in sequential RPN standard file (%s)\n",
             nrec,string);
   else {
     if (! file_entry->fstd_vintage_89) {
       fprintf(stdout,"Number of directory entries \t %d\n",file_entry->header->nrec);
       fprintf(stdout,"Number of valid records     \t %d\n",nrec);
       fprintf(stdout,"File size                   \t %d Words\n",
               W64TOWD(file_entry->header->fsiz));
       fprintf(stdout,"Number of writes            \t %d\n",file_entry->header->nxtn);
       fprintf(stdout,"Number of rewrites          \t %d\n",file_entry->header->nrwr);
       fprintf(stdout,"Number of erasures          \t %d\n",
               file_entry->header->neff - file_entry->header->nrwr);
     }
     fprintf(stdout,"\n%d records in random RPN standard file (%s)\n\n",
             nrec,string);
   }
   return(0);
}

/*splitpoint c_fstweo */
/*****************************************************************************
 *                            C _ F S T W E O                                *
 *                                                                           *
 *Object                                                                     *
 *   Writes a logical end of file on a sequential file.                      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  iun     unit number associated to the file                           *
 *  IN  level   level of logical end of file                                 *
 *                                                                           *
 *****************************************************************************/

int c_fstweo(int iun, int level)
{
  int index_fnom, index;
  file_table_entry *f;
  xdf_record_header header64;

  index_fnom = fnom_index(iun);
  if (index_fnom == -1) {
    sprintf(errmsg,"file (unit=%d) is not connected with fnom",iun);
    return(error_msg("c_fstweo",ERR_NO_FNOM,ERROR));
  }

  if ((index = file_index(iun)) == ERR_NO_FILE) {
    sprintf(errmsg,"file (unit=%d) is not open",iun);
    return(error_msg("c_fstweo",ERR_NO_FILE,ERROR));
  }

  f = file_table[index];

  if (!f->xdf_seq) {
    sprintf(errmsg,"file (unit=%d) is not sequential",iun);
    return(error_msg("c_fstweo",ERR_BAD_FTYPE,WARNING));
  }

  VALID(level,1,15,"level","c_fstweo")
  ;
  if (level < 15) {
    header64.idtyp = 112 + level;
    header64.lng = 1;
    header64.addr = f->cur_addr;
    f->nxtadr = f->cur_addr;
    c_wawrit(iun,&header64,f->nxtadr,W64TOWD(1));
    f->nxtadr += W64TOWD(1);
    f->cur_addr += W64TOWD(1);
  }
  header64.idtyp = 127;
  header64.addr = f->cur_addr;
  c_wawrit(iun,&header64,f->cur_addr,W64TOWD(1));
  f->nxtadr = f->cur_addr;
  return(0);
}

/*splitpoint c_fst_env_var */
/*****************************************************************************
 *                         C _ F S T _ E N V _ V A R                         *
 *                                                                           *
 *Object                                                                     *
 *   Process all sub options from FST_OPTIONS environment variable           *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  cle     suboption name                                               *
 *  IN  index   index of suboption values                                    *
 *  IN  content suboption value                                              *
 *                                                                           *
 *****************************************************************************/

void c_fst_env_var(char *cle, int index, char *content)
{
//  int i;
  char *carac;

/*
  fprintf(stderr, "c_fst_env_var(), fstoption: %s\n", cle);
  fprintf(stderr, "c_fst_env_var(),    index = %d\n", index);
  fprintf(stderr, "c_fst_env_var(), content  : %s\n", content);
*/

  if (strcasecmp(cle,"TURBOCOMP") == 0) {
    carac = content;
    while (*carac) {
      *carac = toupper(*carac);
      carac++;
      }
    c_fstopc("TURBOCOMP",content,0);
    }
  else if (strcasecmp(cle,"DATATYPE_REMAP") == 0) {
    remap_table[(index-1) % 2][(index-1) / 2] = atoi(content);
    /*fprintf(stderr,"Debug+ remap_table[%d][%d]=%d\n",(index-1)%2,(index-1)/2,remap_table[(index-1) % 2][(index-1) / 2]); */
    nb_remap = index / 2;
    }
   else if (strcasecmp(cle,"DEBUGFILE") == 0) {
    debug_filename = malloc(256);
    strncpy(debug_filename,content,256);
/*    fprintf(stderr,"Debug+ debug_filename=%s\n",debug_filename); */
    }
   else if (strcasecmp(cle,"FST_FILTER_FILE") == 0) {
    requetes_filename = malloc(256);
    strncpy(requetes_filename,content,256);
/*    fprintf(stderr,"Debug+ requetes_filename=%s\n",requetes_filename); */
    }
  else
    fprintf(stderr, "c_fst_env_var(), cle %s non reconnue, index=%d valeur=%s\n",cle,index,content);

/*
  for (i=0; i<nb_remap; i++)
    fprintf(stderr,"datatyp %d remap to %d\n",remap_table[0][i],remap_table[1][i]);

  fprintf(stderr, "-----------------------------------------------\n");
*/
}


/*splitpoint c_ip1_all */
/*****************************************************************************
 *                          C _ I P 1 _ A L L                                *
 *                                                                           *
 *Object                                                                     *
 *   Generates all possible coded ip1 values for a given level               *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  level          ip1 level (float value)                               *
 *  IN  kind           level kind as defined in convip                       *
 *                                                                           *
 *****************************************************************************/

int c_ip1_all(float level, int kind)
{
  int ip_old, ip_new, mode, flag;
  char s[128];

  ip1s_flag = 1;
  flag = 0;

  mode = 2;
  ConvipPlus(&ip_new,&level,&kind,&mode);
  ips_tab[0][ip_nb[0]] = ip_new;
  ip_nb[0]++;
  if (ip_nb[0] >= Max_Ipvals) {
    fprintf(stderr,"ip1 table full (ip_nb=%d)\n",ip_nb[0]);
    return(-1);
    }

  mode = 3;
  if (kind < 4)
    ConvipPlus(&ip_old,&level,&kind,&mode);
  else
    ip_old = -9999;     /* no valid value for oldtype */
  ips_tab[0][ip_nb[0]] = ip_old;
  ip_nb[0]++;

  if (ip_nb[0] > Max_Ipvals) {
    fprintf(stderr,"ip1 table full (i1_ind=%d)\n",ip_nb[0]);
    return(-1);
    }
/*  printf("Debug+ c_ip1_all level=%f kind=%d ip_new=%d ip_old=%d\n",level,kind,ip_new,ip_old); */
  return(ip_new);
}


/*splitpoint c_ip2_all */
/*****************************************************************************
 *                          C _ I P 2 _ A L L                                *
 *                                                                           *
 *Object                                                                     *
 *   Generates all possible coded ip2 values for a given level               *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  level          ip2 level (float value)                               *
 *  IN  kind           level kind as defined in convip                       *
 *                                                                           *
 *****************************************************************************/

int c_ip2_all(float level, int kind)
{
  int ip_old, ip_new, mode, flag;
  char s[128];

  ip2s_flag = 1;
  flag = 0;

  mode = 2;
  ConvipPlus(&ip_new,&level,&kind,&mode);
  ips_tab[1][ip_nb[1]] = ip_new;
  ip_nb[1]++;
  if (ip_nb[1] >= Max_Ipvals) {
    fprintf(stderr,"ip2 table full (ip_nb=%d)\n",ip_nb[1]);
    return(-1);
    }

  mode = 3;
  if (kind < 4)
    ConvipPlus(&ip_old,&level,&kind,&mode);
  else
    ip_old = -9999;     /* no valid value for oldtype */
  ips_tab[1][ip_nb[1]] = ip_old;
  ip_nb[1]++;

  if (ip_nb[1] > Max_Ipvals) {
    fprintf(stderr,"ip2 table full (i2_ind=%d)\n",ip_nb[1]);
    return(-1);
    }
  return(ip_new);
}


/*splitpoint c_ip3_all */
/*****************************************************************************
 *                          C _ I P 3 _ A L L                                *
 *                                                                           *
 *Object                                                                     *
 *   Generates all possible coded ip3 values for a given ip3                 *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  level          ip3  (float value)                                    *
 *  IN  kind           level kind as defined in convip                       *
 *                                                                           *
 *****************************************************************************/

int c_ip3_all(float level, int kind)
{
  int ip_old, ip_new, mode, flag;
  char s[128];

  ip3s_flag = 1;
  flag = 0;

  mode = 2;
  ConvipPlus(&ip_new,&level,&kind,&mode);
  ips_tab[2][ip_nb[2]] = ip_new;
  ip_nb[2]++;
  if (ip_nb[2] >= Max_Ipvals) {
    fprintf(stderr,"ip3 table full (ip_nb=%d)\n",ip_nb[2]);
    return(-1);
    }

  mode = 3;
  if (kind < 4)
    ConvipPlus(&ip_old,&level,&kind,&mode);
  else
    ip_old = -9999;     /* no valid value for oldtype */
  ips_tab[2][ip_nb[2]] = ip_old;
  ip_nb[2]++;

  if (ip_nb[2] > Max_Ipvals) {
    fprintf(stderr,"ip3 table full (i3_ind=%d)\n",ip_nb[2]);
    return(-1);
    }
  return(ip_new);
}

/*splitpoint c_ip1_val */
/*****************************************************************************
 *                          C _ I P 1 _ V A L                                *
 *                                                                           *
 *Object                                                                     *
 *   Generates coded ip1 value for a given level (shorthand for convip)      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  level          ip1 level (float value)                               *
 *  IN  kind           level kind as defined in convip                       *
 *                                                                           *
 *****************************************************************************/

int c_ip1_val(float level, int kind)
{
  int ip_new, mode, flag;
  char s[128];

  ip1s_flag = 1;
  flag = 0;

  mode = 2;
  ConvipPlus(&ip_new,&level,&kind,&mode);
  ips_tab[0][ip_nb[0]] = ip_new;
  ip_nb[0]++;
  if (ip_nb[0] >= Max_Ipvals) {
    fprintf(stderr,"ip1 table full (ip_nb=%d)\n",ip_nb[0]);
    return(-1);
    }
  return(ip_new);
}

/*splitpoint c_ip2_val */
/*****************************************************************************
 *                          C _ I P 2 _ V A L                                *
 *                                                                           *
 *Object                                                                     *
 *   Generates coded ip2 value for a given level (shorthand for convip)      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  level          ip2  (float value)                                    *
 *  IN  kind           kind as defined in convip                             *
 *                                                                           *
 *****************************************************************************/

int c_ip2_val(float level, int kind)
{
  int ip_new, mode, flag;
  char s[128];

  ip2s_flag = 1;
  flag = 0;

  mode = 2;
  ConvipPlus(&ip_new,&level,&kind,&mode);
  ips_tab[1][ip_nb[1]] = ip_new;
  ip_nb[1]++;
  if (ip_nb[1] >= Max_Ipvals) {
    fprintf(stderr,"ip2 table full (ip_nb=%d)\n",ip_nb[1]);
    return(-1);
    }
  return(ip_new);
}

/*splitpoint c_ip3_val */
/*****************************************************************************
 *                          C _ I P 3 _ V A L                                *
 *                                                                           *
 *Object                                                                     *
 *   Generates coded ip3 value for a ip3 (shorthand for convip)              *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  level          ip3  (float value)                                    *
 *  IN  kind           kind as defined in convip                             *
 *                                                                           *
 *****************************************************************************/

int c_ip3_val(float level, int kind)
{
  int ip_new, mode, flag;
  char s[128];

  ip3s_flag = 1;
  flag = 0;

  mode = 2;
  ConvipPlus(&ip_new,&level,&kind,&mode);
  ips_tab[2][ip_nb[2]] = ip_new;
  ip_nb[2]++;
  if (ip_nb[2] >= Max_Ipvals) {
    fprintf(stderr,"ip3 table full (ip_nb=%d)\n",ip_nb[2]);
    return(-1);
    }
  return(ip_new);
}


/*splitpoint crack_std_parms */
/*****************************************************************************
 *                      C R A C K _ S T D _ P A R M S                        *
 *                                                                           *
 *Object                                                                     *
 *   Reassembles the splited variables.                                      *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  stdf_entry     directory entry that contains the parameters          *
 * OUT  cracked_parms  reassembled parameters                                *
 *                                                                           *
 *****************************************************************************/

static void crack_std_parms(stdf_dir_keys *stdf_entry,
                    stdf_special_parms *cracked_parms)
{
  int i, run;
  unsigned int datexx;
  double r8_diff;
  ftnword ftn_date;
  int diff;

  for (i=0; i <= 4; i++)
    cracked_parms->etiket[i] = ((stdf_entry->etik15 >> ((4-i)*6)) & 0x3f) + 32;

  for (i=5; i <=9; i++)
    cracked_parms->etiket[i] = ((stdf_entry->etik6a >> ((9-i)*6)) & 0x3f) + 32;

  cracked_parms->etiket[10] = ((stdf_entry->etikbc >> 6) & 0x3f) + 32;
  cracked_parms->etiket[11] = (stdf_entry->etikbc  & 0x3f) + 32;
  cracked_parms->etiket[12] = '\0';

  for (i=0; i <=3; i++)
    cracked_parms->nomvar[i] = ((stdf_entry->nomvar >> ((3-i)*6)) & 0x3f) + 32;

  cracked_parms->nomvar[4] = '\0';

  cracked_parms->typvar[0] = ((stdf_entry->typvar >> 6) & 0x3f) + 32;
  cracked_parms->typvar[1] = ((stdf_entry->typvar  & 0x3f)) + 32;
  cracked_parms->typvar[2] = '\0';

  cracked_parms->gtyp[0] = stdf_entry->gtyp;
  cracked_parms->gtyp[1] = '\0';

  cracked_parms->ig2 = (stdf_entry->ig2a << 16) | (stdf_entry->ig2b << 8) |
                       stdf_entry->ig2c;
  /*
   * de-octalise the date_stamp and convert from valid date to origin date
   */
  run = stdf_entry->date_stamp & 0x7;
  datexx = (stdf_entry->date_stamp >> 3) * 10 + run;
  r8_diff = -((((double)stdf_entry->deet) * ((double)stdf_entry->npas))/3600.0);
  diff = -(((stdf_entry->deet) * (stdf_entry->npas) + 1800)/3600);
  cracked_parms->date_valid = datexx;
  cracked_parms->date_stamp = datexx;
  ftn_date = (ftnword) datexx;
  if ((stdf_entry->deet * stdf_entry->npas) != 0) {           /* compute origin date */
    f77name(incdatr)(&ftn_date,&ftn_date,&r8_diff);
    cracked_parms->date_stamp = (int) ftn_date;
  }
/*    f77name(incdat)(&cracked_parms->date_stamp,&datexx,&diff); */
  cracked_parms->aammjj = 0;
  cracked_parms->hhmmss = 0;
}


/*splitpoint init_ip_vals */
/*****************************************************************************
 *                             I N I T _ I P _ V A L S                       *
 *                                                                           *
 *Object                                                                     *
 *   Reset to zeros ip1-2-3 tables and counters              l               *
 *                                                                           *
 *****************************************************************************/

int init_ip_vals()
{
  int i,j;

  for (i=0; i<Max_Ipvals; i++)
    for (j=0; j<3; j++)
      ips_tab[j][i] = -1;

  for (j=0; j<3; j++)
    ip_nb[j] = 0;

  ip1s_flag = 0;
  ip2s_flag = 0;
  ip3s_flag = 0;
  return(0);
}


/*splitpoint ip_is_equal */
/*****************************************************************************
 *                          I P _ I S _ E Q U A L                            *
 *                                                                           *
 *Object                                                                     *
 *   Compares different coded values of an ip for equality                   *
 *                                                                           *
 *Arguments                                                                  *
 *                                                                           *
 *  IN  target         must be first value in the table of coded value to    *
 *                     compare with                                          *
 *  IN  ip             current ip record value to compare                    *
 *  IN  ind            index (1,2 or 3) representing ip1, ip2 or ip3         *
 *                     comparisons                                           *
 *                                                                           *
 *****************************************************************************/

int ip_is_equal(int target, int ip, int ind)
{
  int kind1, kind2, exp1, exp2, nb, j;
  long long mantis1, mantis2;
//  double level, target_level;
//  double tolr = .0000001;
/*
  double exptab[] = {0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0,
                     1000.0, 10000.0, 100000.0, 1000000.0, 10000000.0,
                     100000000.0, 1000000000.0, 10000000000.0, 100000000000.0};
*/
  /* ipold_0_9: table of old ip1 0-9 (mb) values encoded oldstyle with convip */
  int ipold_0_9[10] = {0,1820,1840,1860,1880,1900,1920,1940,1960,1980};

  ind--;       /* ind is passed in base 1, for ip1 ip2 or ip3 */
  if (target != ips_tab[ind][0]) {
    fprintf(stderr,"ip_is_equal error: target not first element of ip table\n");
    fprintf(stderr,"target = %d ips_tab[%d]=%d\n",target,ind,ips_tab[ind][0]);
    return(0);
    }

  nb = ip_nb[ind];            /* number of elements in ip[1-2-3] table */
  for (j=0; j<nb; j++) {
    target = ips_tab[ind][j] ;
    if (ip == target) return(1);
    if (ip < 10)
      if (ipold_0_9[ip] == target) return(1);
    if ((ip > 32767) && (target > 32767)) {
      kind1 = (ip >> 24) & 0xF;
      kind2 = (target >> 24) & 0xF;
      if (kind1 != kind2) continue;
      exp1 = (ip >> 20) & 0xF;
      exp2 = (target >> 20) & 0xF;
      mantis1 = ip & 0xFFFFF;
      if (mantis1 > 1000000) mantis1 = -(mantis1 - 1000000);
      mantis2 = target & 0xFFFFF;
      if (mantis2 > 1000000) mantis2 = -(mantis2 - 1000000);
      if ((mantis1 ^ mantis2) < 0) continue;                /* mantis1 and mantis2 must be same signed */
/*
      level = mantis1 / exptab[exp1];
      target_level = mantis2 / exptab[exp2];
      printf("Debug ip_is_equal level=%f\n           target_level=%f\n",level,target_level);
      printf("Debug avant kind1=%d kind2=%d exp1=%d exp2=%d\n",kind1,kind2,exp1,exp2);
      printf("Debug avant \n     mantis1=%Ld\n     mantis2=%Ld\n",mantis1,mantis2);
*/
      if (exp1 > exp2)
        while (exp1 > exp2) {
          exp2++;
          mantis2 *= 10;
          }
      else
        while (exp2 > exp1) {
          exp1++;
          mantis1 *= 10;
          }
/*
      printf("Debug kind1=%d kind2=%d exp1=%d exp2=%d\n",kind1,kind2,exp1,exp2);
      printf("Debug\n     mantis1=%Ld\n     mantis2=%Ld\n",mantis1,mantis2);
*/
      if (labs(mantis1-mantis2) <= 1) return(1);
      }
  }
  return(0);
}

#include "if_fstd98.hc"
