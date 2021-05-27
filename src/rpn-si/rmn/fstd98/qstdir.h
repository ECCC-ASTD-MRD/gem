//! @file qstdir.h XDF constants and data structures

#include <stdint.h>
#include <rpnmacros.h>
#include "../INTRALIB_INCLUDES/fnom.h"
#define WRITE_PAGE(a,b,c)

#define MODE3232

/*  DEFINITIONS FOR 32 BIT MACHINES IN 32 BIT MODE */
#ifdef MODE3232

//! ftnword to 64 bit word conversion
#define WDINT64(nwds) (nwds>>1)
//! 64 bit word to ftnword conversion
#define WD64INT(nw64) (nw64<<1)
//! maximum length of primary keys
#define MAX_PRIMARY_LNG 32
//! Maximum length of info keys
#define MAX_SECONDARY_LNG 16
//! Unsigned 32 bit word
typedef unsigned INT_32 word32;
//! Word to 64 bit word conversion
#define WDTO64(nwds) (nwds>>1)
//! 64 bit word to word conversion
#define W64TOWD(nw64) (nw64<<1)
//! 64 bit word to 32 bit word conversion
#define W64TOwd(nw64) (nw64<<1)
//! Fortran word to 32 bit word conversion
#define FWTOwd(nwds) (nwds)

#endif

/*  DEFINITIONS FOR 64 BIT MACHINES IN 64 BIT MODE */
#ifdef MODE6464

//! ftnword to 64 bit word conversion
#define WDINT64(nwds) (nwds)
//! 64 bit word to ftnword conversion
#define WD64INT(nw64) (nw64)
//! Maximum length of primary keys
#define MAX_PRIMARY_LNG 16
//! Maximum length of info keys
#define MAX_SECONDARY_LNG 8
//! Unsigned 32 bit word
typedef unsigned INT_32 word32;
//! Word to 64 bit word conversion
#define WDTO64(nwds) (nwds)
//! 64 bit word to word conversion
#define W64TOWD(nw64) (nw64<<1)
//! 64 bit word to 32 bit word conversion
#define W64TOwd(nw64) (nw64<<1)
//! Fortran word to 32 bit word conversion
#define FWTOwd(nwds) (nwds<<1)

#endif


/* macro definitions */
#define upper_case(car) ( car & (~((car & 64) >> 1)) )
#define upper_case_word(wd) ( wd & (~((wd & 0x40404040) >> 1)) )
#define ascii64(car) ( (car > 95) ? car-32 : car)
#define ascii6(car) ((upper_case(car)-32) & 0x3f)
#define string_copy(dest,src,l) while(--l >= 0) dest[l]=src[l]
#define str_cp_init(strc,lc,strf,lf)                               \
   {                                                               \
     int i;                                                        \
     for (i=0; i < lc-1; i++) strc[i] = (i < lf) ? strf[i] : ' ';  \
     strc[lc-1] = '\0';                                            \
   }                                                               \

#define Min(x,y) ((x < y) ? x : y)
#define Max(x,y) ((x > y) ? x : y)
#define MRBCOV(a,b,c) ((a << 14) | (b << 8) | c)
#define INFOPRINT if (1 >= msg_level)
#define WARNPRINT if (2 >= msg_level)
/* end macro definitions */

#define MAX_XDF_FILES 1024
//! Size of a directory pages in an XDF file
#define ENTRIES_PER_PAGE 256
//! Maximum of 256K records in a random access XDF file
#define MAX_DIR_PAGES 1024
#define MAX_RECORD_LENGTH 33554400
/* Maximum of 128MB per record */
/* Maximum of 8G per file (wa) */

//! Maximum number of statistics for xdfsta
#define MAX_STAT 12
//! Maximum number of primary or info keys
#define MAX_KEYS 100

/* Error levels */
//! Trivial error
#define TRIVIAL 0
//! Informative message
#define INFORM 1
//! Warning message
#define WARNING 2
//! Important error
#define ERROR 3
//! Fatal error
#define ERRFATAL 4
//! Internal software error
#define SYSTEM 5
//! Error not tolerable, automatic abort
#define KAPUT 6

/* Error codes */
/* (codes -1 to -23 taken from previous xdf version when appliable) */
#define ERR_OK 0
//! File is not open or does not exist
#define ERR_NO_FILE -1
//! File table is full
#define ERR_FTAB_FULL -3
//! Short read, truncated record
#define ERR_SHORT_READ -4
//! Invalid unit number
#define ERR_BAD_UNIT -5
//! Src and dest files not compatible
#define ERR_NOT_COMP -6
//! No write permission
#define ERR_NO_WRITE -7
//! Invalid page number
#define ERR_BAD_PAGENO -8
//! Invalid handle
#define ERR_BAD_HNDL -9
//! Special record with idtyp = 0
#define ERR_SPECIAL -10
//! Deleted record
#define ERR_DELETED -11
//! Search target not found
#define ERR_NOT_FOUND -12
//! Error in record initialisation
#define ERR_BAD_INIT -13
//! Invalid datyp
#define ERR_BAD_DATYP -16
//! Addressing error (not a multiple of 64)
#define ERR_BAD_ADDR -18
//! Dimension of buf too small
#define ERR_BAD_DIM -19
//! Invalid option name or value
#define ERR_BAD_OPT -20
//! File already in used in write mode
#define ERR_STILL_OPN -21
//! Read only file
#define ERR_RDONLY -22
//! Invalid header length
#define ERR_BAD_LEN -23

//! No more memory can be allocated
#define ERR_MEM_FULL -24
//! No valid file position
#define ERR_NO_POS -25
//! Bad linked file reference
#define ERR_BAD_LINK -26
//! File is already open
#define ERR_FILE_OPN -27
//! Directory of random file is full
#define ERR_DIR_FULL -28
//! File not connected with fnom
#define ERR_NO_FNOM -29
//! File is not linked
#define ERR_NO_LINK -30
//! No valid search target
#define ERR_NO_TARGET -31
//! Bad checksum
#define ERR_BAD_CHKS -32
//! Incorect number of directory pages
#define ERR_BAD_DIR -33
//! Not an XDF file
#define ERR_NOT_XDF -34
//! Incorect number of stats
#define ERR_BAD_NSTAT -35
//! Value out of valid range
#define ERR_OUT_RANGE -36
//! Incorect file type
#define ERR_BAD_FTYPE -37
//! Routine not implemented in FSTD97
#define ERR_NOT_IMPL -38
//! Invalid or no sequential record postfix
#define ERR_NO_POSTFIX -39
//! Wrong version of standard file software
#define ERR_STDF_VERSION -40
//! Wrong file type, BURP instead of STDF or vice versa
#define ERR_WRONG_FTYPE -41
//! File probably damaged
#define ERR_DAMAGED -45

//! Too many supplementary keys
#define BURP_ERR_CLEF -32
//! Incorect block number
#define BURP_ERR_BNUM -33
//! Value out of range for 32bits and datyp = 2
#define BURP_ERR_CMPR -41
//! Incorect bdesc
#define BURP_ERR_BDESC -43
//! Incorect element code for datyp 7,8,9
#define BURP_ERR_CODE -44

/* File types */
#define STDF89_RND 20
#define STDF89_SEQ 21
#define XDF1_RND 30
#define XDF1_SEQ 31
//! End of file marker for sequential xdf
#define MARK_EOI 0xf0f0f0f0
/* EOI marker for old sequential standard */
#define MARK_EOI_89 31
//! Random standard file signature
#define STDF_RND_SIGN 0x55555555
//! Sequential standard file signature
#define STDF_SEQ_SIGN 0xaaaaaaaa

//! Read only
#define RDMODE 0
//! Write only
#define WMODE  1
//! Create
#define CREATE 2
//! Read and write
#define RWMODE 3
//! Append
#define APPEND 4

//! Starting position of record in buf->data
#define RECADDR 0

/* Position of bit1 */
#define BPBIT1 44
/* Position of key length */
#define BPLCLE 49
/* Position of key type  */
#define BPTCLE 55
/* Number of bits for bit1 */
#define LBIT1 13
//! Number of bits for key length
#define LLCLE 5
//! Number of bits for key type
#define LTCLE 6

/* burp defines */

/*! Burp dim header (32bit word)*/
#define DIMENT 4
//! Number of bits for burp report header
#define NBENTR 320
//! Number of bits for burp block header
#define NBENTB 128
//! Maximum number for nele with 7 bits
#define GROSNELE 128
//! Maximum number for nval and nt (8 bits)
#define GROSDIM 256
//! Maximum number of primary keys
#define NPRIDEF 18
//! Maximum number of auxiliary keys
#define NAUXDEF 5
//! Maximum number of supplementary prim keys
#define NPRISUP 0
//! Maximum number of supplementary aux  keys
#define NAUXSUP 0

//! Fix 64bit fortran buffer for c 32bit functions
#define BUF_C buf[1] = 2 * buf[1] -1
//! Restore to original buffer
#define BUF_F buf[1] = (buf[1] +1) >> 1

#define IUN_INDEX(ifile,unit) { \
   ifile = MAX_XDF_FILES; \
   while (--ifile >= 0) {\
     if (file_table[ifile] != NULL) \
       if (file_table[ifile].iun == unit) break; \
     } \
   }

#define FREE_INDEX(ifile,unit) { \
   ifile = MAX_XDF_FILES; \
   while (--ifile >= 0) {\
     if (file_table[ifile] == NULL) break;\
     } \
   }

#define VALID(val, minval, maxval, what, caller) \
   if ((val < minval) || (val > maxval)) { \
     sprintf(errmsg, "%s = %d must be between %d and %d", what, val, minval, maxval); \
     return(error_msg(caller, ERR_OUT_RANGE, ERROR));\
   }

#define MOREFILES

#ifdef MOREFILES
/* how to make ordinary handles for random and sequential files */
/* handle is limited to 128 MB for a sequential file */
/*       HANDLE RANDOM */
/*       sign #page #record index */
/*         1    12     9     10   */
#define MAKE_RND_HANDLE(pageno,recno,file_index) ((file_index &0x3FF) | ((recno & 0x1FF)<<10) | ((pageno & 0xFFF)<<19))
/*       HANDLE SEQ */
/*       sign cluster addr index */
/*         1    2      22    7   */
#define MAKE_SEQ_HANDLE(cluster,address,file_index) (file_index | ((address & 0x3FFFFF) << 7) | (cluster << 29))
/* how to extract the file index, record number and page number from a handle */
#define INDEX_FROM_HANDLE(handle) ((STDSEQ_opened==1) ? ( 0x7F & handle) : ( 0x3FF & handle))
/* #define INDEX_FROM_HANDLE(handle) ( 0x3FF & handle) */
#define RECORD_FROM_HANDLE(handle) ( 0x1FF & (handle>>10))
#define PAGENO_FROM_HANDLE(handle) ( 0xFFF & (handle>>19))
/* how to extract the record address from a sequential handle */
#define ADDRESS_FROM_HNDL(handle) ( 0x3FFFFF & (handle>>7))
#define CLUSTER_FROM_HANDLE(handle) ( 0x3 & (handle>>29))
#define ADDRESS_FROM_HANDLE(handle) ( ADDRESS_FROM_HNDL(handle) << (2 *CLUSTER_FROM_HANDLE(handle)) )
#else
/* how to make ordinary handles for random and sequential files */
/* handle is limited to 128 MB for a sequential file */
/*       HANDLE RANDOM */
/*       sign #page #record index */
/*         1    12     12     7   */
#define MAKE_RND_HANDLE(pageno,recno,file_index) (file_index | (recno<<7) | (pageno<<19))
/*       HANDLE SEQ */
/*       sign cluster addr index */
/*         1    2      22    7   */
#define MAKE_SEQ_HANDLE(cluster,address,file_index) (file_index | ((address & 0x3FFFFF) << 7) | (cluster << 29))
/* how to extract the file index, record number and page number from a handle */
#define INDEX_FROM_HANDLE(handle) ( 0x7F & handle)
#define RECORD_FROM_HANDLE(handle) ( 0xFFF & (handle>>7))
#define PAGENO_FROM_HANDLE(handle) ( 0xFFF & (handle>>19))
/* how to extract the record address from a sequential handle */
#define ADDRESS_FROM_HNDL(handle) ( 0x3FFFFF & (handle>>7))
#define CLUSTER_FROM_HANDLE(handle) ( 0x3 & (handle>>29))
#define ADDRESS_FROM_HANDLE(handle) ( ADDRESS_FROM_HNDL(handle) << (2 *CLUSTER_FROM_HANDLE(handle)) )
#endif


//! Maximum size of the primary key entry
typedef word max_dir_keys[MAX_PRIMARY_LNG];
typedef word max_info_keys[MAX_SECONDARY_LNG];

//! Directory page record
typedef struct {
#if !defined(Little_Endian)
        //! Type ID (usualy 0)
        word idtyp:8;
        //! Header length (in 64 bit units)
        word lng:24;
        //! Address of directory page (origin 1, 64 bit units)
        addr:32;
#else
        //! Header length (in 64 bit units)
        word lng:24;
        //! Type ID (usualy 0)
        word idtyp:8;
        //! Address of directory page (origin 1, 64 bit units)
        word addr:32;
#endif

        //! idrep (4 ascii char 'DIR0')
        word reserved1:32;
        //! Reserved (0)
        word reserved2:32;

        //! Address of next directory page (origin 1, 64 bit units)
        word nxt_addr:32;
        //! Number of entries in page
        word nent:32;

        //! Checksum (not valid when in core)
        word chksum:32;
        //! page_no, record_no, file_index: handle templage
        word reserved3:32;

        //! (real allocated dimension will be ENTRIES_PER_PAGE * primary_len)
        word entry[];
} xdf_dir_page;



//! Zero size directory page record
typedef struct {
        /* each line describes 64 bits */
#if !defined(Little_Endian)
        word idtyp:8, lng:24,  addr:32;
        word reserved1:32, reserved2:32;
        word nxt_addr:32,  nent:32;
        word chksum:32, page_no:16, record_no:8, file_index:8;
#else
        word lng:24, idtyp:8, addr:32;
        word reserved1:32, reserved2:32;
        word nxt_addr:32,  nent:32;
        word chksum:32, file_index:8, record_no:8, page_no:16;
#endif
} base_dir_page;



/* Directory page + forward/backward chain ptrs */
typedef struct full_dir_page {
        struct full_dir_page* next_page;
        struct full_dir_page* prev_page;
        int modified;
        int true_file_index;
        xdf_dir_page dir;
} full_dir_page;

//! Pointer to a chainable directory page
typedef full_dir_page* page_ptr;


//! XDF record header
typedef struct {
#if !defined(Little_Endian)
        word idtyp:8, lng:24, addr:32;
#else
        word lng:24, idtyp:8, addr:32;
#endif
} xdf_record_header;


//! XDF data record
typedef struct {
#if !defined(Little_Endian)
        word idtyp:8, lng:24,  addr:32;        /* XDF record header */
#else
        word lng:24, idtyp:8,  addr:32; /* XDF record header */
#endif
        word data[2];                        /* primary keys, info keys, data */
} file_record;


typedef struct {
#if !defined(Little_Endian)
        word page_no:16, record_no:8, file_index:8;
#else
        word file_index:8, record_no:8, page_no:16;
#endif
} random_record_handle;


typedef struct {
#if !defined(Little_Endian)
        word address:24, file_index:8;
#else
        word file_index:8, address:24;
#endif
} seq_record_handle;


/*****************************************************************************/
/*                        format of standard files                           */
/*****************************************************************************/

/* define maximum values for stdf descriptors */
#define DEET_MAX    0xFFFFFF
#define NBITS_MAX   64
#define NI_MAX      0xFFFFFF
#define NJ_MAX      0xFFFFFF
#define NK_MAX      0xFFFFF
#define NPAS_MAX    0xFFFFFF
#define IG1_MAX     0xFFFFFF
#define IG2_MAX     0xFFFFFF
#define IG3_MAX     0xFFFFFF
#define IG4_MAX     0xFFFFFF
#define IP1_MAX     0xFFFFFFF
#define IP2_MAX     0xFFFFFFF
#define IP3_MAX     0xFFFFFFF

//! Search tags part of standard file directory entry :  header + 8 x 64 bits
typedef struct {
        /* there is one 64 bit group per line */
#if !defined(Little_Endian)
        word deleted:1, select:7, lng:24, addr:32;
        word deet:24, nbits: 8, ni:   24, gtyp:  8;
        word nj:24,  datyp: 8, nk:   20, ubc:  12;
        word npas: 26, pad7: 6, ig4: 24, ig2a:  8;
        word ig1:  24, ig2b:  8, ig3:  24, ig2c:  8;
        word etik15:30, pad1:2, etik6a:30, pad2:2;
        word etikbc:12, typvar:12, pad3:8, nomvar:24, pad4:8;
        word ip1:28, levtyp:4, ip2:28, pad5:4;
        word ip3:28, pad6:4, date_stamp:32;
#else
        word lng:24, select:7, deleted:1, addr:32;
        word nbits: 8, deet:24, gtyp:  8, ni:   24;
        word datyp: 8, nj:24, ubc:  12, nk:   20;
        word pad7: 6, npas: 26, ig2a:  8, ig4: 24;
        word ig2b:  8, ig1:  24, ig2c:  8, ig3:  24;
        word pad1:2, etik15:30, pad2:2, etik6a:30;
        word pad3:8, typvar:12, etikbc:12, pad4:8, nomvar:24;
        word levtyp:4, ip1:28, pad5:4, ip2:28;
        word pad6:4, ip3:28, date_stamp:32;
#endif
} stdf_dir_keys;

//! Information part of a standard file record header, 1 x 64 bits
typedef struct {
        /* there is one 64 bit group per line */
#if !defined(Little_Endian)
        word nblks: 8, blk1: 24, blk2: 32;
#else
        word blk1: 24, nblks: 8, blk2: 32;
#endif
        /* number of blocks, first 2 block pointers */
} stdf_dir_info;

//! Collection area for cracked record parameter values
typedef struct {
        char etiket[13], nomvar[5], typvar[3], gtyp[2], extra;
        word lng, addr;
        INT_32 aammjj, hhmmss;
        INT_32 ni, nj, nk, nbits, datyp, deet, npas, ip1, ip2, ip3, ig1, ig2, ig3, ig4;

} stdf_rec_parms;

/*! Collection area for some special cracked record parameters that need to
 * be reassembled */
typedef struct {
        char etiket[13], nomvar[5], typvar[3], gtyp[2];
        INT_32 date_stamp, aammjj, hhmmss;
        INT_32 ig2, date_valid;

} stdf_special_parms;

//! Collection area for parameter addresses/values
typedef struct {
        char *etiket, *nomvar, *typvar, *gtyp, *extra;
        ftnword *date;
        int l_etiket, l_nomvar, l_typvar, l_extra;
        INT_32 lng, addr, ni, nj, nk, nbits, datyp, deet, npas, ip1, ip2, ip3, ig1, ig2, ig3, ig4;

} stdf_adr_parms;

typedef struct {
#if !defined(Little_Endian)
   word swa : 32, npas1 : 16, nk : 12, epce1 : 4;
   word ni : 16, nj : 16, nomvar : 16, typvar : 8, nbits : 8;
   word ip1 :16, ip2 : 16, ip3 : 16, epce2 : 7, dltf : 1, npas2 : 8;
   word etiq14 :32, etiq56 : 16, etiq78 : 16;
   word epce3 : 32, epce4 : 16, ig2 : 16;
   word ig3 : 16, ig4 : 16, grtyp : 8, datyp : 8, ig1 : 16;
   word date : 32, ubc : 16, deet : 16;
   word lng : 32, eof : 32;
#else
   word swa : 32, epce1 : 4, nk : 12, npas1 : 16;
   word nj : 16, ni : 16, nbits : 8, typvar : 8, nomvar : 16;
   word ip2 :16, ip1 : 16, npas2 : 8, dltf : 1, epce2 : 7, ip3 : 16;
   word etiq14 :32, etiq78 : 16, etiq56 : 16;
   word epce3 : 32, ig2 : 16, epce4 : 16;
   word ig4 : 16, ig3 : 16, ig1 : 16, datyp : 8, grtyp : 8;
   word date : 32, deet : 16, ubc : 16;
   word lng : 32, eof : 32;
#endif
   /****Format de bits 512-959 SEQ/SQI****/
   /* Le eof defini precedemant appartient a cet interval de nombre de bits */

   word vide1 : 32, swa_last : 32;
   word vide3 : 32, vide4 : 32;
   word epce5 : 32, epce6 : 32;
   word epce7 : 32, epce8 : 32;
   word epce9 : 32, epce10: 32;
   word epce11: 32, epce12: 32;
   word vide5 : 32, vide6 : 32;
} seq_dir_keys;


typedef struct {
   word etiqt1  : 32, etiqt2  : 32;
   word dirsiz  : 32, inuti1  : 32;
   word nutil   : 32, inuti2  : 32;
   word nbecr   : 32, inuti3  : 32;
   word nbrec   : 32, inuti4  : 32;
   word nbext   : 32, inuti5  : 32;
   word nrecup  : 32, inuti6  : 32;
   word nbeff   : 32, inuti7  : 32;
   word nbcorr  : 32, inuti8  : 32;
   word inuti9  : 32, inuti10 : 32;
   word inuti11 : 32, inuti12 : 32;
   word inuti13 : 32, inuti14 : 32;
   word inuti15 : 32, inuti16 : 32;
   word inuti17 : 32, inuti18 : 32;
   word inuti19 : 32, inuti20 : 32;
} stdf_struct_RND;

typedef struct {
#if !defined(Little_Endian)
   word swa : 32, npas1 : 16, nk : 12, epce1 : 4;
   word ni : 16, nj : 16, nomvar : 16, typvar : 8, nbits : 8;
   word ip1 : 16, ip2 : 16, ip3 : 16, epce2 : 7, dltf : 1, npas2 : 8;
   word etiq14 : 32, etiq56 : 16, etiq78 : 16;
   word epce3 : 32, epce4 : 16, ig2 : 16;
   word ig3 : 16, ig4 : 16, grtyp : 8, datyp : 8, ig1 : 16;
   word date : 32, ubc : 16, deet : 16;
   word lng : 32;
#else
   word swa : 32, epce1 : 4, nk : 12, npas1 : 16;
   word nj : 16, ni : 16, nbits : 8, typvar : 8, nomvar : 16;
   word ip2 :16, ip1 : 16, npas2 : 8, dltf : 1, epce2 : 7, ip3 : 16;
   word etiq14 :32, etiq78 : 16, etiq56 : 16;
   word epce3 : 32, ig2 : 16, epce4 : 16;
   word ig4 : 16, ig3 : 16, ig1 : 16, datyp : 8, grtyp : 8;
   word date : 32, deet : 16, ubc : 16;
   word lng : 32;
#endif
} rnd_dir_keys;

typedef struct {
#if !defined(Little_Endian)
   word idtyp : 8, lng : 24, addr : 32;
   word prev_idtyp : 8, prev_lng : 24, prev_addr :32;
#else
   word lng : 24, idtyp : 8, addr : 32;
   word prev_lng : 24, prev_idtyp : 8, prev_addr :32;
#endif
} postfix_seq;


/*****************************************************************************/
/*                           format of BURP files                            */
/*****************************************************************************/


//! Search tags part of burp file directory entry :  header + 3 x 64 bits
typedef struct {
        /* there is one 64 bit group per line */
#if !defined(Little_Endian)
        word idtyp:8,  lng:24,   addr:32;  /* standard XDF record header */
        word sti1:8, sti2:8, sti3:8, sti4:8, sti5:8, sti6:8, sti7:8, sti8:8;
        word sti9:8, flgs:24, lati:16, lon:16;
        word date:20, dx:12, idtp:8, dy:12, heur:6, min:6;
#else
        word lng:24, idtyp:8,   addr:32;  /* standard XDF record header */
        word sti4:8, sti3:8, sti2:8, sti1:8, sti8:8, sti7:8, sti6:8, sti5:8;
        word flgs:24, sti9:8, lon:16, lati:16;
        word dx:12, date:20, min:6, heur:6, dy:12, idtp:8;
#endif
} burp_dir_keys;


//! Information part of a burp file record header, 1 x 64 bits
typedef struct {
        /* there is one 64 bit group per line */
#if !defined(Little_Endian)
        word nblks:16, oars:16, elev:13 ,drcv:11, runn:8;
#else
        word oars:16, nblks:16, runn:8, drcv:11, elev:13;
#endif
} burp_dir_info;


//! Burp file record header + data
typedef struct {
        burp_dir_keys keys;
        burp_dir_info info;
        word data[];
} burp_record;


//! Burp block header
typedef struct {
        /* there is one 64 bit group per line */
#if !defined(Little_Endian)
        word bfamdesc:12, btyp:15, nbit:5, nt:8, datyp:4, bit0:20;
        word flag:1, nele:7, nval:8, elem1:16, elem2:16, elem3:16;
#else
        word nbit:5, btyp:15, bfamdesc:12, bit0:20, datyp:4, nt:8;
        word elem1:16, nval:8, nele:7, flag:1, elem3:16, elem2:16;
#endif
} burp_block_header;

/*****************************************************************************/
/*                      preamble description for XDF files                   */
/*****************************************************************************/

//! Key descriptor structure, 64 bits per key description
typedef struct {
#if !defined(Little_Endian)
        uint32_t ncle:32, bit1:13, lcle:5, tcle:6, reserved:8;
#else
        uint32_t ncle:32, reserved:8, tcle:6, lcle:5, bit1:13;
#endif
} key_descriptor;


//! 2xN array interface from FORTRAN for keys
typedef struct {
        ftnword wd1;
        ftnword wd2;
} ftnword_2;

typedef struct {
        word wd1;
        word wd2;
} word_2;


//! XDF file header template, provision is made for up to 1024 keys
typedef struct {
        /* each line (except last one) describes 64 bits */
#if !defined(Little_Endian)
        uint32_t idtyp:8,  lng:24,   addr:32;  /* standard XDF record header */
        uint32_t vrsn,     sign;               /* char[4] */
        uint32_t fsiz:32,  nrwr:32;
        uint32_t nxtn:32,  nbd:32;
        uint32_t plst:32,  nbig:32;
        uint32_t nprm:16,  lprm:16,  naux:16, laux:16;
#else
        uint32_t lng:24,   idtyp:8,   addr:32; /* standard XDF record header */
        uint32_t vrsn,     sign;               /* char[4] */
        uint32_t fsiz:32,  nrwr:32;
        uint32_t nxtn:32,  nbd:32;
        uint32_t plst:32,  nbig:32;
        uint32_t lprm:16,  nprm:16,  laux:16, naux:16;
#endif
        uint32_t neff:32,  nrec:32;
        uint32_t rwflg:32, reserved:32;
        key_descriptor keys[1024];
/*
 * idtyp:     id type (usualy 0)
 * lng:       header length (in 64 bit units)
 * addr:      address (exception: 0 for a file header)
 * vrsn:      XDF version
 * sign:      application signature
 * fsiz:      file size (in 64 bit units)
 * nrwr:      number of rewrites
 * nxtn:      number of extensions
 * nbd:       number of directory pages
 * plst:      address of last directory page (origin 1, 64 bit units)
 * nbig:      size of biggest record
 * nprm:      number of primary keys
 * lprm:      length of primary keys (in 64 bit units)
 * naux:      number of auxiliary keys
 * laux:      length of auxiliary keys
 * neff:      number of erasures
 * nrec:      number of valid records
 * rwflg:     read/write flag
 * reserved:  reserved
 * keys:      key descriptor table
 */
} file_header;


/*****************************************************************************/
/*                      description of a file table entry                    */
/*****************************************************************************/
//! Pointer to a function returning a word
typedef word* fn_ptr();
//! Pointer to a primary keys building function
typedef word* fn_b_p(word* buf, word* keys, word* mask,
                     word* mskkeys, int index, int mode);
//! Pointer to a info keys building function
typedef word* fn_b_i(word* buf, word* keys, int index, int mode);

typedef struct {
        //! Pointer to directory pages
        page_ptr dir_page[MAX_DIR_PAGES];
        //! Pointer to current directory page
        page_ptr cur_dir_page;
        //! Pointer to primary key building function
        fn_b_p* build_primary;
        //! Pointer to info building function
        fn_ptr* build_info;
        //! Pointer to file scan function
        fn_ptr* scan_file;
        //! Pointer to record filter function
        fn_ptr* file_filter;
        //! Pointer to current directory entry
        word* cur_entry;
        //! Pointer to file header
        file_header* header;
        //! Next write address (in word units)
        INT_32 nxtadr;
        //! Length in 64 bit units of primary keys (including 64 bit header)
        int primary_len;
        //! Length in 64 bit units of info keys
        int info_len;
        //! File index to next linked file,-1 if none
        int link;
        //! Pointer to current general file desc entry
        general_file_info* cur_info;
        //! FORTRAN unit number, -1 if not open, 0 if C file
        int iun;
        //! Index into file table, -1 if not open
        int file_index;
        //! Modified flag
        int modified;
        //! Number of allocated directory pages
        int npages;
        //! Number of records in file
        int nrecords;
        //! Current page number
        int cur_pageno;
        //! Record number within current page
        int page_record;
        //! Number of records in current page
        int page_nrecords;
        //! Version number
        int file_version;
        //! Last search target valid flag
        int valid_target;
        //! File is sequential xdf
        int xdf_seq;
        //! Last position valid flag (seq only)
        int valid_pos;
        //! Current address (WA, sequential xdf)
        int cur_addr;
        //! Address (WA) of first record (seq xdf)
        int seq_bof;
        //! Old standard file flag
        int fstd_vintage_89;
        //! Header & primary keys for last record
        max_dir_keys head_keys;
        //! Info for last read/written record
        max_info_keys info_keys;
        //! Keys for current operation
        max_dir_keys cur_keys;
        //! Current search target
        max_dir_keys target;
        //! Permanent search mask for this file
        max_dir_keys srch_mask;
        //! Current search mask for this file
        max_dir_keys cur_mask;
} file_table_entry;

typedef file_table_entry* file_table_entry_ptr;


//! XDF buffer interface layout
typedef struct {
        //! Pointer to next data block
        struct data_block* next;
        //! Pointer to data itself
        word* data_ptr;
        //! Data length in word units
        int data_lng;
} data_block;
typedef data_block* data_block_ptr;

#define buf7 dummy[0]
#define buf8 dummy[1]
typedef struct {
        //! Dimension of buffer (in words)
        word nwords;
        //! Number of bits for entire record
        word nbits;
        //! Starting position of data (after the keys)
        word data_index;
        //! Starting position of record
        word record_index;
        //! Associated unit number
        word iun;
        //! Starting position of auxiliary keys
        word aux_index;
        /* What the heck!? */
        /* Why is the member nammed dummy? */
        /* Why is it only accessed with the buf8 macro? */
        /* Why is it in a union? */
        /* What do buf7, buf8, and buf9 mean? */
        union {
                data_block_ptr ptr;
                word dummy[2];
        } buf78;
        word buf9;
        //! Record data
        word data[];
} buffer_interface;
typedef buffer_interface* buffer_interface_ptr;

#if defined(XDF_OWNER)
//! File table, exported symbol
file_table_entry_ptr file_table[MAX_XDF_FILES];

//! Buffer to format error messages
char errmsg[1024];
//! Tolerance level before messages are issued
int msg_level = INFORM;
//! Tolerance level before program is aborted
int xdf_toler = ERROR;
//! Stride
int xdf_stride = 1;
//! Double float indicator
int xdf_double = 0;
//! Short integer indicator
int xdf_short = 0;
//! Byte array indicator
int xdf_byte = 0;
//! Enforce 8 char for date specifications
int xdf_enforc8 = 0;
//! Data type of last record read
int xdf_datatyp;
//! Number of splited output files in xdfuse
int xdf_nsplit = 1;
//! Number of bits per FORTRAN word
int FTN_Bitmot = 8 * bytesperword;
//! No pack/unpack, used by editfst
int image_mode_copy = 0;
//! Chekcpoint mode, no closing of the file
int xdf_checkpoint = 0;
//! If one std seq file is opened, the limit of opened files becomes 128
int STDSEQ_opened = 0;
key_descriptor stdfkeys[] = {
#if !defined(Little_Endian)
  {'SF01', 31,31, 0, 0},
  {'SF02', 63,31, 0, 0},
  {'SF03', 95,31, 0, 0},
  {'SF04',127,31, 0, 0},
  {'SF05',159,31, 0, 0},
  {'SF06',191,31, 0, 0},
  {'SF07',223,31, 0, 0},
  {'SF08',255,31, 0, 0},
  {'SF09',287,31, 0, 0},
  {'SF10',319,31, 0, 0},
  {'SF11',351,31, 0, 0},
  {'SF12',383,31, 0, 0},
  {'SF13',415,31, 0, 0},
  {'SF14',447,31, 0, 0},
  {'SF15',479,31, 0, 0},
  {'SF16',511,31, 0, 0}
#else
  {'SF01', 0, 0, 31,  31},
  {'SF02', 0, 0, 31,  63},
  {'SF03', 0, 0, 31,  95},
  {'SF04', 0, 0, 31, 127},
  {'SF05', 0, 0, 31, 159},
  {'SF06', 0, 0, 31, 191},
  {'SF07', 0, 0, 31, 223},
  {'SF08', 0, 0, 31, 255},
  {'SF09', 0, 0, 31, 287},
  {'SF10', 0, 0, 31, 319},
  {'SF11', 0, 0, 31, 351},
  {'SF12', 0, 0, 31, 383},
  {'SF13', 0, 0, 31, 415},
  {'SF14', 0, 0, 31, 447},
  {'SF15', 0, 0, 31, 479},
  {'SF16', 0, 0, 31, 511}
#endif
};
key_descriptor stdf_info_keys[] = {
#if !defined(Little_Endian)
  {'AXI1', 31, 31, 0, 0},
  {'AXI2', 63, 31, 0, 0}
#else
  {'AXI1', 0, 0, 31, 31},
  {'AXI2', 0, 0, 31, 63}
#endif
};
#else
//! File table , exported symbol
extern file_table_entry_ptr file_table[MAX_XDF_FILES];

//! Buffer to format error messages
extern char errmsg[1024];
//! Tolerance level before messages are issued
extern int msg_level;
//! Tolerance level before program is aborted
extern int xdf_toler;
//! !stride
extern int xdf_stride;
//! Double float indicator
extern int xdf_double;
//! Short integer indicator
extern int xdf_short;
//! Byte, char array indicator
extern int xdf_byte;
//! Data type of last record read
extern int xdf_datatyp;
//! Enforce 8 char for date specifications
extern int xdf_enforc8;
//! Number of bits per FORTRAN word
extern int FTN_Bitmot;
//! No pack/unpack, used by editfst
extern int image_mode_copy;
//! Chekcpoint mode, no closing of the file
extern int xdf_checkpoint;
//! If one std seq file is opened, the limit of opened files becomes 128
extern int STDSEQ_opened;
extern key_descriptor stdfkeys[];
extern key_descriptor stdf_info_keys[];
#endif
