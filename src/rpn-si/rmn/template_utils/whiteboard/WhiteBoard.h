#ifndef WHITEBOARD_VERSION

#define WHITEBOARD_VERSION "1.00"   /* Code revision: $Id: WhiteBoard.h 582 2008-10-29 18:57:41Z armnlib $ */

#define WB_FORTRAN_REAL 1
#define WB_FORTRAN_INT 2
#define WB_FORTRAN_CHAR 3
#define WB_FORTRAN_BOOL 4
#define WB_OPTION_SET(options,option) (0 .ne. iand(options,option))

#define WB_IS_ARRAY 4096
#define WB_REWRITE_AT_RESTART 2048
#define WB_REWRITE_MANY 1024
#define WB_REWRITE_UNTIL 512
#define WB_REWRITE_NONE 256
#define WB_DEFAULT WB_REWRITE_NONE
#define WB_READ_ONLY_ON_RESTART 128
#define WB_INITIALIZED 64
#define WB_BADVAL 32
#define WB_HAS_RULES 16
#define WB_IS_LOCAL 8
#define WB_CREATED_BY_RESTART 4
#define WB_NOTINITIALIZED 2
#define WB_CREATE_ONLY 1

#define WB_STRICT_DICTIONARY 2
#define WB_ALLOW_DEFINE 1
#define WB_FORBID_DEFINE 0

#define WB_MSG_DEBUG -1
#define WB_MSG_INFO 0
#define WB_MSG_WARN 1
#define WB_MSG_ERROR 2
#define WB_MSG_SEVERE 3
#define WB_MSG_FATAL 4

#define WB_IS_OK(errcode) (errcode >= 0)
#define WB_IS_ERROR(errcode) (errcode < 0)
#define WB_OK 0
#define WB_ERROR -1
#define WB_ERR_NAMETOOLONG -1000
#define WB_ERR_NOTFOUND -1001
#define WB_ERR_READONLY -1002
#define WB_ERR_WRONGTYPE -1003
#define WB_ERR_WRONGDIMENSION -1004
#define WB_ERR_ALLOC -1005
#define WB_ERR_NOTYPE -1006
#define WB_ERR_NOMEM -1007
#define WB_ERR_NOVAL -1008
#define WB_ERR_BADVAL -1009
#define WB_ERR_WRONGSTRING -1010
#define WB_ERR_CKPT -1011
#define WB_ERR_REDEFINE -1012
#define WB_ERR_BIG -1013
#define WB_ERR_SYNTAX -1014
#define WB_ERR_OPTION -1015
#define WB_ERR_READ -1016

#define WB_MAXSTRINGLENGTH 520
#define WB_MAXNAMELENGTH 27

#endif

#ifdef C_SOURCE_CODE

#define WB_PUT_C(WB,name,value,strglen,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_CHAR,strglen,(unsigned char *)value,0,options,strlen(name))
#define WB_GET_C(WB,name,value,strglen) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_CHAR,strglen,(unsigned char *)value,0,strlen(name))
#define WB_PUT_CV(WB,name,value,strglen,size,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_CHAR,strglen,(unsigned char *)value,size,options,strlen(name))
#define WB_GET_CV(WB,name,value,strglen,size) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_CHAR,strglen,(unsigned char *)value,size,strlen(name))

#define WB_PUT_L1(WB,name,value,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_BOOL,1,(unsigned char *)value,0,options,strlen(name))
#define WB_GET_L1(WB,name,value) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_BOOL,1,(unsigned char *)value,0,strlen(name))
#define WB_PUT_L1V(WB,name,value,size,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_BOOL,1,(unsigned char *)value,size,options,strlen(name))
#define WB_GET_L1V(WB,name,value,size) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_BOOL,1,(unsigned char *)value,size,strlen(name))

#define WB_PUT_I4(WB,name,value,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_INT,4,(unsigned char *)value,0,options,strlen(name))
#define WB_GET_I4(WB,name,value) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_INT,4,(unsigned char *)value,0,strlen(name))
#define WB_PUT_I4V(WB,name,value,size,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_INT,4,(unsigned char *)value,size,options,strlen(name))
#define WB_GET_I4V(WB,name,value,size) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_INT,4,(unsigned char *)value,size,strlen(name))

#define WB_PUT_I8(WB,name,value,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_INT,8,(unsigned char *)value,0,options,strlen(name))
#define WB_GET_I8(WB,name,value) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_INT,8,(unsigned char *)value,0,strlen(name))
#define WB_PUT_I8V(WB,name,value,size,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_INT,8,(unsigned char *)value,size,options,strlen(name))
#define WB_GET_I8V(WB,name,value,size) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_INT,8,(unsigned char *)value,size,strlen(name))

#define WB_PUT_R4(WB,name,value,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_REAL,4,(unsigned char *)value,0,options,strlen(name))
#define WB_GET_R4(WB,name,value) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_REAL,4,(unsigned char *)value,0,strlen(name))
#define WB_PUT_R4V(WB,name,value,size,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_REAL,4,(unsigned char *)value,size,options,strlen(name))
#define WB_GET_R4V(WB,name,value,size) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_REAL,4,(unsigned char *)value,size,strlen(name))

#define WB_PUT_R8(WB,name,value,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_REAL,8,(unsigned char *)value,0,options,strlen(name))
#define WB_GET_R8(WB,name,value) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_REAL,8,(unsigned char *)value,0,strlen(name))
#define WB_PUT_R8V(WB,name,value,size,options) c_wb_put(WB,(unsigned char *)name,WB_FORTRAN_REAL,8,(unsigned char *)value,size,options,strlen(name))
#define WB_GET_R8V(WB,name,value,size) c_wb_get(WB,(unsigned char *)name,WB_FORTRAN_REAL,8,(unsigned char *)value,size,strlen(name))

typedef struct {
   int code;
   char *text;
   }SYMTAB;

typedef struct{
   unsigned int type:3,               /* real/int/logical/char 0 means invalid entry */
                array:1,              /* 1 if array  */
                fromrestart:1,        /* 1 if variable has been created by a restart */
                readonly:1,           /* 1 if variable is now read-only */
                elementsize:10,       /* size of each element 4/8 for real/int, 0-1023 for char, 4 for logical */
                badval:1,             /* 1 if value is dubious due to a bad put */
                islocal:1,            /* 1 if variable is not global (same on all MPI Whiteboards */
                resetuntil:1,         /* 1 if variable can be set until package/variable is locked */
                resetmany:1,          /* 1 if value can be set any number of times */
                noresetafterrestart:1,/* 1 if after restart value becomes read-only */
                initialized:1,        /* variable has been set */
                lines:10;             /* number of datalines occupied by variable. if 0 it means entire page */
   }FLAGS;

typedef union {
   unsigned char      c[WB_MAXNAMELENGTH+1];
   unsigned int       i[(WB_MAXNAMELENGTH+1)/sizeof(int)];
   }NAMETOKEN;  /* large token, full name, integers */

typedef struct {
   int UsedElements;
   int MaxElements;
   }ARRAYDESC;                 /* array metadata container */

typedef union {
   unsigned char      c[8];
   unsigned int       i[2];
   unsigned long long l;
   ARRAYDESC          a;
   }BIGTOKEN;  /* 64 bit token, chars,  ints, long long, array metadata  */

typedef struct {            /* THE SIZE OF THIS STRUCTURE SHOULD BE A MULTIPLE OF 8 BYTES */
   NAMETOKEN name;          /* WB_MAXNAMELENGTH characters name */
   FLAGS     flags;         /* type, unused bytes, size of element, number of lines */
   BIGTOKEN  data;          /* value or array descriptor (this token is 8 byte aligned) */
   }METALINE;               /* whiteboard data line, one line per scalar item <= 8 bytes */

typedef struct {            /* THE SIZE OF THIS STRUCTURE SHOULD BE A MULTIPLE OF 8 BYTES */
   unsigned char data[sizeof(METALINE)];
   }DATALINE;

typedef union {     /* a page line is either a metadata + short data or a long data */
   METALINE m;      /* metadata + data for a simple real integer logical or <=8 characters string */
   DATALINE d;      /* long data container (strings > 8 characters or arrays )  */
   }LINE;

typedef struct {    /* used by wb_read to check declaration/assignation consistency */
   LINE *line;
   int defined;
   int assigned;
   }DEFTBL;

#define WB_MAXLINESPERPAGE 32
#define WB_MAXLINESPERPAGESHIFT 8
typedef struct WB_page{
   struct WB_page *next;      /* really an address only but an 8 byte item should be forced for alignment reasons */
   struct WB_page *not_used;  /* this guarantees a multiple of 8 bytes before line (that needs 8 byte alignment */
   int   NextFreeEntry;
   int   LinesInPage;
   LINE  line[WB_MAXLINESPERPAGE];
   }PAGE;                    /* whiteboard page , size is 16 or 24 bytes + space needed by line */

#else
      type :: whiteboard
        sequence!WARNING: this is requested on AIX but only exists from f95 on
        integer*8 :: wb
      end type whiteboard
#endif
