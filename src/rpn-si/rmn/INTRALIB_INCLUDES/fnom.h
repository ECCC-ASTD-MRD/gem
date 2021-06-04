#define MAXFILES 1024
typedef struct {
  unsigned int 
     stream:1, std:1, burp:1, rnd:1, wa:1, ftn:1,
     unf:1, read_only:1, old:1, scratch:1, notpaged:1,
     pipe:1, write_mode:1, remote:1, padding:18;
} attributs;
typedef struct {
  char * file_name;            /* complete file name */
  char * subname;              /* sub file name for cmcarc files */
  char * file_type;            /* file type and options */
  INT_32 iun;                    /* fnom unit number */
  INT_32 fd;                     /* c file descriptor */
#if defined(NEC)
  INT_32 file_size;         /* file size in words */
  INT_32 eff_file_size;     /* effective file size in words */
#else
  INT_32 file_size;              /* file size in words */
  INT_32 eff_file_size;          /* effective file size in words */
#endif
  INT_32 lrec;                   /* record length when appliable */
  INT_32 open_flag;              /* open/close flag */
  attributs attr;
} general_file_info;
#if defined(FNOM_OWNER)
general_file_info Fnom_General_File_Desc_Table[MAXFILES];
#else
extern general_file_info Fnom_General_File_Desc_Table[MAXFILES];
#endif
#define FGFDT Fnom_General_File_Desc_Table

int c_fretour(int iun);
ftnword f77name(fretour)(ftnword *fiun);
void f77name(d_fgfdt)();
int c_fnom(int *iun,char *nom,char *type,int lrec);
ftnword f77name(fnom)(ftnword *iun,char *nom,char *type,ftnword *flrec,F2Cl l1,F2Cl l2);
int c_fclos(int iun);
ftnword f77name(fclos)(ftnword *fiun);
ftnword f77name(qqqfnom)(ftnword *iun,char *nom,char *type,ftnword *flrec,F2Cl l1,F2Cl l2);
void c_waopen(int iun);
int c_waopen2(int iun);
ftnword f77name(waopen2)(ftnword *fiun);
void f77name(waopen)(ftnword *fiun);
void c_waclos(int iun);
int c_waclos2(int iun);
ftnword f77name(waclos2)(ftnword *fiun);
void f77name(waclos)(ftnword *fiun);
void c_wawrit(int iun,void *buf,unsigned int adr,int nmots);
int c_wawrit2(int iun,void *buf,unsigned int adr,int nmots);
void f77name(wawrit)(ftnword *fiun,void *buf,unsigned ftnword *fadr,ftnword *fnmots);
ftnword f77name(wawrit2)(ftnword *fiun,void *buf,unsigned ftnword *fadr,ftnword *fnmots);
void c_waread(int iun,void *buf,unsigned int adr,int nmots);
int c_waread2(int iun,void *buf,unsigned int adr,int nmots);
void f77name(waread)(ftnword *fiun,void *buf,unsigned ftnword *fadr,ftnword *fnmots);
ftnword f77name(waread2)(ftnword *fiun,void *buf,unsigned ftnword *fadr,ftnword *fnmots);
INT_32 c_wasize(int iun);
ftnword f77name(wasize)(ftnword *fiun);
INT_32 c_numblks(int iun);
ftnword f77name(numblks)(ftnword *fiun);
ftnword f77name(existe)(char *nom,F2Cl lng);
void c_openda(int iun);
void f77name(openda)(ftnword *iun);
void c_closda(int iun);
void f77name(closda)(ftnword *iun);
void c_checda(int iun);
void f77name(checda)(ftnword *iun);
void c_readda(int iun,int *bufptr,int ns,int is);
void f77name(readda)(ftnword *iun,ftnword *bufptr,ftnword *ns,ftnword *is);
void c_writda(int iun,int *bufptr,int ns,int is);
void f77name(writda)(ftnword *iun,ftnword *bufptr,ftnword *ns,ftnword *is);
int c_getfdsc(int iun);
ftnword f77name(getfdsc)( ftnword *iun);
void c_sqopen(int iun);
void f77name(sqopen)(ftnword *iun);
void c_sqclos(int iun);
void f77name(sqclos)(ftnword *iun);
void c_sqrew(int iun);
void f77name(sqrew)(ftnword *iun);
void c_sqeoi(int iun);
void f77name(sqeoi)(ftnword *iun);
int c_sqgetw(int iun, word *bufptr, int nmots);
ftnword f77name(sqgetw)(ftnword *iun, ftnword *bufptr, ftnword *nmots);
int c_sqputw(int iun, word *bufptr, int nmots);
ftnword f77name(sqputw)(ftnword *iun, ftnword *bufptr, ftnword *nmots);
int c_sqgets(int iun, char *bufptr, int nchar);
ftnword f77name(sqgets)(ftnword *iun, char  *bufptr, ftnword *nchar, F2Cl lbuf);
int c_sqputs(int iun, char *bufptr, int nchar);
ftnword f77name(sqputs)(ftnword *iun, char  *bufptr, ftnword *nchar, F2Cl lbuf);
void f77name(d_wafdt)();  
unsigned ftnword f77name(hrjust) (unsigned ftnword *moth, ftnword *ncar);
unsigned ftnword f77name(hljust) (unsigned ftnword *moth, ftnword *ncar);
unsigned INT_32 f77name(check_host_id)();
