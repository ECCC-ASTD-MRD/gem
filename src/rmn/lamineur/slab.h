#define MAX_SLAB_FILES 10      /* maximum number of open slab files */
#define MAX_SLAB_TYPES 50      /* maximum number of open slab types */
#define ERR_NO_FILE -1         /* error opening file */
#define ERR_TAB_FULL -2        /* file table if full */
#define MAX_LEN 257            /* maximum number of char for a file name */
#define MAX_ETIKET 13          /* maximum number of char for etiket */



typedef struct {
                char file_name[MAX_LEN];    /* slab file name */
                int nrows[MAX_SLAB_TYPES];  /* number of rows to write */
                int count[MAX_SLAB_TYPES];  /* counter to check that all columns are written */
                int nio[MAX_SLAB_TYPES];    /* x dimension of full output grid */
                int i1[MAX_SLAB_TYPES];     /* index of 1st point along x */
                int ni[MAX_SLAB_TYPES];     /* x dimension of '#' grid section */
                int njo[MAX_SLAB_TYPES];    /* y dimension of full output grid */
                int j1[MAX_SLAB_TYPES];     /* index of 1st point along y */
                int nj[MAX_SLAB_TYPES];     /* x dimension of '#' grid section */
                                            /* for any grid type but '#', i1=j1=1 */
                                            /*                    ni=nio, nj=njo  */
                unsigned INT_32 *buffer;    /* data buffer */
                int pos;                    /* current write position into buffer */
                } file_table_desc;

typedef struct{
               INT_32 slb0, nBytes, deet, npas, dateo1,dateo2;
               float val15;
               INT_32 Ietiket[3];   /* 4 char variable disguised as 32 bit integer */
               }Id_Block_file;

typedef struct{
               INT_32  slb1, nBytes, slab_id,
	               ig1, ig2, ig3, ig4, Nrows,
                       Niout, Njout, nxgrid, nygrid, 
	               Nextra, ig1_, ig2_, ig3_, ig4_;
               INT_32 Igrtyp ,Igrtyp_;   /* 4 char variable disguised as 32 bit integer */
               }Slab_Descrt_file;

typedef struct{
               INT_32 slb2,nBytes,slab_id, nX, Nrows;
               }Data_Block_file;

typedef struct{
               INT_32 id_end;
               INT_32 nBytes;
	       }Slab_End;

#define SWAP32(temp) ( ((temp>>24) & 0xFF) | ((temp & 0xFF)<<24) | ((temp>>8) & 0xFF00) | ((temp & 0xFF00)<<8) )
