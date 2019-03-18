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

/*****************************************************************************
 *                                                                           *
 *WKOFFIT                                                                    *
 *                                                                           *
 *AUTEUR  M. LEPINE  -  NOV 1994                                             *
 *MODIFICATION : ERIC MICHAUD - JANV 1995                                    *
 *Modification : M. Lepine - Juil 2011 (reconnaissance fichiers NetCDF)      *
 *Modification : M. Lepine - Juil 2011 (bug fix get_mode)                    *
 *                                                                           *
 *OBJET                                                                      *
 *     DETERMINER QUEL EST LE TYPE D'UN FICHIER                              *
 *                                                                           *
 *         EX:   FICHIER STANDARD (89) RANDOM, SEQUENTIEL,                   *
 *               FICHIER FORTRAN ...                                         *
 *                                                                           *
 *     VALEURS DE RETOUR POSSIBLES                                           *
 *                                                                           *
 *        -3     FICHIER INEXISTANT                                          *
 *        -2     FICHIER VIDE                                                *
 *        -1     FICHIER INCONNU                                             *
 *         1     FICHIER STANDARD RANDOM 89                                  *
 *         2     FICHIER STANDARD SEQUENTIEL 89                              *
 *         3     FICHIER STANDARD SEQUENTIEL FORTRAN 89                      *
 *         4     FICHIER CCRN                                                *
 *         5     FICHIER CCRN-RPN                                            *
 *         6     FICHIER BURP                                                *
 *         7     FICHIER GRIB                                                *
 *         8     FICHIER BUFR                                                *
 *         9     FICHIER BLOK                                                *
 *        10     FICHIER FORTRAN                                             *
 *        11     FICHIER COMPRESS                                            *
 *        12     FICHIER GIF89                                               *
 *        13     FICHIER GIF87                                               *
 *        14     FICHIER IRIS                                                *
 *        15     FICHIER JPG                                                 *
 *        16     FICHIER KMW                                                 *
 *        17     FICHIER PBM                                                 *
 *        18     FICHIER PCL                                                 *
 *        19     FICHIER PCX                                                 *
 *        20     FICHIER PDSVICAR                                            *
 *        21     FICHIER PM                                                  *
 *        22     FICHIER PPM                                                 *
 *        23     FICHIER PS                                                  *
 *        24     FICHIER KMW_                                                *
 *        25     FICHIER RRBX                                                *
 *        26     FICHIER SUNRAS                                              *
 *        27     FICHIER TIFF                                                *
 *        28     FICHIER UTAHRLE                                             *
 *        29     FICHIER XBM                                                 *
 *        30     FICHIER XWD                                                 *
 *        31     FICHIER ASCII                                               *
 *        32     FICHIER BMP                                                 *
 *        33     FICHIER STANDARD RANDOM 98                                  *
 *        34     FICHIER STANDARD SEQUENTIEL 98                              *
 *        35     FICHIER NETCDF                                              * 
 *                                                                           *
 *****************************************************************************/

#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64

#define WKF_INEXISTANT            -3
#define WKF_VIDE                  -2
#define WKF_INCONNU               -1
#define WKF_RANDOM89               1
#define WKF_SEQUENTIEL89           2 
#define WKF_SEQUENTIELFORTRAN89    3
#define WKF_CCRN                   4
#define WKF_CCRN_RPN               5
#define WKF_BURP                   6
#define WKF_GRIB                   7
#define WKF_BUFR                   8
#define WKF_BLOK                   9
#define WKF_FORTRAN               10
#define WKF_COMPRESS              11 
#define WKF_GIF89                 12 
#define WKF_GIF87                 13  
#define WKF_IRIS                  14 
#define WKF_JPG                   15
#define WKF_KMW                   16
#define WKF_PBM                   17
#define WKF_PCL                   18
#define WKF_PCX                   19
#define WKF_PDSVICAR              20
#define WKF_PM                    21
#define WKF_PPM                   22
#define WKF_PS                    23
#define WKF_KMW_                  24
#define WKF_RRBX                  25
#define WKF_SUNRAS                26
#define WKF_TIFF                  27
#define WKF_UTAHRLE               28
#define WKF_XBM                   29
#define WKF_XWD                   30
#define WKF_ASCII                 31
#define WKF_BMP                   32
#define WKF_RANDOM98              33
#define WKF_SEQUENTIEL98          34
#define WKF_NETCDF                35

#define SIGN_STD89_RND  012525252525   
#define SIGN_STD89_SEQ  025252525252

#define  FALSE     (0==1)
#define  TRUE      (0==0)

#define  IDGIF87   "GIF87a"
#define  IDGIF89   "GIF89a"

#define  RAS_MAGIC 0x59a66a95


/*  Sun supported ras_type's */

#define  RT_OLD          0      /* Raw pixrect image in 68000 byte order */
#define  RT_STANDARD     1      /* Raw pixrect image in 68000 byte order */
#define  RT_BYTE_ENCODED 2      /* Run-length compression of bytes */
#define  RT_EXPERIMENTAL 0xffff /* Reserved for testing */


/*  Sun registered ras_maptype's */
 

#define  RMT_RAW 2


/*  Sun supported ras_maptype's */
 
#define  RMT_NONE      0 
#define  RMT_EQUAL_RGB 1 



#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#ifdef should_never_be_true
#include <X11/Xmd.h>
#else
  #ifdef WIN32    /*CHC/NRC*/
    typedef unsigned long CARD32;
    #define rindex(a,b) strrchr((a),(b))
  #else
    #include <unistd.h>
    #define CARD32 unsigned int
  #endif
#define B32 :32
#endif

#include <rpnmacros.h>



/*  RRBX stuff */
 
typedef struct {
                  int ftn1;
                  int nbits;
                  int nr;
                  int nw;
                  int ncp;
                  int ftn2[2];
                  char plotid[80];
                  int ftn3[2];
                  char framid[80];
                  int ftn4;
               }  Rrbxfile;


/*  Description of header for files containing raster images */

 struct  rasterfile {
         int ras_magic;     /* magic number */
         int ras_width;     /* width (pixels) of image */
         int ras_height;    /* height (pixels) of image */
         int ras_depth;     /* depth (1, 8, or 24 bits) of pixel */
         int ras_length;    /* length (bytes) of image */
         int ras_type;      /* type of file; see RT_* below */
         int ras_maptype;   /* type of colormap; see RMT_* below */
         int ras_maplength; /* length (bytes) of following map */
         };

/*  XWDFile  stuff */
 
#define XWD_FILE_VERSION 7

#ifdef ALL64
#define sz_XWDheader 104
#else
#define sz_XWDheader 100
#endif
#define sz_XWDColor 12

typedef CARD32 xwdval;      /* for old broken programs */

typedef struct _xwd_file_header {
    CARD32 header_size B32;  /* Size of the entire file header (bytes). */
    CARD32 file_version B32;    /* XWD_FILE_VERSION */
    CARD32 pixmap_format B32;   /* Pixmap format */
    CARD32 pixmap_depth B32;    /* Pixmap depth */
    CARD32 pixmap_width B32;    /* Pixmap width */
    CARD32 pixmap_height B32;   /* Pixmap height */
    CARD32 xoffset B32;     /* Bitmap x offset */
    CARD32 byte_order B32;      /* MSBFirst, LSBFirst */
    CARD32 bitmap_unit B32;     /* Bitmap unit */
    CARD32 bitmap_bit_order B32;    /* MSBFirst, LSBFirst */
    CARD32 bitmap_pad B32;      /* Bitmap scanline pad */
    CARD32 bits_per_pixel B32;  /* Bits per pixel */
    CARD32 bytes_per_line B32;  /* Bytes per scanline */
    CARD32 visual_class B32;    /* Class of colormap */
    CARD32 red_mask B32;        /* Z red mask */
    CARD32 green_mask B32;      /* Z green mask */
    CARD32 blue_mask B32;       /* Z blue mask */
    CARD32 bits_per_rgb B32;    /* Log2 of distinct color values */
    CARD32 colormap_entries B32;    /* Number of entries in colormap */
    CARD32 ncolors B32;     /* Number of Color structures */
    CARD32 window_width B32;    /* Window width */
    CARD32 window_height B32;   /* Window height */
    CARD32 window_x B32;        /* Window upper left X coordinate */
    CARD32 window_y B32;        /* Window upper left Y coordinate */
    CARD32 window_bdrwidth B32; /* Window border width */
#ifdef ALL64
    CARD32 header_end B32;      /* Pad to fill out word */
#endif
} XWDFileHeader;

/*  PPM stuff */
/* Magic constants. */
#define PPM_MAGIC1 'P'
#define PPM_MAGIC2 '3'
#define RPPM_MAGIC2 '6'
#define PPM_FORMAT (PPM_MAGIC1 * 256 + PPM_MAGIC2)
#define RPPM_FORMAT (PPM_MAGIC1 * 256 + RPPM_MAGIC2)

/*  PCL stuff */

static   char *enter_pcl = "\033%-12345X@PJL ENTER LANGUAGE=PCL";
static   char *enter_ps = "\033%-12345X@PJL ENTER LANGUAGE=POSTSCRIPT";
#define  PCL_MAX_LEN    512000

static float    Get_Frac();

static int mutant_kmw=FALSE;

static int ReadFileType(char *fname);
static void Flush_Bytes(unsigned int num, FILE *fp );
static void Flush_To_Term(FILE *fp);
static int ispcl(char *path );
static int isppm(char *path );
static int issun(char *path );
static  int isrrbx(char *path );
static  int isrrbx(char *path );
static int iskmw(char *path);
static int isgif(char *path );
static int isxwd(char *path );
static int isps(char *path );
static int test_fichier (char *nom );

/*****************************************************/

static int retour(pf,code)
FILE *pf;
int code;
{
   fclose(pf);
   return(code);
   }

/******************************************************/

static int
isftnbin(pf,lng)
FILE *pf;
int lng;
{
   int mot;
   INT_32 offset;

   offset = lng +4;
   fseek(pf,offset,0);
   fread32(&mot,sizeof(int),1,pf);
   if (mot == lng) {
      return(1);
   }
   else {
      return(0);
   }
}

/********************************************/

wordint
c_wkoffit(char *nom,int l1) 
{
   FILE *pf;
   char nom2[4096], nom3[4096], *pn2, *pn3;
   int buffer[1024], *ptbuf, lowc=0;
   INT_32 pos,lngf;
   int longnom;

   longnom = ( ( l1 <= 4095 ) ? l1 : 4095 );
   pos = 0;
   ptbuf = &buffer[0];
   if (nom[0] == '+') {     /* garder le nom de fichier tel quel */
     longnom--;
     strncpy(nom2,nom+1,longnom);
     nom2[longnom] = '\0';
     while (nom2[--longnom] == ' ')
       nom2[longnom]='\0';
   }
   else {
     strncpy(nom2,nom,longnom);
     nom2[longnom] = '\0';
     pn2 = &nom2[0];
     pn3 = &nom3[0];
     while ((*pn2 != ' ') && (*pn2 != '\0')) {
       if (islower(*pn2)) {
         *pn3 = *pn2;
         lowc = 1;
       }
       else
         *pn3 = tolower(*pn2);
       pn2++;
       pn3++;
     }
     *pn2 = '\0';
     *pn3 = '\0';
     if (lowc == 0)
       strcpy(nom2,nom3);
   }
   pf = fopen(nom2,"rb");
   if (pf == (FILE *) NULL)
      return(-3);
   else {

     /* positionnement a la fin du fichier */
      fseek(pf,pos,2);       
      lngf=ftell(pf);
      if (lngf == 0) return(retour(pf,-2));

     /* positionnement au debut du fichier */
      fseek(pf,pos,0);     
      fread32(ptbuf,sizeof(int),1024,pf);

     /* RANDOM89 */
      if (*ptbuf == SIGN_STD89_RND && *(ptbuf+1) == SIGN_STD89_RND)
	 return(retour(pf,WKF_RANDOM89));
      else

     /* CCRN */
      if (*(ptbuf) == 64 && *(ptbuf+17) == 64 && *(ptbuf+2) == 0x20202020)
         return(retour(pf,WKF_CCRN));
      else

     /* CCRN-RPN */
      if (*(ptbuf+2) == 0x504b3834 && isftnbin(pf,*ptbuf))  /* PK84 */
         return(retour(pf,WKF_CCRN_RPN));
      else

     /* SEQUENTIEL89 */
      if (*(ptbuf+28) == SIGN_STD89_SEQ && *(ptbuf+29) == SIGN_STD89_SEQ)
         return(retour(pf,WKF_SEQUENTIEL89));
      else

     /* SEQUENTIELFORTRAN89 */ 
      if (*(ptbuf+29) == SIGN_STD89_SEQ && *(ptbuf+30) == SIGN_STD89_SEQ
                       && isftnbin(pf,*ptbuf))
         return(retour(pf,WKF_SEQUENTIELFORTRAN89));
      else
 
    /* STANDARD 98 RANDOM */
      if (*(ptbuf+3) == 'STDR') 
         return(retour(pf,WKF_RANDOM98));
      else

    /* STANDARD 98 SEQUENTIEL */
      if (*(ptbuf+3) == 'STDS') 
         return(retour(pf,WKF_SEQUENTIEL98));
      else

    /* BURP */
      if ((*(ptbuf+3) == 'BRP0') || (*(ptbuf+3) == 'bRp0'))
         return(retour(pf,WKF_BURP));
      else

    /* GRIB */
      if (*(ptbuf) == 0x47524942)   
         return(retour(pf,WKF_GRIB));
      else

    /* BUFR */
      if (*(ptbuf) == 0x42554652)  
         return(retour(pf,WKF_BUFR));
      else

    /* NetCDF classic format */
      if (*(ptbuf) == 'CDF\001')
         return(retour(pf,WKF_NETCDF));
      else
	
    /* NetCDF 64-bit offset format */
      if (*(ptbuf) == 'CDF\002')
         return(retour(pf,WKF_NETCDF));
      else

    /* BLOK */
      if (*(ptbuf) == 0x424c4f4b)   
         return(retour(pf,WKF_BLOK));
      else

    /* FORTRAN */
      if (isftnbin(pf,*ptbuf))
         return(retour(pf,WKF_FORTRAN));
      else {
   
    /* INCONNU  */
	     return(retour(pf,test_fichier (nom2) ));
      }
   }
}

wordint f77name(wkoffit)(char *nom, F2Cl fl1)
{
  int l1=fl1;
  
  return(c_wkoffit(nom,l1));
}

/****************************************************/

static int test_fichier ( nom )
char *nom;
{
  int id;
  int repgif;
 
  repgif = isgif (nom);
  if ( repgif != FALSE ) {
     return repgif;
  } else
    if ( isrrbx( nom ) ) {
     return WKF_RRBX;
  } else
    if ( issun( nom ) ) {
     return WKF_SUNRAS;
  } else
    if ( isxwd( nom ) ) {
     return WKF_XWD;
  } else
   if ( isppm( nom ) ) {
     return WKF_PPM;
  } else
    if ( iskmw( nom ) ) {
     if ( mutant_kmw )
       return WKF_KMW_;
     else
       return WKF_KMW;
  } else
    if ( isps( nom ) ) {
     return WKF_PS;
  } else {

/*  essais la routine de xv-3.00a */
 
     id = ReadFileType( nom );
     if ( id == WKF_INCONNU ) {

/*  dernier espoir */
 
       if ( ispcl( nom ) )
         return WKF_PCL;
       else
         return WKF_INCONNU;
     } else
       return id;
  }
}


/*
 *
 *  module    :  ISPS
 *
 *  auteur    :  VANH SOUVANLASY
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A PostScript
 *
 */

static int isps( path )
char *path;
{
   FILE *fp;
   char buffer[256];
   int  i, j, ps_i;

   if( (fp = fopen( path, "rb")) == NULL ) return(FALSE);

   i = 0;
   while ( fgets( buffer, 256, fp ) != NULL ) {

/*  first must enter Postcript using PJL */
 
     if ( i == 0 ) {
       for ( j = 0 ; (j < 256)&&(buffer[j]!='\0') ; j++ )
         buffer[j] = (char)toupper((int)buffer[j]);
       if ( strncmp( buffer, enter_ps, strlen(enter_ps) ) == 0 ) {
         fclose( fp );
         return (TRUE);
       }
     }

     if ( strncmp( buffer, "%%", 2 ) == 0 ) continue;
     ++i;
     if ( i > 66 ) break;
     if ( strncmp( buffer, "%!", 2 ) == 0 ) {
       while ( fgets( buffer, 256, fp ) != NULL ) {
         ++i;
         if ( i > 66 ) break;
         switch( buffer[0] ) {
         case '%' :
         case '/' :
           fclose( fp );
           return (TRUE);
         break;
         case '\n' :
         break;
         default :
           fclose( fp );
           return (FALSE);
         break;
         }
       }
     }
     if ( strncmp( buffer, "%", 1 ) == 0 ) continue;
     break;
   }
   fclose( fp );
   return (FALSE);
}

/*
 *
 *  module    :  ISXWD
 *
 *  auteur    :  VANH SOUVANLASY
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A XWDFile
 *
 */

static int isxwd( path )
char *path;
{
   FILE   *fp;
   XWDFileHeader xwd;
   int   flen, len;

   if( (fp = fopen( path, "rb")) == NULL ) return(FALSE);

   if (fread32 ((char *) &xwd, sizeof(XWDFileHeader), 1, fp) != 1)
   {
       fclose( fp );
       return (FALSE);
   }

   if (fseek( fp, 0L, SEEK_END ))
   {
       fclose( fp );
       return(FALSE);
   }
   flen = ftell(fp);
   fclose(fp);

   if (xwd.file_version != XWD_FILE_VERSION )
      return (FALSE);

   len = xwd.header_size + xwd.ncolors * sz_XWDColor
         + xwd.bytes_per_line * xwd.pixmap_height;

   if ( len != flen ) return (FALSE);
   return(TRUE);
}



/*
 *
 *  module    :  ISGIF
 *
 *  auteur    :  VANH SOUVANLASY
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A GIF
 *
 */

static int isgif( path )
char *path;
{
   FILE   *fp;
   char    magic[6];
   int     status = FALSE;

   if( (fp = fopen( path, "rb")) == NULL ) return(FALSE);

   if (fread(magic, 6, 1, fp) != 1) return(FALSE);
   if (strncmp( magic, IDGIF87, 6 )==0)
     status = WKF_GIF87;
   if (strncmp( magic, IDGIF89, 6 )==0) 
     status = WKF_GIF89;
   fclose(fp);
   return(status);
}

/*
 *
 *  module    :  GET_MODE
 *
 *  auteur    :  MICHEL GRENIER    CMCOG
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  CE MODULE TENTE D'INDENTIFIER UN FICHIER KMW
 *
 */

 static get_mode ( fp, depth, height, width, length )
 FILE  *fp;
 int  *depth ;
 int  *height ;
 int  *width ;
 INT_32 *length;
     {
     int cpt = 0;
     int kmwndx = 0 ;
     int mode = 0;
     char buf[351];
     static char sigkmw[] = "PLOT$Z";

     register char *data = buf;
     INT_32  curpos;

/*
 *  verifie dans les 350 premiers caracteres si on trouve
 *  la chaine PLOT$Z qui indique un fichier KMW
 */
     while( (*data = getc(fp)) != EOF && cpt++ < 350)
          {
/*	  printf("Debug cpt=%d \n",cpt); */
          if( *data == sigkmw[kmwndx] ) {
/*            printf("Debug kmwndx=%d\n",kmwndx); */
            if( kmwndx == strlen(sigkmw) - 1)
              {
              mode = 1;
              break;
              }
            else 
	      {
              kmwndx++;
	      }
	  }
	  else
	    kmwndx=0;
          data++;
          }
if (mode == 0) return(mode);

*depth = 1;
/*
 *  va lire les dimensions et verifie que la decompression est possible
 */
     *data = '\0';
if (cpt < 10) mode = 0 ;
     if( mode != 0 )
       {
       data = (char *)rindex(buf, '\n');
       sscanf(data - 9, "%4d%4d%1d", width, height, depth );

       if( mode == 1 && (*depth != 1 && *depth != 3) )
        {
                *depth = 777 ;
                *width = 4224 ;
                *height = 6048  ;
        }
       }

    curpos = ftell(fp);
    if (fseek( fp, 0L, SEEK_END ))  {
       fclose( fp );
       perror("fseek");
       exit(0);
    }

    *length = ftell(fp);
    if (fseek( fp, curpos, SEEK_SET ))  {
       fclose( fp );
       perror("fseek");
       exit(0);
    }

     return(mode);
     }

/*
 *
 *  module    :  ISKMW
 *
 *  auteur    :  MICHEL GRENIER
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A KMW
 *
 */

static int iskmw(path)
char *path;
{
   FILE *fp;
   int height,width, depth;
   int ierr=0;
   INT_32 length;

   if( (fp = fopen( path, "rb")) == NULL ) return(ierr);

/*
 *  verifie si c'est un fichier KMW
 */

    ierr = get_mode( fp, &depth, &height, &width, &length );
    if ( depth == 3 ) mutant_kmw = TRUE;

    fclose(fp) ;

    return(ierr);
 }

/*
 *
 *  module    :  ISRRBX
 *
 *  auteur    :  MICHEL GRENIER    CMCOG
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A RRBX FILE
 *
 */

static  int isrrbx( path )
 char *path;
     {
     FILE    *fp;
     Rrbxfile header;

/*
 *  open file
 */

    if( (fp = fopen( path, "rb")) == NULL ) return(FALSE);

/*
 *  read header
 */
    if( fread32( &header, sizeof(Rrbxfile), 1, fp ) == 0 ) return(FALSE);

/*
 *  close file
 */
    fclose(fp);

/*
 *  check magic string
 */
    if( strncmp(&header.plotid[8],"RRBXRRUX",8) != 0 ) return(FALSE);

/*
 *  return result
 */
    return(TRUE);
    }

/*
 *
 *  module    :  ISSUN
 *
 *  auteur    :  MICHEL GRENIER
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A SUNRASTER
 *
 */

 static int issun( path )
 char *path;
     {
     FILE   *fp;
     struct rasterfile header;

/*
 *  open file
 */

    if( (fp = fopen( path, "rb")) == NULL ) return(FALSE);

/*
 *  read sunraster header
 */
    if( fread32( &header, sizeof(struct rasterfile), 1, fp ) == 0 ) return(FALSE
);
/*
 *  check magic word
 */
    if( header.ras_magic != RAS_MAGIC )
      {
      fclose(fp);
      return(FALSE);
      }

/*
 *  close file
 */
    fclose(fp);
    switch( header.ras_maptype ) {
      case RT_OLD :
      case RT_STANDARD :
      case RT_BYTE_ENCODED :
        switch( header.ras_maptype) {
         case RMT_EQUAL_RGB :
         case RMT_NONE :
            return (TRUE);
         break;
         default :
         break;
         }
      break;
      default :
      break;
    }
/*
 *  return result
 */
     return(FALSE);
}

/*
 *
 *  module    :  ISPPM
 *
 *  auteur    :
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPPEMENT
 *
 *  langage   :  C
 *
 *  objet     :  THIS MODULE TEST IF A FILE IS A PPM
 *
 */

static int isppm( path )
char *path;
{
   FILE   *fp;
   int   c0, c1;
   int   magic;

/*
 *  open file
 */

   if( (fp = fopen( path, "rb")) == NULL ) return(FALSE);

   c0 = getc(fp);
   if ( c0 == EOF ) {
     fclose(fp);
     return FALSE;
   }
   c1 = getc(fp);
   if ( c1 == EOF ) {
     fclose(fp);
     return FALSE;
   }

/*
 *  close file
 */
   fclose(fp);

   magic = (c0<<8)+c1 ;
   if ((magic == PPM_FORMAT)||(magic == RPPM_FORMAT)) return (TRUE);

   return FALSE;
}

#define MIN(x,y)        ( ((x) < (y)) ? (x) : (y) )
#define ESC     27





/*
**  These variables are for the new parser.
**  The new parser handles more sequences, and also deals with combined
**  escape sequences better.
*/

static int ispcl( path )
char *path;
{
  FILE *fp;
  int   c,j;
int     parameter;
int     group_char;
int     terminator;
int     value;
float   fvalue;                 /* fractional value */
int     scanf_count;
char    in_sequence = FALSE;
char    pass_seq;
char    plus_sign;              /* for relative values */
char    strip_seq = FALSE;
INT_32    flen;
char    buffer[256];

/*
 *  open file
 */

  if( (fp = fopen( path, "rb")) == NULL )
      return(FALSE);

  if (fseek( fp, 0L, SEEK_END ))  {
       fclose( fp );
       return(FALSE);
  }
  flen = ftell(fp);
  if (fseek( fp, 0L, SEEK_SET ))  {
       fclose( fp );
       return(FALSE);
  }
/*
 * first line must enter PCL if it is not a graphic PCL
 */

   if ( fgets( buffer, 256, fp ) == NULL ) {
     fclose( fp );
     return(FALSE);
   }
   for ( j = 0 ; (j < 256)&&(buffer[j]!='\0') ; j++ )
     buffer[j] = (char)toupper((int)buffer[j]);
   if ( strncmp( buffer, enter_pcl, 13 ) == 0 ) {
     fclose( fp );
     return (TRUE);
   }

/*
 * escape sequence for enter graphic mode must be there
 * or it isnt a PCL
 */
   if (fseek( fp, 0L, SEEK_SET ))  {
       fclose( fp );
       return(FALSE);
   }


/*
 * This is the pcl input parsing loop.
 */
  while( ( c = getc(fp) ) != EOF )
  {

/*  Ignore all chars until an escape char  */

    if ( c != ESC ) continue;
/*
 *  Now we have an escape sequence, get the parameter char.
 */
    parameter = getc(fp);

    if ( parameter == EOF ) {
      fclose( fp );
      return FALSE;
    }

/*
 *  Now check to see if it is a two character sequence.
 */
    if ( parameter >= '0' && parameter <= '~' ) continue;

/*
 *  Now check to make sure that the parameter character is
 *  within range.
 */
     if ( parameter < '!' || parameter > '/' ) {
       fclose( fp );
       return (FALSE);
     }
/*
 *  We are only interested in certain parameters, so pass
 *  the rest of the sequences.
 */
/*
 *  For the moment, we are only interested in '*' (graphics)
 *  '(' and ')' (downloads).  Although we do not do anything
 *  with downloads, we need to pass the binary data thru
 *  untouched.
 *  Oops, '&' is handled too.
 */
    if ( parameter != '*' && parameter != '('
         && parameter != ')' && parameter != '&' ) {
      Flush_To_Term(fp);                /* flush rest of seq. */
      continue;
    }
/*
 *  Parameter character is in range, look for a valid group char
 */
    group_char = getc(fp);
    if ( group_char == EOF )    /* oops, ran out of input */
    {
      fclose( fp );
      return FALSE;
    }
/*
 *  See if in proper range.  If it isn't, it is not an error
 *  because the group character is optional for some sequences.
 *  For the moment, we are not interested in those sequences,
 *  so pass them thru.
 */
    if ( group_char < '`' || group_char > '~' ) {
/*
 *  If the "stripper" is active, we need to suspend
 *  it till graphics are re-started.
 */
      if ( group_char < '@' || group_char > '^' )
        Flush_To_Term(fp);      /* pass rest of seq. */
      continue;
    }
/*
*  Now we have a valid group character, decide if we want
*  to deal with this escape sequence.
*
*  Sequences we want do deal with include:
*
*    <esc>*r    ** graphics
*    <esc>*b    ** graphics
*    <esc>*v    ** graphics
*
*  Sequences we must pass thru binary data:
*
*    <esc>*c    ** pattern
*    <esc>*m    ** download dither
*    <esc>*t    ** obsolete
*    <esc>(f    ** download char set
*    <esc>(s    ** download char
*    <esc>)s    ** download font
*    <esc>&a    ** logical page
*    <esc>&b    ** AppleTalk stuff
*    <esc>&l    ** obsolete
*
*/
     if (  ( parameter == '*'
        && group_char != 'r' && group_char != 'b'
        && group_char != 'v' && group_char != 'c'
        && group_char != 't' && group_char != 'm' )
        || ( parameter == '&'
        && group_char != 'a' && group_char != 'l'
        && group_char != 'b' )
        || ( parameter == '('
        && group_char != 'f' && group_char != 's' )
        || ( parameter == ')' && group_char != 's' ) )
     {
/*
*  Definately not interested in the sequence.
*/
        Flush_To_Term(fp);
        continue;
      }


/*
 *  If the sequence is <esc>&a#H, it will have gotten past
 *  the above, but we need to suspend the "stripper" if
 *  it is active, because the CAP is getting moved.
 *
 *  The <esc>*p#X/Y sequences will have been filtered
 *  thru just above (<esc>*p is not a needed group).
 */

/*
 *  Now set up a pass thru flag so we can ignore the entire
 *  sequences of some of these.
 */
      if ( parameter != '*' )
        pass_seq = TRUE;
      else if ( group_char == 'c' || group_char == 't' || group_char == 'm' )
        pass_seq = TRUE;
      else
        pass_seq = FALSE;

/*
 *  Now we have a sequence that we are definately interested in.
 *
 *  Get the value field and terminator, and loop until final
 *  terminator is found.
 */
/* first see if the value has a plus sign */
        scanf_count = fscanf(fp, " + %d", &value );
        if ( scanf_count == 1 )
          plus_sign = TRUE;
        else
        {
          plus_sign = FALSE;
          scanf_count = fscanf(fp, " %d", &value );
          if ( scanf_count == 0 )
          value = 0;            /* by default */
        }
/*
 *  I wonder if I will get bitten by a trailing
 *  space character right here?
 */
        terminator = getc(fp);
/*
 *  Check for a fractional component.
 */
        fvalue = 0.0;
        if ( terminator == '.' ) {
          fvalue = Get_Frac(fp);
/*
 *  Now get real terminator.
 */
        terminator = getc(fp);
      }

      if ( terminator == EOF ) {
         fclose( fp );
         return (FALSE);
      }
/*
 *  If the pass_seq flag is set, then just pass
 *  it thru to stdout until a 'W' is found.
 */
      if ( pass_seq ) {
        if ( !in_sequence ) in_sequence = TRUE;

/*
 *  See if there was a non-zero fraction.
 */
        if ( fvalue != 0.0 ) {
          if ( value < 0 ) {
            value = -value;
          }
          fvalue += value;
/* if binary data, pass it thru */

          if ( terminator == 'W' ) {
            in_sequence = FALSE;        /* terminates */
            Flush_Bytes ( value, fp );  /* pass data */
          }

          continue;
         }
        }
/*
 * if we had gone so far, a big chance that it is a graphic PCL
 */
        fclose(fp);
        return (TRUE);
    }
    fclose(fp);
    return (FALSE);
}

/*
**  Flush_To_Term() simply passes thru input until a valid terminator
**  character is found.  This is for unwanted escape sequences.
*/

static void Flush_To_Term(fp)
FILE *fp;
{
        int     c;

        do
        {
                c = getc(fp);

                if ( c == EOF )                 /* this is a problem */
                        return;
        } while ( c < '@' || c > '^' );

}


/*
**  Flush_Bytes() simply transfers so many bytes directly from input to outpu
t.
**  This is used to pass thru binary data that we are not interested in so th
at
**  it will not confuse the parser.  I.e. downloads.
*/
static void Flush_Bytes( num, fp )
unsigned int    num;
FILE *fp;
{
   int  bnum;
   char buf[BUFSIZ];

   while ( num > 0 ) {
     bnum = MIN ( BUFSIZ, num );
     fread( buf, 1, bnum, fp );
     num -= bnum;
   }
}

/*
**  Get_Frac() simply gets the fractional part of a value.  This is here
**  because scanf() will consume a trailing 'e' or 'E', which is a problem
**  in PCL.
*/

static float    Get_Frac(fp)
FILE *fp;
{
        float   result = 0.0;
        int     c;
        float   position = 10.0;

        while ( (c = getc(fp)) != EOF )
        {
                /*
                **  Do we have a digit?
                */

                if ( !isdigit(c) )              /* not a digit */
                {
                        ungetc( c, fp );     /* put it back */
                        break;                  /* quit */
                }

                result += ((c - '0') / position);

                position *= 10.0;
        }

        return ( result );
}

/*
 *name: ReadFileType
 *
 *source: XV-3.0a
 */
static int ReadFileType(fname)
     char *fname;
{
  /* opens fname (which *better* be an actual file by this point!) and
     reads the first couple o' bytes.  Figures out what the file's likely
     to be, and returns the appropriate *** code */


  FILE *fp;
  unsigned char  magicno[8];    /* first 8 bytes of file */
  int   rv;


  rv = WKF_INCONNU;
  if (!fname) return rv;   /* shouldn't happen */

  fp = fopen(fname, "rb");

  if (!fp) return rv;

  rv = fread(magicno,8,1,fp);
  fclose(fp);

  if (rv!=1) return WKF_INCONNU;    /* files less than 8 bytes long... */

  rv = WKF_INCONNU;
  if (strncmp((char *) magicno,"GIF87a",6)==0) rv = WKF_GIF87;
 
  else if (strncmp((char *) magicno,"GIF89a",6)==0) rv = WKF_GIF89;

  else if (strncmp((char *) magicno,"VIEW",4)==0 ||
           strncmp((char *) magicno,"WEIV",4)==0) rv = WKF_PM;

  else if (magicno[0] == 'P' && magicno[1]>='1' &&
           magicno[1]<='6') rv = WKF_PBM;

  else if (strncmp((char *) magicno,"#define",7)==0) rv = WKF_XBM;

  else if (magicno[0]==0x59 && (magicno[1]&0x7f)==0x26 &&
           magicno[2]==0x6a && (magicno[3]&0x7f)==0x15) rv = WKF_SUNRAS;

  else if (magicno[0] == 'B' && magicno[1] == 'M') rv = WKF_BMP;

  else if (magicno[0]==0x52 && magicno[1]==0xcc) rv = WKF_UTAHRLE;

  else if ((magicno[0]==0x01 && magicno[1]==0xda) ||
           (magicno[0]==0xda && magicno[1]==0x01)) rv = WKF_IRIS;

  else if (magicno[0]==0x1f && magicno[1]==0x9d) rv = WKF_COMPRESS;

  else if (magicno[0]==0x0a && magicno[1] <= 5) rv = WKF_PCX;

  else if (magicno[0]==0xff && magicno[1]==0xd8 &&
           magicno[2]==0xff) rv = WKF_JPG;

  else if ((magicno[0]=='M' && magicno[1]=='M') ||
           (magicno[0]=='I' && magicno[1]=='I')) rv = WKF_TIFF;

  else if (strncmp((char *) magicno,  "NJPL1I00",8)==0 || /* fixed-len pds */
           strncmp((char *) magicno+2,"NJPL1I",  6)==0 || /* vger+other pds *
/
           strncmp((char *) magicno,  "CCSD3ZF", 7)==0 || /* vikng pds browse
*/
           strncmp((char *) magicno+2,"CCSD3Z",  6)==0 || /* vik. huffman pds
*/
           strncmp((char *) magicno,  "LBLSIZE=",8)==0)   /* vicar */
      rv = WKF_PDSVICAR;

  else if (magicno[0] == '%' && magicno[1] == '!') rv = WKF_PS;

  return rv;
}


