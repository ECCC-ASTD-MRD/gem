#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <rmnlib.h>

/* Cette fonction verifie si le repertoire $TMPDIR existe.
 * Si oui, la fonction retourne sans problemes
 * Si non, on fait un "putenv(TMPDIR=/tmp/)
 * */

int f77name(chk_tmpdir)(void)
  {
  return chk_tmpdir();
  }

int chk_tmpdir(void)
  {
  char tmpdir[512];
  struct stat le_buffer;
  int res;
  
  strcpy(tmpdir, "");
  strcat(tmpdir, (char *)getenv("TMPDIR"));
//   printf("%s\n", tmpdir);
  lstat(tmpdir, &le_buffer);
  res = S_ISDIR(le_buffer.st_mode);
//   printf("%d\n", res);
  if (res != 1)
   {
   fprintf(stderr, "    **************************************************************\n");
   fprintf(stderr, "        Warning : There is a problem with TMPDIR\n");
   fprintf(stderr, "        TMPDIR is %s\n", tmpdir);
   fprintf(stderr, "        Setting TMPDIR=/tmp\n");
   fprintf(stderr, "    **************************************************************\n");
   strcpy(tmpdir, "/tmp/");
   res = putenv("TMPDIR=/tmp");
   if (res != 0)
    {
   fprintf(stderr, "    **************************************************************\n");
    fprintf(stderr, "       Warning : cannot set TMPDIR to /tmp\n");
    fprintf(stderr, "       Program will now exit\n");
   fprintf(stderr, "    **************************************************************\n");
    exit(-13);
    }
   res = strcat(tmpdir, (char *)getenv("TMPDIR"));
//    printf("%s\n", tmpdir);
   }
  return 0;
  }
  
  #include <stdio.h>
#include <rpnmacros.h>

#define NOMVAR 1
#define TYPVAR 2
#define ETIKET 3
#define IP01    4
#define IP02    5
#define IP03    6
#define DATEO  7
#define DATEV  8
#define NI     9
#define NJ    10
#define NK    11
#define LAT   12
#define LON   13

#define NORD   1
#define SUD    2
#define EST    3
#define OUEST  4
#define NONE   5

FILE *ftnUnits[100];

/*
****
****
*/

void strconvdate(char strdate[], int fstdate);
static int c_dateform = 1;

int f77name(pgsmform)(char format[],int *nrepeats, int *lenFormat, F2Cl fortranLenFormat)
{
   char cFormat[32], cNrepeats[8];
   int i, mode, yyyymmhh, hhmmssss;
   
   *nrepeats = 0;

   for (i=0; i < 32; i++)
      {
      cFormat[i] = '\0';
      }

   i = *lenFormat-1;
   while (i >= 0 && format[i] == ' ')
      {
      format[i] = '\0';
      i--;
      }
   
   i = 0;
   while (isdigit((char)format[i]))
     {
     cNrepeats[i] = format[i];
     i++;
     }
   cNrepeats[i] = '\0';

   if (i > 0)
     {
     sscanf(cNrepeats, "%d", nrepeats);
     strcpy(format, &(format[i]));
     }
     
   cFormat[0] = '%';
   strcpy(&cFormat[1],&format[1]);
   
   cFormat[strlen(cFormat)] = (char) tolower(format[0]);
   strcpy(format, cFormat);


     
   
}

/*
****
****
*/

int f77name(pgsmecho)(int *iun, char message[],int *lenMessage, F2Cl fortranLenMessage)
{
   message[*lenMessage] = '\0';
   fprintf(ftnUnits[*iun],"%s\n",message);
}

/*
****
****
*/

int f77name(pgsmof)(int *iun, char *nomFichier,F2Cl lenNomFichier)
{
   int i;

   i = 0;
   while (nomFichier[i] != ' ')
      {
      i++;
      }
   nomFichier[i] = (char)NULL;
      
   ftnUnits[*iun] = fopen(nomFichier, "w");
   if (ftnUnits[*iun] == NULL)
      {
      fprintf(stderr, "Impossible d'ouvrir le fichier: %s\n", nomFichier);
      return -1;
      }
   
   return 0;
   
}

/*
****
****
*/

int f77name(pgsmcf)(int *iun)
{
   fclose(ftnUnits[*iun]);
}

/*
****
****
*/

int f77name(pgsmwr)(int *iun,float *data,int *ni, int *nj, int *nk ,char *format,int *position,int *idents,char *separateur,
             char *nomvar,char *typvar,char *etiket,int *dateo,int *datev,int *dateform, int *ip1,int *ip2,int *ip3,
             float *lat,float *lon, F2Cl len_format, F2Cl len_separateur, F2Cl len_nomvar, F2Cl len_typvar, F2Cl len_etiket)
{
   int i,j;
   char c_etiket[16],c_typvar[4],c_nomvar[8], c_separateur[2];
   char string[256];
   char internalFormat[16],latlonformat[32];
   int latlonflag = 0;
   int longform=16;
   int npts = *ni * *nj * *nk;
   int nrepeats;

   c_dateform = *dateform;
   
   strncpy(c_nomvar,nomvar,4);
   strncpy(c_typvar,typvar,2);
   strncpy(c_etiket,etiket,12);
   strncpy(internalFormat,format,16);
   c_separateur[0] = separateur[0];

   if (c_separateur[0] == 'T')
      {
      c_separateur[0] = '\t';
      }

   c_nomvar[4] = '\0';
   c_typvar[2] = '\0';
   c_etiket[12] = '\0';
   c_separateur[1] = '\0';

   f77name(pgsmform)(internalFormat,&nrepeats, &longform,(F2Cl) 16);
   sprintf(latlonformat,"%s%s%s%s","%s",internalFormat,"%s",internalFormat);
   i = 0;
   while (i < 16)
      {
      if (idents[i] == LAT || idents[i] == LON) 
         {
         latlonflag = 1;
         }
      i++;
      }
   
   for (j=0; j < npts; j++)
      {
      if (*position == NORD)
         {
         if (j == 0)
            {
            ImprimeIdent(string,idents,c_separateur,c_nomvar,c_typvar,c_etiket,*ip1,*ip2,*ip3,*dateo,*datev, *ni, *nj, *nk);
            string[strlen(string)-1] = '\n';
            fprintf(ftnUnits[*iun],"%s",string);
            }
         }

      if (*position == OUEST)
         {
         if (latlonflag)
            {
            ImprimeIdent(string,idents,c_separateur,c_nomvar,c_typvar,c_etiket,*ip1,*ip2,*ip3,*dateo,*datev, *ni, *nj, *nk);
            fprintf(ftnUnits[*iun],"%s%s",string, c_separateur);
            }
         else
            {
            if (j == 0)
               {
               ImprimeIdent(string,idents,c_separateur,c_nomvar,c_typvar,c_etiket,*ip1,*ip2,*ip3,*dateo,*datev, *ni, *nj, *nk);
               fprintf(ftnUnits[*iun],"%s%s",string, c_separateur);
               }
            }
         }

      fprintf(ftnUnits[*iun], internalFormat, data[j]);
      
      if (latlonflag)
         {
         fprintf(ftnUnits[*iun], latlonformat, c_separateur, lat[j],  c_separateur, lon[j]);
         }
      else
         {
         if (j < npts-1)
            {
	    if (nrepeats == 0)
	      {
	      fprintf(ftnUnits[*iun], "%s", c_separateur);
	      }
	    else
	      {
	      if (0 == ((j+1) % nrepeats)) fprintf(ftnUnits[*iun], "\n");
	      else fprintf(ftnUnits[*iun], "%s", c_separateur);
	      }
	    }
         }
      
      if (*position == NORD || *position == OUEST || *position == 0)
        {
        if (latlonflag || j == npts-1)
          fprintf(ftnUnits[*iun], "\n");
        }
      
      
      if (*position == SUD)
         {
         if (latlonflag) 
            {
            if (j < npts-1) fprintf(ftnUnits[*iun],"\n");
            }
         
         
         if (j == npts-1)
            {
            ImprimeIdent(string,idents,c_separateur,c_nomvar,c_typvar,c_etiket,*ip1,*ip2,*ip3,*dateo,*datev, *ni, *nj, *nk);
            string[strlen(string)-1] = '\n';
            fprintf(ftnUnits[*iun],"\n%s",string);
            }
         }

      if (*position == EST)
         {
         if (latlonflag || j == npts-1)
            {
            ImprimeIdent(string,idents,c_separateur,c_nomvar,c_typvar,c_etiket,*ip1,*ip2,*ip3,*dateo,*datev, *ni, *nj, *nk);
            fprintf(ftnUnits[*iun],"%s%s\n",c_separateur,string);
            }
         }

      }
   printf(" ECRIT -- %s \t %s \t %05d \t %04d \t %03d \t %s \t %08d\n", c_nomvar, c_typvar, *ip1, *ip2, *ip3, c_etiket, *dateo);
}


GetIdent(char string[],int item, char *nomvar, char *typvar, char *etiket, 
         int ip1, int ip2, int ip3, int dateo, int datev, int ni, int nj, int nk)
{
  int mode, yyyymmdd, hhmmssss;   
  switch(item)
      {
      case NOMVAR:
        strcpy(string, nomvar);
        break;

      case TYPVAR:
        strcpy(string, typvar);
        break;

      case ETIKET:
        strcpy(string, etiket);
        break;

      case IP01:
        sprintf(string, "%5d", ip1);
        break;

      case IP02:
        sprintf(string, "%4d", ip2);
        break;

      case IP03:
        sprintf(string, "%4d", ip3);
        break;

      case DATEO:
	strconvdate(string, dateo);
        break;

      case DATEV:
	strconvdate(string, datev);
        break;

      case NI:
        sprintf(string, "%03d", ni);
        break;

      case NJ:
        sprintf(string, "%03d", nj);
        break;

      case NK:
        sprintf(string, "%03d", nk);
        break;

      default:
        sprintf(string, "%s", "KABOOM!");
        break;

      }
}


ImprimeIdent(char longString[], int items[],char *separateur, char *nomvar, char *typvar, char *etiket, 
             int ip1, int ip2, int ip3, int dateo, int datev, int ni, int nj, int nk)
{
   int i;
   char string[32];
   
   i = 0;
   strcpy(longString,"");

   while (items[i] != 0 && i < 16)
      {
      switch (items[i])
         {
         case LAT:
         case LON:
           break;
           
         default:
           GetIdent(string,items[i],nomvar,typvar,etiket,ip1,ip2,ip3,dateo,datev,ni,nj,nk);
           strcat(string,separateur);
           strcat(longString,string);
           break;
         }
      i++; 
      }
}

void strconvdate(char strdate[], int fstdate)
{
  int lfstdate, yyyymmdd, hhmmssss, mode;
  int yyyy, month, day, hour, minutes, sec, fracsec;
  
  switch(c_dateform)
    {
    case 0:
      sprintf(strdate, "%08d", fstdate);
      break;

    case 1:
      mode = -3;
      lfstdate = fstdate;
      
      f77name(newdate)(&lfstdate, &yyyymmdd, &hhmmssss, &mode);
      sprintf(strdate, "%08d.%08d", yyyymmdd, hhmmssss);
      break;

    case 2:
      mode = -3;
      lfstdate = fstdate;
      
      f77name(newdate)(&lfstdate, &yyyymmdd, &hhmmssss, &mode);
      
      yyyy = yyyymmdd / 10000;
      month = (yyyymmdd / 100) % 100;
      day   = yyyymmdd % 100;
      
      hour = hhmmssss / 1000000;
      minutes = (hhmmssss / 10000) % 100;
      sec = (hhmmssss % 100);
      
      sprintf(strdate, "%04d-%02d-%02dT%02d:%02d:%02dZ", yyyy, month, day, hour, minutes, sec);
      break;
    }
}
