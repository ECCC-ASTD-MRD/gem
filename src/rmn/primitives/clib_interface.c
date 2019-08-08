/* RMNLIB - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2007  Environnement Canada
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
/*
!-------------------------------------------------------------------
! Modifications: [Date,Who,What]
! 2004-04, Stephane Chamberland
!    Original Code
! 2004-09, Stephane Chamberland
!    Add clib_glob
! 2004-09, Martin Serrer
!    Add WIN32 quirks
! 2007-08, Stephane Chamberland
!    replace dependance on cnf [f77.h] for RMNLIB's new/revamped Fortran interface
!    replace heap-alloc by stack ones, bug fixing
!    added functions of ctype.h
! 2007/09, Michel Valin
!    minor code touchups, replaced basename/dirname to get rid of extra lib needed on IRIX
! 2008/03, Michel Valin
!    added clib_readlink
!    added clib_stat, clib_size, clib_mtime functions to return more stat() information
! 2011/04, Stephane Chamberland
!    added clib_mkdir_r
! 2011/04, Mario Lepine
!    replace 0777 by 0755 for clib_mkdir and clib_mkdir_r function calls
!-------------------------------------------------------------------
! Dependencies:
! This 'package' make use of RPN' FTN2C package
! to link Fortran and C s/r
!
! The interfaces to fortran code are defined in clib_interface.cdk
!-------------------------------------------------------------------
! Description
! This module is an interface for Fortran to some C STD fonctions 
! Description of C STD functions can be found at:
! http://www.opengroup.org/onlinepubs/007908799/headix.html
! 
! stdlib.h: standard library definitions 
!   char *getenv(const char *name);
!   int putenv(char *string);
!   char *realpath(const char *file_name, char *resolved_name);
!   #int system(const char *command); 
!      Use F90 command instead
!   #void qsort(void *base, size_t nel, size_t width, int (*compar)(const void *, const void *));
! stdio.h: standard buffered input/output
!   int remove(const char *path); 
!   int rename(const char *old, const char *new);
! unistd.h: standard symbolic constants and types
!   #int access(const char *path, int amode);
!   #replaced by:
!     clib_fileexist : true if path exist
!     clib_isreadok  : true if path is readable
!     clib_iswriteok : true if path is writable
!     clib_isexecok  : true is path is executable/searchable
!   int chdir(const char *path);
!   char *getcwd(char *buf, size_t size);
!   #use getcwd because getwd is considerd unsafe by the compiler
!   #char *getwd(char *path_name);
!   int getuid();
!   #int getopt(int argc, char * const argv[], const char *optstring);
!   #int link(const char *path1, const char *path2);
!   int rmdir(const char *path);
!   int symlink(const char *path1, const char *path2);
!   int unlink(const char *path);
!
! sys/stat.h - data returned by the stat() function
!   int mkdir(const char *path, mode_t mode);
!   #int stat(const char *path, struct stat *buf);
!   #replaced by:
!     clib_isdir
!     clib_islink
!     clib_ispipe
!     clib_isfile
!     clib_mtime
!     clib_size
!     clib_stat
!   # might develep on top with filetype from rmnlib
!
! libgen.h: http://www.opengroup.org/onlinepubs/007908799/xsh/libgen.h.html
!   char *basename(char *path);  this one is emulated, because of IRIX
!   char *dirname(char *path);  this one is emulated, because of IRIX
! 
! glob.h: http://www.opengroup.org/onlinepubs/007908799/xsh/glob.h.html
!   int glob(const char *pattern, int flags,
!            int(*errfunc)(const char *epath, int errno), 
!            glob_t *pglob);
!   void globfree(glob_t *pglob);
!
! ctype.h
!   int toupper(int c);
!   int tolower(int c);
!   int isalnum(int c);
!   int isalpha(int c);
!   int isblank(int c);
!   int isdigit(int c);
!   int islower(int c);
!   int ispunct(int c);
!   int isspace(int c);
!   int isupper(int c);
!   int isxdigit(int c);
!
! [Pending functions]
! #time.h
! #dirent.h
! #ftw.h - file tree traversal
! #sys/resource.h - definitions for XSI resource operations
! #ulimit.h - ulimit commands
!-------------------------------------------------------------------
! PRIVATE FN: (see below)
! PUBLIC FN: (see below)
!===================================================================*/
/* be ready for files > 2 GB by using long offsets */
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <sys/types.h> /* for mkdir & stat */
#include <libgen.h>
// #include <glob.h>
#include <sys/param.h> /* for MAXPATHLEN = PATH_MAX */
#include <alloca.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include <rmnlib.h>

#define CLIB_OK    1
#define CLIB_ERROR -1
#define CLIB_SUFFIX _schhide

#define CLIB_F77NAME(a) f77_name(a##_schhide)

/* ================================================================
 * Return in value, a string containing the value of env-var name
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_getenv)(F77_CHARACTER *name, F77_CHARACTER *value
                                      HIDDENLEN(name) HIDDENLEN(value) ) {
  char *name_c;
  char *value_c;
  char *defStr = " ";
  int name_c_len;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  name_c_len = 1 + F77STRLEN(name);
  name_c = (char *)alloca((size_t)(name_c_len*sizeof(char)));
  if (!name_c ||
      FTN2C_FSTR2CSTR(name,name_c,F77STRLEN(name),name_c_len) < 0) {
    return(CLIB_ERROR);
  }

  /* Call C function */
  value_c = getenv(name_c);

  /* Translate Back to Fortran strings */
  status = CLIB_ERROR;
  if (value_c &&
      FTN2C_CSTR2FSTR(value_c,value,strlen(value_c)+1,F77STRLEN(value))>=0) {
    status = CLIB_OK;
  } else {
    FTN2C_CSTR2FSTR(defStr,value,2,F77STRLEN(value));
  }
  return(status);
}

/* ================================================================
 * Set an env-var to the name_value value
 * name_value string is of the form "name=value"
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_putenv)(F77_CHARACTER *name_value
                                      HIDDENLEN(name_value)) {
  char *name_value_c;
  int name_c_len;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  name_c_len = 1 + F77STRLEN(name_value);
  name_value_c = (char *)malloc((size_t)(name_c_len*sizeof(char)));
  if (!name_value_c || FTN2C_FSTR2CSTR(name_value,name_value_c,
                                       F77STRLEN(name_value),name_c_len) < 0) {
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (putenv(name_value_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
 /* WARNING: the name_value_c point should not be freed; putenv copy only the pointer, not its content */
  return(status);
}

/* ================================================================
 * The realpath() function derives, from the pathname pointed 
 * to by file_name, an absolute pathname that names the same file, 
 * whose resolution does not involve ".", "..", or symbolic links
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_realpath)(F77_CHARACTER *fnamein,
                                        F77_CHARACTER *fnameout
                                        HIDDENLEN(fnamein) HIDDENLEN(fnameout) ) {
  char fnamein_c[MAXPATHLEN];
  char *fnameout_c;
  char fnameout2_c[MAXPATHLEN];
  char *defStr = " ";
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(fnamein,fnamein_c,F77STRLEN(fnamein),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  fnameout_c = realpath(fnamein_c,fnameout2_c);

  /* Translate Back to Fortran strings */
  status = CLIB_ERROR;
  if (fnameout_c &&
      FTN2C_CSTR2FSTR(fnameout_c,fnameout,
                      strlen(fnameout_c)+1,F77STRLEN(fnameout)) >= 0) {
      status = CLIB_OK;
  } else {
    FTN2C_CSTR2FSTR(defStr,fnameout,2,F77STRLEN(fnameout));
  }
  return(status);
}

/* ================================================================
 * The readlink() function gets, from the pathname pointed 
 * the contents of the link 
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_readlink)(F77_CHARACTER *fnamein,
                                        F77_CHARACTER *fnameout
                                        HIDDENLEN(fnamein) HIDDENLEN(fnameout) ) {
  char fnamein_c[MAXPATHLEN];
  char fnameout2_c[MAXPATHLEN];
  char *defStr = " ";
  F77_INTEGER status;
  ssize_t nc;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(fnamein,fnamein_c,F77STRLEN(fnamein),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  /* Call C function */
  fnameout2_c[0]='\0';
  nc=readlink(fnamein_c,fnameout2_c,MAXPATHLEN-1) ;
  if(nc<0) perror("clib_readlink");
  if(nc>0)fnameout2_c[nc]='\0';
  /* Translate Back to Fortran strings */
  status = CLIB_ERROR;
  if ( (nc > 0) && 
        (FTN2C_CSTR2FSTR(fnameout2_c,fnameout,nc,F77STRLEN(fnameout)) >= 0)  ) {
      status = CLIB_OK;
  } else {
    FTN2C_CSTR2FSTR(defStr,fnameout,2,F77STRLEN(fnameout));
  }
  return(status);
}

/* ================================================================
 * If path does not name a directory, remove(path) is equivalent to unlink(path).*
 * If path names a directory, remove(path) is equivalent to rmdir(path). 
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_remove)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (remove(path_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Rename a file named pathold to pathnew
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_rename)(F77_CHARACTER *pathold, 
                                      F77_CHARACTER *pathnew 
                                      HIDDENLEN(pathold) HIDDENLEN(pathnew) ) {
  char pathold_c[MAXPATHLEN];
  char pathnew_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(pathold,pathold_c,F77STRLEN(pathold),MAXPATHLEN) < 0 ||
      FTN2C_FSTR2CSTR(pathnew,pathnew_c,F77STRLEN(pathnew),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (rename(pathold_c,pathnew_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Checks the file named by the pathname pointed to by the 
 * path argument for accessibility.
 * clib_fileexist : CLIB_OK if path exist
 * clib_isreadok  : CLIB_OK if path is readable
 * clib_iswriteok : CLIB_OK if path is writable
 * clib_isexecok  : CLIB_OK if path is executable/searchable
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_fileexist)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (access(path_c,F_OK)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

F77_INTEGER CLIB_F77NAME(clib_isreadok)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (access(path_c,R_OK)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

F77_INTEGER CLIB_F77NAME(clib_iswriteok)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (access(path_c,W_OK)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

F77_INTEGER CLIB_F77NAME(clib_isexecok)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (access(path_c,X_OK)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Change working directory
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_chdir)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (chdir(path_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Get the current working directory pathname
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_getcwd)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  char *defStr = " ";
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  status = CLIB_ERROR;

  /* Call C function */
  if (getcwd(path_c,(size_t)MAXPATHLEN*sizeof(char)) &&
      FTN2C_CSTR2FSTR(path_c,path,MAXPATHLEN,F77STRLEN(path)) >=0) {
    status = CLIB_OK;
  } else {
    FTN2C_CSTR2FSTR(defStr,path,2,F77STRLEN(path));
  }
  return(status);
}

/* ================================================================
 * Get uid
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_getuid)(F77_INTEGER *uid ) {
  F77_INTEGER status=CLIB_OK;
  /*--------------------------------------------------------------*/
  *uid = getuid();
  return(status);
}

/* ================================================================
 * removes a directory whose name is given by path. 
 * The directory is removed only if it is an empty directory.
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_rmdir)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (rmdir(path_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Creates a symbolic link. Its name is the pathname pointed to 
 * by path2, which must be a pathname that does not name an existing
 * file or symbolic link. The contents of the symbolic link are 
 * the string pointed to by path1.
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_symlink)(F77_CHARACTER *pathold, 
                                       F77_CHARACTER *pathnew 
                                       HIDDENLEN(pathold) HIDDENLEN(pathnew) ) {
  char pathold_c[MAXPATHLEN];
  char pathnew_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(pathold,pathold_c,F77STRLEN(pathold),MAXPATHLEN) < 0 ||
      FTN2C_FSTR2CSTR(pathnew,pathnew_c,F77STRLEN(pathnew),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (symlink(pathold_c,pathnew_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Removes a link to a file. If path names a symbolic link, unlink() 
 * removes the symbolic link named by path and does not affect any 
 * file or directory named by the contents of the symbolic link. 
 * Otherwise, unlink() removes the link named by the pathname 
 * pointed to by path and decrements the link count of the file 
 * referenced by the link.
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_unlink)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (unlink(path_c)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Make a directory
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_mkdir)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  /* mode_t mode; */
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (mkdir(path_c,(mode_t)0755)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}

/* ================================================================
 * Interface to sys/stat.h:stat()
 * Split in multiple fn to avoid defining marcos/parameters
 * clib_isdir
 * clib_islink
 * clib_ispipe
 * clib_isfile
 * clib_size
 * clib_mtime
 * clib_stat (raw output)
 * ================================================================
 * Return CLIB_OK if path is a dir, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_isdir)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if (!stat(path_c,buf) && S_ISDIR(buf->st_mode)) status = CLIB_OK;
  return(status);
}

/* ================================================================
 *
 * Return CLIB_OK if path is a link, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_islink)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if (!lstat(path_c,buf) && S_ISLNK(buf->st_mode)) status = CLIB_OK;
  return(status);
}

/* ================================================================
 *
 * Return CLIB_OK if path is fifo, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_isfifo)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if (!stat(path_c,buf) && S_ISFIFO(buf->st_mode)) status = CLIB_OK;
  return(status);
}

/* ================================================================
 *
 * Return CLIB_OK if path is a file, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_isfile)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if (!stat(path_c,buf) && S_ISREG(buf->st_mode)) status = CLIB_OK;
  return(status);
}

/* ================================================================
 *
 * Return file size if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER8 CLIB_F77NAME(clib_size)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if (!stat(path_c,buf)) status = buf->st_size;
  return(status);
}

/* ================================================================
 *
 * Return time of last modification if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_mtime)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if (!stat(path_c,buf)) status = buf->st_mtime;
  return(status);
}

/* ================================================================
 *
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_stat)(F77_CHARACTER *path, INT_64 *table HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  struct stat buf1;
  struct stat *buf;
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }
  buf = &buf1;

  /* Call C function */
  status = CLIB_ERROR;
  if ( !stat(path_c,buf) ) status = CLIB_OK;
  table[0]=buf->st_dev;
  table[1]=buf->st_ino;
  table[2]=buf->st_mode;
  table[3]=buf->st_nlink;
  table[4]=buf->st_uid;   /* uid of file owner */
  table[5]=buf->st_gid;   /* gid of file owner */
  table[6]=buf->st_rdev;
  table[7]=buf->st_size;  /* file size */
  table[8]=buf->st_blksize;
  table[9]=buf->st_blocks;
  table[10]=buf->st_atime;  /* time of last access */
  table[11]=buf->st_mtime;  /* time of last modification */
  table[12]=buf->st_ctime;  /* time of last status modification */
  return(status);
}


/* ================================================================
 * The basename() function takes the pathname pointed to by path 
 * and returns a pointer to the final component of the pathname, 
 * deleting any trailing '/' characters.
 * If the string consists entirely of the '/' character, basename() 
 * returns a pointer to the string "/" .
 * If path is a null pointer or points to an empty string, 
 * basename() returns a pointer to the string "." . 
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_basename)(F77_CHARACTER *path,
                                        F77_CHARACTER *mybasename
                                        HIDDENLEN(path) HIDDENLEN(mybasename) ) {
  char *defStr = " ";
  F77_INTEGER status;  
  int lpath=F77STRLEN(path)-1;

  FTN2C_CSTR2FSTR(defStr,mybasename,1,F77STRLEN(mybasename));     /* fill destination with blanks, in case ... */
  mybasename[0] = '/';
  if( lpath == 1 && path[0] == '/' ) return(CLIB_OK);        /* in case the path string only contains the / character */
  while(lpath>0 && path[lpath]!='/') lpath--;                     /* scan backwards to find last / in path */
  status = FTN2C_CSTR2FSTR(path+lpath+1,mybasename,F77STRLEN(path)-lpath-1,F77STRLEN(mybasename));  /* copy result */
  return( status<0 ? CLIB_ERROR : CLIB_OK );
}

/* ================================================================
 * The dirname() function takes a pointer to a character string 
 * that contains a pathname, and returns a pointer to a string 
 * that is a pathname of the parent directory of that file. 
 * Trailing '/' characters in the path are not counted as part 
 * of the path.
 * If path does not contain a '/', then dirname() returns a 
 * pointer to the string "." . If path is a null pointer or 
 * points to an empty string, dirname() returns a pointer to 
 * the string "." . 
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_dirname)(F77_CHARACTER *path,
                                   F77_CHARACTER *mydirname
                                   HIDDENLEN(path) HIDDENLEN(mydirname) ) {
  char *defStr = " ";
  F77_INTEGER status;
  int lpath=F77STRLEN(path)-1;

  FTN2C_CSTR2FSTR(defStr,mydirname,1,F77STRLEN(mydirname));     /* fill destination with blanks, in case ... */
  while(lpath>0 && path[lpath]!='/') lpath--;                   /* scan backwards to find last / in path */
  if(lpath==0) {
    if(path[0] != '/') {
      mydirname[0] = '.'; /* no / found, return . as dirname */
      return(CLIB_OK) ;
    } else {              /* path starts with / and has no other / */
      return(CLIB_ERROR);
    }
  }
  status = FTN2C_CSTR2FSTR(path,mydirname,lpath,F77STRLEN(mydirname));  /* copy result */
  return( status<0 ? CLIB_ERROR : CLIB_OK );
}
#if defined(NO_CLIB_INTERFACE2)
/* ================================================================
 * Get the list of files in PWD that match a pattern)
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_glob)(F77_CHARACTER *filelist,
                                    F77_INTEGER *nfiles,
                                    F77_CHARACTER *pattern,
                                    F77_INTEGER *maxnfiles
                                    HIDDENLEN(filelist) HIDDENLEN(pattern) ) {
  glob_t globbuf;
  char pattern_c[MAXPATHLEN];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(pattern,pattern_c,F77STRLEN(pattern),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  *nfiles = 0;
  status = CLIB_ERROR;
  if (!glob(pattern_c, GLOB_NOSORT, NULL, &globbuf)) {
    if ((F77_INTEGER)globbuf.gl_pathc <= *maxnfiles) {
      *nfiles = (F77_INTEGER)globbuf.gl_pathc;
      if (FTN2C_CSTR2FSTR_A(globbuf.gl_pathv,filelist,
                            MAXPATHLEN,(int)F77STRLEN(filelist),
                            (int)*nfiles) >= 0) 
                                 status = CLIB_OK;
    }
  }
  globfree(&globbuf);
  return(status);
}
#endif
/* ================================================================
 * Convert Char array to lowercase/uppercase
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_tolower)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  int ii;
  /*--------------------------------------------------------------*/
  for (ii=0 ; ii<F77STRLEN(mystr) ; ii++) mystr[ii] = (F77_CHARACTER)tolower((int)mystr[ii]);
  return((F77_INTEGER)CLIB_OK);
}

F77_INTEGER CLIB_F77NAME(clib_toupper)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  int ii;
  /*--------------------------------------------------------------*/
  for (ii=0 ; ii<F77STRLEN(mystr) ; ii++) mystr[ii] = (F77_CHARACTER)toupper((int)mystr[ii]);
  return((F77_INTEGER)CLIB_OK);
}

/* ================================================================
 * Check Char type
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_isalnum)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isalnum((int)mystr[0]))? (F77_INTEGER)CLIB_OK : (F77_INTEGER)CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_isalpha)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isalpha((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_isblank)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isblank((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_isdigit)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isdigit((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_islower)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (islower((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_ispunct)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (ispunct((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_isspace)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isspace((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_isupper)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isupper((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

F77_INTEGER CLIB_F77NAME(clib_isxdigit)(F77_CHARACTER *mystr HIDDENLEN(mystr) ) {
  return (isxdigit((int)mystr[0]))? CLIB_OK : CLIB_ERROR;
}

/* Function with behaviour like `mkdir -p' 
Source: http://niallohiggins.com/2009/01/08/mkpath-mkdir-p-alike-in-c-for-unix/
*/
int mkpath(char *s, mode_t mode){
        char *q, *r = NULL, *path = NULL, *up = NULL;
        int rv;

        rv = -1;
        if (strcmp(s, ".") == 0 || strcmp(s, "/") == 0)
                return (0);

        if ((path = strdup(s)) == NULL)
                exit(1);
     
        if ((q = strdup(s)) == NULL)
                exit(1);

        if ((r = dirname(q)) == NULL)
                goto out;
        
        if ((up = strdup(r)) == NULL)
                exit(1);

        if ((mkpath(up, mode) == -1) && (errno != EEXIST))
                goto out;

        if ((mkdir(path, mode) == -1) && (errno != EEXIST))
                rv = -1;
        else
                rv = 0;

out:
        if (up != NULL)
                free(up);
        free(q);
        free(path);

        return (rv);
}

/* ================================================================
 * Make a directory, recursively
 * Return CLIB_OK if success, CLIB_ERROR if not
 * ================================================================*/
F77_INTEGER CLIB_F77NAME(clib_mkdir_r)(F77_CHARACTER *path HIDDENLEN(path) ) {
  char path_c[MAXPATHLEN];
  /* mode_t mode; */
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(path,path_c,F77STRLEN(path),MAXPATHLEN) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  if (mkpath((char *) &path_c,(mode_t)0755)) {
    status = CLIB_ERROR;
  } else {
    status = CLIB_OK;
  }
  return(status);
}
