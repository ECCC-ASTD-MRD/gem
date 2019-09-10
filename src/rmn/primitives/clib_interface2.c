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
!===================================================================*/
#include <stdlib.h>
#include <limits.h>
#include <glob.h>

#include <rmnlib.h>

#define CLIB_OK    1
#define CLIB_ERROR -1
#define CLIB_SUFFIX _schhide

#define CLIB_F77NAME(a) f77_name(a##_schhide)
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
  char pattern_c[PATH_MAX];
  F77_INTEGER status;
  /*--------------------------------------------------------------*/
  /* Translate to C strings */
  if (FTN2C_FSTR2CSTR(pattern,pattern_c,F77STRLEN(pattern),PATH_MAX) < 0){
    return(CLIB_ERROR);
  }

  /* Call C function */
  *nfiles = 0;
  status = CLIB_ERROR;
  if (!glob(pattern_c, GLOB_NOSORT, NULL, &globbuf)) {
    if ((F77_INTEGER)globbuf.gl_pathc <= *maxnfiles) {
      *nfiles = (F77_INTEGER)globbuf.gl_pathc;
      if (FTN2C_CSTR2FSTR_A(globbuf.gl_pathv,filelist,
                            PATH_MAX,(int)F77STRLEN(filelist),
                            (int)*nfiles) >= 0)
                                 status = CLIB_OK;
    }
  }
  globfree(&globbuf);
  return(status);
}

