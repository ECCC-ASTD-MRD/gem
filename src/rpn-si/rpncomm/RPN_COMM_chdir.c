/* RPN_COMM - Library of useful routines for C and FORTRAN programming
 * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
 *                          Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

/* FORTRAN callable code, supports 3 popular name manglings (0,1,2 underscores at end of name) */
/* FORTRAN routine MUST ensure null terminated string ( trim(string)//achar(0) )   */

#pragma weak f_rpn_comm_chdir__=f_rpn_comm_chdir
#pragma weak f_rpn_comm_chdir_=f_rpn_comm_chdir
int f_rpn_comm_chdir__(char *);
int f_rpn_comm_chdir_(char *);
int f_rpn_comm_chdir(char *in_reper)
{
  int ierr;
  fprintf(stderr,"reper='%s' \n",in_reper);
  ierr = chdir(in_reper);

  return(ierr);
}
