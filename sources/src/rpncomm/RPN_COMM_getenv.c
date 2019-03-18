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

/* FORTRAN routine MUST ensure null terminated string ( trim(string)//achar(0) )   */

/* fudged version of strncpy that blank fills destination */
static char*
strncpy_(char *dest, const char *src, size_t n){
    size_t i;

    for (i = 0 ; i < n && src[i] != '\0' ; i++)
	dest[i] = src[i];
    for ( ; i < n ; i++)
	dest[i] = ' ';

    return dest;
}

#pragma weak rpn_comm_getenv__=rpn_comm_getenv
#pragma weak rpn_comm_getenv_=rpn_comm_getenv
int rpn_comm_getenv__(char *,char *, int *);
int rpn_comm_getenv_(char *,char *, int *);
int rpn_comm_getenv(char *name,char *value, int *length)
{
  char *temp=getenv(name);
  if(temp != NULL) strncpy_(value,temp,(size_t) *length);
}
