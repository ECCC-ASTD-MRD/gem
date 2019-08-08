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

#include <stdio.h>
#include <rpnmacros.h>

#ifdef WIN32	/*CHC/NRC*/
#include <Windows.h>
void f77name(micro_sleep) (double *secs)
{
	DWORD dwMilliseconds = 1000 * (*secs);
	Sleep(dwMilliseconds);
}
#else
#include <sys/types.h>
#include <sys/time.h>
#if !defined (HP)
#include <sys/select.h>
#endif
#include <math.h>
void f77name(micro_sleep) (double *secs)
{
  int sec,usec;
  struct timeval timeout;
  sec = *secs / 1;
  usec = fmod(*secs,1.0) * 1000000;
  /*  printf("Debug sec=%d usec=%d\n",sec,usec); */
  timeout.tv_sec=*secs / 1;
  timeout.tv_usec=fmod(*secs,1.0) * 1000000;
  select(0,NULL,NULL,NULL,&timeout);
}
#endif

