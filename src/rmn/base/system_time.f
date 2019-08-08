*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
	subroutine system_time(yyyymmdd,hhmmss00)
	implicit none
	integer yyyymmdd, hhmmss00
	external c_time
	integer c_time
        integer, external :: newdate
        integer status

	integer newstamp,istamp,minutes,secs
	real *8 seconds,hours


*
*	build stamp (new) for jan 1 1980
*
	status = newdate(istamp,[19800101],00000000,3)
*
*	get number of seconds since jan 1 1970
*
	seconds=c_time()
*
*	make it relative to jan 1 1980 (because of new stamp format)
*
	hours=seconds/3600.0-87648.

	call incdatr(newstamp,istamp,hours)

	status = newdate(newstamp,[yyyymmdd],hhmmss00,-3)

	return
	end
