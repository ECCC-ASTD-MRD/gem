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
	integer function gnthread()
	implicit none
	character(len=15) omp_var_S
	character(len=4) result_S
        integer resu2	

	omp_var_S='OMP_NUM_THREADS'

	call getenvc(omp_var_S,result_S)

	if(result_S == '') then
	   gnthread = 1
	   return
	endif
	read(result_S(1:3),10) gnthread
 10	format(I3)
        read(result_S(1:4),20) resu2
 20     format(I4)

        if(resu2 .ne. gnthread) then
	   write(*,*) 'GNTHREAD: Invalid value for OMP_NUM_THREADS'
	   gnthread = -1
	endif
	return

	end

