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
#if !defined (NEC)
      subroutine stkmemw(nw,ptr)
      integer nw,ptr
      integer errcode

      call hpalloc(ptr,nw,errcode,0)
      return
      end
#else
#if !defined (_FLOAT1)
      subroutine stkmemw(nw,ptr)
      integer nw,ptr
      integer errcode

      print *,'ERROR: stkmem2 must be used instead of stkmemw'
      call tracebck()
      return
      end
#endif
#endif
      subroutine unstakw(ptr)
      integer ptr
#if !defined (NEC)
      integer ier
      call hpdeallc(ptr,ier,0)
#endif
      return
      end

