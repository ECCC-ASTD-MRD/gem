!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
      subroutine ez_rgoptc(op, val, flag)
      implicit none
      character *8 op, val
      logical flag

      integer ier,ezsetopt, ezgetopt

!flag = .true.  mode set
!flag = .false. mode get

#include "ez_qqqxtrp0.cdk"

      character*3 localop,localval

      localop =op(1:3)
      localval=val(1:3)

      call up2low(localop,localop)
      call up2low(localval,localval)

      if (flag) then
         if (localop.eq.  'ext') then
            if (localval .eq.  'oui') then
               codxtrap = oui
               ier = ezsetopt('extrap_degree', 'do_nothing')
            else if (localval.eq.  'abo') then
               codxtrap = abort
               ier = ezsetopt('extrap_degree', 'abort')
            else if (localval .eq.  'max') then
               codxtrap = maximum
               ier = ezsetopt('extrap_degree', 'maximum')
            else if (localval .eq.  'min') then
               codxtrap = minimum
               ier = ezsetopt('extrap_degree', 'minimum')
            else if (localval .eq.  'voi') then
               codxtrap = voisin
               ier = ezsetopt('extrap_degree', 'nearest')
            else if (localval .eq.  'val') then
               codxtrap = valeur
               ier = ezsetopt('extrap_degree', 'value')
            else 
               print *, ' <rgoptc>: mauvaise valeur pour val'
               print *, '           val = ', val
               print *, '           val initialisee a ''abort'''
               codxtrap = abort
               ier = ezsetopt('extrap_degree', 'abort')
            endif
         else if (localop.eq.'int') then
            if (localval .eq.  'voi') then
               ordint = voisin
               ier = ezsetopt('interp_degree', 'nearest')
            else if (localval.eq.  'lin') then
               ordint = lineair
               ier = ezsetopt('interp_degree', 'linear')
            else if (localval.eq.  'cub') then
               ordint = cubique
               ier = ezsetopt('interp_degree', 'cubic')
            else 
               print *, ' <rgoptc>: mauvaise valeur pour val'
               print *, '           val = ', val
               print *, '           val initialisee a ''cubique'''
               ordint = cubique
               ier = ezsetopt('interp_degree', 'cubic')
            endif
         else
            print *, ' <rgoptc>: mauvaise valeur pour op'
            print *,             '     op devrait etre egal a ''extrap'' ou ''interp'''
         endif
      else 
         if (localop .eq.  'ext') then
            ier = ezgetopt('extrap_degree', val)
         endif
         if (localop .eq.  'int') then
            ier = ezgetopt('interp_degree', val)
         endif
      endif
      return
      end

      subroutine ez_rgopti(op, val, flag)
      implicit none
      character*8 op
      integer val
      logical flag
      integer ier,ezsetval,ezgetval,ezsetopt,ezgetopt
      
!     flag = .true.  mode set
!     flag = .false. mode get
      
#include "ez_qqqxtrp0.cdk"
      
      data codxtrap / oui /
      data flgxtrap / .false. /
      data ordint   / 3 /
      real rval
      
      character*3 localop
      character*16 local_val
      
      localop=op(1:3)
      call up2low(localop,localop)
      
      if (flag) then
         if (localop .eq.  'ext') then
            if (val .eq.  voisin .or.  val .eq.  lineair             .or.  val .eq. cubique) then
               if (val.eq.100.or.val.eq.0) then
                  ier = ezsetopt('interp_degree','nearest')
               elseif (val.eq.1) then
                  ier = ezsetopt('interp_degree','linear')
               elseif (val.eq.3) then
                  ier = ezsetopt('interp_degree','cubic')
               endif
            else
               valxtrap = real(val)
            endif
            ier = ezsetval('extrap_value',valxtrap)
         else 
            if (localop .eq.  'int') then
               if (val.eq.100.or.val.eq.0) then
                  ier = ezsetopt('interp_degree','nearest')
               elseif (val.eq.1) then
                  ier = ezsetopt('interp_degree','linear')
               elseif (val.eq.3) then
                  ier = ezsetopt('interp_degree','cubic')
               else
                  print *, '<ez_rgopti> Erreur!'
               endif
            endif
         endif
      else
         if (localop.eq.'ext') then
            ier = ezgetval('extrap_value',rval)
            val = nint(rval)
         else if (localop .eq.  'int') then
            ier = ezgetopt('interp_degree', local_val)
            if (local_val.eq.'nearest') then
               val = 0
            elseif (local_val.eq.'linear') then
               val = 1
            else
               val = 3
            endif
         endif
      endif
      return
      end
      
      subroutine ez_rgoptr(op, val, flag)
      implicit none
      character*8 op
      real val
      logical flag
      integer ier, ezsetval,ezgetval
!     flag = .true.  mode set
!     flag = .false. mode get
      
      
#include "ez_qqqxtrp0.cdk"
      
      character*3 localop
      
      localop=op(1:3)
      call up2low(localop,localop)
      
      if (flag) then
         ier = ezsetval('extrap_value',val)
      else 
         ier = ezgetval('extrap_value',val)
      endif
      
      return
      end
