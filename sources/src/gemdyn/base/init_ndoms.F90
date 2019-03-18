!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

      subroutine init_ndoms (ndomains,dom_deb,err)
      use HORgrid_options
      use clib_itf_mod
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer ndomains,dom_deb,err


      character(len=16) :: ndomains_S, last_domain_S
      integer indx1, dom_fin, last
!
!-------------------------------------------------------------------
!
      err  = -1

      if (clib_getenv ('GEM_NDOMAINS',ndomains_S) < 0) then
         write (6,1001) 'GEM_NDOMAINS'
         return
      end if

      indx1= index (ndomains_S,':')

      if (indx1 < 1) then
         write (6,1002) 'GEM_NDOMAINS',ndomains_S
         return
      end if

      read (ndomains_S(1:indx1-1),*,end=33,err=33) dom_deb
      read (ndomains_S(indx1+1: ),*,end=33,err=33) dom_fin

      goto 101
  33  write (6,1002) 'GEM_NDOMAINS',ndomains_S
      return

 101  if (clib_getenv ('DOMAIN_end',last_domain_S) < 0) then
         write (6,1001) 'DOMAIN_end'
         return
      end if

      read (last_domain_S,*,end=43,err=43) last

      goto 201
  43  write (6,1002) 'DOMAIN_end',last_domain_S
      return

 201  Grd_ndomains = dom_fin - dom_deb + 1

      if (Grd_ndomains < 1) then
         write (6,1003) ndomains,ndomains_S
         return
      end if

      Ptopo_last_domain_L = (dom_fin == last)

      Grd_yinyang_L = .false.
      Grd_yinyang_S = ''
      if (clib_getenv ('GEM_YINYANG',ndomains_S) >= 0) then
         Grd_yinyang_L = .true.
      end if

      ndomains = Grd_ndomains
      err      = 0

 1001 format (/' =====> Error in init_ndoms: Env variable ',a,' is undefined'/)
 1002 format (/' =====> Error in init_ndoms: Env variable ',a,' is incorrectly defined ',a/)
 1003 format (/' =====> Error in init_ndoms: ndomains_S= ',i4,' Check Env variable GEM_NDOMAINS ',a/)
!
!-------------------------------------------------------------------
!
      return
      end

