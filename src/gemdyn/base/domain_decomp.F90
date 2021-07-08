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

!**s/r domain_decomp

      integer function domain_decomp ( F_npex, F_npey, F_checkparti_L )
      use gem_options
      use glb_ld
      use HORgrid_options
      use lun
      use ptopo
      implicit none
#include <arch_specific.hf>

      logical, intent(inout) :: F_checkparti_L
      integer, intent(in):: F_npex, F_npey


      logical, external :: decomp
      integer :: errno, g_err, stat
      integer, dimension(F_npex) :: lnis
      integer, dimension(F_npey) :: lnjs

!-------------------------------------------------------------------
!
      errno = -1
      
      if (Lun_out > 0) write (Lun_out,1000) G_ni,F_npex,G_nj,F_npey

      if (decomp (&
          G_ni,l_minx,l_maxx,lnis,G_lnimax,G_halox,l_i0,.true. ,.true.,&
          F_npex, (Grd_extension+1), F_checkparti_L, 0 ) .and.         &
          decomp (&
          G_nj,l_miny,l_maxy,lnjs,G_lnjmax,G_haloy,l_j0,.false.,.true.,&
          F_npey, (Grd_extension+1), F_checkparti_L, 0 ))              &
          errno = 0

      call rpn_comm_Allreduce ( errno,g_err,1,"MPI_INTEGER","MPI_MIN",&
                                "GRID",stat )
      domain_decomp = g_err

      if (domain_decomp < 0)  then
         if  (Lun_out > 0) then
            write(lun_out,*) 'DECOMP: ILLEGAL DOMAIN PARTITIONING'
         end if
         return
      end if

      l_ni  = lnis(1)
      l_nj  = lnjs(1)
      l_nk  = G_nk
      l_njv = l_nj
      l_niu = l_ni
      if (l_north) l_njv = l_nj - 1
      if (l_east ) l_niu = l_ni - 1

      if (.not.F_checkparti_L) call glbpos ()

      if (Lun_debug_L) write(Lun_out,2000) Ptopo_myrow,Ptopo_mycol,&
                          l_ni,l_nj,l_i0,l_i0+l_ni-1,l_j0,l_j0+l_nj-1

 1000 format (/' DOMAIN_DECOMP: checking partitionning of G_ni and G_nj'&
              /2(i6,' in ',i6,' subdomains',5x)/)
 2000 format (' PROCESSOR GRID and DOMAIN_DECOMP: '/&
              ' (myrow,mycolum)= (',i4,',',i4,')'  /&
              ' (L_ni,L_nj)= ('    ,i5,',',i5,')'  /&
              ' (G_I0,G_IN)= ('    ,i5,',',i5,')'  /&
              ' (G_J0,G_JN)= ('    ,i5,',',i5,')'  )
!
!-------------------------------------------------------------------
!
      return
      end


