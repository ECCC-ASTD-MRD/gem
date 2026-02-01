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

!**s/r out_stkecr

      subroutine out_stkecr ( fa,lminx,lmaxx,lminy,lmaxy, &
                              metaf,nplans, g_id,g_if,g_jd,g_jf )
      use iso_c_binding
      use out_collector
      use HORgrid_options
      use out_options
      use glb_ld
      use out_mod
      use out3
      use ptopo
      use omp_timing
      implicit none

      integer lminx,lmaxx,lminy,lmaxy,nplans
      integer g_id,g_if,g_jd,g_jf
      real fa(lminx:lmaxx,lminy:lmaxy,nplans)
      type(fst_record), dimension(:), pointer :: metaf
      type(fst_query) :: query
      
#include <rmnlib_basics.hf>
      include 'mpif.h'
      include "rpn_comm.inc"

      logical iope_L
      integer  nz, err, ni, nis, njs, k, kk, wk_njs,tag,stat,l
      integer, dimension (:)    , pointer :: zlist
      real   , dimension (:,:,:), pointer :: wk_glb
!
!----------------------------------------------------------------------
!
      nis = g_if - g_id + 1
      njs = g_jf - g_jd + 1
      wk_njs = -1 ; nz = 0
      if ( (nis < 1) .or. (njs < 1) ) return

      if (out_type_S == 'REGDYN') then
         call gtmg_start ( 81, 'OUT_DUCOL', 80)
      else
         call gtmg_start ( 92, 'OUT_PUCOL', 48)
      end if

      if (Out3_ezcoll_L) then
         iope_L= (Out3_iome >= 0)
         nz    = (nplans + Out3_npes -1) / Out3_npes
         wk_glb(1:G_ni,1:G_nj,1:nz) => Glb_fld(1:)
         if (Out3_iome >= 0) then
            wk_glb(1:G_ni,1:G_nj,1:nz) => Glb_fld(1:)
            zlist(1:nz) => List_nk(1:nz)
         else
            wk_glb(1:1,1:1,1:1) => Glb_fld(1:)
            zlist(1:1) => List_nk(1:)
         end if
         zlist= -1
         err= RPN_COMM_shuf_ezcoll ( Out3_comm_setno, Out3_comm_id, &
                                     wk_glb, nz, fa, nplans, zlist )
      else
         iope_L= (Bloc_me == 0)
         call block_collect_fullp ( fa, l_minx,l_maxx,l_miny,l_maxy, &
                                    nplans, Glb_fld, nz, List_nk )
      end if

      if ( (iope_L) .and. (nz>0) ) then
         if ((Grd_yinyang_L) .and. (Ptopo_couleur == 0)) then
            wk_njs = njs*2
         else
            wk_njs = njs
         end if
      end if

      if (out_type_S == 'REGDYN') then
         call gtmg_stop (81)
         call gtmg_start ( 82, 'OUT_DUECR', 80)
      else
         call gtmg_stop (92)
         call gtmg_start ( 93, 'OUT_PUECR', 48)
      end if

      IOPE: if (iope_L) then

         do k= nz, 1, -1

            if (List_nk(k) > 0) then
               kk= List_nk(k)
               
               l=G_ni*G_nj*(k-1)+1
               if ( (Grd_yinyang_L) .and. (.not.Out_reduc_l) ) then

                  !Merge from Yang (couleur 1) to Yin (couleur 0)

                  tag=401

                  if (Ptopo_couleur == 0) then
                     call reduc (Reduc_fld,Glb_fld(l),nis,njs,&
                                 G_ni,G_nj,g_id,g_if,g_jd,g_jf)
                     call RPN_COMM_recv ( Reduc_fld(nis*njs+1), nis*njs,&
                             'MPI_REAL', 1, tag, 'GRIDPEERS', stat, err )
                     Out_rec=metaf(kk)
                     Out_rec%typvar=Out_typvar_S
                     Out_rec%etiket=Out_etik_S
                     Out_rec%dateo=Out_dateo
                     Out_rec%deet=Out_deet
                     Out_rec%npas=Out_npas
                     Out_rec%grtyp='U'
                     Out_rec%ni=nis
                     Out_rec%nj=2*njs
                     Out_rec%nk=1
                     Out_rec%ig4=Out_ig4
                     Out_rec%data=c_loc(Reduc_fld(1))

                     success = Out_file%write(Out_rec,rewrite=FST_SKIP)

                  else
                     call reduc (Reduc_fld,Glb_fld(l),nis,njs,&
                                 G_ni,G_nj,g_id,g_if,g_jd,g_jf)
                     call RPN_COMM_send ( Reduc_fld, nis*njs,&
                         'MPI_REAL', 0, tag, 'GRIDPEERS', err )
                  end if

               else

                  call reduc (Reduc_fld,Glb_fld(l),nis,njs,&
                              G_ni,G_nj,g_id,g_if,g_jd,g_jf)
                  Out_rec=metaf(kk)
                  Out_rec%typvar=Out_typvar_S
                  Out_rec%etiket=Out_etik_S
                  Out_rec%dateo=Out_dateo
                  Out_rec%deet=Out_deet
                  Out_rec%npas=Out_npas
                  Out_rec%grtyp='Z'
                  Out_rec%ni=ni
                  Out_rec%nj=njs
                  Out_rec%nk=1
                  Out_rec%ig4=Out_ig4
                  Out_rec%data=c_loc(Reduc_fld(1))
                  
                  success = Out_file%write(Out_rec,rewrite=FST_SKIP)
                  
               end if

            end if

         end do

      end if IOPE

      if (out_type_S == 'REGDYN') then
         call gtmg_stop (82)
      else
         call gtmg_stop (93)
      end if
!
!--------------------------------------------------------------------
!
      return
      end

      subroutine reduc (F_dest,F_src,ni,nj,Gni,Gnj,i0,in,j0,jn)
      implicit none
      integer :: ni,nj,Gni,Gnj,nk,i0,in,j0,jn
      real :: F_dest(ni,nj), F_src(Gni,Gnj)
      F_dest(1:ni,1:nj) = F_src(i0:in,j0:jn)
      return
      end
      
