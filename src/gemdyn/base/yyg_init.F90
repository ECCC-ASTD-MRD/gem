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

!**s/r yyg_init

      subroutine yyg_init
      use ISO_C_BINDING
      use glb_ld
      use gem_options
      use yyg_param
      implicit none
#include <arch_specific.hf>

      include "intrp_bicub_yx.inc"

      integer sendlen, recvlen, maxproc, k
      integer wb,eb,sb,nb,wbu,ebu,sbu,nbu,wbv,ebv,sbv,nbv
!
!-------------------------------------------------------------------
!
!     Initialization for Yin-Yang communications

      call yyg_extend_grid ()
      allocate (YYG_pe_indx(4,Ptopo_numproc))

      YYG_i0 = l_i0
      YYG_j0 = l_j0
      YYG_lni= l_ni
      YYG_lnj= l_nj
      YYG_pe_indx = Ptopo_gindx

      call yyg_initblen ( YYG_BLEN_q2q      , &
                          YYG_xg_8, YYG_yg_8, &
                          YYG_xg_8, YYG_yg_8, &
                          YYG_xg_8, YYG_yg_8, &
                          G_ni, G_nj , &
                          G_ni, G_nj , &
                          G_ni ,G_nj  )

      call yyg_initblen ( YYG_BLEN_uv2u       , &
                          YYG_xgu_8, YYG_yg_8 , &
                          YYG_xgu_8, YYG_yg_8 , &
                          YYG_xg_8 , YYG_ygv_8, &
                          G_niu, G_nj , &
                          G_niu, G_nj , &
                          G_ni , G_njv  )

      call yyg_initblen ( YYG_BLEN_uv2v       , &
                          YYG_xg_8 , YYG_ygv_8, &
                          YYG_xgu_8, YYG_yg_8 , &
                          YYG_xg_8 , YYG_ygv_8, &
                          G_ni , G_njv, &
                          G_niu, G_nj , &
                          G_ni , G_njv  )

      call yyg_initcomm ( YYG_PILT_q2q,  &
                          YYG_xg_8, YYG_yg_8 , &
                          YYG_xg_8, YYG_yg_8 , &
                          YYG_xg_8, YYG_yg_8 , &
                          G_ni, G_nj , 0,0   , &
                          G_ni, G_nj , &
                          G_ni ,G_nj  )

      call yyg_initcomm ( YYG_PILT_uv2u       , &
                          YYG_xgu_8, YYG_yg_8 , &
                          YYG_xgu_8, YYG_yg_8 , &
                          YYG_xg_8 , YYG_ygv_8, &
                          G_niu, G_nj , 0,0   , &
                          G_niu, G_nj , &
                          G_ni , G_njv  )

      call yyg_initcomm ( YYG_PILT_uv2v       , &
                          YYG_xg_8 , YYG_ygv_8, &
                          YYG_xgu_8, YYG_yg_8 , &
                          YYG_xg_8 , YYG_ygv_8, &
                          G_ni , G_njv, 0,0   , &
                          G_niu, G_nj , &
                          G_ni , G_njv  )

      call yyg_initcomm ( YYG_NEAR_q2q       , &
                          YYG_xg_8, YYG_yg_8 , &
                          YYG_xg_8, YYG_yg_8 , &
                          YYG_xg_8, YYG_yg_8 , &
                          G_ni, G_nj , 0,0   , &
                          G_ni, G_nj , &
                          G_ni ,G_nj , F_inttype_S='NEAR' )

      YYG_i0 = l_i0 - west *G_halox
      YYG_j0 = l_j0 - south*G_haloy
      YYG_lni= l_ni + east *G_halox
      YYG_lnj= l_nj + north*G_haloy
      do k= 1, Ptopo_numproc
         YYG_pe_indx(1,k)= Ptopo_gindx(1,k)-G_halox*west
         YYG_pe_indx(2,k)= Ptopo_gindx(2,k)+G_halox*east
         YYG_pe_indx(3,k)= Ptopo_gindx(3,k)-G_haloy*south
         YYG_pe_indx(4,k)= Ptopo_gindx(4,k)+G_halox*north
      end do

      call yyg_initcomm ( YYG_HALO_q2q           ,&
                          YYG_xg_8, YYG_yg_8     ,&
                          YYG_xg_8, YYG_yg_8     ,&
                          YYG_xg_8, YYG_yg_8     ,&
                          G_ni, G_nj , G_halox, G_haloy,&
                          G_ni, G_nj , &
                          G_ni ,G_nj  )

!for yyg_xchng_vec_uv2uv
      sendlen= max( YYG_PILT_uv2u%maxsend,YYG_PILT_uv2v%maxsend,&
                    YYG_BLEN_uv2u%maxsend,YYG_BLEN_uv2v%maxsend )
      maxproc= max( YYG_PILT_uv2u%sendmaxproc,YYG_BLEN_uv2u%sendmaxproc,&
                    YYG_PILT_uv2v%sendmaxproc,YYG_BLEN_uv2v%sendmaxproc )
      if (sendlen > 0) allocate (YYG_usenduv(sendlen*(G_nk+1),maxproc), &
                                 YYG_vsenduv(sendlen*(G_nk+1),maxproc)  )

      recvlen= max( YYG_PILT_uv2u%maxrecv,YYG_PILT_uv2v%maxrecv,&
                    YYG_BLEN_uv2u%maxrecv,YYG_BLEN_uv2v%maxrecv )
      maxproc= max( YYG_PILT_uv2u%recvmaxproc,YYG_BLEN_uv2u%recvmaxproc,&
                    YYG_PILT_uv2v%recvmaxproc,YYG_BLEN_uv2v%recvmaxproc )
      if (recvlen > 0) allocate (YYG_urecvuv(recvlen*(G_nk+1),maxproc), &
                                 YYG_vrecvuv(recvlen*(G_nk+1),maxproc)  )

      YYG_PILT_recv_uv= max( YYG_PILT_uv2u%maxrecv,YYG_PILT_uv2v%maxrecv)
      YYG_BLEN_recv_uv= max( YYG_BLEN_uv2u%maxrecv,YYG_BLEN_uv2v%maxrecv)

!for yyg_xchng_sca_q2q and yyg_xchng_vec_q2q
      sendlen= max( YYG_PILT_q2q%maxsend    , YYG_BLEN_q2q%maxsend    ,&
                    YYG_HALO_q2q%maxsend    , YYG_NEAR_q2q%maxsend     )
      maxproc= max( YYG_PILT_q2q%sendmaxproc, YYG_BLEN_q2q%sendmaxproc,&
                    YYG_HALO_q2q%sendmaxproc, YYG_NEAR_q2q%sendmaxproc )
      if (sendlen > 0) allocate (YYG_uvsend(sendlen*(G_nk+1)*2,maxproc))

      recvlen= max( YYG_PILT_q2q%maxrecv    , YYG_BLEN_q2q%maxrecv    ,&
                    YYG_HALO_q2q%maxrecv    , YYG_NEAR_q2q%maxrecv     )
      maxproc= max( YYG_PILT_q2q%recvmaxproc, YYG_BLEN_q2q%recvmaxproc,&
                    YYG_HALO_q2q%recvmaxproc, YYG_NEAR_q2q%recvmaxproc )
      if (recvlen > 0) allocate (YYG_uvrecv(recvlen*(G_nk+1)*2,maxproc))

      wb = 1-G_halox ; eb = l_ni+G_halox
      sb = 1-G_haloy ; nb = l_nj+G_haloy
      wbu = wb ; ebu = eb ; wbv = wb ; ebv = eb
      sbu = sb ; nbu = nb ; sbv = sb ; nbv = nb

      if (l_west ) then
         wb = pil_w+1 ; wbu = wb ; wbv = wb
      end if
      if (l_east ) then
         eb = l_ni-pil_e-3 ; ebv = eb
         ebu= l_niu-pil_e-3
      end if
      if (l_south) then
         sb = pil_s+1 ; sbu = sb ; sbv = sb
      end if
      if (l_north) then
         nb = l_nj-pil_n-3 ; nbu = nb
         nbv= l_njv-pil_n-3
      end if

      wb  =max(wb  +l_i0-1,    1+Glb_pil_w  )
      eb  =min(eb  +l_i0-1, G_ni-Glb_pil_e-3)
      sb  =max(sb  +l_j0-1,    1+Glb_pil_s  )
      nb  =min(nb  +l_j0-1, G_nj-Glb_pil_n-3)
      wbu =max(wbu +l_i0-1,    1+Glb_pil_w  )
      ebu =min(ebu +l_i0-1, G_ni-Glb_pil_e-4)
      sbu =max(sbu +l_j0-1,    1+Glb_pil_s  )
      nbu =min(nbu +l_j0-1, G_nj-Glb_pil_n-3)
      wbv =max(wbv +l_i0-1,    1+Glb_pil_w  )
      ebv =min(ebv +l_i0-1, G_ni-Glb_pil_e-3)
      sbv =max(sbv +l_j0-1,    1+Glb_pil_s  )
      nbv =min(nbv +l_j0-1, G_nj-Glb_pil_n-4)

      call set_intrp_bicub_off_yx (1-l_i0,1-l_j0)

      call set_intrp_bicub_quv (wb ,eb ,sb ,nb ,&
                                wbu,ebu,sbu,nbu,&
                                wbv,ebv,sbv,nbv)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine yyg_init
