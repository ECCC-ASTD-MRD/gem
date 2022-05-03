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

!**s/r nest_glbstat - 

      subroutine nest_glbstat (F_lframe_S,nf)
      use glb_ld
      use mem_nest
      use tr3d
      use  clib_itf_mod
      implicit none
      
      integer, intent(IN) :: nf
      character(len=*), intent(inout):: F_lframe_S(nf)

      integer n,k,istat
!     
!     ---------------------------------------------------------------
!
!$omp single
      do n=1, nf
         istat = clib_toupper(F_lframe_S(n))
      end do
      
      do n=1, nf
         print*, n,F_lframe_S(n)
      if (F_lframe_S(n) == 'DEB') then
      call glbstat (nest_u_deb,'NUD','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_v_deb,'NVD','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_t_deb,'NTD','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_w_deb,'NWD','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_zd_deb,'NZDD','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_q_deb,'NQD','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk+1,&
                    1,G_ni,1,G_nj,1,l_nk+1)
      call glbstat (nest_s_deb,'NSD','', l_minx,l_maxx,l_miny,l_maxy,1,1,&
                    1,G_ni,1,G_nj,1,1)
      do k=1,Tr3d_ntr
         call glbstat (nest_tr_deb(l_minx,l_miny,(k-1)*l_nk+1),'NTRD',Tr3d_name_S(k), &
              l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
      end do
      endif
      if (F_lframe_S(n) == 'FIN') then
      call glbstat (nest_u_fin,'NUF','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_v_fin,'NVF','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_t_fin,'NTF','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_w_fin,'NWF','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_zd_fin,'NZDF','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_q_fin,'NQF','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk+1,&
                    1,G_ni,1,G_nj,1,l_nk+1)
      call glbstat (nest_s_fin,'NSF','', l_minx,l_maxx,l_miny,l_maxy,1,1,&
                    1,G_ni,1,G_nj,1,1)
      do k=1,Tr3d_ntr
         call glbstat (nest_tr_fin(l_minx,l_miny,(k-1)*l_nk+1),'NTRF',Tr3d_name_S(k), &
              l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
      end do
      endif
      if (F_lframe_S(n) == 'NOW') then
      call glbstat (nest_u,'NUN','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_v,'NVN','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_t,'NTN','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_w,'NWN','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_zd,'NZDN','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk,&
                    1,G_ni,1,G_nj,1,l_nk)
      call glbstat (nest_q,'NQN','', l_minx,l_maxx,l_miny,l_maxy,1,l_nk+1,&
                    1,G_ni,1,G_nj,1,l_nk+1)
      call glbstat (nest_s,'NSN','', l_minx,l_maxx,l_miny,l_maxy,1,1,&
                    1,G_ni,1,G_nj,1,1)
      do k=1,Tr3d_ntr
         call glbstat (nest_tr(l_minx,l_miny,(k-1)*l_nk+1),'NTRN',Tr3d_name_S(k), &
              l_minx,l_maxx,l_miny,l_maxy,1,l_nk,1,G_ni,1,G_nj,1,l_nk)
      end do
      endif
      end do
!$omp end single
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine nest_glbstat
      
