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

!**s/r iau_data  - Reads FST iau increments files

      subroutine iau_data ( F_datev_S )
      use lun
      use mem_iau
      use dyn_fisl_options
      use gem_options
      use init_options
      use inp_base
      use inp_mod
      use svri_mod
      use Openfst
      use gmm_pw
      use glb_ld
      use omp_timing
      implicit none
      
      character(len=*), intent(in):: F_datev_S

      logical :: success
      integer :: err,i,j,k,n1,n2,n3
      integer :: ivar, nka, nka_tt, nka_uv, n123(3), n, &
         nreqs, nplans, indx_tt, indx_uv, nm, deb, dim
      integer, dimension (:), pointer :: ip1_list,TT_ip1_list,UV_ip1_list
      real, dimension(:,:  ), pointer :: Sp0_q,Slsp0_q
      real, dimension(:,:,:), pointer :: p0,pres,tt,uur,vvr,tr
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
      integer, parameter :: nlis = 1024
      integer :: liste_sorted(nlis)
      type(fst_query)  :: query
      type(fst_record) :: recs(nlis) 
!
!-----------------------------------------------------------------------
!
      if (Lun_out > 0) write(lun_out,9000) trim(F_datev_S)

      call gtmg_start(53, 'IAU_read', 50)
      if (INs_server_L) then
         Iau_datev  = F_datev_S
         IAU_recv%nreq  = IAU_recv%iBUF(1)
         IAU_recv%nplans= IAU_recv%iBUF(2)
         nreqs= IAU_recv%nreq ; nplans= IAU_recv%nplans
         if (associated(IAU_recv%nk ) ) deallocate (IAU_recv%nk )
         if (associated(IAU_recv%deb) ) deallocate (IAU_recv%deb)
         allocate (IAU_recv%nk(nreqs),IAU_recv%deb(nreqs))
         nm= 2
         IAU_recv%nk (1:nreqs)= IAU_recv%iBUF(nm+1:nm+nreqs) ; nm=nm+nreqs
         IAU_recv%deb(1:nreqs)= IAU_recv%iBUF(nm+1:nm+nreqs) ; nm=nm+nreqs
         IAU_recv%ip1 =>   IAU_recv%iBUF(nm+1:nm+nplans) ; nm=nm+nplans
         n1= IAU_recv%iBUF(nm+1) ; n2= IAU_recv%iBUF(nm+2) ; n3= IAU_recv%iBUF(nm+3)
         indx_tt= -999 ; indx_uv= -999
         do n= 1, IAU_recv%nreq
            if (IAU_recv%cBUF(n+nreqs)(1:2) == 'TT') indx_tt= n
            if (IAU_recv%cBUF(n+nreqs)(1:2) == 'UV') indx_uv= n
         end do
         nka_tt= IAU_recv%nk(indx_tt)
         nka_uv= IAU_recv%nk(indx_uv)/2
         allocate (TT_ip1_list(nka_tt),UV_ip1_list(nka_uv))
         indx_tt= (indx_tt-1)*nka_tt+1
         indx_uv= (indx_uv-1)*nka_uv+1
         TT_ip1_list(1:nka_tt)= IAU_recv%ip1(indx_tt:indx_tt+nka_tt-1)
         UV_ip1_list(1:nka_uv)= IAU_recv%ip1(indx_uv:indx_uv+nka_uv-1)
         allocate ( vtbl_8(n1,n2,n3) )
         vtbl_8 = reshape(IAU_recv%VGD,(/n1,n2,n3/))
         err= vgd_new ( Iau_vgd_src, vtbl_8 )
         deallocate (vtbl_8)
         
         !##### aiguillage pour inp_read_mt
         SRL(1:nreqs)%nk = IAU_recv%nk(1:nreqs)
         SRL(1:nreqs)%deb= IAU_recv%deb(1:nreqs)
         SRL(1:nreqs)%vname(1) = IAU_recv%cBUF(1:nreqs)
         SRL(1:nreqs)%vname(2) = IAU_recv%cBUF(nreqs+1:2*nreqs)
         INS_ND => IAU_recv%RBUF
         INs_DIP1 => IAU_recv%ip1
         INs_recv_nreqs= IAU_recv%nreq
         !#######################
         
         do n=1,INs_recv_nreqs
            if ((trim(SRL(n)%vname(2)) == 'UVRT1' ) .or. &
                (trim(SRL(n)%vname(2)) == 'UV'    ) .and.&
                (SRL(n)%nk>0) ) then
               nka_uv= SRL(n)%nk/2
               deb= (SRL(n)%deb - 1) * INs_dimgzH + 1
               dim= INs_dimgzH*nka_uv
               allocate ( uur(l_minx:l_maxx,l_miny:l_maxy,nka_uv),&
                          vvr(l_minx:l_maxx,l_miny:l_maxy,nka_uv ))
               call reshapeH ( INS_ND(deb:),uur,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,nka_uv)
               call reshapeH ( INS_ND(deb+dim:),vvr,&
                         l_minx,l_ni+G_halox,l_miny,l_nj+G_haloy,&
                         l_minx,l_maxx,l_miny,l_maxy,nka_uv)
               if (lun_out>0) write(6,'(a,a13,a,i4)') &
               ' I-svr FOUND: ','UV',' at indexe: ',n
            endif
         end do
      else
         call open_fst (F_datev_S, Iau_file, Iau_vgd_src, TT_ip1_list, nka_tt, 'IAUREP',.false.)
         nullify (UV_ip1_list,uur,vvr)
         call inp_read_uv ( uur, vvr, 'UV' , UV_ip1_list, nka_uv, F_datev_S, Iau_file )
      endif
      call gtmg_stop ( 53 )

      call gtmg_start(54, 'IAU_interp', 50)
      nullify(Sp0_q,Slsp0_q,ip1_list,tt,p0)
      err = inp_read ('P0', 'Q', p0,1,ip1_list, nka, F_datev_S, Iau_file, F_type_S='r')
      iau_p0(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) = &
          p0(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1)*100.0
      deallocate (ip1_list,p0) ; nullify(ip1_list,p0)
      
      err = inp_read ('P0', 'Q', p0,1,ip1_list, nka, F_datev_S, Iau_file, F_type_S='a')
      allocate (Sp0_q(l_minx:l_maxx,l_miny:l_maxy),&
                 pres(l_minx:l_maxx,l_miny:l_maxy,nka_tt))
      Sp0_q(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) = &
         p0(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,1)*100.0
      deallocate (ip1_list,p0) ; nullify(ip1_list,p0)

      call inp_3dpres ( Iau_vgd_src,TT_ip1_list,Sp0_q,Slsp0_q,&
                        pres,1,nka_tt,F_inlog_S='in_log' )

      err = inp_read ('TT', 'Q', tt, 1, ip1_list, nka, F_datev_S, Iau_file, F_type_S='r')
      call vertint2 ( iau_t,pw_log_pt,G_nk, tt,pres,nka_tt          ,&
                      l_minx,l_maxx,l_miny,l_maxy                   ,&
                      1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                      varname='', inttype= 'cubic', levtype='P' )
      deallocate (ip1_list,tt) ; nullify(ip1_list,tt)

      err = inp_read ('HU', 'Q', tt, 1, ip1_list, nka, F_datev_S, Iau_file,F_type_S='r')
      call vertint2 ( iau_hu,pw_log_pt,G_nk, tt,pres,nka_tt         ,&
                      l_minx,l_maxx,l_miny,l_maxy                   ,&
                      1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                      varname='', inttype= 'cubic', levtype='P' )
      deallocate (ip1_list,tt) ; nullify(ip1_list,tt)

      allocate (tr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
      do ivar=1, IAU_ntr
         err = inp_read ('TR/'//trim(IAU_trname(ivar)), 'Q', &
                             tt, 1, ip1_list, nka,  F_datev_S, Iau_file,F_type_S='r')
         call vertint2 ( tr ,pw_log_pt,G_nk,&
                         tt,pres,nka_tt,l_minx,l_maxx,l_miny,l_maxy  ,&
                       1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                       varname='', inttype= 'cubic', levtype='P' )
         if (Iau_stats_L) call glbstat (tr,trim(IAU_trname(ivar)),'',&
              l_minx,l_maxx,l_miny,l_maxy,1,G_nk,1,G_ni,1,G_nj,1,G_nk)
         iau_tr(:,:,:,ivar) = tr(:,:,:)
         deallocate (ip1_list,tt) ; nullify(ip1_list,tt)
      end do
      deallocate (pres,tr) ; nullify(pres,tr)
      
      allocate (pres(l_minx:l_maxx,l_miny:l_maxy,nka_uv))
      call inp_3dpres ( Iau_vgd_src,UV_ip1_list,Sp0_q,Slsp0_q,&
                        pres,1,nka_uv,F_inlog_S='in_log' )

      call vertint2 ( iau_u,pw_log_pm,G_nk, uur,pres,nka_uv          ,&
                      l_minx,l_maxx,l_miny,l_maxy                   ,&
                      1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                      varname='', inttype= 'cubic',&
                      levtype='P' )
      call vertint2 ( iau_v,pw_log_pm,G_nk, vvr,pres,nka_uv          ,&
                      l_minx,l_maxx,l_miny,l_maxy                   ,&
                      1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,&
                      varname='', inttype= 'cubic',&
                      levtype='P' )
      deallocate (uur,vvr) ; nullify(uur,vvr)

      if (Iau_stats_L) then
         call glbstat (iau_p0,'p0','', l_minx,l_maxx,l_miny,l_maxy,1,1,1,G_ni,1,G_nj,1,1)
         call glbstat (iau_t ,'tt','', l_minx,l_maxx,l_miny,l_maxy,1,G_nk,1,G_ni,1,G_nj,1,G_nk)
         call glbstat (iau_hu,'hu','', l_minx,l_maxx,l_miny,l_maxy,1,G_nk,1,G_ni,1,G_nj,1,G_nk)
         call glbstat (iau_u ,'uu','', l_minx,l_maxx,l_miny,l_maxy,1,G_nk,1,G_ni,1,G_nj,1,G_nk)
         call glbstat (iau_v ,'vv','', l_minx,l_maxx,l_miny,l_maxy,1,G_nk,1,G_ni,1,G_nj,1,G_nk)
      endif
      call gtmg_stop(54)
      
      err = vgd_free(Iau_vgd_src)
      deallocate (pres,Sp0_q,TT_ip1_list,UV_ip1_list)
      if (.not. INs_server_L) success= IAU_file%close()
      INs_recv_nreqs= 0

 9000 format(/,' TREATING IAU INPUT DATA VALID AT: ',a,&
             /,' ===============================================')
!
!-----------------------------------------------------------------------
!
      return
      end subroutine iau_data
