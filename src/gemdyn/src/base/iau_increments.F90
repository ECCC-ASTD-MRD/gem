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

      subroutine iau_increments (F_kount)
      use iso_c_binding
      use ctrl
      use dyn_fisl_options
      use glb_ld
      use cstv
      use inp_mod
      use gmm_iau
      use mem_iau
      use gmm_pw
      use gmm_vt1
      use init_options
      use step_options
      use lun
      use metric
      use mem_tracers
      use ptopo
      use svri_mod
      use tr3d
      use tdpack
      use omp_timing
      implicit none
      
      integer, intent(in) :: F_kount
      
      include 'mpif.h'
      character(len=16) :: datev, date_n, next_S
      integer :: n,i,j,k,yy,mo,dd,hh,mm,ss,dum,ivar,dim,&
                 i0,in,j0,jn,err
      real, pointer, dimension(:,:,:) :: hu,vt
      real(kind=REAL64) :: dayfrac,tx,a,b,weight
      real(kind=REAL64), parameter :: one=1.0d0, &
                    sid=86400.0d0, rsid=one/sid
      real(kind=REAL64), dimension(l_minx:l_maxx,l_miny:l_maxy) :: delq
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: vtm
!     
!-----------------------------------------------------------------------
!
      if (.not. Ctrl_iau_L) return
      Ctrl_iau_L= (Cstv_dt_8*(F_kount+1) <= Iau_period)

      call gtmg_start(50, 'IAU', 1)
      dayfrac = dble(Step_kount) * Cstv_dt_8 * rsid
      call incdatsd  (date_n, Step_runstrt_S, dayfrac)

      if (Step_kount > IAU_ubstp) then
         if (IAU_ubstp < 0) then
            call iau_fisrt_datev (datev,IAU_ubstp,Step_kount)
            IAU_now= datev
         else
            n       = Iau_interval/Cstv_dt_8
            dayfrac = Iau_interval*rsid
            call incdatsd (datev, IAU_now, dayfrac)
            IAU_now  = datev
            IAU_ubstp= IAU_ubstp+n
         endif
         
         if (INs_server_L) then
            call gtmg_start ( 51, 'INs_wait', 50)
            call gemtime (Lun_out, 'Input-svr: wait for IAU data', .false.)
            call MPI_waitall (size(INs_Iau_irecv),INs_Iau_irecv,&
                              MPI_STATUSES_IGNORE,err)
            call gemtime (Lun_out, 'Input-svr: IAU data received', .false.)
            call gtmg_stop ( 51 )
            dim= size(IAU_recv%cBUF)*len(IAU_recv%cBUF(1))
            call MPI_bcast (IAU_recv%cBUF,dim,MPI_CHARACTER, 0,&
                            COMM_multigrid, err)         
            call MPI_bcast (IAU_recv%iBUF,size(IAU_recv%iBUF),MPI_INTEGER, 0,&
                            COMM_multigrid, err)         
            call MPI_bcast (IAU_recv%VGD,size(IAU_recv%VGD),MPI_DOUBLE_PRECISION, 0,&
                            COMM_multigrid, err)         
            dayfrac = Iau_interval*rsid
            call incdatsd (next_S, IAU_now, dayfrac)
            if ( next_S < IAU_last_S ) then
               call itf_Iserv_request (next_S,INs_Iaulist_S,&
                               INs_Iau_tag,INs_Iau_nrequests)
            endif
         endif
         call iau_data (datev)
      endif

      call gtmg_start(56, 'IAU_increments', 50)
      call tt2virt(vtm, .true., l_minx, l_maxx, l_miny, l_maxy, G_nk) 
      weight= Cstv_dt_8 / Iau_period
      i0= 1+pil_w ; in= l_ni-pil_e
      j0= 1+pil_s ; jn= l_nj-pil_n
      
      if (Lun_out > 0) write(6,'(a,i6,a,1x,a,f15.12/20("#"))') &
         'IAU_increments: ',Step_kount,date_n,IAU_now,weight
      hu=> tracers_P(Tr3d_hu)%pntr
      pw_tt_plus(i0:in,j0:jn,:) = &
                    pw_tt_plus(i0:in,j0:jn,:) + weight * iau_t(i0:in,j0:jn,:)
      hu(i0:in,j0:jn,:) = &
                    hu(i0:in,j0:jn,:) + weight * iau_hu(i0:in,j0:jn,:)
      pw_uu_plus(i0:in,j0:jn,:) = &
                    pw_uu_plus(i0:in,j0:jn,:) + weight * iau_u(i0:in,j0:jn,:)
      pw_vv_plus(i0:in,j0:jn,:) = &
                    pw_vv_plus(i0:in,j0:jn,:) + weight * iau_v(i0:in,j0:jn,:)
!special case here for qt1
      allocate(vt(l_minx:l_maxx, l_miny:l_maxy, G_nk))
      call tt2virt(vt, .true., l_minx, l_maxx, l_miny, l_maxy, G_nk) !# compute VT from incremented TT,HU,...
      qt1(i0:in,j0:jn,l_nk+1)= qt1(i0:in,j0:jn,l_nk+1) + rgasd_8*Cstv_Tstr_8* &
                                 log(1.d0 + weight*iau_p0(i0:in,j0:jn) / &
                                 exp(GVM%lg_pstar_8(i0:in,j0:jn,l_nk+1)+qt1(i0:in,j0:jn,l_nk+1)/ &
                                 (rgasd_8*Cstv_Tstr_8) ) )
      !Initializing delq with the surface value         
      delq(i0:in,j0:jn)= rgasd_8*Cstv_Tstr_8 * log(1.d0 + iau_p0(i0:in,j0:jn) / &
                           exp(GVM%lg_pstar_8(i0:in,j0:jn,l_nk+1)+qt1(i0:in,j0:jn,l_nk+1)/ &
                           (rgasd_8*Cstv_Tstr_8) ) )
      do k = l_nk, 1, -1
         !Computing delq at level k from k+1 for IAU increments
         delq(i0:in,j0:jn)=delq(i0:in,j0:jn)+ grav_8*Cstv_Tstr_8* &
              (1.d0/vt(i0:in,j0:jn,k) - 1.d0/tt1(i0:in,j0:jn,k))/ & 
                           GVM%mc_iJz_8(i0:in,j0:jn,k)
                  
         !Updating qt1 at level k with the IAU increments
         qt1(i0:in,j0:jn,k) = qt1(i0:in,j0:jn,k) + rgasd_8*Cstv_Tstr_8* &
                                log(1.d0 + weight*(exp(delq(i0:in,j0:jn)/ &
                                (rgasd_8*Cstv_Tstr_8)) - 1.d0) )                 
         iau_tv_tend(1:l_ni,1:l_nj,k)=(vt(1:l_ni,1:l_nj,k)-vtm(1:l_ni,1:l_nj,k))/&
                  Cstv_dt_8
               end do
      deallocate (vt)
 888  do ivar=1, IAU_ntr
         tracers_P(IAU_trindx(ivar))%pntr(i0:in,j0:jn,:)= &
         tracers_P(IAU_trindx(ivar))%pntr(i0:in,j0:jn,:)  &
                        + weight * iau_tr(i0:in,j0:jn,:,ivar)
      end do

      call gtmg_stop(56)
      call gtmg_stop(50)
!     
!-----------------------------------------------------------------------
!
      return
      end subroutine iau_increments
