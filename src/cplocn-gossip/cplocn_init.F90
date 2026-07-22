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

      integer function cplocn_init (F_path_S, F_print_L, F_unout, F_dateo, &
                                    F_z0mtype, F_z0lat, F_z0ttype)
use iso_c_binding
use rmn_gmm
use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_path_S,F_z0mtype,F_z0ttype
      logical, intent(in)          :: F_print_L
      integer, intent(in)          :: F_dateo, F_unout
      real   , intent(in)          :: F_z0lat(2)

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version
! v5.0b1 - JM Belanger (dec 2017) - z0ttype

!Purpose
! Initialize coupling with ocean

#include <rmn/WhiteBoard.hf>
      include "thermoconsts.inc"
      include "cpl.cdk"
      include "cplocn.cdk"
      include "mpif.h"
      include "rpn_comm.inc"

      integer, external :: mgi_init,mgi_open,mgi_write,mgi_read

      integer i,j,ier,icpl_ou,icpl_in,isnd,nsend,err,cnt,ivar,ibidon
      integer icplo,icpla
      parameter (nsend = 51)
      character(len=512) :: s_send(nsend)
      logical       l_send(nsend) 
      integer       i_send(nsend)
      real          r_send(nsend), Z0TLAT_W(2), Z0TLAT_R(2)
      character(len=16) :: z0mtype_R,z0ttype_R
      type(gmm_metadata) :: meta3d_ocn_busin, meta2d_surf0
      integer, parameter :: n0=0, n1=1

      logical :: l_yinyang, l_lam
      character(len=3) :: yysubgrid_S
      character(len=2) :: grd_S
      character(len=7) :: cplocn_R_chan_name,cplocn_W_chan_name
      character(len=512) :: weightfile_S
!
!     ---------------------------------------------------------------
!
      cplocn_init=-1

      if (F_print_L) write(F_unout,2000)

      err= RPN_COMM_mype (cplocn_myproc, cplocn_mycol, cplocn_myrow)
      if ( err < 0 ) then
        write(F_unout,9900) 'PROBLEMS WHEN CALLING RPN_COMM_mype'
        call flush(F_unout)
        goto 995
      endif
      call RPN_COMM_size ('GRID', cplocn_nprocs, err)
      if (err < 0) then
        write(F_unout,9900) 'PROBLEMS WHEN CALLING RPN_COMM_size'
        call flush(F_unout)
        goto 995
      endif

      err = 0

      s_send    =''
      s_send(1) = 'GEMatm'
      s_send(2) = F_z0mtype
      s_send(3) = F_z0ttype

      l_send = .false.
      i_send = 0
      r_send = 0.0

      do icpl_ou=1, cplocn_n_fldou
        isnd=3+icpl_ou
        if (isnd.gt.nsend) then
          write(F_unout,9900) ' NOT ENOUGH STORAGE FOR s_send'
          call flush(F_unout)
          err=-1
          goto 995
        endif
        s_send(isnd)=cplocn_cvou_S(icpl_ou)//'_'//cplocn_cvot_S(icpl_ou)
      enddo

      do icpl_in=1, cplocn_n_fldin
        isnd=3+cplocn_n_fldou+icpl_in
        if (isnd.gt.nsend) then
          write(F_unout,9900) ' NOT ENOUGH STORAGE FOR s_send'
          call flush(F_unout)
          err=-1
          goto 995
        endif
        s_send(isnd)=cplocn_cvin_S(icpl_in)//'_'//cplocn_cvit_S(icpl_in)
      enddo

      call datf2p (cplocn_runstrt_S,F_dateo)

      l_yinyang = .false.
      l_lam = .true.  ! LU or yinyang only now
      ier = wb_get('model/Hgrid/is_yinyang',l_yinyang)
      yysubgrid_S=''
      if (l_yinyang) then 
        ier = wb_get('model/Hgrid/yysubgrid',yysubgrid_S)
      endif

      cplocn_dynphy_offset = 0 
      if (l_lam) cplocn_dynphy_offset = 2

      select case (yysubgrid_S)
      case('YIN')
        cplocn_W_chan_name='yin2ocn'
        cplocn_R_chan_name='ocn2yin'
        weightfile_S=trim(F_path_S)//'/MODEL_INPUT/cplwgt.fst_YIN'
      case('YAN')                  
        cplocn_W_chan_name='yan2ocn'
        cplocn_R_chan_name='ocn2yan'
        weightfile_S=trim(F_path_S)//'/MODEL_INPUT/cplwgt.fst_YAN'
      case default
        cplocn_W_chan_name='atm2ocn'
        cplocn_R_chan_name='ocn2atm'
        weightfile_S=trim(F_path_S)//'/MODEL_INPUT/cplwgt.fst'
      end select      

      i_send(1) = F_dateo
      i_send(2) = cplocn_dynphy_offset

      r_send(1) = cpl_drv_delt

      if (F_print_L) then
        write (F_unout,*) 'CPLOCN_INIT: l_yinyang=', l_yinyang
        write (F_unout,*) 'CPLOCN_INIT: yysubgrid_S=', yysubgrid_S
        write (F_unout,*) 'CPLOCN_INIT: cplocn_R_chan_name=', cplocn_R_chan_name
        write (F_unout,*) 'CPLOCN_INIT: cplocn_W_chan_name=', cplocn_W_chan_name

        write (F_unout,*) 'CPLOCN_INIT: cplocn_dynphy_offset=',cplocn_dynphy_offset
      endif

      call cplao_init  ( 'GEMatm', trim(F_path_S)//                        &
                                       '/MODEL_INPUT/cplao_settings.nml',  &
                          cplocn_myproc.eq.0,                              &
                          s_send,l_send,i_send,r_send,nsend,               &
                          cplocn_W_chan,cplocn_W_chan_name,                &
                          cplocn_R_chan,cplocn_R_chan_name,                &
                          cpl_drv_gni,cpl_drv_gnj,                         &
                          cplocn_gni,cplocn_gnj,                           &
                          cplocn_n_fldin,cplocn_n_fldou,                   &
                          ibidon, cplocn_atm_nspread, err )
      if ( err < 0 ) then
         write (F_unout,9900) 'UNABLE TO INITIALIZE COUPLER'
         call flush(F_unout)
         goto 995
      endif
      call RPN_COMM_bcast (cplocn_gni, n1, "MPI_INTEGER", n0,"GRID",ier)
      call RPN_COMM_bcast (cplocn_gnj, n1, "MPI_INTEGER", n0,"GRID",ier)

      err = 0

      if (cplocn_myproc.eq.0) then

         cplocn_oc_dt   = nint(r_send(1))
         cplocn_1st_L   = l_send(1)
         cplocn_off_L   = l_send(2)
         cplocn_bzone_L = l_send(3)
         if (cplocn_1st_L.and.cplocn_off_L) then
            write (F_unout,9900) ' WRONG LOGIC FOR OCEAN COUPLING'
            call flush(F_unout)
            err=-1
            goto 995
         endif

         if (F_print_L) then
           write(F_unout,*)  'cplocn_init: cplocn_oc_dt  =',cplocn_oc_dt
           write(F_unout,*)  'cplocn_init: cplocn_1st_L  =',cplocn_1st_L
           write(F_unout,*)  'cplocn_init: cplocn_off_L  =',cplocn_off_L
           write(F_unout,*)  'cplocn_init: cplocn_bzone_L=',cplocn_bzone_L
         endif
         if (F_print_L) then
            if (cplocn_off_L) then
               write(F_unout,*)  'cplocn_init: WARNING, CPL mode: off ==> ocean model completely bypassed'
            elseif (cplocn_1st_L) then
               write(F_unout,*)  'cplocn_init: WARNING, CPL mode: only 1st ocean receive considered'
            else
               write(F_unout,*)  'cplocn_init: Normal CPL mode'
            endif
         endif

         if (trim(s_send(1)).ne.'NEMoce') then
            write (F_unout,9900) ' WRONG NAME FOR OTHER MODEL: should be NEMoce'
            call flush(F_unout)
            err = -1
            goto 995
         endif

         do icpl_in=1, cplocn_n_fldin
           isnd=3+icpl_in
           if (trim(s_send(isnd)) /=  &
               cplocn_cvin_S(icpl_in)//'_'//cplocn_cvit_S(icpl_in)) then
             write (F_unout,*) 'cplocn_cvin_S(icpl_in)_cplocn_cvit_S(icpl_in)=', &
                                cplocn_cvin_S(icpl_in)//'_'//cplocn_cvit_S(icpl_in)
             write (F_unout,*) 's_send(isnd)=', s_send(isnd)
             write (F_unout,9900) 'ORDER NOT RESPECTED IN COUPLING SEND RECEIVED FIELDS'
             call flush(F_unout)
             err = -1
             goto 995
           endif
         enddo

         do icpl_ou=1, cplocn_n_fldou
           isnd=3+cplocn_n_fldin+icpl_ou
           if (trim(s_send(isnd)) /=  &
               cplocn_cvou_S(icpl_ou)//'_'//cplocn_cvot_S(icpl_ou)) then
             write (F_unout,*) 'cplocn_cvou_S(icpl_ou)_cplocn_cvot_S(icpl_ou)=', &
                                cplocn_cvou_S(icpl_ou)//'_'//cplocn_cvot_S(icpl_ou)
             write (F_unout,*) 's_send(isnd)=', s_send(isnd)
             write (F_unout,9900) 'ORDER NOT RESPECTED IN COUPLING SEND RECEIVED FIELDS'
             call flush(F_unout)
             err = -1
             goto 995
           endif
         enddo

         icpla = 1
         icplo = 0
         err = mgi_read (cplocn_R_chan, icplo, 1, "I" )
         err = min ( mgi_write(cplocn_W_chan, icpla, 1, "I" ), err )

         icpla=icplo+icpla
         if ( icpla /= 2 .or. err < 0 ) then
            write (F_unout,9900) 'UNABLE TO INITIALIZE COUPLING WITH OCEAN'
            call flush(F_unout)
            err = -1
            goto 995 
         endif

         Z0TLAT_W(:) = F_z0lat(:)*180./pi

         err = mgi_write(cplocn_W_chan, Z0TLAT_W, 2, 'R')
         err = min( mgi_read (cplocn_R_chan, Z0TLAT_R, 2, 'R'), err )

         if ( ABS(Z0TLAT_W(1)-Z0TLAT_R(1)) > 0.01 .or. &
              ABS(Z0TLAT_W(2)-Z0TLAT_R(2)) > 0.01 .or. &
              err < 0 ) then
            write (F_unout,*) 'err        =',err
            write (F_unout,*) 'Z0TLAT_W(1)=',Z0TLAT_W(1)
            write (F_unout,*) 'Z0TLAT_W(2)=',Z0TLAT_W(2)
            write (F_unout,*) 'Z0TLAT_R(1)=',Z0TLAT_R(1)
            write (F_unout,*) 'Z0TLAT_R(2)=',Z0TLAT_R(2)
            write (F_unout,9900) 'INCONSISTENT Z0TLAT VALUES GEM VS OCEAN'
            call flush(F_unout)
            err = -1
            goto 995
         endif

         err = 0 
 
         z0mtype_R = trim(s_send(2))
         z0ttype_R = trim(s_send(3))

         write (F_unout,*) 'z0mtype (GEM)   =',F_z0mtype,len(F_z0mtype)
         write (F_unout,*) 'z0mtype_R (NEMO)  =',z0mtype_R,len(z0mtype_R)
         write (F_unout,*) 'z0ttype (GEM)   =',F_z0ttype,len(F_z0ttype)
         write (F_unout,*) 'z0ttype_R (NEMO)  =',z0ttype_R,len(z0ttype_R)

         if (trim(z0mtype_R) /= trim(F_z0mtype)) then
            write (F_unout,9900) 'INCONSISTENT z0mtype VALUES GEM VS OCEAN'
            write (F_unout,*) 'z0mtype    =',F_z0mtype
            write (F_unout,*) 'z0mtype_R  =',z0mtype_R
            call flush(F_unout)
            err = -1
            goto 995
         endif

         if (trim(z0ttype_R) /= trim(F_z0ttype)) then
            write (F_unout,9900) 'INCONSISTENT z0ttype VALUES GEM VS OCEAN'
            write (F_unout,*) 'z0ttype    =',F_z0ttype
            write (F_unout,*) 'z0ttype_R  =',z0ttype_R
            call flush(F_unout)
            err = -1
            goto 995
         endif

         if (i_send(1) /= F_dateo) then
            if (F_print_L) then
               write (F_unout,9900) 'WARNING: MODEL INITIAL TIME INCONSISTENT'
               write (F_unout,*) 'i_send(1) (ocean)  =',i_send(1)
               write (F_unout,*) 'F_dateo   (atmos)  =',F_dateo
            endif
         endif

      endif

      allocate ( ocn_busou(cpl_drv_lni,cpl_drv_lnj,cplocn_n_fldou) )
      ocn_busou = 0.

      gmmk_ocn_busin_s = 'ocn_busin'
      gmmk_gli_0_s     = 'gli_0'
      gmmk_i8i_0_s     = 'i8i_0'
      gmmk_sdi_0_s     = 'sdi_0'
      gmmk_tmo_0_s     = 'tmo_0'

      nullify(ocn_busin,gli_0,i8i_0,sdi_0,tmo_0)
      call gmm_build_meta3D(meta3d_ocn_busin,                    &
                            1,cpl_drv_lni,0,0,cpl_drv_lni,       &
                            1,cpl_drv_lnj,0,0,cpl_drv_lnj,       &
                            1,cplocn_n_fldin,0,0,cplocn_n_fldin, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta2d_surf0,                        &
                            1,cpl_drv_lni,0,0,cpl_drv_lni,       &
                            1,cpl_drv_lnj,0,0,cpl_drv_lnj,       &
                            0,GMM_NULL_FLAGS)

      ier = gmm_create(gmmk_ocn_busin_s , ocn_busin , meta3d_ocn_busin , GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_gli_0_s ,     gli_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_i8i_0_s ,     i8i_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_sdi_0_s ,     sdi_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_tmo_0_s ,     tmo_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)

      ier = gmm_get(gmmk_ocn_busin_s , ocn_busin)
      ier = gmm_get(gmmk_gli_0_s ,     gli_0)
      ier = gmm_get(gmmk_i8i_0_s ,     i8i_0)
      ier = gmm_get(gmmk_sdi_0_s ,     sdi_0)
      ier = gmm_get(gmmk_tmo_0_s ,     tmo_0)

      if ( .NOT. cpl_rstn_L ) then
         ocn_busin(:,:,:) = 0.
         gli_0(:,:) = cplocn_missval
         i8i_0(:,:) = cplocn_missval
         sdi_0(:,:) = cplocn_missval
         tmo_0(:,:) = cplocn_missval
       else
         if (F_print_L) write (F_unout,*) 'CPLOCN_INIT: RESTART MODE'
      endif

      do ivar=1,cplocn_n_fldin

         SELECT CASE ( cplocn_cvin_S(ivar) )

         CASE ( 'MCP' ) ; icvin_MCP=ivar
         CASE ( 'ALO' ) ; icvin_ALO=ivar
         CASE ( 'ALI' ) ; icvin_ALI=ivar
         CASE ( 'T4O' ) ; icvin_T4O=ivar
         CASE ( 'T4I' ) ; icvin_T4I=ivar
         CASE ( 'SHO' ) ; icvin_SHO=ivar
         CASE ( 'SHI' ) ; icvin_SHI=ivar
         CASE ( 'LHO' ) ; icvin_LHO=ivar
         CASE ( 'LHI' ) ; icvin_LHI=ivar
         CASE ( 'TXO' ) ; icvin_TXO=ivar
         CASE ( 'TYO' ) ; icvin_TYO=ivar
         CASE ( 'TXI' ) ; icvin_TXI=ivar
         CASE ( 'TYI' ) ; icvin_TYI=ivar
         CASE ( 'ZTO' ) ; icvin_ZTO=ivar
         CASE ( 'ZQO' ) ; icvin_ZQO=ivar
         CASE ( 'ZUO' ) ; icvin_ZUO=ivar
         CASE ( 'ZVO' ) ; icvin_ZVO=ivar
         CASE ( 'ZTI' ) ; icvin_ZTI=ivar
         CASE ( 'ZQI' ) ; icvin_ZQI=ivar
         CASE ( 'ZUI' ) ; icvin_ZUI=ivar
         CASE ( 'ZVI' ) ; icvin_ZVI=ivar
         CASE ( 'GLI' ) ; icvin_GLI=ivar
         CASE ( 'I8I' ) ; icvin_I8I=ivar
         CASE ( 'SDI' ) ; icvin_SDI=ivar
         CASE ( 'TMO' ) ; icvin_TMO=ivar
         CASE ( 'UUO' ) ; icvin_UUO=ivar
         CASE ( 'VVO' ) ; icvin_VVO=ivar
         CASE ( 'I7I' ) ; icvin_I7I=ivar
         CASE ( 'UUI' ) ; icvin_UUI=ivar
         CASE ( 'VVI' ) ; icvin_VVI=ivar
         CASE ( 'QSO' ) ; icvin_QSO=ivar
         CASE ( 'QSI' ) ; icvin_QSI=ivar
         CASE ( 'ILO' ) ; icvin_ILO=ivar
         CASE ( 'ILI' ) ; icvin_ILI=ivar
         CASE ( 'ZMO' ) ; icvin_ZMO=ivar
         CASE ( 'ZMI' ) ; icvin_ZMI=ivar
         CASE ( 'ZHO' ) ; icvin_ZHO=ivar
         CASE ( 'ZHI' ) ; icvin_ZHI=ivar

         CASE DEFAULT

           ier=-1
           call handle_error(ier,'cplocn_init','WRONG NAME TAG')

         END SELECT

      enddo

      call cplocn_distwgt (weightfile_S,F_print_L,F_unout)

 995  call handle_error(err,'cplocn_init','Problems with cplocn_init')

      call RPN_COMM_bcast (cplocn_1st_L, n1, "MPI_LOGICAL", n0,"GRID",ier)
      call RPN_COMM_bcast (cplocn_off_L, n1, "MPI_LOGICAL", n0,"GRID",ier)
      call RPN_COMM_bcast (cplocn_oc_dt, n1, "MPI_INTEGER", n0,"GRID",ier)

      do i = 1, cplocn_n_fldin
        cnt = len(cplocn_cvin_S(i)) 
        call RPN_COMM_bcastc(cplocn_cvin_S(i), cnt, "MPI_CHARACTER", n0,"GRID",ier)
        cnt = len(cplocn_cvit_S(i))
        call RPN_COMM_bcastc(cplocn_cvit_S(i), cnt, "MPI_CHARACTER", n0,"GRID",ier)
      enddo

      do i = 1, cplocn_n_fldou
        cnt = len(cplocn_cvou_S(i))
        call RPN_COMM_bcastc(cplocn_cvou_S(i), cnt, "MPI_CHARACTER", n0,"GRID",ier)
        cnt = len(cplocn_cvot_S(i))
        call RPN_COMM_bcastc(cplocn_cvot_S(i), cnt, "MPI_CHARACTER", n0,"GRID",ier)
      enddo

      if (cplocn_off_L) then
         if (F_print_L) then
            write(F_unout,*) 'cplocn_init: WARNING - cplocn_off_L !!!!'
            write(F_unout,*) 'resetting CPLOCN to FALSE'
         endif
         cplocn_init = 0
      else
         cplocn_init = 1
      endif
!
      if (F_print_L) write(F_unout,2001)

 2000 format( &
      /,'INITIALIZATION OF COUPLING INTERFACE S/R cplocn_init', &
      /,'=====================================================')
 2001 format( &
      /,'INITIALIZATION OF COUPLING INTERFACE ENDED S/R cplocn_init', &
      /,'===========================================================')
 9900 format (/,1x,a)
!
!     ---------------------------------------------------------------
!
      return
      end function cplocn_init
