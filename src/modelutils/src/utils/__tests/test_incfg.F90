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

#include <rmn/msg.h>

!/@
subroutine test_incfg()
   use, intrinsic :: iso_fortran_env, only: INT64
   use testutils
   use incfg_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>
   integer,parameter :: NSTEPMAX = 4
   integer,parameter :: NVARMAX = 12

   character(len=512) :: filename_S
   integer :: istat
   ! ---------------------------------------------------------------------
   call testutils_verbosity()

   filename_S = './__to-delete-test_incfg_table__'
   call testutils_set_name('test_incfg1')
   call test_incfg_1(trim(filename_S)//'1')
   call testutils_set_name('test_incfg2')
   call test_incfg_2(trim(filename_S)//'2')
   ! ---------------------------------------------------------------------
   return

contains
   
   
   !/@
   subroutine test_incfg_1(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer,parameter :: NSTEP = 4
      integer,parameter :: NVAR = 12
      character(len=8) :: result_S(NSTEP+1,NVAR)

      integer :: istat,fileid
      character(len=512) :: filename_S
      !invalid/ignored lines:
      ! '#IN=SD;	SEARCH=ANAL+CLIM;	 INTERP=NEAREST;	FREQ=1,MONTH,12'
      ! ''
      ! 'SEARCH=ANAL+CLIM;	INTERP=NEAREST;FREQ=0; timeint=near'
      ! 'IN=tw;INTERP=NEAREST;FREQ=0; timeint=near'
      !valid lines:
      ! 'IN=V0; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
      ! 'IN=V0b; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
      ! 'IN=V1;SEARCH=CLIM;	FREQ=1;'
      ! 'IN=V01;SEARCH=CLIM;	FREQ=0,1;'
      ! 'IN=V213;SEARCH=CLIM;	FREQ=step,2,1,3;levels=1,26'
      ! 'IN=b123;SEARCH=CLIM;	FREQ=1,2,3'
      ! 'IN=vh01;SEARCH=ANAL,CLIM;	INTERP=NEAREST;FREQ=HOUR,0,1; timeint=near;toignore=2'
      ! 'IN=vh1;SEARCH=ANAL,CLIM;	INTERP=NEAREST;FREQ=HOUR,1'
      ! 'IN=vmi1;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,1'
      ! 'IN=vm60;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,60'
      ! 'IN=m030;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,0,30'
      ! 'IN=m015;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,0,15;levels=-1'
      ! ---------------------------------------------------------------------
!!$      skip_list_S(:) = (/'v0b','q1'/)

      fileid = 0
      istat = fnom(fileid ,F_filename_S,'SEQ/FMT+R/W',0)
      write(fileid,'(a)',iostat=istat) &
           '#IN=SNODP;	SEARCH=ANAL+CLIM;	 INTERP=NEAREST;	FREQ=1,MONTH,12'
      write(fileid,'(a)',iostat=istat) &
           ''
      write(fileid,'(a)',iostat=istat) &
           'SEARCH=ANAL+CLIM;	INTERP=NEAREST;FREQ=0; timeint=nearest'
      write(fileid,'(a)',iostat=istat) &
           'IN=twater;INTERP=NEAREST;FREQ=0; timeint=nearest'
      write(fileid,'(a)',iostat=istat) &
           'IN=V0; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
!!$      write(fileid,'(a)',iostat=istat) &
!!$           'IN=V0b; SEARCH=ANAL;	 interp=NEAREST;	FREQ=0'
      write(fileid,'(a)',iostat=istat) &
           'IN=V1;SEARCH=CLIM;	FREQ=1;'
      write(fileid,'(a)',iostat=istat) &
           'IN=V01;SEARCH=CLIM;	FREQ=0,1;'
      write(fileid,'(a)',iostat=istat) &
           'IN=V213;SEARCH=CLIM;	FREQ=step,2,1,3;levels=1,26'
      write(fileid,'(a)',iostat=istat) &
           'IN=b123;SEARCH=CLIM;	FREQ=1,2,3'
      write(fileid,'(a)',iostat=istat) &
           'IN=vh01;SEARCH=ANAL,CLIM;	INTERP=NEAREST;FREQ=HOUR,0,1; timeint=nearest;toignore=2'
      write(fileid,'(a)',iostat=istat) &
           'IN=vh1;SEARCH=ANAL,CLIM;	INTERP=NEAREST;FREQ=HOUR,1'
      write(fileid,'(a)',iostat=istat) &
           'IN=vmi1;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,1'
      write(fileid,'(a)',iostat=istat) &
           'IN=vm60;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,60'
      write(fileid,'(a)',iostat=istat) &
           'IN=m030;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,0,30;levels=-1'
      write(fileid,'(a)',iostat=istat) &
           'IN=m015;SEARCH=ANAL;	INTERP=NEAREST;FREQ=MINUTe,0,15'
      istat = fclos(fileid)

      result_S = ' '
      result_S(1,1:4) = (/'v0  ','v01 ','vh01','m030'/)
      result_S(2,1:4) = (/'v1  ','v01 ','b123','m030'/)
      result_S(3,1:6) = (/'v01 ','v213','vh01','vh1 ','vm60','m030'/)
      result_S(4,1:4) = (/'v01 ','v213','b123','m030'/)
      result_S(5,1:3) = (/'v01 ','vh01','m030'/)

      call test_incfg_chk(F_filename_S,result_S,9,NSTEP,'1')
      istat = clib_unlink(trim(F_filename_S))
      ! ---------------------------------------------------------------------
      return
   end subroutine test_incfg_1

   !/@
   subroutine test_incfg_2(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer,parameter :: NSTEP = 4
      integer,parameter :: NVAR = 9
      character(len=8) :: result_S(NSTEP+1,NVAR)

      integer :: istat,fileid
      character(len=512) :: filename_S
      !'IN=LA ; FREQ=0; SEARCH=ANAL ; timeint=near'
      !'IN=LO ; FREQ=1; SEARCH=ANAL ; timeint=near'
      !'IN=AL ; FREQ=0; SEARCH=CLIM'
      !'IN=MG ; FREQ=0; SEARCH=GEOP'
      !'IN=VG ; FREQ=0; SEARCH=ANAL,GEOP'
      !'IN=SD ; FREQ=0; SEARCH=CLIM ; timeint=linear'
      !'IN=TM ; FREQ=0,1; SEARCH=CLIM ; timeint=linear'
      !'IN=VF ; FREQ=0; SEARCH=GEOP ; lvl=1,26'
      ! ---------------------------------------------------------------------
      fileid = 0
      istat = fnom(fileid ,F_filename_S,'SEQ/FMT+R/W',0)
      write(fileid,'(a)',iostat=istat) &
           'IN=LA ; FREQ=0; SEARCH=ANAL ; timeint=near'
      write(fileid,'(a)',iostat=istat) &
           'IN=LO ; FREQ=2; SEARCH=ANAL ; timeint=step'
      write(fileid,'(a)',iostat=istat) &
           'IN=AL ; FREQ=0; SEARCH=CLIM'
      write(fileid,'(a)',iostat=istat) &
           'IN=MG ; FREQ=0; SEARCH=GEOP'
      write(fileid,'(a)',iostat=istat) &
           'IN=VG ; FREQ=0; SEARCH=ANAL,GEOP'
      write(fileid,'(a)',iostat=istat) &
           'IN=SD ; FREQ=0; SEARCH=CLIM ; timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'IN=TM ; FREQ=0,2; SEARCH=CLIM ; timeint=next'
      write(fileid,'(a)',iostat=istat) &
           'IN=VF ; FREQ=0; SEARCH=GEOP ; levels=1,26'
      write(fileid,'(a)',iostat=istat) &
           'IN=UU; IN2=VV ; FREQ=0; SEARCH=ANAL '
      istat = fclos(fileid)

      result_S = ' '
      result_S(1,1:9) = (/'la','al','mg','vg','sd','tm','vf','uu','vv'/)
!!$      result_S(0,1:7) = (/'la','al','mg','vg','sd','tm','vf'/)
      result_S(3,1:2) = (/'lo','tm'/)
      result_S(5,1) = 'tm'

      call test_incfg_chk(F_filename_S,result_S,9,NSTEP,'2')
      istat = clib_unlink(trim(filename_S))
      ! ---------------------------------------------------------------------
      return
   end subroutine test_incfg_2


   !/@
   subroutine test_incfg_chk(F_filename_S,result_S,F_nvar,F_nstep,F_msg_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S,F_msg_S
      integer,intent(in) :: F_nvar,F_nstep
      character(len=8),intent(in) :: result_S(F_nstep+1,F_nvar)
      !@/
      character(len=8) :: result2_S(F_nstep+1,F_nvar*2)
      character(len=512) :: dateo_S,varname_S,varname2_S,skip_list_S(2),string_S
      integer :: istat,incfgid,dateo,istep,index,dt,ii,nbvar,ivar,lvl0,lvl1,ip1list(99),ip1list01(99),ip1list02(99),nip1,ipkind
      real :: zp1
      ! ---------------------------------------------------------------------
      result2_S = ' '

      skip_list_S(:) = (/'v0b','q1 '/)
      dateo_S = '20090427.000000'
      call datp2f(dateo,dateo_S)
      dt = 1800
      ipkind = RMN_CONV_ARBITRARY
      do ii=1,26
         zp1 = real(ii)
         call convip_plus(ip1list01(ii),zp1,ipkind,RMN_CONV_P2IPNEW,' ',.not.RMN_CONV_USEFORMAT_L)
      enddo
      ip1list02(1:3) = (/123,234,345/)
      incfgid = incfg_new(dateo,dt,F_filename_S,ip1list02(1:3))
      nbvar = incfg_nbvar(incfgid)
      call testutils_assert_eq(nbvar,F_nvar,'nbvar')
      write(string_S,*) 'incfg file opened: ',trim(F_filename_S),incfgid,nbvar,F_nstep
      call msg(MSG_DEBUG,string_S)
      STEPLOOP: do istep = 0,F_nstep
         ii = 0
         VARLOOP: do ivar=1,nbvar
            istat = incfg_isvarstep(incfgid,ivar,istep)
            if (.not.RMN_IS_OK(istat)) cycle VARLOOP
            ip1list = -1
            istat = incfg_meta(incfgid,ivar,varname_S,varname2_S,F_lvl0=lvl0,F_lvl1=lvl1,F_ip1list=ip1list,F_nip1=nip1)
            if (.not.RMN_IS_OK(istat)) cycle VARLOOP
            if (any(varname_S == skip_list_S)) then
               print *,'skipping '//trim(varname_S)
               cycle VARLOOP
            endif
            ii = ii + 1
            result2_S(istep+1,ii) = varname_S
            if (result2_S(istep+1,ii) /= result_S(istep+1,ii)) &
                 print *,'ERROR v1',istep,ii,varname_S(1:8),index
            if (varname2_S/=' ') then
               ii = ii + 1
               result2_S(istep+1,ii) = varname2_S
               if (result2_S(istep+1,ii) /= result_S(istep+1,ii)) &
                    print *,'ERROR v2',istep,ii,varname2_S(1:8),index
            endif
            if (varname_S=='v213'.or.varname_S=='vf') then
               call testutils_assert_eq((/lvl0,lvl1/),(/1,26/),trim(varname_S)//' lvl0/1')
               call testutils_assert_eq(nip1,26,trim(varname_S)//' nip1')
               call testutils_assert_eq(ip1list(1:26),ip1list01(1:26),trim(varname_S)//' ip1list')
            elseif (varname_S=='m030') then
               call testutils_assert_eq((/lvl0,lvl1/),(/-1,-1/),trim(varname_S)//' lvl0/1')
               call testutils_assert_eq(nip1,3,trim(varname_S)//' nip1')
               call testutils_assert_eq(ip1list(1:3),ip1list02(1:3),trim(varname_S)//' ip1list')
            else 
               call testutils_assert_eq((/lvl0,lvl1/),(/0,0/),trim(varname_S)//' lvl0/1')
               call testutils_assert_eq(nip1,1,trim(varname_S)//' nip1')
               call testutils_assert_eq(ip1list(1),0,trim(varname_S)//' ip1list')
            endif
         enddo VARLOOP
      enddo STEPLOOP
      call testutils_assert_ok(all(result_S(1:F_nstep+1,1:F_nvar)==result2_S(1:F_nstep+1,1:F_nvar)),F_msg_S)
      !TODO: test incfg_varindex, incfg_add
      ! ---------------------------------------------------------------------
      return
   end subroutine test_incfg_chk

end subroutine test_incfg
