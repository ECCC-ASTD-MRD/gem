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
subroutine test_input()
   use, intrinsic :: iso_fortran_env, only: INT64
   use testutils
   use input_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>
   integer,parameter :: NDIGITS = 4
   logical,parameter :: FSTOPC_SET = .false.
   character(len=512) :: filename_S,dfiles_S,bcmk_S, &
        anal_S,anal2_S,clim_S,clim2_S,geop_S,geop2_S,inrep_S,inrep2_S
   integer :: istat,myproc
   ! ---------------------------------------------------------------------
   call testutils_verbosity()
   call testutils_set_name('test_input')
   istat = fstopc('MSGLVL','SYSTEM',FSTOPC_SET)

   istat = clib_getenv('ATM_MODEL_DFILES',dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_ERROR,'ATM_MODEL_DFILES not defined')
      return
   endif
   bcmk_S = trim(dfiles_S)//'/bcmk/'
   anal_S = trim(dfiles_S)//'/bcmk/2009042700_000'
   anal2_S = './ANALYSIS' !'./analysis'
   clim_S = trim(dfiles_S)//'/bcmk/climato'
   clim2_S = './CLIMATO'
   geop_S = trim(dfiles_S)//'/bcmk/geophy/Gem_geophy.fst'
   geop2_S = './GEOPHY'
   inrep_S = trim(dfiles_S)//'/bcmk'
   inrep2_S = './INREP'
   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))
   istat = clib_symlink(trim(anal_S),trim(anal2_S))
   istat = min(clib_symlink(trim(clim_S),trim(clim2_S)),istat)
   istat = min(clib_symlink(trim(geop_S),trim(geop2_S)),istat)
   istat = min(clib_symlink(trim(inrep_S),trim(inrep2_S)),istat)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_WARNING,'problem creating symlink')
   endif

   filename_S = './__to-delete-test_input_table__'
   call test_input_write2(trim(filename_S))
   call test_input_get_data(trim(filename_S))
   call test_input_write3(trim(filename_S)//'2')
   call test_input_get_data3(trim(filename_S)//'2')

   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))
   istat = clib_unlink(trim(filename_S))
   istat = clib_unlink(trim(filename_S)//'2')
   ! ---------------------------------------------------------------------
   return

contains

   !/@
   subroutine test_input_write2(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer :: istat,fileid
      character(len=512) :: filename_S
      ! ---------------------------------------------------------------------
      fileid = 0
      istat = fnom(fileid ,F_filename_S,'SEQ/FMT+R/W',0)

      write(fileid,'(a)',iostat=istat) &
           'IN=ME ; FREQ=0; SEARCH=GEOP'
      write(fileid,'(a)',iostat=istat) &
           'IN=MG ; FREQ=0; SEARCH=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'IN=VF ; FREQ=0; SEARCH=GEOP ; levels=1,26'
      write(fileid,'(a)',iostat=istat) &
           'IN=VG ; FREQ=1; SEARCH=ANAL,GEOP;interp=nearest'

      write(fileid,'(a)',iostat=istat) &
           'IN=ts ; FREQ=0; SEARCH=CLIM'
      write(fileid,'(a)',iostat=istat) &
           'IN=al ; FREQ=1; SEARCH=CLIM ; timeint=any'
      write(fileid,'(a)',iostat=istat) &
           'IN=TM ; FREQ=0,1; SEARCH=ANAL,CLIM ; timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'IN=HS ; FREQ=1; SEARCH=ANAL,CLIM ; timeint=near'

      write(fileid,'(a)',iostat=istat) &
           'IN=P0 ; FREQ=0; SEARCH=ANAL'
      write(fileid,'(a)',iostat=istat) &
           'IN=LA ; FREQ=0; SEARCH=ANAL ; timeint=near'
      write(fileid,'(a)',iostat=istat) &
           'IN=LO ; FREQ=1; SEARCH=ANAL ; timeint=any'
      write(fileid,'(a)',iostat=istat) &
           'IN=al ; FREQ=1; SEARCH=INREP ; timeint=any'
      istat = fclos(fileid)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input_write2


   !/@
   subroutine test_input_write3(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer :: istat,fileid
      ! ---------------------------------------------------------------------
      fileid = 0
      istat = fnom(fileid ,F_filename_S,'SEQ/FMT+R/W',0)
      write(fileid,'(a)',iostat=istat) &
           'in=MG;search=GEOP;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=TM;search=ANAL;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=LH;search=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=Y7;search=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=Y8;search=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=Y9;search=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=GA;search=GEOP;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=I8;search=ANAL,CLIM;interp=nearest;timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=I9;search=ANAL;interp=linear;levels= 1, 2'
      write(fileid,'(a)',iostat=istat) &
           'in=I7;search=ANAL;interp=linear;levels= 1, 3'
      write(fileid,'(a)',iostat=istat) &
           'in=ZP;search=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=SD;search=ANAL,CLIM;interp=nearest;levels= 1, 5;timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=I0;search=ANAL,CLIM;interp=linear;levels= 1, 2;timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=VF;search=GEOP;interp=nearest;levels= 1,26'
      write(fileid,'(a)',iostat=istat) &
           'in=LG;search=ANAL,CLIM;interp=nearest;timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=ME;search=GEOP;interp=cubic'
      write(fileid,'(a)',iostat=istat) &
           'in=AL;search=GEOP,CLIM;interp=nearest ; timeint=any'
!!$           'in=AL;search=GEOP,CLIM;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=I1;search=ANAL;interp=nearest;levels= 1, 2'
      write(fileid,'(a)',iostat=istat) &
           'in=I2;search=ANAL;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=I3;search=ANAL;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=I4;search=ANAL;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=I6;search=ANAL;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=XA;search=ANAL;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=J1;search=GEOP;interp=nearest;levels= 1, 3'
      write(fileid,'(a)',iostat=istat) &
           'in=J2;search=GEOP;interp=nearest;levels= 1, 3'
      write(fileid,'(a)',iostat=istat) &
           'in=DN;search=ANAL;interp=nearest'
      write(fileid,'(a)',iostat=istat) &
           'in=ICEL;search=ANAL,CLIM;interp=nearest;timeint=linear'
      write(fileid,'(a)',iostat=istat) &
           'in=FSA;search=GEOP;interp=linear'
      write(fileid,'(a)',iostat=istat) &
           'IN=UU; SEARCH=ANAL '
      write(fileid,'(a)',iostat=istat) &
           'IN=LO; IN2=LA ; FREQ=0; SEARCH=ANAL '
      write(fileid,'(a)',iostat=istat) &
           'IN=TM ; IN2=HS ; FREQ=1; SEARCH=CLIM ; timeint=linear'
!!$      write(fileid,'(a)',iostat=istat) &
!!$           'IN=UU; IN2=VV ; FREQ=0; SEARCH=ANAL '
!!$      write(fileid,'(a)',iostat=istat) &
!!$           'IN=UU; IN2=VV ; FREQ=1; SEARCH=ANAL;timeint=step'
      istat = fclos(fileid)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input_write3


   !/@
   subroutine test_input_get_data(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer,parameter :: NSTEP = 1
      integer,parameter :: NVAR = 12
      integer,parameter :: NI = 8
      integer,parameter :: NJ = 4
      character(len=512) :: dateo_S,varname_S,skip_list_S(2),dummy_S
      character(len=8) :: result_S(0:NSTEP,NVAR),result2_S(0:NSTEP,NVAR)
      integer :: istat,inputid,dateo,istep,index,dt,ii,gridid,nbvar,ivar
      logical :: ok_L(0:NSTEP,NVAR)
      real, pointer :: data(:,:,:)
      integer :: vals(8),me(8,NSTEP+1),mg(8,NSTEP+1),vf(8,NSTEP+1),ts(8,NSTEP+1),tm(8,NSTEP+1),p0(8,NSTEP+1),la(8,NSTEP+1),vg(8,NSTEP+1),al(8,NSTEP+1),hs(8,NSTEP+1),lo(8,NSTEP+1),allvals(8,NSTEP+1,NVAR),nn
      real :: fact
      ! ---------------------------------------------------------------------
      result_S = ' '
      result2_S = ' '
      allvals = -999

!!$      me(:,1) = (/     -52352,     932156, 4, 3, 3, 7, 1, 1/)
!!$      mg(:,1) = (/          0,       1000, 1, 1, 4, 1, 1, 1/)
!!$      vf(:,1) = (/        -45,       1045, 8, 8, 2, 2, 4, 1/)
!!$      ts(:,1) = (/     -11159,      30049, 1, 8, 2, 7, 1, 4/)
!!$      tm(:,1) = (/     271278,     302299, 1, 7, 6, 7, 1, 1/)
!!$      p0(:,1) = (/      91663,     102534, 8, 8, 2, 6, 3, 3/)
!!$      la(:,1) = (/     -59470,      59383, 1, 1, 1, 1, 1, 1/)
!!$
!!$      vg(:,2) = (/       1000,      26000, 1, 1, 3, 1, 1, 1/)
!!$      al(:,2) = (/         52,        400, 1, 1, 3, 8, 2, 1/)
!!$      tm(:,2) = (/      -1830,      29602, 1, 1, 6, 7, 1, 1/)
!!$      hs(:,2) = (/          6,       1027, 1, 4, 1, 2, 3, 3/)
!!$      lo(:,2) = (/    -180004,     134871, 5, 5, 5, 5, 1, 1/)
      allvals(:,1, 1) = (/       -294,      92370, 0, 0, 0, 0, 0, 0/) !me
      allvals(:,1, 2) = (/          0,      10000, 0, 0, 0, 0, 0, 0/) !mg
      allvals(:,1, 3) = (/       -279,      10495, 0, 0, 0, 0, 0, 0/) !vf
      allvals(:,1, 4) = (/     -11116,      30155, 0, 0, 0, 0, 0, 0/) !ts
!!$      allvals(:,1, 5) = (/      27134,      30234, 0, 0, 0, 0, 0, 0/) !tm
!!$      allvals(:,1, 6) = (/       8853,      10243, 0, 0, 0, 0, 0, 0/) !p0
!!$      allvals(:,1, 7) = (/     -60016,      59922, 0, 0, 0, 0, 0, 0/) !la
      allvals(:,1, 5) = (/      27128,      30230, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,1, 6) = (/       9166,      10253, 0, 0, 0, 0, 0, 0/) !p0
      allvals(:,1, 7) = (/     -59470,      59383, 0, 0, 0, 0, 0, 0/) !la
      allvals(:,2, 1) = (/       1000,      26000, 0, 0, 0, 0, 0, 0/) !vg
      allvals(:,2, 2) = (/        600,       4010, 0, 0, 0, 0, 0, 0/) !al
!!$      allvals(:,2, 3) = (/      27134,      30234, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,2, 3) = (/      27128,      30230, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,2, 4) = (/         11,      10001, 0, 0, 0, 0, 0, 0/) !hs
      allvals(:,2, 5) = (/     -18000,      13487, 0, 0, 0, 0, 0, 0/) !lo
      allvals(:,2, 6) = (/        600,       4010, 0, 0, 0, 0, 0, 0/) !al

      result_S(0,1:7) = (/'me','mg','vf','ts','tm','p0','la'/)
      result_S(1,1:6) = (/'vg','al','tm','hs','lo','al'/)

      ok_L = .false.
      gridid = ezqkdef(NI,NJ, 'G', 0,0,0,0,0)
      dateo_S = '20090427.000000'
      call datp2f(dateo,dateo_S)
      dt = 1800
      inputid = input_new(dateo,dt,F_filename_S)
      call testutils_assert_ok(RMN_IS_OK(inputid),'1 input_new status')
      nbvar = input_nbvar(inputid)
      call testutils_assert_ok(nbvar==NVAR,'1 input_nbvar')

      STEPLOOP: do istep = 0,NSTEP
         ii = 0
         VARLOOP: do ivar=1,nbvar
!!$            istat = input_meta(inputid,ivar,varname_S)
!!$            call testutils_assert_ok(RMN_IS_OK(istat),'1 input_meta status')
            istat = input_isvarstep(inputid,ivar,istep)
            if (.not.RMN_IS_OK(istat)) then
!!$               print *,istep,ivar,' ',trim(varname_S),' Not Needed'
               cycle VARLOOP
!!$            else
!!$               print *,istep,ivar,' ',trim(varname_S),' Needed'
            endif
            nullify(data)
            istat = input_get(inputid,ivar,istep,gridid,data,varname_S)
            if (.not.(RMN_IS_OK(istat) .and. associated(data))) then
               call msg(MSG_WARNING,'test_input - var not found: '//trim(varname_S))
               ok_L = .false.
               cycle VARLOOP !exit STEPLOOP
            endif
            ii = ii + 1
            result2_S(istep,ii) = varname_S
!!$            vals = (/nint(1000.*minval(data)),nint(1000.*maxval(data)),minloc(data,1),minloc(data,2),minloc(data,2),maxloc(data,1),maxloc(data,2),maxloc(data,3)/)
            nn = log10(abs(maxval(data))+1e-5)
            fact = 10.**(NDIGITS-nn)
            vals(1:8) = (/nint(fact*minval(data)),nint(fact*maxval(data)),0,0,0,0,0,0/)
            if (all(vals(1:2) == allvals(1:2,istep+1,ii))) ok_L(istep,ii) = .true.
            !TODO: check minloc and maxloc on AIX... 
!!$            print '(a,i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)',varname_S(1:2)//"= (/",vals,"/)"
!!$            select case(varname_S(1:2))
!!$            case('me')
!!$               if (all(vals(1:2) == me(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('mg')
!!$               if (all(vals(1:2) == mg(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('vf')
!!$               if (all(vals(1:2) == vf(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('ts')
!!$               if (all(vals(1:2) == ts(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('tm')
!!$               if (all(vals(1:2) == tm(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('p0')
!!$               vals(1:2) = vals(1:2) / 10 !AIX v/s Linux last digit difference
!!$               if (all(vals(1:2) == p0(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('la')
!!$               if (all(vals(1:2) == la(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('lo')
!!$               if (all(vals(1:2) == lo(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('vg')
!!$               if (all(vals(1:2) == vg(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('al')
!!$               if (all(vals(1:2) == al(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            case('hs')
!!$               if (all(vals(1:2) == hs(1:2,istep+1))) ok_L(istep,ii) = .true.
!!$            end select
            deallocate(data,stat=istat)
            if (.not.ok_L(istep,ii)) then
               print '(a,i1,a,i2,a,i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)',"allvals(:,",istep+1,",",ii,") = (/",vals,"/) !"//varname_S(1:2)
            endif

         enddo VARLOOP
      enddo STEPLOOP
      call testutils_assert_ok(all(result_S==result2_S),'1 var list')
      if (.not.all(result_S==result2_S)) then
         do istep=0,ubound(result_S,1)
            do ivar=1,size(result_S,2)
               if (result_S(istep,ivar)/=''.or.result2_S(istep,ivar)/='') print *,istep,ivar,result_S(istep,ivar)==result2_S(istep,ivar),result_S(istep,ivar),result2_S(istep,ivar)
            enddo
         enddo
      endif
      do istep = 0,NSTEP
         do ii=1,NVAR
            if (result_S(istep,ii) == ' ') cycle
            write(dummy_S,'(i4)') istep
            call testutils_assert_ok(ok_L(istep,ii),'1 values for "'//trim(result_S(istep,ii))//'" at step='//trim(dummy_S))
         enddo
      enddo
      istat = input_close_files()
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input_get_data


   !/@
   subroutine test_input_get_data3(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer,parameter :: NSTEP = 1
      integer,parameter :: NVAR = 29
      integer,parameter :: NI = 8
      integer,parameter :: NJ = 4
      character(len=512) :: dateo_S,varname_S,varname2_S,skip_list_S(4),dummy_S
      character(len=8) :: result_S(0:NSTEP,NVAR),result2_S(0:NSTEP,NVAR)
      integer :: istat,inputid,dateo,istep,index,dt,ii,gridid,nbvar,ivar
      logical :: ok_L(0:NSTEP,NVAR)
      real, pointer :: data(:,:,:),data2(:,:,:)
      integer :: vals(8),allvals(8,NSTEP+1,NVAR),nn
      real :: fact
      ! ---------------------------------------------------------------------
      result_S = ' '
      result2_S = ' '
      allvals = -999
      
      result_S(0,1:27) = (/'mg','tm','lh','y7','y8','y9','ga','i8','i9','i7','zp','sd','i0','vf','lg','me','al','i1','i2','i3','i4','i6','j1','j2','dn','lo','la'/)
      result_S(1,1:2) = (/'tm','hs'/)
!!$      result_S(0,1:27) = (/'mg','tm','lh','y7','y8','y9','ga','i8','i9','i7','zp','sd','i0','vf','lg','me','al','i1','i2','i3','i4','i6','j1','j2','dn','uu','vv'/)
!!$      result_S(1,1:2) = (/'uu','vv'/)
      skip_list_S(1:3) = (/'xa  ','icel','fsa '/)
      
!!$      allvals(:,1, 1) = (/          0,       1000, 1, 1, 4, 1, 1, 1/) !mg
!!$      allvals(:,1, 2) = (/     271372,     302295, 1, 7, 6, 7, 1, 1/) !tm
!!$      allvals(:,1, 3) = (/          0,     932722, 1, 1, 4, 1, 1, 1/) !lh
!!$      allvals(:,1, 4) = (/          0,          2, 1, 1, 4, 7, 1, 1/) !y7
!!$      allvals(:,1, 5) = (/          0,          1, 1, 1, 4, 7, 1, 1/) !y8
!!$      allvals(:,1, 6) = (/          0,          0, 1, 8, 1, 8, 3, 1/) !y9
!!$      allvals(:,1, 7) = (/          0,        602, 1, 1, 1, 1, 1, 1/) !ga
!!$      allvals(:,1, 8) = (/         -4,       1977, 1, 1, 1, 5, 1, 1/) !i8
!!$      allvals(:,1, 9) = (/     253906,     273126, 1, 1, 1, 7, 8, 1/) !i9
!!$      allvals(:,1,10) = (/     258042,     273122, 1, 1, 1, 7, 4, 1/) !i7
!!$      allvals(:,1,11) = (/      -6908,       1426, 1, 1, 5, 1, 1, 1/) !zp
!!$      allvals(:,1,12) = (/      -3096,      73904, 1, 1, 1, 1, 1, 1/) !sd
!!$      allvals(:,1,13) = (/     234906,     301969, 3, 7, 6, 7, 3, 7/) !i0
!!$      allvals(:,1,14) = (/          0,       1000, 8, 4, 1, 2, 1, 1/) !vf
!!$      allvals(:,1,15) = (/          0,        948, 1, 1, 1, 1, 1, 1/) !lg
!!$      allvals(:,1,16) = (/     -52352,     932156, 4, 3, 3, 7, 1, 1/) !me
!!$      allvals(:,1,17) = (/         60,        401, 1, 1, 3, 8, 2, 1/) !al
!!$      allvals(:,1,18) = (/          0,       1000, 1, 4, 2, 6, 1, 4/) !i1
!!$      allvals(:,1,19) = (/          0,        181, 1, 1, 1, 1, 1, 1/) !i2
!!$      allvals(:,1,20) = (/         -1,        554, 1, 8, 1, 3, 1, 1/) !i3
!!$      allvals(:,1,21) = (/        -28,        159, 1, 1, 1, 1, 1, 1/) !i4
!!$      allvals(:,1,22) = (/        509,        800, 1, 1, 1, 3, 4, 4/) !i6
!!$      allvals(:,1,23) = (/          0,      62904, 1, 1, 4, 4, 1, 1/) !j1
!!$      allvals(:,1,24) = (/          0,      37208, 1, 1, 4, 7, 1, 1/) !j2
!!$      allvals(:,1,25) = (/      99992,     302992, 1, 1, 3, 1, 1, 1/) !dn
!!$      allvals(:,1,26) = (/      0,     0, 1, 1, 3, 1, 1, 1/) !uu
!!$      allvals(:,1,27) = (/      0,     0, 1, 1, 3, 1, 1, 1/) !vv
!!$      allvals(:,2, 1) = (/      0,     0, 1, 1, 3, 1, 1, 1/) !uu
!!$      allvals(:,2, 2) = (/      0,     0, 1, 1, 3, 1, 1, 1/) !vv
      allvals(:,1, 1) = (/          0,      10000, 0, 0, 0, 0, 0, 0/) !mg
!!$      allvals(:,1, 2) = (/      27134,      30234, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,1, 2) = (/      27137,      30229, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,1, 3) = (/          0,      70995, 0, 0, 0, 0, 0, 0/) !lh
      allvals(:,1, 4) = (/          0,       2449, 0, 0, 0, 0, 0, 0/) !y7
      allvals(:,1, 5) = (/          0,       8466, 0, 0, 0, 0, 0, 0/) !y8
      allvals(:,1, 6) = (/       -457,       2925, 0, 0, 0, 0, 0, 0/) !y9
      allvals(:,1, 7) = (/          0,       2082, 0, 0, 0, 0, 0, 0/) !ga
      allvals(:,1, 8) = (/        -38,      19767, 0, 0, 0, 0, 0, 0/) !i8
!!$      allvals(:,1, 9) = (/      25383,      27313, 0, 0, 0, 0, 0, 0/) !i9
!!$      allvals(:,1,10) = (/      25784,      27312, 0, 0, 0, 0, 0, 0/) !i7
      allvals(:,1, 9) = (/      25391,      27313, 0, 0, 0, 0, 0, 0/) !i9
      allvals(:,1,10) = (/      25804,      27312, 0, 0, 0, 0, 0, 0/) !i7
      allvals(:,1,11) = (/     -69078,      14006, 0, 0, 0, 0, 0, 0/) !zp
      allvals(:,1,12) = (/      -3096,      73904, 0, 0, 0, 0, 0, 0/) !sd
!!$      allvals(:,1,13) = (/      23340,      30109, 0, 0, 0, 0, 0, 0/) !i0
      allvals(:,1,13) = (/      23491,      30197, 0, 0, 0, 0, 0, 0/) !i0
      allvals(:,1,14) = (/          0,      10000, 0, 0, 0, 0, 0, 0/) !vf
      allvals(:,1,15) = (/          0,       9529, 0, 0, 0, 0, 0, 0/) !lg
      allvals(:,1,16) = (/       -294,      92370, 0, 0, 0, 0, 0, 0/) !me
      allvals(:,1,17) = (/        600,       4010, 0, 0, 0, 0, 0, 0/) !al
      allvals(:,1,18) = (/         -5,       9998, 0, 0, 0, 0, 0, 0/) !i1
      allvals(:,1,19) = (/         -1,       1813, 0, 0, 0, 0, 0, 0/) !i2
      allvals(:,1,20) = (/        -11,       5536, 0, 0, 0, 0, 0, 0/) !i3
      allvals(:,1,21) = (/       -282,       1593, 0, 0, 0, 0, 0, 0/) !i4
      allvals(:,1,22) = (/       5088,       7998, 0, 0, 0, 0, 0, 0/) !i6
      allvals(:,1,23) = (/          0,      65914, 0, 0, 0, 0, 0, 0/) !j1
      allvals(:,1,24) = (/          0,      36881, 0, 0, 0, 0, 0, 0/) !j2
      allvals(:,1,25) = (/       9999,      30299, 0, 0, 0, 0, 0, 0/) !dn
      allvals(:,1,26) = (/     -18000,      13487, 0, 0, 0, 0, 0, 0/) !lo
!!$      allvals(:,1,27) = (/     -60016,      59922, 0, 0, 0, 0, 0, 0/) !la
!!$      allvals(:,2, 1) = (/      27134,      30234, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,1,27) = (/     -59470,      59383, 0, 0, 0, 0, 0, 0/) !la
      allvals(:,2, 1) = (/      27128,      30230, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,2, 2) = (/         33,      10000, 0, 0, 0, 0, 0, 0/) !hs

      ok_L = .false.
      gridid = ezqkdef(NI,NJ, 'G', 0,0,0,0,0)
      dateo_S = '20090427.000000'
      call datp2f(dateo,dateo_S)
      dt = 1800
      inputid = input_new(dateo,dt,F_filename_S)
      if (.not.RMN_IS_OK(index)) then
         call msg(MSG_ERROR,'test_input3 - cannot init input_mod for: '//trim(F_filename_S))
         return
      endif
      nbvar = input_nbvar(inputid)

      STEPLOOP: do istep = 0,NSTEP
         ii = 0
         VARLOOP: do ivar=1,nbvar
            istat = input_isvarstep(inputid,ivar,istep)
            if (.not.RMN_IS_OK(istat)) cycle VARLOOP
            istat = input_meta(inputid,ivar,varname_S,varname2_S)
            if (any(varname_S == skip_list_S)) cycle VARLOOP
            nullify(data,data2)
            istat = input_get(inputid,ivar,istep,gridid,data,data2)
            if (.not.(RMN_IS_OK(istat) .and. associated(data))) then
               call msg(MSG_WARNING,'test_input3 - var not found: '//trim(varname_S))
               cycle VARLOOP !exit STEPLOOP
            endif
            ii = ii + 1
            result2_S(istep,ii) = varname_S
            
!!$            vals = (/nint(1000.*minval(data)),nint(1000.*maxval(data)),minloc(data,1),minloc(data,2),minloc(data,2),maxloc(data,1),maxloc(data,2),maxloc(data,3)/)
!!$            if (all(vals(:) == allvals(:,istep+1,ii))) ok_L(istep,ii) = .true.
!!$            if (.not.ok_L(istep,ii)) then
!!$               print '(a,i1,a,i2,a,i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)',"allvals(:,",istep+1,",",ii,") = (/",vals,"/) !"//varname_S(1:2)
!!$            endif

            nn = log10(abs(maxval(data))+1e-5)
            fact = 10.**(NDIGITS-nn)
            vals(1:8) = (/nint(fact*minval(data)),nint(fact*maxval(data)),0,0,0,0,0,0/)
            if (all(vals(1:2) == allvals(1:2,istep+1,ii))) ok_L(istep,ii) = .true.
            if (.not.ok_L(istep,ii)) then
               print '(a,i1,a,i2,a,i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)',"allvals(:,",istep+1,",",ii,") = (/",vals,"/) !"//varname_S(1:2)
            endif
!!$            print *,'READ: ',istep,ivar,ii,ok_L(istep,ii),trim(varname_S)//':'//trim(varname2_S),associated(data),associated(data2)

            if (varname2_S/=' ') then
               ii = ii + 1
               result2_S(istep,ii) = varname2_S
!!$               if (result2_S(istep,ii) /= result_S(istep,ii)) &
!!$                    print *,'ERROR v2',istep,ii,varname2_S(1:8),index
               nn = log10(abs(maxval(data2))+1.e-5)
               fact = 10.**(NDIGITS-nn)
               vals(1:8) = (/nint(fact*minval(data2)),nint(fact*maxval(data2)),0,0,0,0,0,0/)
               if (all(vals(1:2) == allvals(1:2,istep+1,ii))) ok_L(istep,ii) = .true.
               if (.not.ok_L(istep,ii)) then
                  print '(a,i1,a,i2,a,i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)',"allvals(:,",istep+1,",",ii,") = (/",vals,"/) !"//varname2_S(1:2)
               endif
            endif
            if (associated(data)) deallocate(data,stat=istat)
            if (associated(data2)) deallocate(data2,stat=istat)
         enddo VARLOOP
      enddo STEPLOOP

      call testutils_assert_ok(all(result_S==result2_S),'2 var list')
      if (.not.all(result_S==result2_S)) then
         do istep=0,ubound(result_S,1)
            do ivar=1,size(result_S,2)
               if (result_S(istep,ivar)/=''.or.result2_S(istep,ivar)/='') print *,istep,ivar,result_S(istep,ivar)==result2_S(istep,ivar),result_S(istep,ivar),result2_S(istep,ivar)
            enddo
         enddo
      endif
      do istep = 0,NSTEP
         do ii=1,NVAR
            if (result_S(istep,ii) == ' ') cycle
            write(dummy_S,'(i4)') istep
            call testutils_assert_ok(ok_L(istep,ii),'2 values for "'//trim(result_S(istep,ii))//'" at step='//trim(dummy_S))   
         enddo
      enddo
      istat = input_close_files()
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input_get_data3

end subroutine test_input
