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
subroutine test_input2_mpi()
use iso_c_binding
   use testutils
   use clib_itf_mod
   use input_mod
   use vGrid_Descriptors
   use vgrid_wb
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"
   integer,parameter :: NDIGITS = 4
   character(len=512) :: dir_S,filename_S,dfiles_S,bcmk_S, &
        anal_S,anal2_S,clim_S,clim2_S,geop_S,geop2_S,inrep_S,inrep2_S
   integer :: istat,myproc
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   call testutils_verbosity()
   call testutils_set_name('test_input2_mpi')
   istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
   call msg_set_p0only(0)

   istat = clib_getenv('ATM_MODEL_DFILES',dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_ERROR,'ATM_MODEL_DFILES not defined')
      return
   endif

   write(dir_S,'(a,i4.4)') '__to-delete_test_input2_mpi-',myproc
   istat = clib_mkdir(trim(dir_S))
   istat = clib_chdir(trim(dir_S))

   bcmk_S = trim(dfiles_S)//'/bcmk/'
   anal_S = trim(dfiles_S)//'/bcmk/2009042700_000'
!!$!   anal_S = '/users/dor/armn/sch/storage_model/zeta/ANALYSIS-mid'
!!$   anal_S = '/users/dor/armn/sch/storage_model/zeta/ANALYSIS-mid2'
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

   filename_S='../test_input2_table'
   call test_input2_get_data(trim(filename_S))

   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))

   istat = clib_chdir('./..')
   istat = clib_unlink(trim(dir_S))

   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return

contains

   !/@
   subroutine test_input2_get_data(F_filename_S)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      !@/
      integer,parameter :: NSTEP = 1
      integer,parameter :: NVAR = 32
      integer,parameter :: NI = 8
      integer,parameter :: NJ = 4
      integer,parameter :: NK = 2
      real,parameter :: hyb(4) = (/0.5492443,0.7299818,0.8791828,0.9950425/)
      character(len=512) :: dateo_S,varname_S,varname2_S,skip_list_S(4),dummy_S
      character(len=8) :: result_S(0:NSTEP,NVAR),result2_S(0:NSTEP,NVAR)
      integer :: istat,inputid,dateo,istep,index,dt,ii,gridid,nbvar,ivar,nhyb,k,ikind,vals(8),allvals(8,NSTEP+1,NVAR),ip1list(NK),nn
      logical :: ok_L(0:NSTEP,NVAR)
      real :: hyblist(NK),fact
      real(8) :: ptop_8
      real, pointer :: data(:,:,:),data2(:,:,:),p0data(:,:)
      integer,pointer :: ip1listt(:)
      type(vgrid_descriptor) :: vgrid1
      type(gmm_metadata) :: p0meta
      ! ---------------------------------------------------------------------
      result_S = ' '
      result2_S = ' '
      allvals = 0
      
      result_S(0,1:NVAR) = (/ &
           'mg  ','tm  ','lh  ','y7  ','y8  ', &
           'y9  ','ga  ','i8  ','i9  ','i7  ', &
           'zp  ','sd  ','i0  ','vf  ','lg  ', &
           'hs  ','vg  ','me  ','al  ','i1  ', &
           'i2  ','i3  ','i4  ','i6  ','j1  ', &
           'j2  ','dn  ','icel','al  ','tt  ', &
           'uu  ','vv  '/)

!!$      skip_list_S(1:3) = (/'xa  ','icel','fsa '/)
      
      allvals(:,1, 1) = (/          0,      10000, 0, 0, 0, 0, 0, 0/) !mg
      allvals(:,1, 2) = (/      27137,      30229, 0, 0, 0, 0, 0, 0/) !tm
      allvals(:,1, 3) = (/          0,      70995, 0, 0, 0, 0, 0, 0/) !lh
      allvals(:,1, 4) = (/          0,       2449, 0, 0, 0, 0, 0, 0/) !y7
      allvals(:,1, 5) = (/          0,       8466, 0, 0, 0, 0, 0, 0/) !y8
      allvals(:,1, 6) = (/       -457,       2925, 0, 0, 0, 0, 0, 0/) !y9
      allvals(:,1, 7) = (/          0,       2082, 0, 0, 0, 0, 0, 0/) !ga
      allvals(:,1, 8) = (/        -38,      19767, 0, 0, 0, 0, 0, 0/) !i8
      allvals(:,1, 9) = (/      25391,      27313, 0, 0, 0, 0, 0, 0/) !i9
      allvals(:,1,10) = (/      25804,      27312, 0, 0, 0, 0, 0, 0/) !i7
      allvals(:,1,11) = (/     -69078,      14006, 0, 0, 0, 0, 0, 0/) !zp
      allvals(:,1,12) = (/      -3096,      73904, 0, 0, 0, 0, 0, 0/) !sd
      allvals(:,1,13) = (/      23491,      30197, 0, 0, 0, 0, 0, 0/) !i0
      allvals(:,1,14) = (/          0,      10000, 0, 0, 0, 0, 0, 0/) !vf
      allvals(:,1,15) = (/          0,       9529, 0, 0, 0, 0, 0, 0/) !lg
      allvals(:,1,16) = (/         33,      10000, 0, 0, 0, 0, 0, 0/) !hs
      allvals(:,1,17) = (/       1000,      26000, 0, 0, 0, 0, 0, 0/) !vg
      allvals(:,1,18) = (/       -294,      92369, 0, 0, 0, 0, 0, 0/) !me
      allvals(:,1,19) = (/        600,       4010, 0, 0, 0, 0, 0, 0/) !al
      allvals(:,1,20) = (/         -5,       9998, 0, 0, 0, 0, 0, 0/) !i1
      allvals(:,1,21) = (/         -1,       1813, 0, 0, 0, 0, 0, 0/) !i2
      allvals(:,1,22) = (/        -11,       5536, 0, 0, 0, 0, 0, 0/) !i3
      allvals(:,1,23) = (/       -282,       1593, 0, 0, 0, 0, 0, 0/) !i4
      allvals(:,1,24) = (/       5088,       7998, 0, 0, 0, 0, 0, 0/) !i6
      allvals(:,1,25) = (/          0,      65914, 0, 0, 0, 0, 0, 0/) !j1
      allvals(:,1,26) = (/          0,      36881, 0, 0, 0, 0, 0, 0/) !j2
      allvals(:,1,27) = (/       9999,      30299, 0, 0, 0, 0, 0, 0/) !dn
      allvals(:,1,28) = (/          0,      10000, 0, 0, 0, 0, 0, 0/) !ic
      allvals(:,1,29) = (/        600,       4010, 0, 0, 0, 0, 0, 0/) !al
      allvals(:,1,30) = (/     -14073,      32563, 0, 0, 0, 0, 0, 0/) !tt
      allvals(:,1,31) = (/     -29433,      40536, 0, 0, 0, 0, 0, 0/) !uu
      allvals(:,1,32) = (/     -23352,      24180, 0, 0, 0, 0, 0, 0/) !vv

      ptop_8 = 9575.0d0
      istat = vgd_new(vgrid1, &
           kind     = VGRID_HYBS_KIND, &
           version  = VGRID_HYBS_VER, &
           hyb      = hyb, &
           ptop_8   = ptop_8, &
           pref_8   = 100000.d0, &
           rcoef1   = 1., &
           rcoef2   = 1.)
      istat = vgd_get(vgrid1,key='VIPT',value=ip1listt)
      istat = vgrid_wb_put('gemthlvl',vgrid1,ip1listt,'P0')
      call gmm_build_meta2D(p0meta, 1,NI,0,0,NI, 1,NJ,0,0,NJ, 0,GMM_NULL_FLAGS)
      nullify(p0data)
      istat = gmm_create('P0',p0data,p0meta)
      p0data = 99900.

      ok_L = .false.
      gridid = ezqkdef(NI,NJ, 'G', 0,0,0,0,0)
!      dateo_S = '20100107.120000'
!!$      dateo_S = '20081215.000000'
      dateo_S = '20090427.000000'
      call datp2f(dateo,dateo_S)
      nhyb = 2
      hyblist(1:2) = (/0.985,0.995/)
      do k=1,nhyb
         ikind = RMN_CONV_HY
         call convip_plus(ip1list(k),hyblist(k),ikind,RMN_CONV_P2IPNEW,'',.false.)
      enddo
      dt = 1800
      inputid = input_new(dateo,dt,F_filename_S,ip1list(1:nhyb))
      if (.not.RMN_IS_OK(index)) then
         call msg(MSG_ERROR,'test_input3 - cannot init input_mod for: '//trim(F_filename_S))
         return
      endif
      nbvar = input_nbvar(inputid)
      call testutils_assert_eq(nbvar,34,'input_nbvar')

      STEPLOOP: do istep = 0,NSTEP
         ii = 0
         VARLOOP: do ivar=1,nbvar
            istat = input_isvarstep(inputid,ivar,istep)
            if (.not.RMN_IS_OK(istat)) cycle VARLOOP
            istat = input_meta(inputid,ivar,varname_S,varname2_S)
            if (any(varname_S == skip_list_S)) cycle VARLOOP
            nullify(data,data2)
            istat = input_get(inputid,ivar,istep,gridid,'gemthlvl',data,data2)
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
               nn = log10(abs(maxval(data2))+1e-5)
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

      call testutils_assert_ok(all(result_S==result2_S),'test_input2_get_data','var list')
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
            call testutils_assert_ok(ok_L(istep,ii),'test_input2_get_data','values for "'//trim(result_S(istep,ii))//'" at step='//trim(dummy_S))   
         enddo
      enddo
      istat = input_close_files()
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input2_get_data

end subroutine test_input2_mpi
