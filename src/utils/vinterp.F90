!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

!/@*
module vinterp_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use vGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use vgrid_wb
   use tdpack
   use samevgrid_mod
   use rmn_gmm
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2011-04
   !@description
   ! Public functions
   public :: vte_intvertx, vte_intvertx3, vte_intvertx4, vte_intvertx_isodst
   public :: vinterp,vinterp0,vinterp1,vinterp10,vinterp01
   ! Public constants
   character(len=*), parameter, public :: VINTERP_DIAG_CUBIC = 'cubic'
   character(len=*), parameter, public :: VINTERP_DIAG_LINEAR = 'linear'
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface vte_intvertx
      module procedure vte_intvertx3
      module procedure vte_intvertx4
   end interface vte_intvertx
   
   interface vinterp
      module procedure vinterp0
      module procedure vinterp0ls
      module procedure vinterp1
      module procedure vinterp10
      module procedure vinterp01
   end interface

   real, parameter :: EPSILON_R4 = 1.e-6
   real, parameter :: P0MIN = 10000.
   real, parameter :: P0MAX = 200000.
   real, parameter :: H0MIN = -1000.  !# AdHoc negative value
   real, parameter :: H0MAX = 10000.


contains
   
   !/@*
   subroutine vte_intvertx3(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,F_type)
      implicit none
      !@ojective cubic vertical interpolation from eta/sigma to hybrid
      !@arguments
      ! F_dch                   Interpolated field
      ! F_sch                   source       field
      ! F_srclev                source levels
      ! F_dstlev                destination levels
      ! F_var_S                 name of the variable
      ! F_type                  order of interpolation ('linear' or 'cubic')
      integer,intent(in) :: n, nkd, nks
      real :: F_dch(n,nkd), F_sch(n,nks), F_srclev(n,nks),F_dstlev(n,nkd)
      character(len=*) :: F_var_S,F_type
      !---------------------------------------------------------------
      if (any(F_type(1:1) == (/'l','L'/))) then
         call vte_intvertx4(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,nkd)
      else
         call vte_intvertx4(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,0)
      endif
      !---------------------------------------------------------------
      return
   end subroutine vte_intvertx3


   !/@*
   subroutine vte_intvertx4(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,F_nlinbot)
      implicit none
      !@ojective cubic vertical interpolation from eta/sigma to hybrid
      !@arguments
      ! F_dch                   Interpolated field
      ! F_sch                   source       field
      ! F_srclev                source levels
      ! F_dstlev                destination levels
      ! F_var_S                 name of the variable
      ! F_nlinbot               number of levels at bottom linearly interpolated
      integer,intent(in) :: n, nkd, nks,F_nlinbot
      real :: F_dch(n,nkd), F_sch(n,nks), F_srclev(n,nks),F_dstlev(n,nkd)
      character(len=*) :: F_var_S
      !@author - Methot/Laroche - April 97 - v1_01
      !@revision
      ! v2_10 - L. Corbeil          - rewrited for optimization,
      ! v2_10                         removed e_vcubique
      ! v2_30 - L. Corbeil          - renamed vte_ (no more called in gemntr)
      ! v2_31 - Lee V.              - F_srclev,F_dstlev is calculated outside
      ! v3_32 - Plante A.           - correction of F_dstlev(n,nks) to F_dstlev(n,nkd)
      !@description
      !     General case, we'll care about VT later
      !
      !     First, we find the level we are by squeezing the destination
      !     between increasing bot() and decreasing top(). We need log_2_(nks)
      !     iteration to squeeze it completely (each integer between 1 and
      !     nks can be expressed as a sum of power of 2 such as sum c_i 2^i
      !     where c_i = 0 or 1 and i lower or equal than log_2_(nks)
      !
      !     WARNING
      !  niter calculation is ok for nks.lt. 2097152: should be ok for now...
      !  (Maybe the grid will be that precise in 2010!) (I don't even bother
      !  to add an if statement, for performance purpose...)
      !*@/
      integer :: i,k,iter,niter,lev,lev_lin,k_surf,k_ciel
      integer :: top(n),bot(n),topcub(n),botcub(n),ref(n)
      real(REAL64) :: deltalev,prxd,prda,prdb,prsaf,prsbf,prsad,prsbd
      logical :: ascending_L
      !---------------------------------------------------------------
      ascending_L = (F_srclev(1,1) < F_srclev(1,nks))
      k_surf = 1
      k_ciel = nks
      if (.not.ascending_L) then
         k_surf = nks
         k_ciel = 1
      endif
      if (real(int(log(real(nks))/log(2.0))).eq. &
           log(real(nks))/log(2.0)) then
         niter=int(log(real(nks))/log(2.0))
      else
         niter=int(log(real(nks))/log(2.0))+1
      endif

      !- squeeze...

      do k=1,nkd
         top(1:n)=nks
         bot(1:n)=1
         if (ascending_L) then
            do iter=1,niter
               do i=1,n
                  !     divide by two (the old fashioned way...)
                  ref(i)=ishft(top(i)+bot(i),-1)
                  !     adjust top or bot
                  if(F_dstlev(i,k).lt.F_srclev(i,ref(i))) then
                     top(i)=ref(i)
                  else
                     bot(i)=ref(i)
                  endif
               enddo
            enddo
         else
            do iter=1,niter
               do i=1,n
                  !     divide by two (the old fashioned way...)
                  ref(i)=ishft(top(i)+bot(i),-1)
                  !     adjust top or bot
                  if(F_dstlev(i,k).gt.F_srclev(i,ref(i))) then
                     top(i)=ref(i)
                  else
                     bot(i)=ref(i)
                  endif
               enddo
            enddo
         endif
         !- adjusting top and bot to ensure we can perform cubic interpolation
         do i=1,n
            botcub(i)=max(2,bot(i))
            topcub(i)=min(nks-1,top(i))
         enddo
         !- cubic or linear interpolation
         do i=1,n
            lev=botcub(i)
            lev_lin=bot(i)
            deltalev=(F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin))

            !- Interpolation: either if not enough points to perform cubic interp
            !                 or if linear interpolation is requested use linear
            !                 interpolation and constant extrapolation
            if((lev.ne.lev_lin).or.(topcub(i).ne.top(i)).or.k<=F_nlinbot) then
               !- persistancy of this interval
               if(F_dstlev(i,k).le.F_srclev(i,k_surf)) then
                  F_dch(i,k) = F_sch(i,k_surf)
               else if(F_dstlev(i,k).ge.F_srclev(i,k_ciel)) then
                  F_dch(i,k) = F_sch(i,k_ciel)
               else
                  !- linear interpolation
                  prxd=(F_dstlev(i,k)-F_srclev(i,lev_lin))/deltalev
                  F_dch(i,k) = (1.0-prxd)*F_sch(i,lev_lin) &
                       +prxd*F_sch(i,lev_lin+1)
               endif
            else
               !- cubic interpolation
               prxd = (F_dstlev(i,k)-F_srclev(i,lev_lin))/ &
                    deltalev
               prda = ((F_sch(i,lev_lin+1)-F_sch(i,lev_lin-1))/ &
                    (F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin-1))* &
                    deltalev)
               prdb = ((F_sch(i,lev_lin+2)-F_sch(i,lev_lin))/ &
                    (F_srclev(i,lev_lin+2)-F_srclev(i,lev_lin))* &
                    deltalev)
               prsaf= (1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
               prsbf= (3.0-2.0*prxd)*prxd*prxd
               prsad= prxd*(1.0-prxd)*(1.0-prxd)
               prsbd= (1.0-prxd)*prxd*prxd
               F_dch(i,k) = F_sch(i,lev_lin  )*prsaf &
                    +F_sch(i,lev_lin+1)*prsbf+prda*prsad &
                    -prdb*prsbd
            endif
         enddo
      enddo

      !- special case for VT
      if(F_var_S == 'VT') then
         do k=1,nkd
            do i=1,n
               if(F_srclev(i,k_ciel).lt.F_dstlev(i,k)) then
                  F_dch(i,k) = F_sch(i,k_ciel) * exp (  &
                       rgasd_8*stlo_8*(F_dstlev(i,k)-F_srclev(i,k_ciel)) )
               endif
            enddo
         enddo
      endif
      !     ---------------------------------------------------------------
      return
   end subroutine vte_intvertx4

   
   !/@*
   subroutine vte_intvertx_isodst(F_dch,F_sch,F_srclev,F_dstlev,n,nks,nkd,F_var_S,F_nlinbot)
      implicit none
      !@ojective cubic vertical interpolation from eta/sigma to hybrid
      !@arguments
      ! F_dch                   Interpolated field
      ! F_sch                   source       field
      ! F_srclev                source levels
      ! F_dstlev                destination levels
      ! F_var_S                 name of the variable
      ! F_nlinbot               number of levels at bottom linearly interpolated
      integer,intent(in) :: n, nkd, nks,F_nlinbot
      real :: F_dch(n,nkd), F_sch(n,nks), F_srclev(n,nks),F_dstlev(nkd)
      character(len=*) :: F_var_S
      !@author - Methot/Laroche - April 97 - v1_01
      !@revision
      ! v2_10 - L. Corbeil          - rewrited for optimization,
      ! v2_10                         removed e_vcubique
      ! v2_30 - L. Corbeil          - renamed vte_ (no more called in gemntr)
      ! v2_31 - Lee V.              - F_srclev,F_dstlev is calculated outside
      ! v3_32 - Plante A.           - correction of F_dstlev(n,nks) to F_dstlev(n,nkd)
      !@description
      !     General case, we'll care about VT later
      !
      !     First, we find the level we are by squeezing the destination
      !     between increasing bot() and decreasing top(). We need log_2_(nks)
      !     iteration to squeeze it completely (each integer between 1 and
      !     nks can be expressed as a sum of power of 2 such as sum c_i 2^i
      !     where c_i = 0 or 1 and i lower or equal than log_2_(nks)
      !
      !     WARNING
      !  niter calculation is ok for nks.lt. 2097152: should be ok for now...
      !  (Maybe the grid will be that precise in 2010!) (I don't even bother
      !  to add an if statement, for performance purpose...)
      !*@/
      integer :: i,k,iter,niter,lev,lev_lin,k_surf,k_ciel
      integer :: top(n),bot(n),topcub(n),botcub(n),ref(n)
      real(REAL64) :: deltalev,prxd,prda,prdb,prsaf,prsbf,prsad,prsbd
      logical :: ascending_L
      !---------------------------------------------------------------
      ascending_L = (F_srclev(1,1) < F_srclev(1,nks))
      k_surf = 1
      k_ciel = nks
      if (.not.ascending_L) then
         k_surf = nks
         k_ciel = 1
      endif
      if (real(int(log(real(nks))/log(2.0))).eq. &
           log(real(nks))/log(2.0)) then
         niter=int(log(real(nks))/log(2.0))
      else
         niter=int(log(real(nks))/log(2.0))+1
      endif

      !- squeeze...

      do k=1,nkd
         top(1:n)=nks
         bot(1:n)=1
         if (ascending_L) then
            do iter=1,niter
               do i=1,n
                  !     divide by two (the old fashioned way...)
                  ref(i)=ishft(top(i)+bot(i),-1)
                  !     adjust top or bot
                  if(F_dstlev(k).lt.F_srclev(i,ref(i))) then
                     top(i)=ref(i)
                  else
                     bot(i)=ref(i)
                  endif
               enddo
            enddo
         else
            do iter=1,niter
               do i=1,n
                  !     divide by two (the old fashioned way...)
                  ref(i)=ishft(top(i)+bot(i),-1)
                  !     adjust top or bot
                  if(F_dstlev(k).gt.F_srclev(i,ref(i))) then
                     top(i)=ref(i)
                  else
                     bot(i)=ref(i)
                  endif
               enddo
            enddo
         endif
         !- adjusting top and bot to ensure we can perform cubic interpolation
         do i=1,n
            botcub(i)=max(2,bot(i))
            topcub(i)=min(nks-1,top(i))
         enddo
         !- cubic or linear interpolation
         do i=1,n
            lev=botcub(i)
            lev_lin=bot(i)
            deltalev=(F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin))

            !- Interpolation: either if not enough points to perform cubic interp
            !                 or if linear interpolation is requested use linear
            !                 interpolation and constant extrapolation
            if((lev.ne.lev_lin).or.(topcub(i).ne.top(i)).or.k<=F_nlinbot) then
               !- persistancy of this interval
               if(F_dstlev(k).le.F_srclev(i,k_surf)) then
                  F_dch(i,k) = F_sch(i,k_surf)
               else if(F_dstlev(k).ge.F_srclev(i,k_ciel)) then
                  F_dch(i,k) = F_sch(i,k_ciel)
               else
                  !- linear interpolation
                  prxd=(F_dstlev(k)-F_srclev(i,lev_lin))/deltalev
                  F_dch(i,k) = (1.0-prxd)*F_sch(i,lev_lin) &
                       +prxd*F_sch(i,lev_lin+1)
               endif
            else
               !- cubic interpolation
               prxd = (F_dstlev(k)-F_srclev(i,lev_lin))/ &
                    deltalev
               prda = ((F_sch(i,lev_lin+1)-F_sch(i,lev_lin-1))/ &
                    (F_srclev(i,lev_lin+1)-F_srclev(i,lev_lin-1))* &
                    deltalev)
               prdb = ((F_sch(i,lev_lin+2)-F_sch(i,lev_lin))/ &
                    (F_srclev(i,lev_lin+2)-F_srclev(i,lev_lin))* &
                    deltalev)
               prsaf= (1.0+2.0*prxd)*(1.0-prxd)*(1.0-prxd)
               prsbf= (3.0-2.0*prxd)*prxd*prxd
               prsad= prxd*(1.0-prxd)*(1.0-prxd)
               prsbd= (1.0-prxd)*prxd*prxd
               F_dch(i,k) = F_sch(i,lev_lin  )*prsaf &
                    +F_sch(i,lev_lin+1)*prsbf+prda*prsad &
                    -prdb*prsbd
            endif
         enddo
      enddo

      !- special case for VT
      if(F_var_S == 'VT') then
         do k=1,nkd
            do i=1,n
               if(F_srclev(i,k_ciel).lt.F_dstlev(k)) then
                  F_dch(i,k) = F_sch(i,k_ciel) * exp (  &
                       rgasd_8*stlo_8*(F_dstlev(k)-F_srclev(i,k_ciel)) )
               endif
            enddo
         enddo
      endif
      !     ---------------------------------------------------------------
      return
   end subroutine vte_intvertx_isodst

   
   !/@*
   function vinterp0( &
        F_dataout, F_vgridout, F_ip1listout, &
        F_datain, F_vgridin, F_ip1listin, &
        F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgridout,F_vgridin
      integer,intent(in) :: F_ip1listout(:),F_ip1listin(:)
      real,pointer :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S, tmp_S
      integer :: nlinbot
      real,pointer :: sfcfldout2(:,:), sfcfldin2(:,:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp0 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      if (.not.associated(F_datain)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, null pointer: '//trim(msg_S))
         return
      endif
      nlinbot = 0
      if (present(F_nlinbot)) nlinbot = F_nlinbot
      nullify(sfcfldout2, sfcfldin2)
      F_istat = vinterp0ls( &
           F_dataout, F_vgridout, F_ip1listout, F_sfcfldout, sfcfldout2, &
           F_datain, F_vgridin, F_ip1listin, F_sfcfldin, sfcfldin2, &
           F_nlinbot, F_msg_S)
      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp0 [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp0


   !/@*
   function vinterp0ls( &
        F_dataout, F_vgridout, F_ip1listout, F_sfcfldout, F_sfcfldout2,&
        F_datain, F_vgridin, F_ip1listin, F_sfcfldin, F_sfcfldin2, &
        F_nlinbot, F_msg_S, F_levelsout, F_levelsin, F_altin_S, F_altout_S) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgridout,F_vgridin
      integer,intent(in) :: F_ip1listout(:),F_ip1listin(:)
      real,pointer :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S, F_altin_S, F_altout_S
      real,pointer,optional :: F_levelsout(:,:,:), F_levelsin(:,:,:)
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S, tmp_S, altin_S, altout_S
      integer,pointer :: samevert_list(:)
      integer :: istat, nlinbot, nijkout(3), nijkin(3), lijk(3), uijk(3), k
      integer, dimension(size(F_datain,dim=3),1) :: slist
      real,pointer :: levelsin(:,:,:), levelsout(:,:,:), sfcfldin(:,:), sfcfldin2(:,:)
      real,dimension(size(F_datain,dim=3)) :: scol
      real,dimension(size(F_datain,dim=1),size(F_datain,dim=2),size(F_datain,dim=3)) :: sdatain, slevelsin

      logical :: use_same_sfcfld_L, samevert_L, ispressin_L, ispressout_L, usealt_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp0ls [BEGIN]')
      F_istat = RMN_ERR
      
      nlinbot = 0
      if (present(F_nlinbot)) nlinbot = F_nlinbot
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nullify(levelsout, levelsin)
      if (present(F_levelsout)) then
         if (associated(F_levelsout)) levelsout => F_levelsout
      endif
      if (present(F_levelsin)) then
         if (associated(F_levelsin)) levelsin => F_levelsin
      endif
      altin_S = ''
      if (present(F_altin_S)) altin_S = F_altin_S
      altout_S = ''
      if (present(F_altout_S)) altout_S = F_altout_S

      if (.not.associated(F_datain)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, null input data pointer: '//trim(msg_S))
         return
      endif
      
      ispressin_L = vgrid_wb_is_press_kind(F_vgridin)
      ispressout_L = vgrid_wb_is_press_kind(F_vgridout)
      usealt_L = (.not.(ispressin_L .eqv. ispressout_L))
      IF_SAME_VTYPE: if (usealt_L) then
         if (altout_S /= ' ' .and. .not.associated(levelsout)) then
            istat = gmm_get(altout_S, levelsout)
            if (.not.RMN_IS_OK(istat)) nullify(levelsout)
            if (associated(levelsout)) &
                 call msg(MSG_INFOPLUS, '(vinterp) Using output alt vcoor desc: '//trim(altout_S))
         else if (altin_S /= ' ' .and. .not.associated(levelsin)) then
            istat = gmm_get(altin_S, levelsin)
            if (.not.RMN_IS_OK(istat)) nullify(levelsin)
            if (associated(levelsin)) &
                 call msg(MSG_INFOPLUS, '(vinterp) Using input alt vcoor desc: '//trim(altin_S))
         endif
         if (.not.(associated(levelsin) .or. associated(levelsout))) then
            if (ispressout_L) then
               call msg(MSG_WARNING,'(vinterp) Cannot Interpolate from Height based to Pressure based vert. coor.')            
            else
               call msg(MSG_WARNING,'(vinterp) Cannot Interpolate from Pressure based to Height based vert. coor.')
            endif
            return
         endif
      endif IF_SAME_VTYPE

      if (.not.associated(F_dataout)) then
         lijk = lbound(F_datain)
         uijk = ubound(F_datain)
         if (associated(levelsout)) then
            uijk(3) = size(levelsout, 3)
         else
            uijk(3) = size(F_ip1listout)
         endif
         nijkout = shape(F_datain)
         nijkout(3) = uijk(3) - lijk(3) + 1
      else
         nijkout = shape(F_dataout)
      endif

      nijkin  = shape(F_datain)
      if (.not.all(nijkout(1:2) == nijkin(1:2))) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, not same shape: '//trim(msg_S))
         print *,'(vinterp) nijkin(1:2) =', nijkin(1:2)
         print *,'(vinterp) nijkout(1:2)=', nijkout(1:2)
         call flush(6)
         return
      endif

      !#TODO: what about PRESS coor where sfcfld is not avail?
      if (.not.associated(levelsin)) then
         nullify(sfcfldin2)
         if (associated(F_sfcfldin) .and. .not.associated(F_sfcfldin,F_sfcfldout)) then
            use_same_sfcfld_L = .false.
            sfcfldin => F_sfcfldin
            if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
         else if (associated(F_sfcfldout)) then
            use_same_sfcfld_L = .true.
            sfcfldin => F_sfcfldout
            if (associated(F_sfcfldout2)) sfcfldin2 => F_sfcfldout2
!!$         else
!!$            call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, missing input sfc field')
!!$            return
         endif
      endif
      
      IF_SAME_VTYPE2: if (.not.usealt_L) then

         allocate(samevert_list(nijkout(3)))
         samevert_L = samevgrid(F_vgridin, F_ip1listin, F_vgridout, &
              F_ip1listout, samevert_list)
         if (samevert_L .and. (.not.use_same_sfcfld_L) .and. &
              associated(F_sfcfldout) .and. associated(sfcfldin)) then
            samevert_L = (.not.(any(abs(F_sfcfldout-sfcfldin) > EPSILON_R4)))
         endif

         if (samevert_L) then
            call msg(MSG_INFO,'(vinterp) No Interpolation for: '//trim(msg_S))
            if (.not.associated(F_dataout)) &
                 allocate(F_dataout(lijk(1):uijk(1),lijk(2):uijk(2),uijk(3)))
            do k=1,size(samevert_list)
               F_dataout(:,:,k) = F_datain(:,:,samevert_list(k))
            enddo
            deallocate(samevert_list, stat=istat)
            F_istat = RMN_OK
            call msg(MSG_DEBUG,'(vinterp) vinterp0ls [END] 0')
            return
         endif
         deallocate(samevert_list, stat=istat)

      endif IF_SAME_VTYPE2

      if (.not.associated(levelsin)) then
         F_istat = priv_calc_vcoor_cube(levelsin, F_vgridin, F_ip1listin, &
              sfcfldin, sfcfldin2, shape(F_datain), ispressin_L, &
              msg_S//' [inputFld]')
         if (.not.RMN_IS_OK(F_istat)) return
      endif

      if (.not.associated(levelsout)) then
         F_istat = priv_calc_vcoor_cube(levelsout, F_vgridout, F_ip1listout, &
              F_sfcfldout, F_sfcfldout2, nijkout, ispressout_L, &
              msg_S//' [outputFld]')
         if (.not.RMN_IS_OK(F_istat)) return
      endif 

      call msg(MSG_INFO,'(vinterp) Interpolating: '//trim(msg_S))

      ! Sort input levels into increasing pressures
      !#TODO: sorting should be done in vgrid from file, could be ascending or
      ! descending depending on model (ascending for GEM), only check monotonicity in here
      if (size(scol) /= size(levelsin,3)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, scol size missmach')
         return         
      endif
      scol = levelsin(1,1,:)
      do k=1,size(scol)
         slist(k,:) = minloc(scol)
         scol(slist(k,1)) = maxval(scol)+1.
      enddo

      do k=1,size(scol)
         sdatain(:,:,k) = F_datain(:,:,slist(k,1))
         slevelsin(:,:,k) = levelsin(:,:,slist(k,1))
      enddo

!!$      do k=1,size(slevelsin)
!!$         print *,'i',k,slevelsin(1,1,k)
!!$      enddo
!!$      do k=1,size(levelsout)
!!$         print *,'o',k,levelsout(1,1,k)
!!$      enddo
      
      if (.not.associated(F_dataout)) &
           allocate(F_dataout(lijk(1):uijk(1),lijk(2):uijk(2),uijk(3)))
      call vte_intvertx4(F_dataout,sdatain,slevelsin,levelsout, &
           nijkout(1)*nijkout(2),nijkin(3),nijkout(3), &
           msg_S,nlinbot)

      if (associated(levelsout)) then
         if (present(F_levelsout)) then
            if (.not.associated(levelsout, F_levelsout)) F_levelsout => levelsout
         else
            if (.not.(usealt_L .and. altout_S /= '')) deallocate(levelsout, stat=istat)
         endif
      endif
      if (associated(levelsin)) then
         if (present(F_levelsin)) then
            if (.not.associated(levelsin, F_levelsin)) F_levelsin => levelsin
         else
            if (.not.(usealt_L .and. altin_S /= '')) deallocate(levelsin, stat=istat)
         endif
      endif
      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp0ls [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp0ls


   !/@*
   function vinterp1(F_dataout, F_vgridout_S, F_datain, F_vgridin_S, &
        F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S, F_use_same_sfcfld_L, &
        F_sfcfldout2, F_sfcfldin2) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      character(len=*),intent(in) :: F_vgridout_S,F_vgridin_S
      real,pointer,optional :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer,optional :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      logical,intent(in),optional :: F_use_same_sfcfld_L
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S,tmp_S,altin_S,altout_S
      integer :: nlinbot, istat
      type(vgrid_descriptor) :: vgridout,vgridin
      real,pointer :: sfcfldout(:,:),sfcfldin(:,:)
      real,pointer :: sfcfldout2(:,:),sfcfldin2(:,:)
      integer,pointer :: ip1listout(:),ip1listin(:)
      logical :: use_same_sfcfld_L
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp1 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nlinbot = 0 
      if (present(F_nlinbot)) nlinbot = F_nlinbot
      use_same_sfcfld_L = .false.
      if (present(F_use_same_sfcfld_L)) use_same_sfcfld_L = F_use_same_sfcfld_L

      nullify(ip1listout, ip1listin, sfcfldout, sfcfldin, sfcfldout2, sfcfldin2)
      if (present(F_sfcfldout)) then
         if (associated(F_sfcfldout)) sfcfldout => F_sfcfldout
      endif
      if (present(F_sfcfldin)) then
         if (associated(F_sfcfldin)) sfcfldin => F_sfcfldin
      endif
      if (present(F_sfcfldout2)) then
         if (associated(F_sfcfldout2)) sfcfldout2 => F_sfcfldout2
      endif
      if (present(F_sfcfldin2)) then
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
      endif

      F_istat = priv_vgrid_details(F_vgridout_S, vgridout, ip1listout, &
           sfcfldout, sfcfldout2, F_alt_S=altout_S)
      if (use_same_sfcfld_L) sfcfldin => sfcfldout
      if (use_same_sfcfld_L) sfcfldin2 => sfcfldout2
      if (RMN_IS_OK(F_istat)) &
           F_istat = priv_vgrid_details(F_vgridin_S, vgridin, ip1listin, &
                                        sfcfldin, sfcfldin2, F_alt_S=altin_S)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, missing vgrid info: '//trim(msg_S))
         return
      endif

      F_istat = vinterp0ls( &
           F_dataout, vgridout, ip1listout, sfcfldout, sfcfldout2, &
           F_datain, vgridin, ip1listin, sfcfldin, sfcfldin2, nlinbot, msg_S, &
           F_altin_S=altin_S, F_altout_S=altout_S)

      istat = vgd_free(vgridout)
      if (associated(ip1listout)) deallocate(ip1listout,stat=istat)

      istat = vgd_free(vgridin)
      if (associated(ip1listin)) deallocate(ip1listin,stat=istat)

      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp1 [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp1


   !/@*
   function vinterp01(F_dataout, F_vgridout, F_ip1listout, &
        F_datain, F_vgridin_S, F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S, &
        F_sfcfldout2, F_sfcfldin2) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgridout
      integer,intent(in) :: F_ip1listout(:)
      character(len=*),intent(in) :: F_vgridin_S
      real,pointer,optional :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer,optional :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S,tmp_S,alt_S
      integer :: nlinbot, istat
      type(vgrid_descriptor) :: vgridin
      real,pointer :: sfcfldout(:,:),sfcfldin(:,:)
      real,pointer :: sfcfldout2(:,:),sfcfldin2(:,:)
      integer,pointer :: ip1listin(:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp01 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nlinbot = 0
      if (present(F_nlinbot)) nlinbot = F_nlinbot

      nullify(ip1listin, sfcfldout, sfcfldin, sfcfldout2, sfcfldin2)
      if (present(F_sfcfldout)) then
         if (associated(F_sfcfldout)) sfcfldout => F_sfcfldout
      endif
      if (present(F_sfcfldin)) then
         if (associated(F_sfcfldin)) sfcfldin => F_sfcfldin
      endif
      if (present(F_sfcfldout2)) then
         if (associated(F_sfcfldout2)) sfcfldout2 => F_sfcfldout2
      endif
      if (present(F_sfcfldin2)) then
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
      endif
      F_istat = priv_vgrid_details(F_vgridin_S, vgridin, ip1listin, &
           sfcfldin, sfcfldin2, F_alt_S=alt_S)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, vgrid info: '//trim(msg_S))
         return
      endif

      F_istat = vinterp0ls( &
           F_dataout, F_vgridout, F_ip1listout, sfcfldout, sfcfldout2, &
           F_datain, vgridin, ip1listin, sfcfldin, sfcfldin2, nlinbot, msg_S, &
           F_altin_S=alt_S)

      istat = vgd_free(vgridin)
      if (associated(ip1listin)) deallocate(ip1listin,stat=istat)

      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp01 [END] '//trim(tmp_S))
     !------------------------------------------------------------------
      return
   end function vinterp01


   !/@*
   function vinterp10(F_dataout, F_vgridout_S, F_datain, F_vgridin, &
        F_ip1listin, F_sfcfldout, F_sfcfldin, F_nlinbot, F_msg_S, &
        F_sfcfldout2, F_sfcfldin2) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real,pointer :: F_dataout(:,:,:),F_datain(:,:,:)
      character(len=*),intent(in) :: F_vgridout_S
      type(vgrid_descriptor),intent(in) :: F_vgridin
      integer,intent(in) :: F_ip1listin(:)
      real,pointer,optional :: F_sfcfldout(:,:),F_sfcfldin(:,:)
      real,pointer,optional :: F_sfcfldout2(:,:),F_sfcfldin2(:,:)
      integer,intent(in),optional :: F_nlinbot
      character(len=*),intent(in),optional :: F_msg_S
      !@return
      integer :: F_istat
      !*@/
      character(len=256) :: msg_S,tmp_S,alt_S
      integer :: nlinbot, istat
      type(vgrid_descriptor) :: vgridout
      real,pointer :: sfcfldout(:,:),sfcfldin(:,:)
      real,pointer :: sfcfldout2(:,:),sfcfldin2(:,:)
      integer,pointer :: ip1listout(:)
      !------------------------------------------------------------------
      call msg(MSG_DEBUG,'(vinterp) vinterp10 [BEGIN]')
      F_istat = RMN_ERR
      msg_S = ''
      if (present(F_msg_S)) msg_S = F_msg_S
      nlinbot = 0 
      if (present(F_nlinbot)) nlinbot = F_nlinbot

      nullify(ip1listout, sfcfldout, sfcfldin, sfcfldout2, sfcfldin2)
      if (present(F_sfcfldout)) then
         if (associated(F_sfcfldout)) sfcfldout => F_sfcfldout
      endif
      if (present(F_sfcfldin)) then
         if (associated(F_sfcfldin)) sfcfldin => F_sfcfldin
      endif
      if (present(F_sfcfldout2)) then
         if (associated(F_sfcfldout2)) sfcfldout2 => F_sfcfldout2
      endif
      if (present(F_sfcfldin2)) then
         if (associated(F_sfcfldin2)) sfcfldin2 => F_sfcfldin2
      endif

      F_istat = priv_vgrid_details(F_vgridout_S, vgridout, ip1listout, &
           sfcfldout, sfcfldout2, F_alt_S=alt_S)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, vgrid info: '//trim(msg_S))
         return
      endif

      F_istat = vinterp0ls( &
           F_dataout, vgridout, ip1listout, sfcfldout, sfcfldout2, &
           F_datain, F_vgridin, F_ip1listin, sfcfldin, sfcfldin2, &
           nlinbot, msg_S, F_altout_S=alt_S)

      istat = vgd_free(vgridout)
      if (associated(ip1listout)) deallocate(ip1listout,stat=istat)

      write(tmp_S,'(i6)') F_istat
      call msg(MSG_DEBUG,'(vinterp) vinterp10 [END] '//trim(tmp_S))
      !------------------------------------------------------------------
      return
   end function vinterp10


   !==== Private Functions =================================================


   function priv_vgrid_details(F_vgrid_S, F_vgrid, F_ip1list, &
        F_sfcfld, F_sfcfld2, F_alt_S) result(F_istat)
      implicit none
      character(len=*),intent(in) :: F_vgrid_S
      type(vgrid_descriptor),intent(out) :: F_vgrid
      integer,pointer :: F_ip1list(:)
      real,pointer :: F_sfcfld(:,:)
      real,pointer :: F_sfcfld2(:,:)
      character(len=*),intent(out),optional :: F_alt_S

      integer :: F_istat, istat, vtype
      character(len=32) :: sfcfld_S, sfcfld2_S, alt_S
      !------------------------------------------------------------------
      sfcfld2_S = ' '
      call vgrid_nullify(F_vgrid)
      F_istat = vgrid_wb_get(F_vgrid_S, F_vgrid, F_ip1list, vtype, &
           sfcfld_S, sfcfld2_S, F_altfld_S=alt_S)
      if (RMN_IS_OK(F_istat)) then
         if (sfcfld_S /= ' ' .and. .not.associated(F_sfcfld)) &
              istat = gmm_get(sfcfld_S, F_sfcfld)
         if (sfcfld2_S /= ' ' .and. .not.associated(F_sfcfld2)) &
              istat = gmm_get(sfcfld2_S, F_sfcfld2)
      endif
      if (.not.(RMN_IS_OK(F_istat) .and. &
           (associated(F_sfcfld) .or. sfcfld_S == ' '))) then
         istat = vgd_free(F_vgrid)
         if (associated(F_ip1list)) deallocate(F_ip1list, stat=istat)
         F_istat = RMN_ERR
         return
      endif
      if (sfcfld2_S /= ' ' .and. .not.associated(F_sfcfld2)) then
         istat = vgd_free(F_vgrid)
         if (associated(F_ip1list)) deallocate(F_ip1list, stat=istat)
         F_istat = RMN_ERR
      endif
      if (present(F_alt_S)) F_alt_S = alt_S
      !------------------------------------------------------------------
      return
   end function priv_vgrid_details


   function priv_calc_vcoor_cube(F_levels, F_vgrid, F_ip1list, F_sfcfld, F_sfcfld2, F_nijk, F_ispress_L, F_msg_S) result(F_istat)
      implicit none
      real, pointer :: F_levels(:,:,:)
      type(vgrid_descriptor),intent(in) :: F_vgrid
      integer,intent(in) :: F_ip1list(:)
      real, pointer :: F_sfcfld(:,:), F_sfcfld2(:,:)
      integer, intent(in) :: F_nijk(3)
      logical, intent(in) :: F_ispress_L
      character(len=*),intent(in) :: F_msg_S
      integer :: F_istat
      
      logical :: ok_L, need_rfld_L
      integer :: istat, k
      real :: sfc_min_val, sfc_max_val
      real, pointer :: plevels(:)
      real, target :: levels(size(F_ip1list))
      character(len=32) :: rfld_S, rfls_S
      !------------------------------------------------------------------
      F_istat = RMN_ERR

      rfld_S = ''
      istat = vgd_get(F_vgrid, key='RFLD', value=rfld_S)
      need_rfld_L = (rfld_S /= '')

      IF_SFCFLD: if (associated(F_sfcfld) .and. need_rfld_L) then

         ok_L = all(F_nijk(1:2) == shape(F_sfcfld))
         if (ok_L .and. associated(F_sfcfld2)) ok_L = all(F_nijk(1:2) == shape(F_sfcfld2))
         if (.not.ok_L) then
            call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, not same shape: '//trim(F_msg_S))
            print *,'(vinterp) ok_L=', ok_L
            print *,'(vinterp) shape(sfcfld) =', shape(F_sfcfld)
            if (associated(F_sfcfld2)) &
                 print *,'(vinterp) shape(sfcfld2)=', shape(F_sfcfld2)
            print *,'(vinterp) nijk(1:2)=', F_nijk(1:2)
            call flush(6)
            return
         endif

         sfc_min_val = P0MIN
         sfc_max_val = P0MAX
         if (.not.F_ispress_L) then
            sfc_min_val = H0MIN
            sfc_max_val = H0MAX
         endif
         ok_L = .true.
         if (associated(F_sfcfld2)) &
              ok_L = (minval(F_sfcfld2) >= sfc_min_val .and. maxval(F_sfcfld2) <= sfc_max_val)
         if (minval(F_sfcfld) < sfc_min_val .or. maxval(F_sfcfld) > sfc_max_val .or. .not.ok_L) then
            call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, provided sfcfld has wrong values: '//trim(F_msg_S))
            rfls_S = ''
            if (associated(F_sfcfld2)) &
               istat = vgd_get(F_vgrid, key='RFLS', value=rfls_S)
            print *,'(vinterp) rfld_S =', trim(rfld_S), '  ', trim(rfls_S), F_ispress_L
            print *,'(vinterp) ref minmax =', sfc_min_val, sfc_max_val
            print *,'(vinterp) minmax(F_sfcfld) =', minval(F_sfcfld), maxval(F_sfcfld)
            if (associated(F_sfcfld2)) &
                 print *,'(vinterp) minmax(F_sfcfld2)=', minval(F_sfcfld2), maxval(F_sfcfld2)
            call flush(6)
            return
         endif

      endif IF_SFCFLD

      if (F_nijk(3) /= size(F_ip1list)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, not same nk: '//trim(F_msg_S))
         print *,'(vinterp) nijk(3)=', F_nijk(3)
         print *,'(vinterp) size(F_ip1list)=', size(F_ip1list)
         call flush(6)
         return
      endif

      if (need_rfld_L .and. associated(F_sfcfld)) then
         if (associated(F_sfcfld2)) then
            F_istat = vgd_levels(F_vgrid, F_ip1list, F_levels, &
                 F_sfcfld, in_log=.false., sfc_field_ls=F_sfcfld2)
         else 
            F_istat = vgd_levels(F_vgrid, F_ip1list, F_levels, &
                 F_sfcfld, in_log=.false.)
         endif
      else
         plevels => levels
         F_istat = vgd_levels(F_vgrid, F_ip1list, plevels, in_log=.false.)
         if (.not.associated(F_levels)) &
              allocate(F_levels(F_nijk(1),F_nijk(2),F_nijk(3)))
         do k=1,size(levels)
            F_levels(:,:,k) = levels(k)
         enddo
      endif
!!$      print *,'(vinterp) F_istat=', F_istat
!!$      print *,'(vinterp) size(F_ip1list)=', size(F_ip1list)
!!$      print *,'(vinterp) shape(F_levels)=', shape(F_levels)
!!$      do k=1,size(F_levels,3)
!!$         print *,'(vinterp) minmax(F_levels)=', k, minval(F_levels(:,:,k)),maxval(F_levels(:,:,k))
!!$      enddo

      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING,'(vinterp) Cannot Interpolate, problem getting levels: '//trim(F_msg_S))
         return
      endif
      !------------------------------------------------------------------
      return
   end function priv_calc_vcoor_cube
   
end module vinterp_mod
