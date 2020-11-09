module my_sedi_mod
   use, intrinsic :: iso_fortran_env, only: REAL64

!================================================================================!
!  The following subroutines are used by the schemes in the multimoment package. !
!                                                                                !
!  Package version:  2.18.0     (internal bookkeeping)                           !
!  Last modified  :  2011-01-07                                                  !
!================================================================================!

   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

  private

  public :: SEDI_main_1,SEDI_main_1b,SEDI_main_2,SEDI_ISGH_V33,countColumns_v33, &
            blg4sedi,blg5sedi,countColumns_a,countColumns

   contains

!=====================================================================================!
 subroutine SEDI_main_2(QX,NX,cat,Q,T,DE,iDE,gamfact,epsQ,epsN,afx,bfx,cmx,dmx,      &
                        ckQx1,ckQx2,ckQx4,LXP,ni,nk,VxMax,DxMax,dt,DZ,massFlux,      &
                        ktop_sedi,GRAV,massFlux3D)

!-------------------------------------------------------------------------------------!
!  DOUBLE-MOMENT version of sedimentation subroutine for categories whose
!  fall velocity equation is V(D) = gamfact * afx * D^bfx
!-------------------------------------------------------------------------------------!

! Passing parameters:
!
!  VAR   Description
!  ---   ------------
!  QX    mass mixing ratio of category x
!  NX    number concentration of category x
!  cat:  hydrometeor category:
!   1     rain
!   2     ice
!   3     snow
!   4     graupel
!   5     hail
!-------------------------------------------------------------------------------------!

  use my_fncs_mod

  implicit none

! PASSING PARAMETERS:
  real, dimension(:,:), intent(inout) :: QX,NX,Q,T
  real, dimension(:),    intent(out)  :: massFlux
  real, optional, dimension(:,:), intent(out) :: massFlux3D
  real, dimension(:,:), intent(in)    :: DE,iDE,DZ
  real, intent(in) :: epsQ,epsN,VxMax,LXP,afx,bfx,cmx,dmx,ckQx1,ckQx2,ckQx4,DxMax,dt,GRAV
  integer, dimension(:), intent(in)   :: ktop_sedi
  integer, intent(in)                 :: ni,nk,cat

! LOCAL PARAMETERS:
  logical                :: slabHASmass,locallim,QxPresent
  integer                :: nnn,a,i,k,counter,l,km1,kp1,ks,kw,idzmin
  integer, dimension(nk) :: flim_Q,flim_N
  integer, dimension(ni) :: activeColumn,npassx,ke
  real                   :: VqMax,VnMax,iLAMx,iLAMxB0,tmp1,tmp2,tmp3,Dx,iDxMax,icmx,     &
                            VincFact,ratio_Vn2Vq,zmax_Q,zmax_N,tempo,idmx,Nos_Thompson,  &
                            No_s,iLAMs
  real, dimension(ni,nk) :: VVQ,VVN,RHOQX,gamfact
  real, dimension(ni)    :: dzMIN,dtx,VxMaxx
  real, dimension(nk)    :: vp_Q,vp_N,zt_Q,zt_N,zb_Q,zb_N,dzi,Q_star,N_star
  real, dimension(0:nk)  :: zz
  real, parameter        :: epsilon = 1.e-2
  real, parameter        :: thrd    = 1./3.
  real, parameter        :: sxth    = 1./6.
  real, parameter        :: CoMAX   = 2.0

!-------------------------------------------------------------------------------------!

   massFlux = 0.

  !Factor to estimate increased V from size-sorting:
  ! - this factor should be higher for categories with more time-splitting, since Vmax
  !   increases after each sedimentation split step (to be tuned)
   VincFact = 1.
   if (present(massFlux3D)) massFlux3D= 0.  !(for use in MAIN for diagnostics)

  !Determine for which slabs and columns sedimentation should be computes:
   call countColumns(QX,ni,nk,epsQ,counter,activeColumn,ktop_sedi)

   ratio_Vn2Vq= ckQx2/ckQx1
   iDxMax= 1./DxMax
   icmx  = 1./cmx
   idmx  = 1./dmx
   ks    = nk
   ke    = ktop_sedi  !(i-array) - formerly ke=1; now depends on max. level with hydrometeor
   kw    = -1         !direction of vertical leveling; -1 implies nk is bottom

   VVQ  = 0.
   VVN  = 0.
   VqMax= 0.
   VnMax= 0.

   do a= 1,counter
      i= activeColumn(a)

      VVQ(i,:) = 0.
      do k= ktop_sedi(i),nk  !formerly do k= 1,nk
         QxPresent =  (QX(i,k)>epsQ .and. NX(i,k)>epsN)
         if (QxPresent) VVQ(i,k)= calcVV()*ckQx1
         if (present(massFlux3D)) massFlux3D(i,k)= VVQ(i,k)*DE(i,k)*QX(i,k)  !(for use in MAIN)
      enddo  !k-loop
      Vxmaxx(i)= min( VxMax, maxval(VVQ(i,:))*VincFact )

     !note: dzMIN is min. value in column (not necessarily lowest layer in general)
      dzMIN(i) = minval(DZ(i,:))
      npassx(i)= max(1, nint( dt*Vxmaxx(i)/(CoMAX*dzMIN(i)) ))
      dtx(i)   = dt/float(npassx(i))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      do nnn= 1,npassx(i)

         locallim = (nnn==1)

         do k= ktop_sedi(i),nk  !formerly do k= 1,nk
           RHOQX(i,k) = DE(i,k)*QX(i,k)
           QxPresent  = (QX(i,k)>epsQ .and. NX(i,k)>epsN)
           if (QxPresent) then
              if (locallim) then     !to avoid re-computing VVQ on first pass
                 VVQ(i,k)= -VVQ(i,k)
              else
                 VVQ(i,k)= -calcVV()*ckQx1
              endif
              VVN(i,k)= VVQ(i,k)*ratio_Vn2Vq
              VqMax   = max(VxMAX,-VVQ(i,k))
              VnMax   = max(VxMAX,-VVN(i,k))
           else
              VVQ(i,k)= 0.
              VVN(i,k)= 0.
              VqMax   = 0.
              VnMax   = 0.
           endif
         enddo  !k-loop

        !sum instantaneous surface mass flux at each split step: (for division later)
         massFlux(i)= massFlux(i) - VVQ(i,nk)*DE(i,nk)*QX(i,nk)

     !-- Perform single split sedimentation step:
     !   (formerly by calls to s/r 'blg4sedi', a modified [JM] version of 'blg2.ftn')
         zz(ks)= 0.
         do k= ks,ke(i),kw
            zz(k+kw)= zz(k)+dz(i,k)
            dzi(k)  = 1./dz(i,k)
            vp_Q(k) = 0.
            vp_N(k) = 0.
         enddo

         do k=ks,ke(i),kw
            zb_Q(k)= zz(k) + VVQ(i,k)*dtx(i)
            zb_N(k)= zz(k) + VVN(i,k)*dtx(i)
         enddo

         zt_Q(ke(i))= zb_Q(ke(i)) + dz(i,ke(i))
         zt_N(ke(i))= zb_N(ke(i)) + dz(i,ke(i))
         do k= ks,ke(i)-kw,kw
            zb_Q(k)= min(zb_Q(k+kw)-epsilon*dz(i,k), zz(k)+VVQ(i,k)*dtx(i))
            zb_N(k)= min(zb_N(k+kw)-epsilon*dz(i,k), zz(k)+VVN(i,k)*dtx(i))
            zt_Q(k)= zb_Q(k+kw)
            zt_N(k)= zb_N(k+kw)
         enddo

         do k=ks,ke(i),kw    !formerly k=1,nk
            Q_star(k)= RHOQX(i,k)*dz(i,k)/(zt_Q(k)-zb_Q(k))
            N_star(k)=    NX(i,k)*dz(i,k)/(zt_N(k)-zb_N(k))
         enddo

         if (locallim) then
            zmax_Q= abs(VqMax*dtx(i))
            zmax_N= abs(VnMax*dtx(i))
            do l=ks,ke(i),kw
               flim_Q(l)= l
               flim_N(l)= l
               do k= l,ke(i),kw
                  if (zmax_Q.ge.zz(k)-zz(l+kw)) flim_Q(l)= k
                  if (zmax_N.ge.zz(k)-zz(l+kw)) flim_N(l)= k
               enddo
            enddo
         endif

         do l=ks,ke(i),kw
            do k=l,flim_Q(l),kw
               vp_Q(l)= vp_Q(l) + Q_star(k)*max(0.,min(zz(l+kw),zt_Q(k))-max(zz(l),zb_Q(k)))
            enddo
            do k=l,flim_N(l),kw
               vp_N(l)= vp_N(l) + N_star(k)*max(0.,min(zz(l+kw),zt_N(k))-max(zz(l),zb_N(k)))
            enddo
         enddo

         do k=ks,ke(i),kw
            RHOQX(i,k)= vp_Q(k)*dzi(k)
               NX(i,k)= vp_N(k)*dzi(k)
         enddo
     !--

         do k= ktop_sedi(i),nk  !formerly do k= 1,nk
           QX(i,k)= RHOQX(i,k)*iDE(i,k)

         !Prevent levels with zero N and nonzero Q and size-limiter:
           QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN)
           if (QxPresent) then    !size limiter
              Dx= (DE(i,k)*QX(i,k)/(NX(i,k)*cmx))**idmx
              if (cat==1 .and. Dx>3.e-3) then
                 tmp1   =  Dx-3.e-3;   tmp1= tmp1*tmp1
                 tmp2   = (Dx/DxMAX);  tmp2= tmp2*tmp2*tmp2
                 NX(i,k)= NX(i,k)*max((1.+2.e4*tmp1),tmp2)
              else
                 NX(i,k)= NX(i,k)*(max(Dx,DxMAX)*iDxMAX)**dmx   !impose Dx_max
              endif
           else   !here, "QxPresent" implies correlated QX and NX
              Q(i,k) = Q(i,k) + QX(i,k)
              T(i,k) = T(i,k) - LXP*QX(i,k)   !LCP for rain; LSP for i,s,g,h
              QX(i,k)= 0.
              NX(i,k)= 0.
           endif

         enddo

       enddo  !nnn-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !compute average mass flux during the full time step: (used to compute the
      !instantaneous sedimentation rate [liq. equiv. volume flux] in the main s/r)
       massFlux(i)= massFlux(i)/float(npassx(i))

    enddo  !a(i)-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

contains

   real function calcVV()
   !Calculates portion of moment-weighted fall velocities
      iLAMx   = ((QX(i,k)*DE(i,k)/NX(i,k))*ckQx4)**idmx
      iLAMxB0 = iLAMx**bfx
      calcVV  = gamfact(i,k)*iLAMxB0
   end function calcVV

 end subroutine SEDI_main_2

!=====================================================================================!
 subroutine SEDI_main_1(QX,cat,T,DE,gamfact,epsQ,afx,bfx,icmx,dmx,dtx,cx6,ckQx1,    &
                      ckQx2,ckQx4,npassx,ni,nk,VxMax,DxMax,DZ,PR,No_x,ktop_sedi,    &
                      massFlux3D)

!                  ** Used by 'my_smom' (from 2.7.1) only **                          !

!-------------------------------------------------------------------------------------!
!  SINGLE-MOMENT version of sedimentation subroutine for categories whose
!  fall velocity equation is V(D) = gamfact * afx * D^bfx
!
!     ** ASSUMES INVERSE-EXPONENTIAL DISTRIBTIONS (alpha_x=0) **
!-------------------------------------------------------------------------------------!

  use my_fncs_mod
  implicit none

! PASSING ARGUMENTS:
  real, dimension(ni,nk), intent(inout) :: QX,T
  real, optional, dimension(:,:), intent(out) :: massFlux3D
  real, dimension(ni),    intent(inout) :: PR
  real, dimension(ni,nk), intent(in)    :: DE,DZ
  real, intent(in)       :: dtx,epsQ,cx6,VxMax,afx,bfx,icmx,dmx,ckQx1,ckQx2,ckQx4,  &
                            DxMax,No_x
  integer, intent(in)    :: npassx,ni,nk,cat,ktop_sedi


! LOCAL ARGUMENTS:
  logical                :: slabHASmass,LOCALLIM,QxPresent
  real, dimension(ni,nk) :: VVQ,VVN,VVZ,RHOQX,gamfact
  integer, dimension(nk) :: FLIM
  integer, dimension(ni) :: activeColumn
  real                   :: VqMax,iLAMx,tmp1,tmp2,Dx,iDxMax,No_s,NX,iNo_x
  integer                :: nnn,a,i,k,counter
  real, parameter        :: thrd    = 1./3.
  real, parameter        :: sxth    = 1./6.

!-------------------------------------------------------------------------------------!

   if (present(massFlux3D)) massFlux3D= 0.  !(for use in MAIN for diagnostics)

  !Determine for which slabs and columns sedimentation should be computes:
   call countColumns_a(QX,ni,nk,epsQ,counter,activeColumn,ktop_sedi)

    iNo_x = 1./No_x
    iDxMax= 1./DxMax
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    do nnn= 1,npassx

       RHOQX= DE*QX
       VVQ  = 0.;  VqMax= 0.
       do a= 1,counter
         i= activeColumn(a)
         do k= ktop_sedi,nk   !do k=1,nk

           QxPresent=  (QX(i,k)>epsQ)
           if (QxPresent) then
            !ice:
              if (cat==2) then
                 NX    = 5.*exp(0.304*(273.15-max(233.,T(i,k))))
                 iLAMx = (ckQx4*QX(i,k)*DE(i,k)/NX)**thrd
            !snow:
              else if (cat==3) then
                 No_s  = min(2.e+8, 2.e+6*exp(-0.12*min(-0.001,T(i,k)-273.15))) !T2004
                 iNo_x = 1./No_s
                 iLAMx = sqrt(sqrt(QX(i,k)*DE(i,k)*icmx*sxth*iNo_x))
            !rain, graupel, hail:
              else
                 iLAMx = sqrt(sqrt(QX(i,k)*DE(i,k)*icmx*sxth*iNo_x))
              endif
              VVQ(i,k) = -gamfact(i,k)*ckQx1*iLAMx**bfx
              VqMax    = max(VxMAX,-VVQ(i,k))

              if (present(massFlux3D)) massFlux3D(i,k)= -VVQ(i,k)*DE(i,k)*QX(i,k)  !(for use in MAIN)

           endif

         enddo  !k-loop
       enddo    !i(a)-loop

       locallim= (nnn==1)
       call blg4sedi(RHOQX,DZ,VVQ,nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)
       QX= RHOQX/DE

       PR(:)= PR(:) - cx6*VVQ(:,nk)*DE(:,nk)*QX(:,nk)

    enddo  !nnn-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

 end subroutine SEDI_main_1

!=====================================================================================!
 subroutine SEDI_main_1b(QX,cat,T,DE,iDE,gamfact,epsQ,afx,bfx,icmx,dmx,ckQx1,ckQx4, &
                         ni,nk,VxMax,DxMax,dt,DZ,massFlux,No_x,ktop_sedi,GRAV,      &
                         massFlux3D)

!-------------------------------------------------------------------------------------!
!  SINGLE-MOMENT version of sedimentation subroutine for categories whose
!  fall velocity equation is V(D) = gamfact * afx * D^bfx
!-------------------------------------------------------------------------------------!

! Passing parameters:
!
!  VAR   Description
!  ---   ------------
!  QX    mass mixing ratio of category x
!  cat:  hydrometeor category:
!   1     rain
!   2     ice
!   3     snow
!   4     graupel
!   5     hail
!-------------------------------------------------------------------------------------!

  use my_fncs_mod

  implicit none

! PASSING PARAMETERS:
  real, dimension(:,:), intent(inout) :: QX,T
  real, dimension(:),    intent(out)   :: massFlux
  real, optional, dimension(:,:), intent(out) :: massFlux3D
  real, dimension(:,:), intent(in)    :: DE,iDE,DZ
  real,    intent(in)    :: epsQ,VxMax,afx,bfx,icmx,dmx,ckQx1,ckQx4,DxMax,dt,GRAV,No_x
  integer, dimension(:), intent(in) :: ktop_sedi
  integer, intent(in)    :: ni,nk,cat !,ktop_sedi

! LOCAL PARAMETERS:
  logical                :: slabHASmass,locallim,QxPresent
  integer                :: nnn,a,i,k,counter,l,km1,kp1,ks,kw,idzmin !,ke
  integer, dimension(nk) :: flim_Q
  integer, dimension(ni) :: activeColumn,npassx,ke
  real                   :: VqMax,iLAMx,iLAMxB0,tmp1,tmp2,Dx,iDxMax,VincFact,NX,iNo_x,   &
                            zmax_Q,zmax_N,tempo
  real, dimension(ni,nk) :: VVQ,RHOQX,gamfact
  real, dimension(ni)    :: dzMIN,dtx,VxMaxx
  real, dimension(nk)    :: vp_Q,zt_Q,zb_Q,dzi,Q_star
  real, dimension(0:nk)  :: zz
  real, parameter        :: epsilon = 1.e-2
  real, parameter        :: thrd  = 1./3.
  real, parameter        :: sxth  = 1./6.
  real, parameter        :: CoMAX = 2.0

!-------------------------------------------------------------------------------------!

   massFlux= 0.
  !Factor to estimate increased V from size-sorting:
  ! - this factor should be higher for categories with more time-splitting, since Vmax
  !   increases after each sedimentation split step (to be tuned)
   VincFact= 1.
   if (present(massFlux3D)) massFlux3D= 0.  !(for use in MAIN for diagnostics)

  !Determine for which slabs and columns sedimentation should be computes:
   call countColumns(QX,ni,nk,epsQ,counter,activeColumn,ktop_sedi)
   iNo_x = 1./No_x
   iDxMax= 1./DxMax
   ks    = nk
   ke    = ktop_sedi  !(i-array) - old: ke=1
   kw    = -1         !direction of vertical leveling

   VVQ  = 0.
   VqMax= 0.

   do a= 1,counter
      i= activeColumn(a)

      VVQ(i,:) = 0.
      do k= ktop_sedi(i),nk  !do k= 1,nk
         QxPresent =  (QX(i,k)>epsQ)
!        if (QxPresent) VVQ(i,k)= calcVV()*ckQx1

         if (QxPresent) then
            !ice:
              if (cat==2) then
                 NX    = 5.*exp(0.304*(273.15-max(233.,T(i,k))))
                 iLAMx = (ckQx4*QX(i,k)*DE(i,k)/NX)**thrd
            !snow:
              else if (cat==3) then
                 iNo_x = 1./min(2.e+8, 2.e+6*exp(-0.12*min(-0.001,T(i,k)-273.15)))
                 iLAMx = sqrt(sqrt(QX(i,k)*DE(i,k)*icmx*sxth*iNo_x))
            !rain, graupel, hail:
              else
                 iLAMx = sqrt(sqrt(QX(i,k)*DE(i,k)*icmx*sxth*iNo_x))
              endif
              VVQ(i,k) = -gamfact(i,k)*ckQx1*iLAMx**bfx
         !    VqMax    = max(VxMAX,-VVQ(i,k))
         endif
         if (present(massFlux3D)) massFlux3D(i,k)= -VVQ(i,k)*DE(i,k)*QX(i,k)  !(for use in MAIN)

      enddo  !k-loop
      Vxmaxx(i)= min( VxMax, maxval(VVQ(i,:))*VincFact )

     !note: dzMIN is min. value in column (not necessarily lowest layer in general)
      dzMIN(i) = minval(DZ(i,:))
      npassx(i)= max(1, nint( dt*Vxmaxx(i)/(CoMAX*dzMIN(i)) ))
      dtx(i)   = dt/float(npassx(i))

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      do nnn= 1,npassx(i)

         locallim = (nnn==1)

         do k= ktop_sedi(i),nk  !do k= 1,nk
           RHOQX(i,k) = DE(i,k)*QX(i,k)
           QxPresent  = (QX(i,k)>epsQ)
            if (QxPresent) then
             !ice:
               if (cat==2) then
                  NX    = 5.*exp(0.304*(273.15-max(233.,T(i,k))))
                  iLAMx = (ckQx4*QX(i,k)*DE(i,k)/NX)**thrd
             !snow:
               else if (cat==3) then
                  iNo_x = 1./min(2.e+8, 2.e+6*exp(-0.12*min(-0.001,T(i,k)-273.15)))
                  iLAMx = sqrt(sqrt(QX(i,k)*DE(i,k)*icmx*sxth*iNo_x))
             !rain, graupel, hail:
               else
                  iLAMx = sqrt(sqrt(QX(i,k)*DE(i,k)*icmx*sxth*iNo_x))
               endif
               VVQ(i,k) = -gamfact(i,k)*ckQx1*iLAMx**bfx
               VqMax    = max(VxMAX,-VVQ(i,k))
            endif

         enddo  !k-loop

     !-- Perform single split sedimentation step:  (formerly by calls to s/r 'blg4sedi')
         zz(ks)= 0.
         do k= ks,ke(i),kw
            zz(k+kw)= zz(k)+dz(i,k)
            dzi(k)  = 1./dz(i,k)
            vp_Q(k) = 0.
         enddo

         do k=ks,ke(i),kw
            zb_Q(k)= zz(k) + VVQ(i,k)*dtx(i)
         enddo

         zt_Q(ke(i))= zb_Q(ke(i)) + dz(i,ke(i))
         do k= ks,ke(i)-kw,kw
            zb_Q(k)= min(zb_Q(k+kw)-epsilon*dz(i,k), zz(k)+VVQ(i,k)*dtx(i))
            zt_Q(k)= zb_Q(k+kw)
         enddo

         do k=ks,ke(i),kw    !k=1,nk
            Q_star(k)= RHOQX(i,k)*dz(i,k)/(zt_Q(k)-zb_Q(k))
         enddo

         if (locallim) then
            zmax_Q= abs(VqMax*dtx(i))
            do l=ks,ke(i),kw
               flim_Q(l)= l
               do k= l,ke(i),kw
                  if (zmax_Q.ge.zz(k)-zz(l+kw)) flim_Q(l)= k
               enddo
            enddo
         endif

         do l=ks,ke(i),kw
            do k=l,flim_Q(l),kw
               vp_Q(l)= vp_Q(l) + Q_star(k)*max(0.,min(zz(l+kw),zt_Q(k))-max(zz(l),zb_Q(k)))
            enddo
         enddo

         do k=ks,ke(i),kw
            RHOQX(i,k)= vp_Q(k)*dzi(k)
         enddo
     !--

         do k= ktop_sedi(i),nk  ! do k= 1,nk
           QX(i,k)= RHOQX(i,k)*iDE(i,k)
         enddo

        !sum instantaneous flux at each split step: (for division later)
         massFlux(i)= massFlux(i) - VVQ(i,nk)*DE(i,nk)*QX(i,nk)

       enddo  !nnn-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
      !compute average flux during the full time step: (this will be used to compute
      ! the instantaneous sedimentation rate [volume flux] in the main s/r)
       massFlux(i)= massFlux(i)/float(npassx(i))

    enddo  !a(i)-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

 end subroutine SEDI_main_1b

!=====================================================================================!
 subroutine countColumns_a(QX,ni,nk,minQX,counter,activeColumn,ktop_sedi)

!    ***  Used by SEDI_main_1, which is ued my my_smom (2.5.1) only ***      !

! Searches the hydrometeor array QX(ni,nk) for non-zero (>minQX) values.
! Returns the array if i-indices (activeColumn) for the columns (i)
! which contain at least one non-zero value, as well as the number of such
! columns (counter).

 implicit none

!PASSING PARAMETERS:
  integer, intent(in)                   :: ni,nk,ktop_sedi
  integer, intent(out)                  :: counter
  integer, dimension(ni), intent(out)   :: activeColumn
  real,    dimension(ni,nk), intent(in) :: QX
  real,    intent(in)                   :: minQX

!LOCAL PARAMETERS:
  integer                               :: i,k

   k=ktop_sedi-1; counter=0; activeColumn=0
!  k=0; counter=0; activeColumn=0
   do i=1,ni

      do
         k=k+1
         if (QX(i,k)>minQX) then
            counter=counter+1
            activeColumn(counter)=i
!             slabHASmass=.true.     <-- results in i-dependency
            k=0
            exit
         else
            if (k==nk) then
               k=0
               exit
            endif
         endif
      enddo

   enddo  !i-loop

 end subroutine countColumns_a

!=====================================================================================!
 subroutine countColumns(QX,ni,nk,minQX,counter,activeColumn,ktop_sedi)

! Searches the hydrometeor array QX(ni,nk) for non-zero (>minQX) values.
! Returns the array if i-indices (activeColumn) for the columns (i)
! which contain at least one non-zero value, as well as the number of such
! columns (counter).

  implicit none

!PASSING PARAMETERS:
  integer, intent(in)                   :: ni,nk !,ktop_sedi
  integer, dimension(:), intent(in)     :: ktop_sedi
  integer,                 intent(out)  :: counter
  integer, dimension(:), intent(out)    :: activeColumn
  real,    dimension(:,:), intent(in)   :: QX
  real,    intent(in)                   :: minQX

!LOCAL PARAMETERS:
  integer                               :: i !,k
  integer, dimension(ni)                :: k

!    k= ktop_sedi-1  !  k=0
   counter     = 0
   activeColumn= 0

   do i=1,ni
      k(i)= ktop_sedi(i)-1  !  k=0
      do
         k(i)=k(i)+1
         if (QX(i,k(i))>minQX) then
            counter=counter+1
            activeColumn(counter)=i
            k(i)=0
            exit
         else
            if (k(i)==nk) then
               k(i)=0
               exit
            endif
         endif
      enddo
   enddo  !i-loop

 end subroutine countColumns

!=====================================================================================!

! The following subroutines are used only by 'my_tmom_mod.cdk90'.  They are somewhat
! redundant from above routines, though there are small differences.  Eventually,
! all versions will use the same routines.
! 2008-04-15

!=====================================================================================!
 subroutine SEDI_ISGH_V33(QX,NX,ZX,cat,Q,T,DE,gamfact,epsQ,epsN,epsZ,afx,bfx,cmx,dmx,    &
                          dtx,cx6,ALFxfix,Noxfix,LXP,npassx,ni,nk,VxMax,DxMax,DZ,SR,     &
                          scheme,ktop_sedi)

!-------------------------------------------------------------------------------------!
! Sedimentation subroutine for categories whose fall velocity equation is
! V(D) = gamfact * afx * D^bfx
!
!  ***  for my_main_full.ftn90  ***
!-------------------------------------------------------------------------------------!

  use my_fncs_mod

  implicit none

! PASSING PARAMETERS:
  real, dimension(ni,nk), intent(inout) :: QX,NX,ZX,Q,T
  real, dimension(ni),    intent(inout) :: SR
  real, dimension(ni,nk), intent(in)    :: DE,DZ
  real,    intent(in)    :: dtx,epsQ,epsN,epsZ,cx6,VxMax,LXP
  real(REAL64),  intent(in)    :: afx,bfx,cmx,dmx,ALFxfix,Noxfix,DxMax
  integer, intent(in)    :: npassx,ni,nk,scheme,cat !,ktop_sedi
  integer, dimension(ni), intent(in) :: ktop_sedi

! LOCAL PARAMETERS:
  real, dimension(ni,nk) :: VVQ,VVN,VVZ,RHOQX,gamfact
  integer, dimension(nk) :: FLIM
  logical                :: slabHASmass,LOCALLIM,QxPresent
  integer, dimension(ni) :: activeColumn
  real                   :: VqMax,VnMax,Vzmax,cmxSP
  real(REAL64)                 :: ALFx,GX2,GX5,ckQx1,ckQx2,ckQx3,iLAMx,iLAMxB0,tmpdp1,tmpdp2,Dx
  integer                :: nnn,a,i,k,counter
  real(REAL64), parameter      :: thrd  = 1.d0/3.d0

!-------------------------------------------------------------------------------------!

  cmxSP= sngl(cmx)
!
  !Determine for which slabs and columns sedimentation should be computes:
   call countColumns_v33(QX,ni,nk,epsQ,counter,activeColumn,slabHASmass,ktop_sedi)

   if (slabHASmass) then

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    do nnn= 1,npassx

       RHOQX= DE*QX
       VVQ= 0.;  VVN= 0.;  VVZ= 0.;  VqMax= 0.;  VnMax= 0.;  VzMax= 0.
       do a= 1,counter
         i=activeColumn(a)
!        do k= 1,nk
         do k= ktop_sedi(i),nk

           if (scheme==1)                    then
              QxPresent = (QX(i,k)>epsQ)
              ALFx      = ALFxfix
              if (cat==2) then  ![ice]
                NX(i,k) = 5.*exp(0.304*(273.15-max(233.,T(i,k))))  !Cooper eqn.
              else if (cat>=3.and.cat<=5) then  ![snow, grpl, or hail]
                tmpdp1  = gammaDP(1.d0+ALFx)
                tmpdp2  = gammaDP(4.d0+ALFx)
                NX(i,k) = (Noxfix*tmpdp1)**(3./(4.+ALFx))*(tmpdp1/tmpdp2*DE(i,k)*     &
                           QX(i,k)/cmx)**((1.+ALFx)/(4.+ALFx))
              endif
           else if (scheme==2.or.scheme==3)  then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN)
              if (QxPresent) then
                 Dx  = (dble(DE(i,k)*QX(i,k)/NX(i,k))/cmx)**thrd
                 ALFx= diagAlpha_v33(Dx,cat)
              endif
           else if (scheme==4)               then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN .and. ZX(i,k)>epsZ)
              if (QxPresent) ALFx= solveAlpha_v33(QX(i,k),NX(i,k),ZX(i,k),cmxSP,DE(i,k))
           endif

           if (QxPresent) then
              GX2      = 1.d0/gammaDP(1.d0+ALFx)
              GX5      = gammaDP(1.d0+ALFx+dmx)
              ckQx1    = afx*gammaDP(1.d0+ALFx+dmx+bfx)/GX5
              ckQx2    = afx*gammaDP(1.d0+ALFx+bfx)*GX2
              ckQx3    = afx*gammaDP(7.d0+ALFx+bfx)/gammaDP(7.d0+ALFx)
              iLAMx    = (dble(QX(i,k)*DE(i,k)/NX(i,k))/(cmx*GX5*GX2))**thrd
              iLAMxB0  = iLAMx**bfx
              tmpdp1   = -gamfact(i,k)*iLAMxB0
              VVQ(i,k) = tmpdp1*ckQx1;   VqMax= max(VxMAX,-VVQ(i,k))
              VVN(i,k) = tmpdp1*ckQx2;   VnMax= max(VxMAX,-VVN(i,k))
              VVZ(i,k) = tmpdp1*ckQx3;   VzMax= max(VxMAX,-VVZ(i,k))
           endif

         enddo  !k-loop
       enddo    !i-loop
       locallim= (nnn==1)
       call blg5sedi(RHOQX,DZ,VVQ,   nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)
       if (scheme >1)  &
          call blg5sedi(NX,DZ,VVN,   nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)
      if (scheme==4)  &
          call blg5sedi(ZX,DZ,VVZ,   nk,dtx,locallim,VqMax,FLIM,counter,activeColumn,ktop_sedi)
      QX= RHOQX/DE

    ! Prevent levels with zero N and nonzero Q and size-limiter:
       if (scheme>1) then
       do a= 1,counter
         i=activeColumn(a)
!        do k= 1,nk
         do k= ktop_sedi(i),nk
           if (scheme==2.or.scheme==3)  then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN)
           elseif (scheme==4)           then
              QxPresent=  (QX(i,k)>epsQ .and. NX(i,k)>epsN .and. ZX(i,k)>epsZ)
           endif
           if (.not. QxPresent) then
              Q(i,k) = Q(i,k)+QX(i,k)
              T(i,k) = T(i,k) - LXP*QX(i,k)   !LCP for rain; LSP for i,s,g,h
              QX(i,k)= 0.;  NX(i,k)= 0.;  ZX(i,k)= 0.
           else  ! size limiter:
              Dx     = (dble(DE(i,k)*QX(i,k)/NX(i,k))/cmx)**thrd
              tmpdp1 = sngl(max(Dx,DxMAX)/DxMAX)
              NX(i,k)= NX(i,k)*tmpdp1*tmpdp1*tmpdp1
           endif
         enddo
       enddo
       endif !(if scheme>1)

       SR(:)= SR(:) - cx6*VVQ(:,nk)*DE(:,nk)*QX(:,nk)

    enddo  !nnn-loop
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!

   endif  !slabHASmass

 end subroutine SEDI_ISGH_V33

!=====================================================================================!

 subroutine countColumns_v33(QX,ni,nk,minQX,counter,activeColumn,slabHASmass,ktop_sedi)

! Searches the hydrometeor array QX(ni,nk) for non-zero (>minQX) values.
! Returns slabHASmass=TRUE if there is a single non-zero value in the array
! and returns the array if i-indices (activeColumn) for the columns (i)
! which contain at least one non-zero value, as well as the number of such
! columns (counter).

  implicit none

 ! PASSING PARAMETERS:
  integer, intent(in)                   :: ni,nk !,ktop_sedi
  integer, intent(out)                  :: counter
  integer, dimension(ni), intent(in)    :: ktop_sedi
  integer, dimension(ni), intent(out)   :: activeColumn
  real,    dimension(ni,nk), intent(in) :: QX
  real,    intent(in)                   :: minQX
  logical, intent(out)                  :: slabHASmass

! LOCAL PARAMETERS:
   integer                              :: i !,k
   integer, dimension(ni)               :: k

!  k=0; counter=0; activeColumn=0; slabHASmass=.false.
   k=ktop_sedi-1; counter=0; activeColumn=0; slabHASmass=.false.
   do i=1,ni

      do
         k(i)=k(i)+1
         if (QX(i,k(i))>minQX) then
            counter=counter+1
            activeColumn(counter)=i
            slabHASmass=.true.
            k(i)=0
            exit
         else
            if (k(i)==nk) then
               k(i)=0
               exit
            endif
         endif
      enddo

   enddo  !i-loop


 end subroutine countColumns_v33

!=====================================================================================!
!This subroutine is modified from S/P BLG.FTN

 subroutine blg4sedi (RO,DZ,WW,NK,DT,COMPLIM,WMAX,FLIM,counter,activeColumn,ktop_sedi)

 implicit none

!PASSING ARGUMENTS:

 integer, intent(in)                             :: NK,counter
 real, dimension(:,:), intent(inout)             :: ro
 real, dimension(:,:), intent(in)                :: dz,ww
 real, intent(in)                                :: dt,wmax
 integer, dimension(nk)                          :: FLIM
 integer, dimension(size(RO,dim=1)), intent(in)  :: activeColumn
 logical, intent(inout)                          :: COMPLIM
 integer, intent(in)                             :: ktop_sedi

! Author
!          C. Girard and A. Plante (2001)
!
! Revisions
!
! 001      A. Glazer and A. Plante (Oct 2001)
!             - introduced complim and ind_limit in computation of precipitation
! 002      A. Plante (June 2003)
!             - IBM conversion, added ro_star. This imporved sedimentaion
!               in mixphase4 by 38%.
!             - change computation limit to insure reentrance in OpenMP
!             - blg2.ftn validates perfectly with blg.ftn of PHY 3.8.
! 003      J. Milbrandt (Nov 2004)
!             - Modified for use in Milbrandt-Yau microphysics scheme
!               + add condition for mass-in-column (activeColumn); for all i-loops, the
!                 line 'do i=1,ni' was replaced with 'do a=1,counter' & 'i=activeColumn(a)'
!               + remove RT calculation (and pass)
!               + hard-wired various options for multimoment.ftn90 (removed unnecessary IF statements)
! 004      J. Milbrandt (Jan 2007)
!             - corrected 'idzmin' initial value (for use in activeColumn i-loops)
!             - added option to exclude upper levels (ktop_sedi)
! 005      J. Milbrandt and R. McTaggart-Cowan (March 2007)
!             - removed i-dependency; cleaned up arrays, etc.  RO, DZ, WW are (ni,nk) arrays but there is
!               only one outer i (a)-loop (which could be moved outside of this subroutine, leaving it
!               as a single column subroutine).
!
! Object
!
!
!   Version 1.0
!
!   CALCULATES
!                sedimentation of a quantity of matter ro(kg/m3)
!                falling with negative downward velocity ww(m/s)
!
!   ACCORDING TO
!                the BOX-LAGRANGIAN SCHEME
!                         (Kato,1995:JMSJ,73,241-245)
!             OR
!                the ADJACENT-LAGRANGIAN SCHEME
!                Girard and Plante
!
!   PLUS
!                a conservative two-grid-length filter
!                in order to control noise if needed
!
!Arguments
!
!          - Input/Output -
!
! RO       density in cell before and after sedimentation.
! COMPLIM  logical switch to force recomputing the table FLIM
!
!          - Input -
!
! DZ       cell thickness
! WW       vertical velocity (negative downward)
! NI       number of profiles
! NK       number of grid point in profiles
! DT       timestep
! KF       number of filtering pass
! FLIM     index of hier grid point where to look for mass in sedimentation.
! WMAX     maximum vertical velocity (positive or negative value).
!          This is used to save computation time. For each
!          level, the index of the heighest level from which
!          mass can be received will be computed if COMPLIM is .true.
! IDIR           direction of index:
!                idir=0 if index one is at top of model
!                idir=1 if index one is at surface of model
! ktop_sedi      uppermost level below which sedimentation is computed
! COUNTER        number of columns in (ni,nk) slab with non-zero RO
! ACTIVECOLUMN   array of i-indices with non-zero columns in RO(ni,nk)

!#TODO: get LEVMAX from phy_options module
!     NOMBRE MAXIMUM DE NIVEAUX POUR LA PHYSIQUE
      integer LEVMAX
      parameter (LEVMAX = 200)
!--

! LOCAL VARIABLES AND PARAMETERS:

 integer               :: i,k,l,km1,kp1,ks,ke,kw,a
 real, dimension(nk)   :: vp,zt,zb,dzi,ro_star
 real, dimension(0:nk) :: zz
 real                  :: zmax,tempo
 integer               :: idzmin
 real, parameter       :: epsilon = 1.e-2

!---------------------------------------------------------------------------
!     Set parameters related to index direction.

!       ks=1
!       ke=nk
!       kw=1
!       if(idir.eq.0)then
!          ks=nk
!          ke=1
!          kw=-1
!       endif

  !For nk=bottom:  (hard-wired to remove pass of 'idir' parameter)
     ks=  nk
     ke=  ktop_sedi   !ke=1 (old)
     kw= -1

    !------------------------------------------------------------------------
     do a= 1,counter   !i=1,ni
        i=activeColumn(a)

!     Compute cell height and final position of the top (zt) and bottom (zb)
!     of moving boxes:

        zz(ks)=0.
        do k=ks,ke,kw
            zz (k+kw)=zz(k)+dz(i,k)
            dzi(k)=1./dz(i,k)
            vp (k)=0.
        enddo

        do k=ks,ke,kw
            zb(k)=zz(k)+ww(i,k)*dt
        enddo

!     Scheme='Girard' (not 'Kato')
        zt(ke)=zb(ke)+dz(i,ke)
        do k=ks,ke-kw,kw
               zb(k)=min(zb(k+kw)-epsilon*dz(i,k),zz(k)+ww(i,k)*dt)
               zt(k)=zb(k+kw)
        enddo

        do k=ks,ke,kw    !k=1,nk
            ro_star(k)= ro(i,k)*dz(i,k)/(zt(k)-zb(k))
        enddo

!      Compute limit index where to look for mass:
        if (complim) then
           zmax=abs(wmax*dt)
           do l=ks,ke,kw
              flim(l)=l
              do k=l,ke,kw
                  if (zmax.ge.zz(k)-zz(l+kw)) flim(l)=k
              enddo
           enddo
        endif

        do l=ks,ke,kw
          do k=l,flim(l),kw
               vp(l)= vp(l) + ro_star(k)*max( 0.,min(zz(l+kw),zt(k)) &
                         - max(zz(l),zb(k)) )
          enddo
        enddo

        do k=ks,ke,kw
            ro(i,k)=vp(k)*dzi(k)
        enddo


     enddo !i-loop
    !------------------------------------------------------------------------

 return
 end subroutine blg4sedi

!=====================================================================================!
!This subroutine is modified from S/P BLG.FTN

 subroutine blg5sedi (RO,DZ,WW,NK,DT,COMPLIM,WMAX,FLIM,counter,activeColumn,ktop_sedi)

 implicit none

!PASSING ARGUMENTS:

 integer, intent(in)                             :: NK,counter
 real, dimension(:,:), intent(inout)             :: ro
 real, dimension(:,:), intent(in)                :: dz,ww
 real, intent(in)                                :: dt,wmax
 integer, dimension(nk)                          :: FLIM
 integer, dimension(size(RO,dim=1)), intent(in)  :: activeColumn
 logical, intent(inout)                          :: COMPLIM
! integer, intent(in)                             :: ktop_sedi
 integer, dimension(size(RO,dim=1))              :: ktop_sedi

! Author
!          C. Girard and A. Plante (2001)
!
! Revisions
!
! 001      A. Glazer and A. Plante (Oct 2001)
!             - introduced complim and ind_limit in computation of precipitation
! 002      A. Plante (June 2003)
!             - IBM conversion, added ro_star. This imporved sedimentaion
!               in mixphase4 by 38%.
!             - change computation limit to insure reentrance in OpenMP
!             - blg2.ftn validates perfectly with blg.ftn of PHY 3.8.
! 003      J. Milbrandt (Nov 2004)
!             - Modified for use in Milbrandt-Yau microphysics scheme
!               + add condition for mass-in-column (activeColumn); for all i-loops, the
!                 line 'do i=1,ni' was replaced with 'do a=1,counter' & 'i=activeColumn(a)'
!               + remove RT calculation (and pass)
!               + hard-wired various options for multimoment.ftn90 (removed unnecessary IF statements)
! 004      J. Milbrandt (Jan 2007)
!             - corrected 'idzmin' initial value (for use in activeColumn i-loops)
!             - added option to exclude upper levels (ktop_sedi)
! 005      J. Milbrandt and R. McTaggart-Cowan (March 2007)
!             - removed i-dependency; cleaned up arrays, etc.  RO, DZ, WW are (ni,nk) arrays but there is
!               only one outer i (a)-loop (which could be moved outside of this subroutine, leaving it
!               as a single column subroutine).
! 006      J. Milbrandt (Nov 2008)
!             - changed 'ktop_sedi' from scalar to i-array
!
! Object
!
!
!   Version 1.0
!
!   CALCULATES
!                sedimentation of a quantity of matter ro(kg/m3)
!                falling with negative downward velocity ww(m/s)
!
!   ACCORDING TO
!                the BOX-LAGRANGIAN SCHEME
!                         (Kato,1995:JMSJ,73,241-245)
!             OR
!                the ADJACENT-LAGRANGIAN SCHEME
!                Girard and Plante
!
!   PLUS
!                a conservative two-grid-length filter
!                in order to control noise if needed
!
!Arguments
!
!          - Input/Output -
!
! RO       density in cell before and after sedimentation.
! COMPLIM  logical switch to force recomputing the table FLIM
!
!          - Input -
!
! DZ       cell thickness
! WW       vertical velocity (negative downward)
! NI       number of profiles
! NK       number of grid point in profiles
! DT       timestep
! KF       number of filtering pass
! FLIM     index of hier grid point where to look for mass in sedimentation.
! WMAX     maximum vertical velocity (positive or negative value).
!          This is used to save computation time. For each
!          level, the index of the heighest level from which
!          mass can be received will be computed if COMPLIM is .true.
! IDIR           direction of index:
!                idir=0 if index one is at top of model
!                idir=1 if index one is at surface of model
! ktop_sedi      uppermost level below which sedimentation is computed
! COUNTER        number of columns in (ni,nk) slab with non-zero RO
! ACTIVECOLUMN   array of i-indices with non-zero columns in RO(ni,nk)


! LOCAL VARIABLES AND PARAMETERS:

 integer                :: i,k,l,km1,kp1,ks,kw,a !,ke
 integer, dimension(size(RO,dim=1)) :: ke
 real, dimension(nk)    :: vp,zt,zb,dzi,ro_star
 real, dimension(0:nk)  :: zz
 real                   :: zmax,tempo
 integer                :: idzmin
 real, parameter        :: epsilon = 1.e-2

!---------------------------------------------------------------------------
!     Set parameters related to index direction.

!       ks=1
!       ke=nk
!       kw=1
!       if(idir.eq.0)then
!          ks=nk
!          ke=1
!          kw=-1
!       endif

  !For nk=bottom:  (hard-wired to remove pass of 'idir' parameter)
     ks=  nk
     ke=  ktop_sedi   !ke=1 (old)
     kw= -1

    !------------------------------------------------------------------------
     do a= 1,counter   !i=1,ni
        i=activeColumn(a)

!     Compute cell height and final position of the top (zt) and bottom (zb)
!     of moving boxes:

        zz(ks)=0.
        do k=ks,ke(i),kw
            zz (k+kw)=zz(k)+dz(i,k)
            dzi(k)=1./dz(i,k)
            vp (k)=0.
        enddo

        do k=ks,ke(i),kw
            zb(k)=zz(k)+ww(i,k)*dt
        enddo

!     Scheme='Girard' (not 'Kato')
        zt(ke(i))=zb(ke(i))+dz(i,ke(i))
        do k=ks,ke(i)-kw,kw
               zb(k)=min(zb(k+kw)-epsilon*dz(i,k),zz(k)+ww(i,k)*dt)
               zt(k)=zb(k+kw)
        enddo

        do k=ks,ke(i),kw    !k=1,nk
            ro_star(k)= ro(i,k)*dz(i,k)/(zt(k)-zb(k))
        enddo

!      Compute limit index where to look for mass:
        if (complim) then
           zmax=abs(wmax*dt)
           do l=ks,ke(i),kw
              flim(l)=l
              do k=l,ke(i),kw
                  if (zmax.ge.zz(k)-zz(l+kw)) flim(l)=k
              enddo
           enddo
        endif

        do l=ks,ke(i),kw
          do k=l,flim(l),kw
               vp(l)= vp(l) + ro_star(k)*max( 0.,min(zz(l+kw),zt(k)) &
                         - max(zz(l),zb(k)) )
          enddo
        enddo

        do k=ks,ke(i),kw
            ro(i,k)=vp(k)*dzi(k)
        enddo


     enddo !i-loop
    !------------------------------------------------------------------------

 return
 end subroutine blg5sedi
!=====================================================================================!

end module my_sedi_mod
