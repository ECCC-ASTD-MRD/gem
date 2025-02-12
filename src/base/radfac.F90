
!*@/
subroutine radfac4(qoz,ozotoit,sig,nlev,nk,lref, &
     dlat,press,np,npmax, &
     o3f,x1,x2,x3,x4,i1,i2, &
     nlat,fozo,clat,pref)
   use phy_status, only: phy_error_L
   implicit none
!!!#include <arch_specific.hf>
   integer nlev,lref,np,npmax,nlat,nk

   real qoz(npmax,nlev), &
        sig(np,nk+1),dlat(np),press(np),o3f(npmax,lref), &
        x1(np),x2(np),x3(np),x4(np)
   real fozo(nlat,lref),clat(nlat),pref(lref)
   real ozotoit(npmax)
   integer i1(np),i2(np)

   !@Author  L.Garand RPN (June 1989)
   !@Revision
   ! 001      see version 5.5.0 for previous history
   !@Object
   !          to calculate ozone mixing ratio at model levels
   !@Arguments
   !          - Output -
   ! qoz      ozone mixing ratio (kg o3/kg air) for each
   !          sigma level
   ! ozotoit  total ozone (cm stp) above model roof
   !          - input -
   ! sig      sigma levels at the centre of the layers
   ! nlev     number of flux levels
   ! nk       number of layers
   ! lref     number of ozone climatological levels
   ! dlat     latitude of np points to process in radians
   ! press    np points of surface pressure
   ! np       number of points to process
   ! npmax    maximum number of points allowed
   ! maxlev   number of maximum flux levels in the model
   !          - output -
   ! o3f      ozone (kg o3/kg air) at  each reference
   !          climatological level for dlat latitudes
   !          - input -
   ! x1       work field
   ! x2       work field
   ! x3       work field
   ! x4       work field
   ! i1       work field
   ! i2       work field
   ! s2       work field
   ! s3       work field
   ! nlat     number of climatological latitudes
   ! fozo     ozone climatological field in ppmv
   ! clat     ozone climatological latitudes
   ! pref     ozone climatological pressures
   !@notes
   !          this routine calls:
   !          ozoref2 (kg o3/kg air at climatological levels)
   !          qozon3  (kg o3/kg air at desired sigma  levels)
   !*@/

   call ozoref3(o3f,lref,dlat,np,npmax,i1,nlat,clat,fozo)
   if (phy_error_L) return

   call qozon3(qoz,ozotoit,o3f,press,sig,nlev,nk,np,npmax, &
        lref,pref,x1,x2,x3,x4,i1,i2)

   return
end subroutine radfac4
