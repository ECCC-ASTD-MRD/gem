subroutine test_sfclayer
   use sfclayer, only:sl_put,sl_get,sl_prelim,sl_sfclayer,SL_OK
   implicit none
 
   integer, parameter :: NI=1,STDERR=0,STDOUT=6
   real :: as_value
   real, dimension(NI) :: t=300.,q=0.01,u=5.,v=10.,p0=1.e5,zm=40.,zt=20.,ts=298.,qs=0.012, &
        z0m=0.1,z0t=0.1,lat=.7,fcor=1e-4,spd,dir,tv,rho,ilmo,frv,tdiag
 
   ! Begin messaging
   write(STDOUT,*) 'Testing Surface Layer Module'

   ! Change the value of the inverse Prandtl number from 1 (default) to 0.85
   if ( sl_put('BETA',0.85) /= SL_OK ) write(STDERR,*) 'Error setting BETA'
 
   ! Retrieve the stability function coefficient for the stable case
   if ( sl_get('AS',as_value) /= SL_OK ) write(STDERR,*) 'Error retrieving AS'
 
   ! Preliminary calculations for inputs to the surface layer parameterization using helper function
   if ( sl_prelim(t,q,u,v,p0,zm,tv_air=tv,rho_air=rho,spd_air=spd,dir_air=dir) /= SL_OK ) &
        write(STDERR,*) 'Error returned by sl_prelim'
 
   ! Estimate selected surface layer properties using the surface layer parameterization
   if ( sl_sfclayer(t,q,spd,dir,zm,zt,ts,qs,z0m,z0t,lat,fcor,hghtt_diag=2.,ilmo=ilmo,ue=frv,t_diag=tdiag) /= SL_OK ) &
        write(STDERR,*) 'Error returned by sl_sfclayer'
 
   ! Test for expecte result
   if (close_enough(as_value,12.) .and. close_enough(1./ilmo(1),231.2455597) .and. &
        close_enough(frv(1),0.6263254) .and. close_enough(tdiag(1),299.0290833)) then
      write(STDOUT,*) ' * Pass'
   else
      write(STDOUT,*) ' * Fail'
   endif

contains
   
   function close_enough(val,expected) result(is_close)
      real, intent(in) :: val,expected
      logical :: is_close
      is_close = .true.
      if (abs(val-expected)/expected > 10.*epsilon(val)) is_close = .false.
   end function close_enough
 
end subroutine test_sfclayer
