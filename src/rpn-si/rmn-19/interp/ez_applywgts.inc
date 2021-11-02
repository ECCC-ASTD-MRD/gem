   subroutine ez_applywgts(outfld, wts, idxs, infld, x, y, masque, ni_src, nj_src, ni_dst, nj_dst, n_wts)

   implicit none
   integer ni_src, nj_src, ni_dst, nj_dst, n_wts, n
   real :: x(ni_src,nj_src),y(ni_src,nj_src), dist, total_wgt
   real :: infld(ni_src*nj_src), outfld(ni_dst*nj_dst)
   real rmin, rmax
   integer i,j,k

   real :: wts(ni_dst, nj_dst, n_wts)
   integer :: idxs(ni_dst, nj_dst, n_wts), masque(ni_dst*nj_dst)
   integer :: ezgetval, ezgetopt, ier

   character(len=32) :: interopt, xtrapopt
   real xtrapval

   ier = ezgetopt('INTERP_DEGREE', interopt)
   ier = ezgetopt('EXTRAP_DEGREE', xtrapopt)
   ier = ezgetval('EXTRAP_VALUE', xtrapval)
!   print *, interopt, xtrapopt, xtrapval

   if (xtrapopt(1:5) == 'value') then
      ier = ezgetval('EXTRAP_VALUE', xtrapval)
      outfld = xtrapval
   else
      rmin = minval(infld)
      rmax = maxval(infld)
      rmin  = rmin - 0.1*(rmax-rmin)
      outfld = rmin
   endif


   do k=1,ni_dst*nj_dst
     i = mod((k-1),ni_dst) + 1
     j = 1+k/ni_dst
     if (masque(k) == 1) then
        outfld(k) = 0.0
        do n=1,n_wts
           if (idxs(i,j,n) < 1) exit
           outfld(k) = outfld(k)+wts(i,j,n)*infld(idxs(i,j,n))
        enddo
     endif
   enddo

   end subroutine ez_applywgts
